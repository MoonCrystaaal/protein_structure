"""Biohub Platform을 통한 ESMFold2 신규 예측."""

import os
from math import isfinite
from pathlib import Path
from statistics import fmean
from typing import Any, Sequence

from .config import ESMFOLD2_MODEL
from .errors import PredictionUnavailableError
from .sequence import sequence_sha256
from .storage import atomic_write_text, utc_now, validate_mmcif


THREE_TO_ONE = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def _flatten_numbers(values: Any) -> list[float]:
    if hasattr(values, "detach"):
        values = values.detach().cpu()
    if hasattr(values, "tolist"):
        values = values.tolist()
    if isinstance(values, (list, tuple)):
        result: list[float] = []
        for value in values:
            result.extend(_flatten_numbers(value))
        return result
    try:
        return [float(values)]
    except (TypeError, ValueError):
        return []


def normalize_plddt(values: Any) -> list[float]:
    numbers = _flatten_numbers(values)
    if numbers and max(numbers) <= 1.0:
        numbers = [value * 100.0 for value in numbers]
    return numbers


def summarize_plddt(values: Sequence[float]) -> dict[str, Any]:
    if not values:
        return {
            "mean_plddt": None,
            "min_plddt": None,
            "fraction_very_high": None,
            "fraction_confident_or_higher": None,
            "fraction_very_low": None,
        }
    count = len(values)
    return {
        "mean_plddt": round(fmean(values), 3),
        "min_plddt": round(min(values), 3),
        "fraction_very_high": round(
            sum(value >= 90 for value in values) / count, 6
        ),
        "fraction_confident_or_higher": round(
            sum(value >= 70 for value in values) / count, 6
        ),
        "fraction_very_low": round(
            sum(value < 50 for value in values) / count, 6
        ),
    }


def find_low_confidence_regions(
    values: Sequence[float], threshold: float = 70.0
) -> list[dict[str, Any]]:
    regions: list[dict[str, Any]] = []
    start: int | None = None
    for index, value in enumerate([*values, threshold], start=1):
        if value < threshold and start is None:
            start = index
        elif value >= threshold and start is not None:
            region_values = values[start - 1 : index - 1]
            regions.append(
                {
                    "start": start,
                    "end": index - 1,
                    "length": index - start,
                    "mean_plddt": round(fmean(region_values), 3),
                }
            )
            start = None
    return regions


def validate_prediction_result(
    result: Any, sequence: str
) -> tuple[list[float], list[float], dict[str, Any]]:
    complex_result = getattr(result, "complex", None)
    if complex_result is None:
        raise ValueError("API 응답에 molecular complex가 없습니다.")

    residue_tokens = list(getattr(complex_result, "sequence", []) or [])
    predicted_sequence = "".join(
        THREE_TO_ONE.get(str(token).upper(), "X")
        for token in residue_tokens
    )
    if predicted_sequence != sequence:
        raise ValueError(
            "예측 결과의 잔기 서열이 입력 서열과 일치하지 않습니다."
        )

    plddt = normalize_plddt(result.plddt)
    if len(plddt) != len(sequence):
        raise ValueError(
            "pLDDT 개수가 입력 서열 길이와 일치하지 않습니다: "
            f"{len(plddt)} != {len(sequence)}"
        )
    if not all(isfinite(value) and 0.0 <= value <= 100.0 for value in plddt):
        raise ValueError("pLDDT에 유효 범위를 벗어난 값이 있습니다.")

    ptm_values = _flatten_numbers(result.ptm)
    if not ptm_values or not isfinite(ptm_values[0]):
        raise ValueError("유효한 pTM 값이 없습니다.")
    if not 0.0 <= ptm_values[0] <= 1.0:
        raise ValueError("pTM 값이 0~1 범위를 벗어났습니다.")

    coordinates = _flatten_numbers(
        getattr(complex_result, "atom_positions", [])
    )
    if not coordinates or not all(isfinite(value) for value in coordinates):
        raise ValueError("원자 좌표가 비어 있거나 NaN/무한대가 포함되어 있습니다.")

    return plddt, ptm_values, {
        "sequence_exact": True,
        "residue_count": len(residue_tokens),
        "plddt_count": len(plddt),
        "coordinate_value_count": len(coordinates),
        "coordinates_finite": True,
    }


def run_esmfold2(sequence: str, output: Path) -> dict[str, Any]:
    token = os.environ.get("BIOHUB_API_TOKEN")
    if not token:
        raise PredictionUnavailableError(
            "BIOHUB_API_TOKEN 환경변수가 없어 ESMFold2를 실행할 수 없습니다."
        )
    try:
        from esm.sdk.api import FoldingConfig
        from esm.sdk.forge import SequenceStructureForgeInferenceClient
        from esm.utils.structure.input_builder import (
            ProteinInput,
            StructurePredictionInput,
        )
    except ImportError as exc:
        raise PredictionUnavailableError(
            "ESMFold2 SDK가 설치되지 않았습니다. "
            "requirements-esmfold2.txt의 설치 명령을 확인하세요."
        ) from exc

    try:
        client = SequenceStructureForgeInferenceClient(
            model=ESMFOLD2_MODEL,
            url="https://biohub.ai",
            token=token,
        )
        prediction_input = StructurePredictionInput(
            sequences=[ProteinInput(id="A", sequence=sequence)]
        )
        result = client.fold_all_atom(
            prediction_input,
            config=FoldingConfig(num_loops=3, num_sampling_steps=32),
        )
        plddt, ptm_values, prediction_validation = (
            validate_prediction_result(result, sequence)
        )
        cif_text = result.complex.to_mmcif()
        validate_mmcif(cif_text)
    except Exception as exc:
        raise PredictionUnavailableError(
            f"ESMFold2 예측에 실패했습니다: {exc}"
        ) from exc

    structure_path = output / "structure.cif"
    atomic_write_text(structure_path, cif_text)
    return {
        "status": "success",
        "source": "ESMFold2",
        "experimental": False,
        "precomputed": False,
        "input": {
            "sequence_length": len(sequence),
            "sequence_sha256": sequence_sha256(sequence),
        },
        "model": ESMFOLD2_MODEL,
        **summarize_plddt(plddt),
        "ptm": round(ptm_values[0], 6) if ptm_values else None,
        "low_confidence_regions": find_low_confidence_regions(plddt),
        "prediction_validation": prediction_validation,
        "config": {"num_loops": 3, "num_sampling_steps": 32},
        "structure_file": structure_path.name,
        "retrieved_at": utc_now(),
    }
