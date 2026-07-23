"""RCSB에 색인된 AlphaFold DB 모델의 exact 검색과 저장."""

from math import isfinite
from pathlib import Path
from typing import Any

import requests

from .config import (
    ALPHAFOLD_API_URL,
    ALPHAFOLD_RANKING_RULE_VERSION,
    COMPUTATIONAL_ENTITY_QUERY,
)
from .errors import DatabaseUnavailableError
from .http_client import request_with_retry
from .rcsb import fetch_entities_graphql, search_rcsb_entities
from .sequence import is_global_exact, sequence_sha256
from .storage import atomic_write_text, download_text, utc_now


def _uniprot_accessions(entity: dict[str, Any]) -> list[str]:
    identifiers = (
        entity.get("rcsb_polymer_entity_container_identifiers") or {}
    ).get("reference_sequence_identifiers") or []
    return [
        item["database_accession"]
        for item in identifiers
        if item.get("database_name") == "UniProt"
        and item.get("database_accession")
    ]


def fetch_alphafold_prediction(
    session: requests.Session, accession: str
) -> list[dict[str, Any]]:
    response = request_with_retry(
        session, "GET", f"{ALPHAFOLD_API_URL}/{accession}"
    )
    if response.status_code == 404:
        return []
    try:
        response.raise_for_status()
        body = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise DatabaseUnavailableError(
            f"AlphaFold DB 응답을 해석할 수 없습니다: {accession}"
        ) from exc
    return body if isinstance(body, list) else [body]


def _model_id(model: dict[str, Any]) -> str:
    return str(model.get("modelEntityId") or model.get("entryId") or "")


def _model_metric(model: dict[str, Any]) -> float:
    try:
        value = float(model.get("globalMetricValue"))
    except (TypeError, ValueError):
        return -1.0
    return value if isfinite(value) else -1.0


def _model_version(model: dict[str, Any]) -> int:
    try:
        return int(model.get("latestVersion"))
    except (TypeError, ValueError):
        return -1


def rank_alphafold_models(
    models: list[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not models:
        raise ValueError("선택할 AlphaFold DB 모델이 없습니다.")

    indexed = list(enumerate(models))
    decisions: dict[int, dict[str, str]] = {}

    best_metric = max(_model_metric(model) for _, model in indexed)
    metric_winners = [
        (index, model)
        for index, model in indexed
        if _model_metric(model) == best_metric
    ]
    for index, model in indexed:
        if _model_metric(model) < best_metric:
            decisions[index] = {
                "stage": "mean_plddt",
                "reason": "lower_mean_plddt",
            }

    latest_version = max(
        _model_version(model) for _, model in metric_winners
    )
    version_winners = [
        (index, model)
        for index, model in metric_winners
        if _model_version(model) == latest_version
    ]
    for index, model in metric_winners:
        if _model_version(model) < latest_version:
            decisions[index] = {
                "stage": "model_version",
                "reason": "older_model_version",
            }

    selected_index, selected = min(
        version_winners,
        key=lambda item: (_model_id(item[1]), item[0]),
    )
    for index, _ in version_winners:
        if index != selected_index:
            decisions[index] = {
                "stage": "deterministic_tiebreak",
                "reason": "higher_or_duplicate_model_identifier",
            }
    decisions[selected_index] = {
        "stage": "selected",
        "reason": "selected_by_alphafold_ranking_rule",
    }

    records = []
    for index, model in indexed:
        records.append(
            {
                "uniprot_accession": model.get("uniprotAccession"),
                "model_id": _model_id(model),
                "model_version": model.get("latestVersion"),
                "model_created_date": model.get("modelCreatedDate"),
                "mean_plddt": model.get("globalMetricValue"),
                "selected": index == selected_index,
                "decision": decisions[index],
            }
        )
    return selected, records


def find_exact_alphafold_model(
    session: requests.Session, sequence: str
) -> tuple[
    dict[str, Any] | None,
    dict[str, int],
    list[dict[str, Any]],
]:
    identifiers, raw_count = search_rcsb_entities(
        session, sequence, "computational"
    )
    entities = fetch_entities_graphql(
        session, identifiers, COMPUTATIONAL_ENTITY_QUERY
    )
    exact_entities = [
        entity
        for entity in entities
        if is_global_exact(
            sequence,
            (entity.get("entity_poly") or {}).get(
                "pdbx_seq_one_letter_code_can"
            ),
        )
        and (
            (
                (entity.get("entry") or {}).get(
                    "rcsb_comp_model_provenance"
                )
                or {}
            ).get("source_db")
            == "AlphaFoldDB"
        )
    ]

    models: list[dict[str, Any]] = []
    seen_accessions: set[str] = set()
    for entity in exact_entities:
        for accession in _uniprot_accessions(entity):
            if accession in seen_accessions:
                continue
            seen_accessions.add(accession)
            for model in fetch_alphafold_prediction(session, accession):
                if (
                    is_global_exact(sequence, model.get("sequence"))
                    and model.get("cifUrl")
                ):
                    models.append(model)

    selected, candidate_records = (
        rank_alphafold_models(models) if models else (None, [])
    )
    return selected, {
        "raw_candidate_count": raw_count,
        "prefiltered_candidate_count": len(identifiers),
        "exact_alphafold_entity_count": len(exact_entities),
        "verified_alphafold_model_count": len(models),
    }, candidate_records


def save_alphafold_result(
    session: requests.Session,
    model: dict[str, Any],
    search_summary: dict[str, int],
    candidate_records: list[dict[str, Any]],
    sequence: str,
    output: Path,
) -> dict[str, Any]:
    structure_path = output / "structure.cif"
    atomic_write_text(
        structure_path, download_text(session, str(model["cifUrl"]))
    )
    return {
        "status": "success",
        "source": "AlphaFold_DB",
        "experimental": False,
        "precomputed": True,
        "input": {
            "sequence_length": len(sequence),
            "sequence_sha256": sequence_sha256(sequence),
        },
        "search": search_summary,
        "uniprot_accession": model.get("uniprotAccession"),
        "model_id": model.get("modelEntityId") or model.get("entryId"),
        "model_version": model.get("latestVersion"),
        "model_created_date": model.get("modelCreatedDate"),
        "mean_plddt": model.get("globalMetricValue"),
        "plddt_url": model.get("plddtDocUrl"),
        "pae_url": model.get("paeDocUrl"),
        "structure_file": structure_path.name,
        "selection": {
            "ranking_rule_version": ALPHAFOLD_RANKING_RULE_VERSION,
            "rule": (
                "global exact -> maximum mean pLDDT -> latest model "
                "version -> ascending model ID"
            ),
            "selection_reason": "selected_by_alphafold_ranking_rule",
        },
        "candidates": candidate_records,
        "retrieved_at": utc_now(),
    }
