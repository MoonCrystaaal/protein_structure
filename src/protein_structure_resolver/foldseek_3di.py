"""Foldseek를 이용해 구조 파일에서 아미노산 및 3Di 서열을 추출합니다."""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import gemmi

from .errors import FoldseekExtractionError
from .storage import atomic_write_text, utc_now


FOLDSEEK_TIMEOUT_SECONDS = 300
AA_FASTA_NAME = "structure_aa.fasta"
THREE_DI_FASTA_NAME = "structure_3di.fasta"


@dataclass(frozen=True)
class FoldseekRunner:
    command: tuple[str, ...]
    mode: str
    executable: str

    def path_argument(self, path: Path) -> str:
        resolved = path.resolve()
        if self.mode != "wsl":
            return str(resolved)
        drive = resolved.drive
        if not drive or len(drive) != 2 or drive[1] != ":":
            raise FoldseekExtractionError(
                "WSL Foldseek에서는 Windows 로컬 드라이브 경로만 지원합니다: "
                f"{resolved}"
            )
        relative = resolved.as_posix()[3:]
        return f"/mnt/{drive[0].lower()}/{relative}"


def _run(
    command: list[str],
    *,
    allow_failure: bool = False,
) -> subprocess.CompletedProcess[str]:
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=FOLDSEEK_TIMEOUT_SECONDS,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise FoldseekExtractionError(
            f"Foldseek 실행을 시작하거나 완료하지 못했습니다: {exc}"
        ) from exc
    if completed.returncode != 0 and not allow_failure:
        detail = (completed.stderr or completed.stdout).strip()
        if len(detail) > 1200:
            detail = detail[-1200:]
        raise FoldseekExtractionError(
            "Foldseek 명령이 실패했습니다"
            + (f": {detail}" if detail else ".")
        )
    return completed


def _find_wsl_foldseek(requested: str | None) -> FoldseekRunner | None:
    if os.name != "nt" or shutil.which("wsl.exe") is None:
        return None
    candidate = requested or "foldseek"
    if candidate.startswith("/"):
        discovered = candidate
    else:
        # 로그인 셸의 PATH에는 사용자가 설치한 ~/.local/bin도 포함됩니다.
        probe = _run(
            [
                "wsl.exe",
                "sh",
                "-lc",
                f"command -v -- {shlex.quote(candidate)}",
            ],
            allow_failure=True,
        )
        discovered = probe.stdout.strip().splitlines()
        if probe.returncode != 0 or not discovered:
            return None
        discovered = discovered[-1]
    return FoldseekRunner(
        command=("wsl.exe", str(discovered)),
        mode="wsl",
        executable=str(discovered),
    )


def find_foldseek_runner(requested: str | None = None) -> FoldseekRunner:
    """네이티브 Foldseek를 우선하고 Windows에서는 WSL도 자동 탐색합니다."""
    if requested:
        native = shutil.which(requested)
        if native is None and Path(requested).is_file():
            native = str(Path(requested).resolve())
    else:
        native = shutil.which("foldseek")
    if native:
        return FoldseekRunner(
            command=(native,), mode="native", executable=native
        )
    wsl = _find_wsl_foldseek(requested)
    if wsl is not None:
        return wsl
    suffix = f" ({requested})" if requested else ""
    raise FoldseekExtractionError(
        "Foldseek 실행 파일을 찾지 못했습니다"
        f"{suffix}. Foldseek를 설치하거나 --foldseek-bin으로 지정하세요."
    )


def _foldseek_version(runner: FoldseekRunner) -> str | None:
    # Foldseek는 별도의 version 옵션 대신 인자 없는 도움말에 빌드 ID를
    # 표시하고 종료 코드 1을 반환하는 버전이 있습니다.
    completed = _run([*runner.command], allow_failure=True)
    banner = f"{completed.stdout}\n{completed.stderr}"
    match = re.search(r"(?:Foldseek|MMseqs) Version:\s*([^\s]+)", banner)
    if match:
        return match.group(1)
    match = re.search(r"Version:\s*([^\s]+)", banner)
    return match.group(1) if match else None


def _normalize_esmfold2_structure(source: Path, destination: Path) -> str:
    """ESMFold2 mmCIF를 Foldseek가 안정적으로 읽는 임시 PDB로 변환합니다."""
    try:
        structure = gemmi.read_structure(str(source))
        structure.setup_entities()
        structure.assign_label_seq_id()
        pdb_text = structure.make_pdb_string()
    except (RuntimeError, ValueError) as exc:
        raise FoldseekExtractionError(
            f"Gemmi가 ESMFold2 구조를 정규화하지 못했습니다: {exc}"
        ) from exc
    if not any(
        line.startswith(("ATOM  ", "HETATM"))
        for line in pdb_text.splitlines()
    ):
        raise FoldseekExtractionError(
            "Gemmi 정규화 결과에 Foldseek가 읽을 원자 좌표가 없습니다."
        )
    atomic_write_text(destination, pdb_text)
    return gemmi.__version__


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        elif header is None:
            raise FoldseekExtractionError(
                f"Foldseek FASTA 출력의 첫 레코드에 헤더가 없습니다: {path.name}"
            )
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    if not records or any(not header or not sequence for header, sequence in records):
        raise FoldseekExtractionError(
            f"Foldseek가 유효한 FASTA 레코드를 만들지 못했습니다: {path.name}"
        )
    return records


def _validate_paired_records(
    aa_records: list[tuple[str, str]],
    three_di_records: list[tuple[str, str]],
) -> list[dict[str, Any]]:
    if len(aa_records) != len(three_di_records):
        raise FoldseekExtractionError(
            "아미노산 서열과 3Di 서열의 레코드 수가 다릅니다."
        )
    summaries: list[dict[str, Any]] = []
    for (aa_id, aa_sequence), (three_di_id, three_di_sequence) in zip(
        aa_records, three_di_records, strict=True
    ):
        if aa_id != three_di_id:
            raise FoldseekExtractionError(
                "아미노산 서열과 3Di 서열의 레코드 ID가 다릅니다: "
                f"{aa_id!r} != {three_di_id!r}"
            )
        if len(aa_sequence) != len(three_di_sequence):
            raise FoldseekExtractionError(
                f"{aa_id!r}의 아미노산/3Di 길이가 다릅니다: "
                f"{len(aa_sequence)} != {len(three_di_sequence)}"
            )
        summaries.append(
            {
                "id": aa_id,
                "residue_count": len(aa_sequence),
                "length_match": True,
            }
        )
    return summaries


def extract_3di(
    structure_path: Path,
    output: Path,
    *,
    source: str | None,
    foldseek_bin: str | None = None,
) -> dict[str, Any]:
    """구조를 Foldseek DB로 변환한 뒤 AA와 3Di FASTA만 보존합니다."""
    runner = find_foldseek_runner(foldseek_bin)
    output.mkdir(parents=True, exist_ok=True)
    aa_output = output / AA_FASTA_NAME
    three_di_output = output / THREE_DI_FASTA_NAME

    normalization: dict[str, Any] = {
        "applied": False,
        "reason": "source_mmcif_used_directly",
    }
    with tempfile.TemporaryDirectory(prefix=".foldseek-3di-", dir=output) as tmp:
        workspace = Path(tmp)
        foldseek_input = structure_path
        if source == "ESMFold2":
            foldseek_input = workspace / "structure.normalized.pdb"
            gemmi_version = _normalize_esmfold2_structure(
                structure_path, foldseek_input
            )
            normalization = {
                "applied": True,
                "reason": "esmfold2_mmcif_compatibility",
                "tool": "Gemmi",
                "tool_version": gemmi_version,
                "temporary_format": "PDB",
            }

        database = workspace / "structure_db"
        aa_temporary = workspace / AA_FASTA_NAME
        three_di_temporary = workspace / THREE_DI_FASTA_NAME
        _run(
            [
                *runner.command,
                "createdb",
                runner.path_argument(foldseek_input),
                runner.path_argument(database),
            ]
        )
        _run(
            [
                *runner.command,
                "convert2fasta",
                runner.path_argument(database),
                runner.path_argument(aa_temporary),
            ]
        )
        # 3Di DB에는 자체 header DB가 없으므로 AA header DB를 연결합니다.
        _run(
            [
                *runner.command,
                "lndb",
                runner.path_argument(Path(f"{database}_h")),
                runner.path_argument(Path(f"{database}_ss_h")),
            ]
        )
        _run(
            [
                *runner.command,
                "convert2fasta",
                runner.path_argument(Path(f"{database}_ss")),
                runner.path_argument(three_di_temporary),
            ]
        )

        aa_records = _parse_fasta(aa_temporary)
        three_di_records = _parse_fasta(three_di_temporary)
        records = _validate_paired_records(aa_records, three_di_records)
        atomic_write_text(
            aa_output, aa_temporary.read_text(encoding="utf-8")
        )
        atomic_write_text(
            three_di_output, three_di_temporary.read_text(encoding="utf-8")
        )

    return {
        "status": "success",
        "tool": "Foldseek",
        "tool_version": _foldseek_version(runner),
        "runner": runner.mode,
        "alphabet": "3Di",
        "similarity_search_performed": False,
        "input_structure_file": structure_path.name,
        "normalization": normalization,
        "amino_acid_fasta_file": aa_output.name,
        "three_di_fasta_file": three_di_output.name,
        "record_count": len(records),
        "records": records,
        "created_at_utc": utc_now(),
    }
