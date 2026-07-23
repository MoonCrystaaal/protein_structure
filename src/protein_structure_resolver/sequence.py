"""입력 서열 정규화, 검증 및 식별자 생성."""

import hashlib
from pathlib import Path
from typing import Protocol

from .config import STANDARD_AA
from .errors import SequenceInputError


class SequenceArguments(Protocol):
    sequence: str | None
    fasta: Path | None


def normalize_sequence(raw: str) -> str:
    """원시 서열 또는 단일 레코드 FASTA 문자열을 정규화합니다."""
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    headers = [index for index, line in enumerate(lines) if line.startswith(">")]
    if len(headers) > 1:
        raise SequenceInputError("한 번에는 단일 FASTA 레코드만 지원합니다.")
    if headers and headers[0] != 0:
        raise SequenceInputError("FASTA 헤더는 첫 번째 줄에 있어야 합니다.")

    sequence_lines = lines[1:] if headers else lines
    if any(line.startswith(">") for line in sequence_lines):
        raise SequenceInputError("한 번에는 단일 FASTA 레코드만 지원합니다.")

    sequence = "".join("".join(sequence_lines).split()).upper()
    if sequence.endswith("*"):
        sequence = sequence[:-1]
    return sequence


def validate_sequence(sequence: str) -> None:
    if not sequence:
        raise SequenceInputError("입력 서열이 비어 있습니다.")
    invalid = [
        (index + 1, residue)
        for index, residue in enumerate(sequence)
        if residue not in STANDARD_AA
    ]
    if invalid:
        position, residue = invalid[0]
        raise SequenceInputError(
            f"지원하지 않는 잔기 '{residue}'가 {position}번 위치에 있습니다."
        )


def load_sequence(args: SequenceArguments) -> str:
    if args.sequence is not None:
        raw = args.sequence
    else:
        if args.fasta is None:
            raise SequenceInputError("FASTA 파일 경로가 없습니다.")
        try:
            raw = args.fasta.read_text(encoding="utf-8-sig")
        except OSError as exc:
            raise SequenceInputError(
                f"FASTA 파일을 읽을 수 없습니다: {args.fasta}"
            ) from exc
    sequence = normalize_sequence(raw)
    validate_sequence(sequence)
    return sequence


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def normalize_rcsb_sequence(sequence: str | None) -> str:
    return "".join((sequence or "").split()).upper()


def is_global_exact(query: str, target: str | None) -> bool:
    normalized = normalize_rcsb_sequence(target)
    return len(query) == len(normalized) and query == normalized
