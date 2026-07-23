"""구조 탐색 단계의 순서, 캐시 및 fallback을 조정합니다."""

import argparse
import logging
from pathlib import Path
from typing import Any

from .alphafold import (
    find_exact_alphafold_model,
    save_alphafold_result,
)
from .cache import (
    cache_entry_dir,
    default_cache_dir,
    lookup_cache,
    materialize_cache_hit,
    store_cache,
)
from .config import MIN_RCSB_SEQUENCE_LENGTH
from .esmfold2 import run_esmfold2
from .http_client import make_session
from .metadata import (
    complete_search_trace,
    decorate_result,
    stage_record,
)
from .rcsb import (
    find_exact_experimental_structures,
    rank_experimental_candidates,
    save_experimental_result,
)
from .sequence import load_sequence, sequence_sha256
from .storage import utc_now


LOGGER = logging.getLogger(__name__)


def _cache_settings(args: argparse.Namespace) -> tuple[bool, Path, bool]:
    enabled = not bool(getattr(args, "no_cache", False))
    configured = getattr(args, "cache_dir", None)
    cache_dir = (
        configured.expanduser().resolve()
        if configured is not None
        else default_cache_dir()
    )
    refresh = bool(getattr(args, "refresh_cache", False))
    return enabled, cache_dir, refresh


def _decorate_and_cache(
    metadata: dict[str, Any],
    sequence: str,
    *,
    source_class: str | None,
    search_trace: list[dict[str, Any]],
    output: Path,
    cache_enabled: bool,
    cache_dir: Path,
    cache_miss_reason: str | None,
) -> dict[str, Any]:
    cache_key = sequence_sha256(sequence)
    result = decorate_result(
        metadata,
        sequence,
        source_class=source_class,
        search_trace=search_trace,
        cache_enabled=cache_enabled,
        cache_key=cache_key,
        cache_miss_reason=cache_miss_reason,
    )
    if cache_enabled and result.get("status") == "success":
        try:
            store_cache(
                cache_dir,
                sequence,
                result,
                output / str(result["structure_file"]),
            )
            result["cache"]["write_status"] = "stored"
            result["cache"]["entry"] = str(
                cache_entry_dir(cache_dir, sequence)
            )
        except (OSError, ValueError) as exc:
            result["cache"]["write_status"] = "failed"
            result["cache"]["write_error"] = str(exc)
            LOGGER.warning("캐시 저장 실패: %s", exc)
    return result


def resolve(
    args: argparse.Namespace,
    *,
    output: Path | None = None,
) -> dict[str, Any]:
    sequence = load_sequence(args)
    output = output or args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    cache_enabled, cache_dir, refresh_cache = _cache_settings(args)
    cache_miss_reason: str | None = None

    if cache_enabled and not refresh_cache:
        LOGGER.info("SHA-256 캐시 확인")
        lookup = lookup_cache(cache_dir, sequence, args)
        if lookup.metadata is not None:
            LOGGER.info("캐시 적중: %s", sequence_sha256(sequence))
            return materialize_cache_hit(lookup, output, sequence)
        cache_miss_reason = lookup.reason
        LOGGER.info("캐시 미적중: %s", lookup.reason)
    elif refresh_cache:
        cache_miss_reason = "refresh_requested"

    trace: list[dict[str, Any]] = []
    sequence_search_supported = len(sequence) >= MIN_RCSB_SEQUENCE_LENGTH

    with make_session() as session:
        if sequence_search_supported:
            LOGGER.info("RCSB PDB exact 실험 구조 검색")
            experimental, experimental_search = (
                find_exact_experimental_structures(session, sequence)
            )
            trace.append(
                stage_record(
                    "rcsb_experimental",
                    attempted=True,
                    status="found" if experimental else "not_found",
                    summary=experimental_search,
                )
            )
        else:
            experimental = []
            trace.append(
                stage_record(
                    "rcsb_experimental",
                    attempted=False,
                    status="skipped",
                    reason=(
                        "RCSB sequence search requires at least "
                        f"{MIN_RCSB_SEQUENCE_LENGTH} residues"
                    ),
                )
            )

        if experimental:
            selected, candidate_records = rank_experimental_candidates(
                experimental
            )
            trace = complete_search_trace(
                trace, stop_reason="experimental_exact_structure_selected"
            )
            metadata = save_experimental_result(
                session,
                selected,
                candidate_records,
                experimental_search,
                sequence,
                output,
            )
            return _decorate_and_cache(
                metadata,
                sequence,
                source_class="experimental",
                search_trace=trace,
                output=output,
                cache_enabled=cache_enabled,
                cache_dir=cache_dir,
                cache_miss_reason=cache_miss_reason,
            )

        computational_search: dict[str, int] | None = None
        if not sequence_search_supported:
            trace.append(
                stage_record(
                    "alphafold_db",
                    attempted=False,
                    status="skipped",
                    reason=(
                        "RCSB computational sequence search requires at "
                        f"least {MIN_RCSB_SEQUENCE_LENGTH} residues"
                    ),
                )
            )
            model = None
            computational_candidates: list[dict[str, Any]] = []
        elif args.skip_computational:
            trace.append(
                stage_record(
                    "alphafold_db",
                    attempted=False,
                    status="skipped",
                    reason="disabled_by_cli",
                )
            )
            model = None
            computational_candidates = []
        else:
            LOGGER.info("AlphaFold DB exact 계산 구조 검색")
            (
                model,
                computational_search,
                computational_candidates,
            ) = find_exact_alphafold_model(session, sequence)
            trace.append(
                stage_record(
                    "alphafold_db",
                    attempted=True,
                    status="found" if model else "not_found",
                    summary=computational_search,
                )
            )

        if model:
            trace = complete_search_trace(
                trace, stop_reason="alphafold_exact_model_selected"
            )
            metadata = save_alphafold_result(
                session,
                model,
                computational_search or {},
                computational_candidates,
                sequence,
                output,
            )
            return _decorate_and_cache(
                metadata,
                sequence,
                source_class="precomputed_prediction",
                search_trace=trace,
                output=output,
                cache_enabled=cache_enabled,
                cache_dir=cache_dir,
                cache_miss_reason=cache_miss_reason,
            )

        if args.skip_prediction:
            trace.append(
                stage_record(
                    "esmfold2",
                    attempted=False,
                    status="skipped",
                    reason="disabled_by_cli",
                )
            )
            metadata = {
                "status": "not_found",
                "source": None,
                "message": (
                    "global exact 구조를 찾지 못했고 예측은 생략했습니다."
                ),
                "retrieved_at": utc_now(),
            }
            return _decorate_and_cache(
                metadata,
                sequence,
                source_class=None,
                search_trace=trace,
                output=output,
                cache_enabled=cache_enabled,
                cache_dir=cache_dir,
                cache_miss_reason=cache_miss_reason,
            )

        LOGGER.info("ESMFold2 신규 예측")
        metadata = run_esmfold2(sequence, output)
        trace.append(
            stage_record(
                "esmfold2",
                attempted=True,
                status="success",
                summary={
                    "model": metadata.get("model"),
                    "mean_plddt": metadata.get("mean_plddt"),
                    "ptm": metadata.get("ptm"),
                },
            )
        )
        return _decorate_and_cache(
            metadata,
            sequence,
            source_class="new_prediction",
            search_trace=trace,
            output=output,
            cache_enabled=cache_enabled,
            cache_dir=cache_dir,
            cache_miss_reason=cache_miss_reason,
        )
