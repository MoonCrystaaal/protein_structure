"""모든 결과 유형이 공유하는 metadata 스키마."""

from copy import deepcopy
from typing import Any, Sequence

from .config import (
    CACHE_FORMAT_VERSION,
    METADATA_SCHEMA_VERSION,
    RESOLUTION_ORDER,
    RESOLVER_VERSION,
)
from .sequence import sequence_sha256
from .storage import utc_now


def stage_record(
    stage: str,
    *,
    attempted: bool,
    status: str,
    summary: dict[str, Any] | None = None,
    reason: str | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "stage": stage,
        "attempted": attempted,
        "status": status,
    }
    if summary is not None:
        record["summary"] = summary
    if reason is not None:
        record["reason"] = reason
    return record


def complete_search_trace(
    records: Sequence[dict[str, Any]],
    *,
    stop_reason: str,
) -> list[dict[str, Any]]:
    result = [deepcopy(record) for record in records]
    recorded = {record["stage"] for record in result}
    for stage in RESOLUTION_ORDER:
        if stage not in recorded:
            result.append(
                stage_record(
                    stage,
                    attempted=False,
                    status="not_needed",
                    reason=stop_reason,
                )
            )
    result.sort(key=lambda item: RESOLUTION_ORDER.index(item["stage"]))
    return result


def decorate_result(
    metadata: dict[str, Any],
    sequence: str,
    *,
    source_class: str | None,
    search_trace: Sequence[dict[str, Any]],
    cache_enabled: bool,
    cache_key: str,
    cache_miss_reason: str | None = None,
) -> dict[str, Any]:
    cache: dict[str, Any] = {
        "enabled": cache_enabled,
        "hit": False,
        "key": cache_key,
        "format_version": CACHE_FORMAT_VERSION,
    }
    if cache_miss_reason:
        cache["miss_reason"] = cache_miss_reason
    return {
        "schema_version": METADATA_SCHEMA_VERSION,
        "resolver_version": RESOLVER_VERSION,
        **metadata,
        "source_class": source_class,
        "input": {
            "sequence_length": len(sequence),
            "sequence_sha256": sequence_sha256(sequence),
        },
        "pipeline": {
            "resolution_order": list(RESOLUTION_ORDER),
            "search_trace": [deepcopy(record) for record in search_trace],
        },
        "cache": cache,
        "created_at_utc": utc_now(),
    }


def mark_cache_hit(
    metadata: dict[str, Any],
    *,
    cache_key: str,
    stored_at_utc: str | None,
) -> dict[str, Any]:
    result = deepcopy(metadata)
    # 구조는 캐시에서 재사용하더라도 metadata를 내보내는 실행기는 현재
    # 버전이므로 과거 캐시의 resolver 버전을 그대로 노출하지 않습니다.
    result["resolver_version"] = RESOLVER_VERSION
    result["cache"] = {
        "enabled": True,
        "hit": True,
        "key": cache_key,
        "format_version": CACHE_FORMAT_VERSION,
        "stored_at_utc": stored_at_utc,
        "reused_at_utc": utc_now(),
    }
    return result
