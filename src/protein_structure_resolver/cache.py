"""서열 SHA-256을 키로 사용하는 검증형 로컬 성공 결과 캐시."""

import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import CACHE_FORMAT_VERSION, METADATA_SCHEMA_VERSION
from .metadata import mark_cache_hit
from .sequence import sequence_sha256
from .storage import atomic_write_text, utc_now, validate_mmcif


@dataclass(frozen=True)
class CacheLookup:
    metadata: dict[str, Any] | None
    structure_path: Path | None
    reason: str


def default_cache_dir() -> Path:
    if os.name == "nt" and os.environ.get("LOCALAPPDATA"):
        root = Path(os.environ["LOCALAPPDATA"])
    elif os.environ.get("XDG_CACHE_HOME"):
        root = Path(os.environ["XDG_CACHE_HOME"])
    else:
        root = Path.home() / ".cache"
    return root / "protein_structure_resolver"


def cache_entry_dir(cache_dir: Path, sequence: str) -> Path:
    return cache_dir.expanduser().resolve() / sequence_sha256(sequence)


def _cache_source_allowed(metadata: dict[str, Any], args: Any) -> bool:
    source_class = metadata.get("source_class")
    if source_class == "precomputed_prediction":
        return not args.skip_computational
    if source_class == "new_prediction":
        return not args.skip_prediction
    return source_class == "experimental"


def lookup_cache(
    cache_dir: Path, sequence: str, args: Any
) -> CacheLookup:
    entry = cache_entry_dir(cache_dir, sequence)
    metadata_path = entry / "metadata.json"
    structure_path = entry / "structure.cif"
    if not metadata_path.is_file() or not structure_path.is_file():
        return CacheLookup(None, None, "not_found")
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not isinstance(metadata, dict):
            return CacheLookup(None, None, "invalid_metadata_type")
        if metadata.get("status") != "success":
            return CacheLookup(None, None, "cached_status_not_success")
        if metadata.get("schema_version") != METADATA_SCHEMA_VERSION:
            return CacheLookup(None, None, "schema_version_mismatch")
        cache = metadata.get("cache") or {}
        if cache.get("format_version") != CACHE_FORMAT_VERSION:
            return CacheLookup(None, None, "cache_format_version_mismatch")
        input_data = metadata.get("input") or {}
        if input_data.get("sequence_sha256") != sequence_sha256(sequence):
            return CacheLookup(None, None, "sequence_hash_mismatch")
        if not _cache_source_allowed(metadata, args):
            return CacheLookup(None, None, "source_disabled_by_cli")
        validate_mmcif(structure_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError, UnicodeError):
        return CacheLookup(None, None, "invalid_cache_entry")
    return CacheLookup(metadata, structure_path, "hit")


def materialize_cache_hit(
    lookup: CacheLookup,
    output: Path,
    sequence: str,
) -> dict[str, Any]:
    if lookup.metadata is None or lookup.structure_path is None:
        raise ValueError("캐시 적중 결과가 없습니다.")
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(lookup.structure_path, output / "structure.cif")
    stored_at = (lookup.metadata.get("cache") or {}).get("stored_at_utc")
    return mark_cache_hit(
        lookup.metadata,
        cache_key=sequence_sha256(sequence),
        stored_at_utc=stored_at,
    )


def store_cache(
    cache_dir: Path,
    sequence: str,
    metadata: dict[str, Any],
    structure_path: Path,
) -> None:
    if metadata.get("status") != "success" or not structure_path.is_file():
        return
    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    entry = cache_entry_dir(cache_dir, sequence)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{entry.name}.tmp-", dir=cache_dir)
    )
    try:
        shutil.copy2(structure_path, stage / "structure.cif")
        cached_metadata = deepcopy_metadata(metadata)
        cached_metadata["cache"] = {
            **(cached_metadata.get("cache") or {}),
            "enabled": True,
            "hit": False,
            "key": sequence_sha256(sequence),
            "format_version": CACHE_FORMAT_VERSION,
            "stored_at_utc": utc_now(),
        }
        atomic_write_text(
            stage / "metadata.json",
            json.dumps(cached_metadata, ensure_ascii=False, indent=2) + "\n",
        )
        if entry.is_dir():
            shutil.rmtree(entry)
        elif entry.exists():
            entry.unlink()
        stage.replace(entry)
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def deepcopy_metadata(metadata: dict[str, Any]) -> dict[str, Any]:
    return json.loads(json.dumps(metadata))
