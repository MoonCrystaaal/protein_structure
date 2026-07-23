"""RCSB PDB exact 검색, 실험 구조 평가 및 선택."""

import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable, Sequence

import requests

from .config import (
    EXPERIMENTAL_ENTITY_QUERY,
    GRAPHQL_BATCH_SIZE,
    PAGE_SIZE,
    RANKING_RULE_VERSION,
    RCSB_FILE_URL,
    RCSB_GRAPHQL_URL,
    RCSB_SEARCH_URL,
    SELECTED_ENTITY_FEATURE_QUERY,
)
from .errors import DatabaseUnavailableError
from .http_client import request_with_retry
from .models import ExperimentalCandidate
from .sequence import is_global_exact, sequence_sha256
from .storage import atomic_write_text, download_text, utc_now


LOGGER = logging.getLogger(__name__)


def build_rcsb_query(
    sequence: str, content_type: str, start: int
) -> dict[str, Any]:
    if content_type not in {"experimental", "computational"}:
        raise ValueError(f"지원하지 않는 RCSB content type: {content_type}")
    return {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "value": sequence,
                "sequence_type": "protein",
                "identity_cutoff": 1.0,
                "evalue_cutoff": 1.0,
            },
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "results_verbosity": "verbose",
            "results_content_type": [content_type],
            "paginate": {"start": start, "rows": PAGE_SIZE},
        },
        "return_type": "polymer_entity",
    }


def _match_contexts(hit: dict[str, Any]) -> Iterable[dict[str, Any]]:
    for service in hit.get("services", []):
        if service.get("service_type") != "sequence":
            continue
        for node in service.get("nodes", []):
            yield from node.get("match_context", [])


def hit_can_be_global_exact(hit: dict[str, Any], query_length: int) -> bool:
    """저비용 사전 필터이며 canonical sequence는 이후 다시 비교합니다."""
    required_fields = {
        "sequence_identity",
        "mismatches",
        "gaps_opened",
        "query_length",
        "subject_length",
        "alignment_length",
        "query_beg",
        "query_end",
        "subject_beg",
        "subject_end",
    }
    contexts = list(_match_contexts(hit))
    if not contexts:
        return True
    uncertain = False
    for context in contexts:
        if not required_fields.issubset(context):
            uncertain = True
            continue
        if (
            context.get("sequence_identity") == 1.0
            and context.get("mismatches") == 0
            and context.get("gaps_opened") == 0
            and context.get("query_length") == query_length
            and context.get("subject_length") == query_length
            and context.get("alignment_length") == query_length
            and context.get("query_beg") == 1
            and context.get("query_end") == query_length
            and context.get("subject_beg") == 1
            and context.get("subject_end") == query_length
        ):
            return True
    return uncertain


def search_rcsb_entities(
    session: requests.Session, sequence: str, content_type: str
) -> tuple[list[str], int]:
    identifiers: list[str] = []
    start = 0
    total_count = 0

    while True:
        response = request_with_retry(
            session,
            "POST",
            RCSB_SEARCH_URL,
            json=build_rcsb_query(sequence, content_type, start),
        )
        if response.status_code == 204:
            return identifiers, total_count
        try:
            response.raise_for_status()
            body = response.json()
        except (requests.RequestException, ValueError) as exc:
            raise DatabaseUnavailableError(
                f"RCSB {content_type} 검색 응답을 해석할 수 없습니다."
            ) from exc

        page = body.get("result_set", [])
        total_count = int(body.get("total_count", len(page)))
        identifiers.extend(
            hit["identifier"]
            for hit in page
            if "identifier" in hit
            and hit_can_be_global_exact(hit, len(sequence))
        )
        start += len(page)
        if not page or start >= total_count:
            return list(dict.fromkeys(identifiers)), total_count


def _chunks(values: Sequence[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def fetch_entities_graphql(
    session: requests.Session, identifiers: Sequence[str], query: str
) -> list[dict[str, Any]]:
    entities: list[dict[str, Any]] = []
    for batch in _chunks(identifiers, GRAPHQL_BATCH_SIZE):
        response = request_with_retry(
            session,
            "POST",
            RCSB_GRAPHQL_URL,
            json={"query": query, "variables": {"ids": batch}},
        )
        try:
            response.raise_for_status()
            body = response.json()
        except (requests.RequestException, ValueError) as exc:
            raise DatabaseUnavailableError(
                "RCSB Data API 응답을 해석할 수 없습니다."
            ) from exc
        if body.get("errors"):
            message = body["errors"][0].get("message", "unknown")
            raise DatabaseUnavailableError(f"RCSB Data API 오류: {message}")
        entities.extend(body.get("data", {}).get("polymer_entities") or [])
    return entities


def _first_number(values: Any) -> float | None:
    if values is None:
        return None
    if not isinstance(values, list):
        values = [values]
    numbers: list[float] = []
    for value in values:
        try:
            numbers.append(float(value))
        except (TypeError, ValueError):
            continue
    return min(numbers) if numbers else None


def _best_instance(
    entity: dict[str, Any],
) -> tuple[str, str | None, int] | None:
    instances = entity.get("polymer_entity_instances") or []
    parsed: list[tuple[str, str | None, int]] = []
    for instance in instances:
        ids = (
            instance.get(
                "rcsb_polymer_entity_instance_container_identifiers"
            )
            or {}
        )
        asym_id = ids.get("asym_id")
        if not asym_id:
            rcsb_id = instance.get("rcsb_id", "")
            asym_id = (
                rcsb_id.rsplit(".", 1)[-1] if "." in rcsb_id else None
            )
        if not asym_id:
            continue
        count = int(
            (instance.get("rcsb_polymer_instance_info") or {}).get(
                "modeled_residue_count", 0
            )
            or 0
        )
        parsed.append((asym_id, ids.get("auth_asym_id"), count))
    if not parsed:
        return None
    max_modeled = max(value[2] for value in parsed)
    return min(
        (value for value in parsed if value[2] == max_modeled),
        key=lambda value: value[0],
    )


def build_experimental_candidate(
    entity: dict[str, Any], sequence: str
) -> ExperimentalCandidate | None:
    canonical = (entity.get("entity_poly") or {}).get(
        "pdbx_seq_one_letter_code_can"
    )
    if not is_global_exact(sequence, canonical):
        return None
    best_instance = _best_instance(entity)
    entry = entity.get("entry") or {}
    if best_instance is None or not entry:
        return None

    asym_id, auth_asym_id, modeled_count = best_instance
    rcsb_entity_id = entity.get("rcsb_id", "")
    if "_" not in rcsb_entity_id:
        return None
    pdb_id, entity_id = rcsb_entity_id.rsplit("_", 1)
    entry_info = entry.get("rcsb_entry_info") or {}
    methods = [
        item.get("method")
        for item in (entry.get("exptl") or [])
        if item.get("method")
    ]
    method = " + ".join(sorted(set(methods)))
    if not method:
        method = entry_info.get("experimental_method") or "UNKNOWN"

    refine = entry.get("refine") or []
    geometry = entry.get("pdbx_vrpt_summary_geometry") or []
    release = (entry.get("rcsb_accession_info") or {}).get(
        "initial_release_date", ""
    )
    return ExperimentalCandidate(
        pdb_id=pdb_id,
        entity_id=entity_id,
        asym_id=asym_id,
        auth_asym_id=auth_asym_id,
        modeled_residue_count=modeled_count,
        sequence_length=len(sequence),
        coordinate_coverage=round(modeled_count / len(sequence), 6),
        experimental_method=method,
        resolution_angstrom=_first_number(
            entry_info.get("resolution_combined")
        ),
        initial_release_date=release,
        r_free=_first_number(
            [item.get("ls_R_factor_R_free") for item in refine]
        ),
        clashscore=_first_number(
            [item.get("clashscore") for item in geometry]
        ),
        ramachandran_outlier_percent=_first_number(
            [
                item.get("percent_ramachandran_outliers")
                for item in geometry
            ]
        ),
    )


def find_exact_experimental_structures(
    session: requests.Session, sequence: str
) -> tuple[list[ExperimentalCandidate], dict[str, int]]:
    identifiers, raw_count = search_rcsb_entities(
        session, sequence, "experimental"
    )
    entities = fetch_entities_graphql(
        session, identifiers, EXPERIMENTAL_ENTITY_QUERY
    )
    candidates = [
        candidate
        for entity in entities
        if (
            candidate := build_experimental_candidate(entity, sequence)
        )
        is not None
    ]
    candidates.sort(
        key=lambda item: (item.pdb_id, item.entity_id, item.asym_id)
    )
    return candidates, {
        "raw_candidate_count": raw_count,
        "prefiltered_candidate_count": len(identifiers),
        "exact_candidate_count": len(candidates),
    }


def select_experimental_candidate(
    candidates: Sequence[ExperimentalCandidate],
) -> ExperimentalCandidate:
    selected, _ = rank_experimental_candidates(candidates)
    return selected


def _candidate_key(
    candidate: ExperimentalCandidate,
) -> tuple[str, str, str]:
    return candidate.pdb_id, candidate.entity_id, candidate.asym_id


def rank_experimental_candidates(
    candidates: Sequence[ExperimentalCandidate],
) -> tuple[ExperimentalCandidate, list[dict[str, Any]]]:
    if not candidates:
        raise ValueError("선택할 실험 구조 후보가 없습니다.")

    decisions: dict[tuple[str, str, str], dict[str, str]] = {}
    max_modeled = max(item.modeled_residue_count for item in candidates)
    coverage_winners = [
        item
        for item in candidates
        if item.modeled_residue_count == max_modeled
    ]
    for item in candidates:
        if item.modeled_residue_count < max_modeled:
            decisions[_candidate_key(item)] = {
                "stage": "coordinate_coverage",
                "reason": "lower_coordinate_coverage",
            }

    by_method: dict[str, list[ExperimentalCandidate]] = {}
    for candidate in coverage_winners:
        by_method.setdefault(candidate.experimental_method, []).append(
            candidate
        )

    method_winners: list[ExperimentalCandidate] = []
    for group in by_method.values():
        with_resolution = [
            item for item in group if item.resolution_angstrom is not None
        ]
        if with_resolution:
            best_resolution = min(
                item.resolution_angstrom for item in with_resolution
            )
            group_winners = [
                item
                for item in with_resolution
                if item.resolution_angstrom == best_resolution
            ]
            method_winners.extend(group_winners)
            for item in group:
                if item not in group_winners:
                    reason = (
                        "missing_resolution_within_method"
                        if item.resolution_angstrom is None
                        else "lower_resolution_within_method"
                    )
                    decisions[_candidate_key(item)] = {
                        "stage": "resolution_within_method",
                        "reason": reason,
                    }
        else:
            method_winners.extend(group)

    newest_release = max(
        item.initial_release_date for item in method_winners
    )
    newest = [
        item
        for item in method_winners
        if item.initial_release_date == newest_release
    ]
    for item in method_winners:
        if item not in newest:
            decisions[_candidate_key(item)] = {
                "stage": "initial_release_date",
                "reason": "older_initial_release_date",
            }

    selected = min(
        newest, key=lambda item: (item.pdb_id, item.entity_id, item.asym_id)
    )
    for item in newest:
        if item != selected:
            decisions[_candidate_key(item)] = {
                "stage": "deterministic_tiebreak",
                "reason": "higher_pdb_or_chain_identifier",
            }
    decisions[_candidate_key(selected)] = {
        "stage": "selected",
        "reason": "selected_by_ranking_rule",
    }

    records: list[dict[str, Any]] = []
    for item in candidates:
        decision = decisions[_candidate_key(item)]
        records.append(
            {
                **asdict(item),
                "selected": item == selected,
                "decision": decision,
            }
        )
    return selected, records


def _feature_ranges(feature: dict[str, Any]) -> list[dict[str, int]]:
    ranges: list[dict[str, int]] = []
    for position in feature.get("feature_positions") or []:
        start = position.get("beg_seq_id")
        if start is None:
            continue
        end = position.get("end_seq_id")
        ranges.append(
            {
                "start": int(start),
                "end": int(end if end is not None else start),
            }
        )
    return ranges


def fetch_selected_structure_annotations(
    session: requests.Session,
    selected: ExperimentalCandidate,
) -> dict[str, Any]:
    response = request_with_retry(
        session,
        "POST",
        RCSB_GRAPHQL_URL,
        json={
            "query": SELECTED_ENTITY_FEATURE_QUERY,
            "variables": {
                "entry_id": selected.pdb_id,
                "entity_id": selected.entity_id,
            },
        },
    )
    try:
        response.raise_for_status()
        body = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise DatabaseUnavailableError(
            "선택 구조의 잔기 annotation을 해석할 수 없습니다."
        ) from exc
    if body.get("errors"):
        message = body["errors"][0].get("message", "unknown")
        raise DatabaseUnavailableError(
            f"선택 구조 annotation API 오류: {message}"
        )

    entity = (body.get("data") or {}).get("polymer_entity") or {}
    modified_monomers = []
    for feature in entity.get("rcsb_polymer_entity_feature") or []:
        if feature.get("type") != "modified_monomer":
            continue
        modified_monomers.append(
            {
                "feature_id": feature.get("feature_id"),
                "name": feature.get("name"),
                "description": feature.get("description"),
                "ranges": _feature_ranges(feature),
            }
        )

    missing_ranges: list[dict[str, int]] = []
    for instance in entity.get("polymer_entity_instances") or []:
        identifiers = (
            instance.get(
                "rcsb_polymer_entity_instance_container_identifiers"
            )
            or {}
        )
        asym_id = identifiers.get("asym_id")
        if not asym_id:
            rcsb_id = instance.get("rcsb_id", "")
            asym_id = (
                rcsb_id.rsplit(".", 1)[-1] if "." in rcsb_id else None
            )
        if asym_id != selected.asym_id:
            continue
        for feature in instance.get("rcsb_polymer_instance_feature") or []:
            if feature.get("type") == "UNOBSERVED_RESIDUE_XYZ":
                missing_ranges.extend(_feature_ranges(feature))

    missing_ranges = [
        {"start": start, "end": end}
        for start, end in sorted(
            {(item["start"], item["end"]) for item in missing_ranges}
        )
    ]
    missing_positions = {
        position
        for item in missing_ranges
        for position in range(item["start"], item["end"] + 1)
    }
    return {
        "retrieval_status": "success",
        "modified_monomers": modified_monomers,
        "has_modified_monomers": bool(modified_monomers),
        "missing_residue_ranges": missing_ranges,
        "missing_residue_count": len(missing_positions),
    }


def quality_warnings(candidate: ExperimentalCandidate) -> list[str]:
    warnings: list[str] = []
    if (
        candidate.ramachandran_outlier_percent is not None
        and candidate.ramachandran_outlier_percent > 1.0
    ):
        warnings.append("Ramachandran outlier 비율이 1%를 초과합니다.")
    if candidate.clashscore is not None and candidate.clashscore > 20:
        warnings.append("clashscore가 20을 초과합니다.")
    if candidate.r_free is not None and candidate.r_free > 0.35:
        warnings.append("R-free가 0.35를 초과합니다.")
    return warnings


def save_experimental_result(
    session: requests.Session,
    selected: ExperimentalCandidate,
    candidate_records: Sequence[dict[str, Any]],
    search_summary: dict[str, int],
    sequence: str,
    output: Path,
) -> dict[str, Any]:
    structure_path = output / "structure.cif"
    text = download_text(
        session, f"{RCSB_FILE_URL}/{selected.pdb_id.upper()}.cif"
    )
    atomic_write_text(structure_path, text)
    annotation_warning: str | None = None
    try:
        annotations = fetch_selected_structure_annotations(session, selected)
    except (DatabaseUnavailableError, TypeError, ValueError) as exc:
        annotation_warning = (
            "선택 구조의 좌표 annotation을 가져오지 못했습니다: "
            f"{exc}"
        )
        LOGGER.warning("%s", annotation_warning)
        annotations = {
            "retrieval_status": "unavailable",
            "retrieval_warning": annotation_warning,
            "modified_monomers": None,
            "has_modified_monomers": None,
            "missing_residue_ranges": None,
            "missing_residue_count": None,
        }
    warnings = quality_warnings(selected)
    if annotation_warning is not None:
        warnings.append(annotation_warning)
    return {
        "status": "success",
        "source": "RCSB_PDB",
        "experimental": True,
        "input": {
            "sequence_length": len(sequence),
            "sequence_sha256": sequence_sha256(sequence),
        },
        "search": {
            **search_summary,
            "global_identity": 1.0,
            "query_coverage": 1.0,
            "target_coverage": 1.0,
        },
        "selected_structure": {
            **asdict(selected),
            "structure_file": structure_path.name,
            "selection_reason": "selected_by_ranking_rule",
            "coordinate_annotations": annotations,
        },
        "structure_file": structure_path.name,
        "selection": {
            "ranking_rule_version": RANKING_RULE_VERSION,
            "rule": (
                "global exact -> maximum coordinate coverage -> best "
                "resolution within each experimental method -> newest "
                "initial release -> ascending PDB ID"
            ),
        },
        "validation": {
            "r_free": selected.r_free,
            "clashscore": selected.clashscore,
            "ramachandran_outlier_percent": (
                selected.ramachandran_outlier_percent
            ),
            "warnings": warnings,
        },
        "candidates": list(candidate_records),
        "retrieved_at": utc_now(),
    }
