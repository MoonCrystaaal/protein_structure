"""파이프라인 내부 데이터 모델."""

from dataclasses import dataclass


@dataclass(frozen=True)
class ExperimentalCandidate:
    pdb_id: str
    entity_id: str
    asym_id: str
    auth_asym_id: str | None
    modeled_residue_count: int
    sequence_length: int
    coordinate_coverage: float
    experimental_method: str
    resolution_angstrom: float | None
    initial_release_date: str
    r_free: float | None
    clashscore: float | None
    ramachandran_outlier_percent: float | None
    global_exact: bool = True
