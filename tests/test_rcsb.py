import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from protein_structure_resolver.errors import DatabaseUnavailableError
from protein_structure_resolver.models import ExperimentalCandidate
from protein_structure_resolver.rcsb import (
    _best_instance,
    fetch_selected_structure_annotations,
    hit_can_be_global_exact,
    rank_experimental_candidates,
    save_experimental_result,
    select_experimental_candidate,
)


def candidate(
    pdb_id: str,
    *,
    modeled: int = 100,
    method: str = "X-RAY DIFFRACTION",
    resolution: float | None = 2.0,
    release: str = "2020-01-01T00:00:00Z",
) -> ExperimentalCandidate:
    return ExperimentalCandidate(
        pdb_id=pdb_id,
        entity_id="1",
        asym_id="A",
        auth_asym_id="A",
        modeled_residue_count=modeled,
        sequence_length=100,
        coordinate_coverage=modeled / 100,
        experimental_method=method,
        resolution_angstrom=resolution,
        initial_release_date=release,
        r_free=None,
        clashscore=None,
        ramachandran_outlier_percent=None,
    )


class SearchTests(unittest.TestCase):
    def test_verbose_match_context_prefilter(self) -> None:
        hit = {
            "services": [
                {
                    "service_type": "sequence",
                    "nodes": [
                        {
                            "match_context": [
                                {
                                    "sequence_identity": 1.0,
                                    "mismatches": 0,
                                    "gaps_opened": 0,
                                    "query_length": 25,
                                    "subject_length": 25,
                                    "alignment_length": 25,
                                    "query_beg": 1,
                                    "query_end": 25,
                                    "subject_beg": 1,
                                    "subject_end": 25,
                                }
                            ]
                        }
                    ],
                }
            ]
        }
        self.assertTrue(hit_can_be_global_exact(hit, 25))
        context = hit["services"][0]["nodes"][0]["match_context"][0]
        context["subject_length"] = 26
        self.assertFalse(hit_can_be_global_exact(hit, 25))

    def test_missing_match_context_is_kept_for_canonical_check(self) -> None:
        self.assertTrue(hit_can_be_global_exact({"identifier": "1ABC_1"}, 25))

    def test_selected_structure_annotations_are_filtered(self) -> None:
        response = Mock()
        response.raise_for_status.return_value = None
        response.json.return_value = {
            "data": {
                "polymer_entity": {
                    "rcsb_polymer_entity_feature": [
                        {
                            "type": "modified_monomer",
                            "feature_id": "MSE",
                            "name": "SELENOMETHIONINE",
                            "description": None,
                            "feature_positions": [
                                {"beg_seq_id": 10, "end_seq_id": None}
                            ],
                        },
                        {
                            "type": "Pfam",
                            "feature_positions": [
                                {"beg_seq_id": 1, "end_seq_id": 50}
                            ],
                        },
                    ],
                    "polymer_entity_instances": [
                        {
                            "rcsb_id": "1AAA.A",
                            "rcsb_polymer_entity_instance_container_identifiers": {
                                "asym_id": "A",
                                "auth_asym_id": "A",
                            },
                            "rcsb_polymer_instance_feature": [
                                {
                                    "type": "UNOBSERVED_RESIDUE_XYZ",
                                    "feature_positions": [
                                        {"beg_seq_id": 20, "end_seq_id": 22}
                                    ],
                                }
                            ],
                        }
                    ],
                }
            }
        }
        with patch(
            "protein_structure_resolver.rcsb.request_with_retry",
            return_value=response,
        ):
            annotations = fetch_selected_structure_annotations(
                Mock(), candidate("1AAA")
            )
        self.assertTrue(annotations["has_modified_monomers"])
        self.assertEqual(
            annotations["modified_monomers"][0]["ranges"],
            [{"start": 10, "end": 10}],
        )
        self.assertEqual(
            annotations["missing_residue_ranges"],
            [{"start": 20, "end": 22}],
        )
        self.assertEqual(annotations["missing_residue_count"], 3)
        self.assertEqual(annotations["retrieval_status"], "success")

    def test_annotation_failure_does_not_discard_structure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with (
                patch(
                    "protein_structure_resolver.rcsb.download_text",
                    return_value="data_test\n",
                ),
                patch(
                    "protein_structure_resolver.rcsb."
                    "fetch_selected_structure_annotations",
                    side_effect=DatabaseUnavailableError("temporary failure"),
                ),
            ):
                metadata = save_experimental_result(
                    Mock(),
                    candidate("1AAA"),
                    [],
                    {"exact_candidate_count": 1},
                    "A" * 100,
                    Path(directory),
                )

            self.assertTrue((Path(directory) / "structure.cif").is_file())
        annotations = metadata["selected_structure"][
            "coordinate_annotations"
        ]
        self.assertEqual(annotations["retrieval_status"], "unavailable")
        self.assertIsNone(annotations["missing_residue_count"])
        self.assertTrue(metadata["validation"]["warnings"])


class RankingTests(unittest.TestCase):
    def test_equal_coverage_instance_uses_ascending_asym_id(self) -> None:
        entity = {
            "polymer_entity_instances": [
                {
                    "rcsb_polymer_entity_instance_container_identifiers": {
                        "asym_id": "B"
                    },
                    "rcsb_polymer_instance_info": {
                        "modeled_residue_count": 100
                    },
                },
                {
                    "rcsb_polymer_entity_instance_container_identifiers": {
                        "asym_id": "A"
                    },
                    "rcsb_polymer_instance_info": {
                        "modeled_residue_count": 100
                    },
                },
            ]
        }
        self.assertEqual(_best_instance(entity), ("A", None, 100))

    def test_coverage_wins_before_resolution(self) -> None:
        selected = select_experimental_candidate(
            [
                candidate("1AAA", modeled=99, resolution=1.0),
                candidate("2BBB", modeled=100, resolution=3.0),
            ]
        )
        self.assertEqual(selected.pdb_id, "2BBB")

    def test_resolution_compared_within_method_then_latest(self) -> None:
        selected = select_experimental_candidate(
            [
                candidate(
                    "1AAA",
                    method="X-RAY DIFFRACTION",
                    resolution=1.5,
                    release="2019-01-01T00:00:00Z",
                ),
                candidate(
                    "2BBB",
                    method="X-RAY DIFFRACTION",
                    resolution=2.0,
                    release="2025-01-01T00:00:00Z",
                ),
                candidate(
                    "3CCC",
                    method="ELECTRON MICROSCOPY",
                    resolution=3.0,
                    release="2024-01-01T00:00:00Z",
                ),
            ]
        )
        self.assertEqual(selected.pdb_id, "3CCC")

    def test_pdb_id_breaks_complete_tie(self) -> None:
        selected = select_experimental_candidate(
            [candidate("2BBB"), candidate("1AAA")]
        )
        self.assertEqual(selected.pdb_id, "1AAA")

    def test_candidate_records_include_rejection_reason(self) -> None:
        selected, records = rank_experimental_candidates(
            [
                candidate("1AAA", modeled=99, resolution=1.0),
                candidate("2BBB", modeled=100, resolution=3.0),
            ]
        )
        self.assertEqual(selected.pdb_id, "2BBB")
        by_id = {record["pdb_id"]: record for record in records}
        self.assertFalse(by_id["1AAA"]["selected"])
        self.assertEqual(
            by_id["1AAA"]["decision"]["reason"],
            "lower_coordinate_coverage",
        )
        self.assertTrue(by_id["2BBB"]["selected"])


if __name__ == "__main__":
    unittest.main()
