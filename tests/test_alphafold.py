import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from protein_structure_resolver.alphafold import (
    rank_alphafold_models,
    save_alphafold_result,
)


def model(
    model_id: str,
    *,
    plddt: float,
    version: int,
) -> dict[str, object]:
    return {
        "modelEntityId": model_id,
        "uniprotAccession": "P00001",
        "globalMetricValue": plddt,
        "latestVersion": version,
        "modelCreatedDate": "2026-01-01",
        "cifUrl": "https://example.test/model.cif",
    }


class AlphaFoldRankingTests(unittest.TestCase):
    def test_records_selection_and_rejection_reasons(self) -> None:
        selected, records = rank_alphafold_models(
            [
                model("AF-Z", plddt=90.0, version=4),
                model("AF-B", plddt=92.0, version=3),
                model("AF-A", plddt=92.0, version=4),
            ]
        )

        self.assertEqual(selected["modelEntityId"], "AF-A")
        by_id = {record["model_id"]: record for record in records}
        self.assertEqual(
            by_id["AF-Z"]["decision"]["reason"], "lower_mean_plddt"
        )
        self.assertEqual(
            by_id["AF-B"]["decision"]["reason"], "older_model_version"
        )
        self.assertTrue(by_id["AF-A"]["selected"])

    def test_model_id_uses_ascending_tiebreak(self) -> None:
        selected, _ = rank_alphafold_models(
            [
                model("AF-B", plddt=92.0, version=4),
                model("AF-A", plddt=92.0, version=4),
            ]
        )
        self.assertEqual(selected["modelEntityId"], "AF-A")

    def test_saved_metadata_contains_selection_rule_and_candidates(self) -> None:
        selected, records = rank_alphafold_models(
            [model("AF-A", plddt=92.0, version=4)]
        )
        with tempfile.TemporaryDirectory() as directory:
            with patch(
                "protein_structure_resolver.alphafold.download_text",
                return_value="data_test\n",
            ):
                metadata = save_alphafold_result(
                    Mock(),
                    selected,
                    {"verified_alphafold_model_count": 1},
                    records,
                    "ACDEFGHIKLMNPQRSTVWYACDEF",
                    Path(directory),
                )

        self.assertEqual(
            metadata["selection"]["ranking_rule_version"], "1.0"
        )
        self.assertEqual(
            metadata["selection"]["selection_reason"],
            "selected_by_alphafold_ranking_rule",
        )
        self.assertTrue(metadata["candidates"][0]["selected"])


if __name__ == "__main__":
    unittest.main()
