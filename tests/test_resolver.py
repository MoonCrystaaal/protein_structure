import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from protein_structure_resolver.resolver import resolve


class ResolverTests(unittest.TestCase):
    def test_short_sequence_skips_database_and_reaches_prediction(self) -> None:
        args = SimpleNamespace(
            sequence="ACDEF",
            fasta=None,
            output=Path("unused"),
            skip_computational=False,
            skip_prediction=False,
            no_cache=True,
            cache_dir=None,
            refresh_cache=False,
        )
        prediction = {
            "status": "success",
            "source": "ESMFold2",
            "experimental": False,
            "precomputed": False,
            "model": "test-model",
            "mean_plddt": 80.0,
            "ptm": 0.8,
            "structure_file": "structure.cif",
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with (
                patch(
                    "protein_structure_resolver.resolver."
                    "find_exact_experimental_structures"
                ) as experimental_search,
                patch(
                    "protein_structure_resolver.resolver."
                    "find_exact_alphafold_model"
                ) as computational_search,
                patch(
                    "protein_structure_resolver.resolver.run_esmfold2",
                    return_value=prediction,
                ) as prediction_call,
            ):
                metadata = resolve(args, output=output)

            experimental_search.assert_not_called()
            computational_search.assert_not_called()
            prediction_call.assert_called_once()
            statuses = {
                item["stage"]: item["status"]
                for item in metadata["pipeline"]["search_trace"]
            }
            self.assertEqual(statuses["rcsb_experimental"], "skipped")
            self.assertEqual(statuses["alphafold_db"], "skipped")
            self.assertEqual(statuses["esmfold2"], "success")
            self.assertEqual(metadata["source_class"], "new_prediction")


if __name__ == "__main__":
    unittest.main()
