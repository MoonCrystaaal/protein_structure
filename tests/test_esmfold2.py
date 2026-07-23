import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from protein_structure_resolver.errors import PredictionUnavailableError
from protein_structure_resolver.esmfold2 import (
    find_low_confidence_regions,
    normalize_plddt,
    run_esmfold2,
    summarize_plddt,
    validate_prediction_result,
)


class ConfidenceTests(unittest.TestCase):
    def test_normalizes_zero_to_one_plddt(self) -> None:
        self.assertEqual(normalize_plddt([0.5, 0.9]), [50.0, 90.0])

    def test_summary_and_regions(self) -> None:
        values = [80.0, 60.0, 50.0, 90.0]
        summary = summarize_plddt(values)
        self.assertEqual(summary["mean_plddt"], 70.0)
        self.assertEqual(
            find_low_confidence_regions(values),
            [
                {
                    "start": 2,
                    "end": 3,
                    "length": 2,
                    "mean_plddt": 55.0,
                }
            ],
        )

    def test_validates_sequence_confidence_and_coordinates(self) -> None:
        result = SimpleNamespace(
            complex=SimpleNamespace(
                sequence=["ALA", "CYS", "ASP"],
                atom_positions=[
                    [[0.0, 1.0, 2.0]],
                    [[1.0, 2.0, 3.0]],
                    [[2.0, 3.0, 4.0]],
                ],
            ),
            plddt=[0.8, 0.9, 0.7],
            ptm=0.6,
        )
        plddt, ptm, validation = validate_prediction_result(result, "ACD")
        self.assertEqual(plddt, [80.0, 90.0, 70.0])
        self.assertEqual(ptm, [0.6])
        self.assertTrue(validation["sequence_exact"])
        self.assertTrue(validation["coordinates_finite"])

    def test_rejects_mismatched_prediction_sequence(self) -> None:
        result = SimpleNamespace(
            complex=SimpleNamespace(
                sequence=["ALA", "CYS"],
                atom_positions=[[[0.0, 1.0, 2.0]]],
            ),
            plddt=[0.8, 0.9],
            ptm=0.6,
        )
        with self.assertRaisesRegex(ValueError, "일치하지 않습니다"):
            validate_prediction_result(result, "ACD")


class PredictionTests(unittest.TestCase):
    def test_requires_biohub_token_before_importing_sdk(self) -> None:
        with (
            tempfile.TemporaryDirectory() as directory,
            patch.dict("os.environ", {}, clear=True),
        ):
            with self.assertRaisesRegex(
                PredictionUnavailableError, "BIOHUB_API_TOKEN"
            ):
                run_esmfold2("A" * 25, Path(directory))


if __name__ == "__main__":
    unittest.main()
