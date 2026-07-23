import unittest

from protein_structure_resolver.errors import SequenceInputError
from protein_structure_resolver.sequence import (
    is_global_exact,
    normalize_sequence,
    validate_sequence,
)


class SequenceTests(unittest.TestCase):
    def test_normalizes_raw_and_fasta(self) -> None:
        self.assertEqual(normalize_sequence(" acd\nef* "), "ACDEF")
        self.assertEqual(normalize_sequence(">sample\nACD\nEF\n"), "ACDEF")

    def test_rejects_multiple_fasta_records(self) -> None:
        with self.assertRaises(SequenceInputError):
            normalize_sequence(">a\nACD\n>b\nEFG\n")

    def test_rejects_nonstandard_residue(self) -> None:
        with self.assertRaisesRegex(SequenceInputError, "'X'.*3번"):
            validate_sequence("ACX" + "A" * 22)

    def test_exact_is_length_sensitive(self) -> None:
        self.assertTrue(is_global_exact("ACDE", "A C D E\n"))
        self.assertFalse(is_global_exact("ACDE", "ACDEF"))

    def test_short_standard_sequence_is_valid_input(self) -> None:
        validate_sequence("ACDEF")


if __name__ == "__main__":
    unittest.main()
