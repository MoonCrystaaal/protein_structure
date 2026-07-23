import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from protein_structure_resolver.cli import main


class CliTests(unittest.TestCase):
    def test_success_message_accepts_experimental_result(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            metadata = {
                "status": "success",
                "source": "RCSB_PDB",
                "structure_file": "structure.cif",
            }
            args = SimpleNamespace(
                output=Path(directory) / "result",
                overwrite=False,
                verbose=False,
            )
            with (
                patch(
                    "protein_structure_resolver.cli.parse_args",
                    return_value=args,
                ),
                patch(
                    "protein_structure_resolver.cli.resolve",
                    return_value=metadata,
                ),
                patch("protein_structure_resolver.cli.write_metadata"),
            ):
                self.assertEqual(main([]), 0)


if __name__ == "__main__":
    unittest.main()
