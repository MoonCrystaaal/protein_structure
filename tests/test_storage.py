import tempfile
import unittest
from pathlib import Path

from protein_structure_resolver.errors import OutputExistsError
from protein_structure_resolver.storage import OutputTransaction


class OutputTransactionTests(unittest.TestCase):
    def test_refuses_existing_output_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "result"
            target.mkdir()
            (target / "old.txt").write_text("old", encoding="utf-8")
            with self.assertRaises(OutputExistsError):
                with OutputTransaction(target, overwrite=False):
                    pass
            self.assertEqual(
                (target / "old.txt").read_text(encoding="utf-8"), "old"
            )

    def test_overwrite_replaces_only_after_commit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "result"
            target.mkdir()
            (target / "old.txt").write_text("old", encoding="utf-8")
            with OutputTransaction(target, overwrite=True) as transaction:
                assert transaction.stage is not None
                (transaction.stage / "new.txt").write_text(
                    "new", encoding="utf-8"
                )
                self.assertTrue((target / "old.txt").exists())
                transaction.commit()
            self.assertFalse((target / "old.txt").exists())
            self.assertEqual(
                (target / "new.txt").read_text(encoding="utf-8"), "new"
            )


if __name__ == "__main__":
    unittest.main()
