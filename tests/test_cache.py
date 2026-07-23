import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from protein_structure_resolver.cache import (
    lookup_cache,
    materialize_cache_hit,
    store_cache,
)
from protein_structure_resolver.config import (
    CACHE_FORMAT_VERSION,
    METADATA_SCHEMA_VERSION,
)
from protein_structure_resolver.sequence import sequence_sha256


class CacheTests(unittest.TestCase):
    def test_round_trip_validated_success_result(self) -> None:
        sequence = "A" * 25
        args = SimpleNamespace(
            skip_computational=False,
            skip_prediction=False,
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.cif"
            source.write_text("data_test\n#\n", encoding="utf-8")
            metadata = {
                "schema_version": METADATA_SCHEMA_VERSION,
                "status": "success",
                "source": "RCSB_PDB",
                "source_class": "experimental",
                "structure_file": "structure.cif",
                "input": {
                    "sequence_length": len(sequence),
                    "sequence_sha256": sequence_sha256(sequence),
                },
                "cache": {
                    "format_version": CACHE_FORMAT_VERSION,
                    "enabled": True,
                    "hit": False,
                },
            }
            cache_dir = root / "cache"
            store_cache(cache_dir, sequence, metadata, source)
            lookup = lookup_cache(cache_dir, sequence, args)
            self.assertEqual(lookup.reason, "hit")

            output = root / "output"
            cached = materialize_cache_hit(lookup, output, sequence)
            self.assertTrue(cached["cache"]["hit"])
            self.assertTrue((output / "structure.cif").is_file())

    def test_disabled_prediction_source_is_not_reused(self) -> None:
        sequence = "A" * 25
        args = SimpleNamespace(
            skip_computational=False,
            skip_prediction=True,
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.cif"
            source.write_text("data_test\n#\n", encoding="utf-8")
            metadata = {
                "schema_version": METADATA_SCHEMA_VERSION,
                "status": "success",
                "source": "ESMFold2",
                "source_class": "new_prediction",
                "structure_file": "structure.cif",
                "input": {
                    "sequence_length": len(sequence),
                    "sequence_sha256": sequence_sha256(sequence),
                },
                "cache": {
                    "format_version": CACHE_FORMAT_VERSION,
                    "enabled": True,
                    "hit": False,
                },
            }
            cache_dir = root / "cache"
            store_cache(cache_dir, sequence, metadata, source)
            lookup = lookup_cache(cache_dir, sequence, args)
            self.assertEqual(lookup.reason, "source_disabled_by_cli")


if __name__ == "__main__":
    unittest.main()
