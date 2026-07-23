#!/usr/bin/env python
"""설치 전에도 실행할 수 있는 호환용 CLI 진입점."""

from pathlib import Path
import sys


SOURCE_ROOT = Path(__file__).resolve().parent / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from protein_structure_resolver.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
