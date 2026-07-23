"""구조 및 metadata 파일 저장."""

import json
import shutil
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from .config import DOWNLOAD_TIMEOUT
from .errors import StructureDownloadError
from .errors import OutputExistsError
from .http_client import request_with_retry


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.part")
    try:
        temporary.write_text(text, encoding="utf-8", newline="\n")
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def validate_mmcif(text: str) -> None:
    if not text.lstrip().startswith("data_"):
        raise StructureDownloadError(
            "다운로드 결과가 유효한 mmCIF로 보이지 않습니다."
        )


def download_text(
    session: requests.Session, url: str, *, expected_mmcif: bool = True
) -> str:
    response = request_with_retry(
        session, "GET", url, timeout=DOWNLOAD_TIMEOUT, allow_redirects=True
    )
    try:
        response.raise_for_status()
    except requests.RequestException as exc:
        raise StructureDownloadError(
            f"구조 파일 다운로드에 실패했습니다: {url}"
        ) from exc
    text = response.text
    if expected_mmcif:
        validate_mmcif(text)
    return text


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_metadata(output: Path, metadata: dict[str, Any]) -> None:
    atomic_write_text(
        output / "metadata.json",
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
    )


class OutputTransaction:
    """완성된 결과 디렉터리만 최종 출력 경로에 노출합니다."""

    def __init__(self, target: Path, *, overwrite: bool) -> None:
        self.target = target.expanduser().resolve()
        self.overwrite = overwrite
        self.stage: Path | None = None
        self._committed = False

    def __enter__(self) -> "OutputTransaction":
        if self.target.exists() and not self.overwrite:
            raise OutputExistsError(
                f"출력 경로가 이미 존재합니다: {self.target}. "
                "--overwrite를 지정해야 교체할 수 있습니다."
            )
        self.target.parent.mkdir(parents=True, exist_ok=True)
        self.stage = Path(
            tempfile.mkdtemp(
                prefix=f".{self.target.name}.tmp-",
                dir=self.target.parent,
            )
        )
        return self

    def commit(self) -> None:
        if self.stage is None or not self.stage.exists():
            raise OSError("커밋할 임시 출력 디렉터리가 없습니다.")
        backup: Path | None = None
        try:
            if self.target.exists():
                backup = self.target.with_name(
                    f".{self.target.name}.backup-{uuid.uuid4().hex}"
                )
                self.target.replace(backup)
            self.stage.replace(self.target)
            self._committed = True
        except Exception:
            if backup is not None and backup.exists() and not self.target.exists():
                backup.replace(self.target)
            raise
        finally:
            if backup is not None and backup.exists() and self._committed:
                if backup.is_dir():
                    shutil.rmtree(backup)
                else:
                    backup.unlink()

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        if self.stage is not None and self.stage.exists():
            shutil.rmtree(self.stage)
