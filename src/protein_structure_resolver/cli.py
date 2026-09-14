"""명령행 인터페이스."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence

from .config import METADATA_SCHEMA_VERSION, RESOLVER_VERSION
from .errors import OutputExistsError, ResolverError
from .foldseek_3di import extract_3di
from .resolver import resolve
from .storage import OutputTransaction, utc_now, write_metadata


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "단백질 서열의 exact 실험 구조, 기존 AlphaFold DB 모델, "
            "또는 ESMFold2 예측 구조를 저장합니다."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence", help="단백질 아미노산 서열")
    source.add_argument("--fasta", type=Path, help="단일 레코드 FASTA 파일")
    parser.add_argument("--output", required=True, type=Path, help="결과 디렉터리")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="기존 출력 디렉터리를 완성된 새 결과로 교체합니다.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="캐시와 검색 단계 진행 상황을 표시합니다.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="SHA-256 결과 캐시 경로(기본값: OS 사용자 캐시 디렉터리)",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="캐시 조회와 저장을 모두 사용하지 않습니다.",
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="기존 캐시를 무시하고 성공 결과로 갱신합니다.",
    )
    parser.add_argument(
        "--skip-computational",
        action="store_true",
        help="AlphaFold DB 계산 구조 검색을 생략합니다.",
    )
    parser.add_argument(
        "--skip-prediction",
        action="store_true",
        help="기존 exact 구조가 없을 때 ESMFold2 실행을 생략합니다.",
    )
    parser.add_argument(
        "--with-3di",
        action="store_true",
        help=(
            "구조 해결 성공 후 Foldseek로 아미노산 및 3Di FASTA를 "
            "추가 생성합니다. 유사도 검색은 수행하지 않습니다."
        ),
    )
    parser.add_argument(
        "--foldseek-bin",
        help=(
            "Foldseek 실행 파일 이름 또는 경로. 생략하면 네이티브 환경을 "
            "먼저 찾고 Windows에서는 WSL도 확인합니다."
        ),
    )
    args = parser.parse_args(argv)
    if args.no_cache and args.refresh_cache:
        parser.error("--no-cache와 --refresh-cache는 함께 사용할 수 없습니다.")
    return args


def configure_console_utf8() -> None:
    """Windows에서 출력이 리디렉션되어도 한글 진단을 유지합니다."""
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(encoding="utf-8", errors="replace")


def main(argv: Sequence[str] | None = None) -> int:
    configure_console_utf8()
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )
    try:
        with OutputTransaction(
            args.output, overwrite=args.overwrite
        ) as transaction:
            if transaction.stage is None:
                raise OSError("임시 출력 디렉터리를 만들지 못했습니다.")
            try:
                metadata = resolve(args, output=transaction.stage)
                if (
                    metadata.get("status") == "success"
                    and bool(getattr(args, "with_3di", False))
                ):
                    LOGGER = logging.getLogger(__name__)
                    LOGGER.info("Foldseek 3Di 구조 알파벳 추출")
                    structure_path = transaction.stage / str(
                        metadata["structure_file"]
                    )
                    metadata["structure_alphabet"] = extract_3di(
                        structure_path,
                        transaction.stage,
                        source=metadata.get("source"),
                        foldseek_bin=getattr(args, "foldseek_bin", None),
                    )
            except ResolverError as exc:
                error = {
                    "schema_version": METADATA_SCHEMA_VERSION,
                    "resolver_version": RESOLVER_VERSION,
                    "status": "error",
                    "stage": exc.stage,
                    "error_code": exc.error_code,
                    "message": str(exc),
                    "failed_at": utc_now(),
                }
                write_metadata(transaction.stage, error)
                transaction.commit()
                print(f"오류: {exc}", file=sys.stderr)
                return 1
            write_metadata(transaction.stage, metadata)
            transaction.commit()
    except OutputExistsError as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"파일 오류: {exc}", file=sys.stderr)
        return 1

    if metadata["status"] == "success":
        destination = args.output.resolve() / metadata["structure_file"]
        print(f"완료: {metadata['source']} 구조를 {destination}에 저장했습니다.")
        return 0
    print(metadata["message"])
    return 2
