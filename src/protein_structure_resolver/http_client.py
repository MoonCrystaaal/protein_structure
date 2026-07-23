"""외부 HTTP 요청의 세션, 재시도 및 오류 변환."""

import time
from typing import Any

import requests

from .config import (
    MAX_ATTEMPTS,
    REQUEST_TIMEOUT,
    RETRYABLE_STATUS,
    USER_AGENT,
)
from .errors import DatabaseUnavailableError


def make_session() -> requests.Session:
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    return session


def _retry_delay(response: requests.Response | None, attempt: int) -> float:
    if response is not None:
        retry_after = response.headers.get("Retry-After")
        if retry_after:
            try:
                return max(0.0, float(retry_after))
            except ValueError:
                pass
    return float(2 ** (attempt + 1))


def request_with_retry(
    session: requests.Session,
    method: str,
    url: str,
    *,
    timeout: int = REQUEST_TIMEOUT,
    **kwargs: Any,
) -> requests.Response:
    last_error: Exception | None = None
    last_response: requests.Response | None = None

    for attempt in range(MAX_ATTEMPTS):
        try:
            response = session.request(method, url, timeout=timeout, **kwargs)
            last_response = response
            if response.status_code not in RETRYABLE_STATUS:
                return response
        except requests.RequestException as exc:
            last_error = exc

        if attempt < MAX_ATTEMPTS - 1:
            time.sleep(_retry_delay(last_response, attempt))

    status = (
        f" (마지막 HTTP 상태: {last_response.status_code})"
        if last_response is not None
        else ""
    )
    raise DatabaseUnavailableError(
        f"API 요청에 실패했습니다: {url}{status}"
    ) from last_error
