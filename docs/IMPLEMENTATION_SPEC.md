# Protein Structure Resolver 최종 구현 명세서

## 1. 문서 정보

- 문서 기준일: 2026-07-24
- 구현 버전: `0.3.0`
- metadata 스키마: `1.0`
- 실험 구조 순위 규칙: `1.0`
- 캐시 형식: `2`
- 구현 경로: `C:\Users\dayou\immune\protein_structure_resolver`
- 상태: 현재 코드 기준 구현 명세 및 검토 결과

이 문서는 현재 구현이 수행하는 동작을 정의한다. 마지막의
`알려진 문제와 수정 권고`에는 의도한 명세와 현재 구현 사이의 차이도
기록한다.

## 2. 목적

단일 단백질 아미노산 서열을 입력받아 다음 순서로 가장 적절한 구조를
반환한다.

1. RCSB PDB의 전역 100% exact 실험 구조
2. RCSB에 색인된 전역 100% exact AlphaFold DB 구조
3. Biohub Platform ESMFold2 신규 예측

부분 정렬만 100%인 결과, 치환·삽입·결실이 있는 결과 및 길이가 다른
결과는 exact로 인정하지 않는다.

## 3. 범위

### 3.1 지원

- 한 번에 단일 단백질 서열
- 직접 문자열 또는 단일 레코드 FASTA 입력
- 20개 표준 아미노산 `ACDEFGHIKLMNPQRSTVWY`
- mmCIF 구조 출력
- 실험 구조, 사전 계산 예측 구조, 신규 예측 구조의 출처 구분
- 성공 결과의 로컬 캐시
- 검색·선택·검증 근거를 담은 JSON metadata

### 3.2 지원하지 않음

- 다중 FASTA 및 복합체 서열 입력
- 비표준·모호 잔기(`X`, `B`, `Z`, `U`, `O` 등)
- 유사 서열 구조를 exact 결과 대신 사용하는 기능
- CDR 또는 항체 전용 분석
- 일반적인 protein interface 분석
- 선택 chain만 자동 추출
- biological assembly 자동 선택
- ESMFold v1 fallback

RCSB 실험 구조는 선택된 chain이 포함된 **전체 PDB entry**를 저장한다.
선택 chain은 metadata의 `selected_structure.asym_id`와
`auth_asym_id`로 식별한다.

## 4. 프로젝트 구조

```text
protein_structure_resolver/
├── docs/
│   └── IMPLEMENTATION_SPEC.md
├── src/protein_structure_resolver/
│   ├── __init__.py
│   ├── cli.py
│   ├── resolver.py
│   ├── sequence.py
│   ├── rcsb.py
│   ├── alphafold.py
│   ├── esmfold2.py
│   ├── http_client.py
│   ├── storage.py
│   ├── cache.py
│   ├── metadata.py
│   ├── models.py
│   ├── config.py
│   └── errors.py
├── tests/
├── outputs/
├── sequence_to_structure.py
├── pyproject.toml
├── requirements.txt
├── requirements-esmfold2.txt
└── README.md
```

`sequence_to_structure.py`는 설치 전 실행을 위한 호환 진입점이고, 실제
구현은 `src/protein_structure_resolver`에 둔다.

## 5. 실행 환경과 의존성

- Python 3.11 이상
- 기본 의존성: `requests>=2.31,<3`
- ESMFold2 의존성:
  `esm @ git+https://github.com/Biohub/esm.git@c94ed8d`
- ESMFold2 인증 환경변수: `BIOHUB_API_TOKEN`

토큰은 CLI 인자, 소스 코드, metadata 또는 캐시 파일에 기록하지 않는다.

권장 설치:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
python -m pip install -r requirements-esmfold2.txt
```

현재 점검한 로컬 기본 Python 3.12 환경에는 ESM SDK가 설치되어 있지만,
프로젝트 자체는 editable 설치되지 않은 상태다. 따라서 현재는
`python sequence_to_structure.py ...` 방식이 확실하며, 콘솔 명령
`protein-structure-resolver`를 사용하려면 `pip install -e .`가
필요하다.

## 6. 입력 처리

### 6.1 CLI

다음 둘 중 하나를 반드시 지정한다.

- `--sequence <SEQUENCE>`
- `--fasta <PATH>`

공통 인자:

- `--output <DIR>`: 필수 결과 디렉터리
- `--overwrite`: 기존 결과 디렉터리 교체
- `--verbose`: 단계별 로그 출력
- `--cache-dir <DIR>`: 캐시 위치 직접 지정
- `--no-cache`: 캐시 조회와 저장 비활성화
- `--refresh-cache`: 기존 캐시를 무시하고 재검색 후 갱신
- `--skip-computational`: AlphaFold DB 검색 생략
- `--skip-prediction`: ESMFold2 신규 예측 생략

`--no-cache`와 `--refresh-cache`는 함께 사용할 수 없다.

### 6.2 정규화

1. 빈 줄과 줄 앞뒤 공백을 제거한다.
2. FASTA 헤더는 첫 줄의 하나만 허용한다.
3. 서열 줄을 연결하고 대문자로 변환한다.
4. 맨 끝의 stop 표기 `*` 하나는 제거한다.
5. 비어 있거나 표준 20개 아미노산 외 문자가 있으면 거부한다.

### 6.3 exact 정의

정규화된 입력 `Q`와 데이터베이스 canonical sequence `T`에 대해 다음을
모두 만족해야 한다.

```text
len(Q) == len(T)
Q == T
```

따라서 global identity, query coverage, target coverage가 모두 100%다.

## 7. 전체 처리 흐름

```text
입력 정규화·검증
  → 성공 캐시 확인
  → RCSB 실험 구조 exact 검색
  → 없으면 AlphaFold DB exact 검색
  → 없으면 ESMFold2 신규 예측
  → structure.cif 및 metadata.json 저장
```

외부 DB 장애는 `not_found`로 간주하지 않는다. 검색 API가 실패하면
파이프라인을 중단하고 오류 metadata를 반환한다. 이는 DB 장애 때문에
실험 구조를 놓친 채 신규 예측을 생성하는 상황을 방지한다.

25 residues 미만은 RCSB sequence search의 입력 제한 때문에 실험 구조와
AlphaFold DB 검색을 건너뛴다. 기본 동작은 ESMFold2로 진행하며,
`--skip-prediction`도 지정하면 `not_found`다.

## 8. 성공 캐시

### 8.1 키

정규화된 서열의 ASCII 바이트에 SHA-256을 적용한 64자리 16진수 값을
캐시 디렉터리 이름으로 사용한다.

```text
cache/<sequence_sha256>/
├── structure.cif
└── metadata.json
```

SHA-256은 서열 자체를 파일명으로 노출하지 않으면서 같은 서열을
안정적으로 식별하기 위한 해시다. 암호화가 아니므로 짧거나 알려진
서열의 완전한 비밀성을 보장하는 수단은 아니다.

### 8.2 조회 검증

- `status == "success"`
- metadata schema version 일치
- cache format version 일치
- 입력 서열 SHA-256 일치
- 현재 CLI에서 해당 source class가 허용됨
- `structure.cif`가 `data_`로 시작하는 최소 mmCIF 형식

성공 결과만 저장한다. 오류와 `not_found`는 캐시하지 않는다.

### 8.3 현재 갱신 정책

기본 캐시는 만료 시간이 없다. RCSB 또는 AlphaFold DB의 최신 상태를
다시 확인하려면 사용자가 `--refresh-cache`를 지정해야 한다.

## 9. RCSB PDB 실험 구조

### 9.1 후보 검색과 exact 재검증

1. RCSB Search API sequence service에 `identity_cutoff=1.0`을 사용한다.
2. `results_content_type=["experimental"]`로 실험 구조만 조회한다.
3. 결과를 100개 단위로 끝까지 pagination한다.
4. verbose alignment가 완전 정렬임이 명확하면 후보로 유지한다.
5. alignment 정보가 불완전한 결과도 오탐 방지를 위해 바로 버리지 않고
   canonical sequence 재검증 단계로 넘긴다.
6. RCSB Data API의 canonical sequence를 입력과 길이·문자 단위로 다시
   비교해 최종 exact 후보를 결정한다.

### 9.2 좌표 coverage

각 exact polymer entity에서 `modeled_residue_count`가 가장 큰 instance를
대표 chain으로 선택한다.

```text
coordinate_coverage = modeled_residue_count / input_sequence_length
```

coverage는 mmCIF 파일 자체를 다시 파싱해서 계산한 값이 아니라 RCSB Data
API가 제공하는 modeled residue count에 기반한다.

### 9.3 최종 선택 규칙

1. global exact 후보만 유지
2. modeled residue count가 최대인 후보만 유지
3. 각 동일 실험법 그룹 안에서 가장 작은 Å 해상도 선택
4. 해상도 정보가 있는 후보가 하나라도 있으면 같은 그룹에서 해상도
   정보가 없는 후보는 탈락
5. 서로 다른 실험법의 우승 후보들 중 initial release date가 가장
   최신인 후보 선택
6. 완전 동률이면 `(PDB ID, entity ID, asym ID)` 오름차순 선택

서로 다른 실험법의 해상도 숫자는 직접 비교하지 않는다. 후보마다
선택 여부, 탈락 단계와 이유를 metadata의 `candidates`에 기록한다.

### 9.4 저장 및 품질 정보

- `https://files.rcsb.org/download/<PDB_ID>.cif`의 전체 entry 저장
- 선택 chain의 좌표가 없는 잔기 구간과 개수
- entity의 modified monomer 정보
- `R-free`, clashscore, Ramachandran outlier percentage
- 단순 경고 기준:
  - Ramachandran outlier > 1%
  - clashscore > 20
  - R-free > 0.35

이 경고는 구조를 탈락시키는 필터가 아니라 사용자가 품질을 해석하기
위한 표시다.

구조 파일 다운로드를 먼저 완료한 뒤 annotation을 조회한다. annotation
API만 실패한 경우 구조는 정상 반환하고
`coordinate_annotations.retrieval_status="unavailable"` 및
`validation.warnings`에 실패 사유를 기록한다.

## 10. AlphaFold DB 사전 계산 구조

RCSB 실험 구조가 없고 `--skip-computational`이 아닐 때 수행한다.

1. RCSB sequence service에서 computational entity 검색
2. canonical sequence를 다시 global exact 비교
3. provenance의 `source_db == "AlphaFoldDB"` 확인
4. entity의 UniProt accession 수집 및 중복 제거
5. AlphaFold DB prediction API 조회
6. API 모델의 sequence도 입력과 global exact인지 재검증
7. 유효한 `cifUrl`이 있는 모델만 유지

현재 여러 모델이 있으면 다음 값의 내림차순 최대값을 선택한다.

1. `globalMetricValue`(mean pLDDT)
2. `latestVersion`
3. `modelEntityId` 또는 `entryId` 문자열 오름차순

출력에는 UniProt accession, model ID/version/생성일, mean pLDDT,
pLDDT 문서 URL, PAE URL과 precomputed prediction임을 기록한다.
`selection`에는 규칙 버전과 선택 근거를, `candidates`에는 후보별 선택
여부 및 탈락 단계·이유를 기록한다.

## 11. ESMFold2 신규 예측

실험 구조와 AlphaFold DB exact 구조가 모두 없고
`--skip-prediction`이 아닐 때 수행한다.

### 11.1 호출 설정

- endpoint: `https://biohub.ai`
- model: `esmfold2-fast-2026-05`
- chain ID: `A`
- `num_loops=3`
- `num_sampling_steps=32`

### 11.2 응답 검증

저장 전에 다음을 모두 확인한다.

- 반환 residue sequence가 입력과 exact
- residue별 pLDDT 개수가 입력 길이와 같음
- 모든 pLDDT가 유한한 0~100
- pTM이 유한한 0~1
- 원자 좌표가 비어 있지 않음
- 모든 좌표가 유한함
- 직렬화된 mmCIF가 `data_`로 시작함

검증 실패는 성공 구조로 저장하지 않고 `PREDICTION_UNAVAILABLE` 오류로
처리한다.

### 11.3 신뢰도 metadata

- mean/min pLDDT
- pLDDT ≥ 90 비율
- pLDDT ≥ 70 비율
- pLDDT < 50 비율
- pLDDT < 70인 연속 구간
- pTM
- 서열·pLDDT 개수·좌표 유한성 검증 결과

이 값은 모든 단백질에 적용되는 일반 예측 신뢰도다. CDR 신뢰도는
계산하지 않는다.

## 12. 출력과 원자성

성공 시:

```text
<output>/
├── structure.cif
└── metadata.json
```

`not_found` 또는 처리된 오류에는 `metadata.json`만 있을 수 있다.

출력은 목표 경로 옆 임시 디렉터리에서 완성한 후 최종 경로로
이동한다. 기존 출력은 기본적으로 거부하며 `--overwrite`일 때만
완성된 새 결과로 교체한다. 이 방식은 이전 구조와 새 metadata가
섞이는 것을 방지한다.

## 13. 공통 metadata 계약

모든 정상 처리 결과가 공통으로 포함할 필드:

- `schema_version`
- `resolver_version`
- `status`: `success`, `not_found`, `error`
- `source`: `RCSB_PDB`, `AlphaFold_DB`, `ESMFold2`, 또는 `null`
- `source_class`: `experimental`, `precomputed_prediction`,
  `new_prediction`, 또는 `null`
- `input.sequence_length`
- `input.sequence_sha256`
- `pipeline.resolution_order`
- `pipeline.search_trace`
- `cache`
- UTC 시각

성공 결과는 `structure_file`을 포함한다. source별 품질·선택 필드는
각각 9~11절을 따른다.

오류 metadata:

- `status="error"`
- `stage`
- `error_code`
- 사용자용 `message`
- `failed_at`

## 14. 오류 및 종료 코드

| 종료 코드 | 의미 |
|---:|---|
| 0 | 구조 해결 성공 |
| 1 | 입력, DB, 다운로드, ESMFold2, 파일 또는 출력 오류 |
| 2 | exact 구조가 없고 신규 예측을 생략함 |

주요 오류 코드:

- `INVALID_SEQUENCE`
- `DATABASE_UNAVAILABLE`
- `STRUCTURE_DOWNLOAD_FAILED`
- `PREDICTION_UNAVAILABLE`
- `OUTPUT_EXISTS`

RCSB/AlphaFold HTTP 요청은 기본 30초, 구조 다운로드는 120초 timeout을
사용한다. 408, 429, 500, 502, 503, 504 및 네트워크 예외는 최대 4회
시도하며 지수 backoff를 사용한다.

## 15. 보안과 개인정보

- 입력 서열은 RCSB, AlphaFold DB 또는 Biohub API로 전송될 수 있다.
- 로컬 캐시와 metadata에는 원문 서열 대신 길이와 SHA-256을 저장한다.
- API 토큰은 환경변수에서만 읽고 출력·metadata에 쓰지 않는다.
- 민감하거나 미공개인 서열은 각 외부 서비스의 데이터 처리 정책을
  확인한 뒤 전송해야 한다.

## 16. 검증 결과

2026-07-24 로컬 점검:

- Python source `compileall`: 통과
- 단위 테스트: 28개 모두 통과
- `pip check`: 손상된 의존성 없음
- CLI help 및 root launcher: 정상
- 기존 출력 metadata 4개 JSON 파싱: 정상
- 성공 출력 3개(RCSB, AlphaFold DB, ESMFold2)의 최소 mmCIF 검증: 정상
- 실제 ESMFold2 API end-to-end 호출: HTTP 200 및 파일 저장 성공
- 실제 호출에서 RCSB not found → AlphaFold DB not found → ESMFold2
  success 검색 이력 확인
- 실제 ESMFold2 결과의 서열 exact, pLDDT 수, pTM 범위, 유한 좌표 확인

실제 ESMFold2 smoke test 결과는
`outputs/esmfold2_api_test/metadata.json`에 보존되어 있다. 이 결과의
낮은 pLDDT/pTM은 테스트에 사용한 짧은 합성 서열의 예측 신뢰도가
낮다는 뜻이며 API 실패를 의미하지 않는다.

## 17. 알려진 문제와 수정 권고

### 17.1 후속 수정

1. **캐시에 자동 만료가 없음**

   캐시는 DB 검색보다 먼저 반환된다. 그 사이 새 실험 구조가 공개되어도
   `--refresh-cache` 없이는 과거 AlphaFold DB 또는 ESMFold2 결과를 계속
   반환할 수 있다. “현재 가장 좋은 구조”가 요구사항이면 experimental과
   precomputed 결과에 configurable TTL을 적용하고, TTL 경과 시 RCSB부터
   재검색하는 정책이 필요하다.

### 17.2 일관성 개선

2. **캐시 정책 fingerprint 부족**

   현재 캐시 디렉터리 키는 서열 SHA-256 하나다. resolver/model/순위
   정책 변경은 수동 cache format version 변경에 의존한다. 키 또는
   검증 metadata에 resolver version, source model과 중요 설정을 포함한
   policy fingerprint를 추가하는 것이 안전하다.

3. **캐시 entry 교체 중 짧은 공백이 생길 수 있음**

   기존 entry를 삭제한 다음 staged entry를 이동한다. 동시 실행 또는
   그 사이의 프로세스 종료에서는 캐시 miss나 캐시 손실이 가능하다.
   캐시는 재생성 가능하므로 결과 손상 위험은 낮지만, lock과
   backup-rename 교체를 쓰면 더 견고하다.

### 17.3 테스트 및 운영 개선

4. API pagination/retry, 손상 캐시와 동시 캐시 쓰기에 대한 자동
   테스트를 추가한다.
5. 실제 ESMFold2 smoke test는 과금·토큰 사용 때문에 기본 단위 테스트와
   분리하고 명시적 opt-in integration test로 유지한다.
6. 배포 시점의 dependency snapshot 또는 lock 파일을 관리하는 것이
   좋다.

## 18. 완료 판정

다음 조건을 만족하면 한 실행을 성공으로 판정한다.

1. 입력이 정규화·검증됨
2. 선택된 DB 구조는 canonical sequence가 global exact임
3. 신규 예측은 응답 서열·신뢰도·좌표 검증을 모두 통과함
4. `structure.cif`가 존재하고 최소 mmCIF 검증을 통과함
5. `metadata.json`이 구조 출처, 검색 이력, 캐시 상태와 검증 정보를
   포함함
6. 출력 transaction이 최종 경로에 commit됨
7. CLI 종료 코드가 결과 상태와 일치함
