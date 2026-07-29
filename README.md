# Protein Structure Resolver

단일 단백질 서열을 다음 우선순위로 구조 파일에 연결합니다.

1. RCSB PDB의 전역(global) 100% exact 실험 구조
2. RCSB에 색인된 전역 exact AlphaFold DB 모델
3. Biohub Platform의 ESMFold2 신규 예측

부분 정렬 구간에서만 identity가 100%인 결과는 exact 구조로 인정하지
않습니다. RCSB 검색 결과를 사전 필터링한 뒤 Data API의 canonical
sequence를 입력 서열과 길이·문자 단위로 다시 비교합니다.

## 프로젝트 구조

```text
protein_structure_resolver/
├── src/protein_structure_resolver/
│   ├── cli.py           # CLI와 종료 코드
│   ├── resolver.py      # 검색 순서와 fallback
│   ├── sequence.py      # 입력 정규화·검증·SHA-256
│   ├── rcsb.py          # 실험 구조 검색·평가·선택
│   ├── alphafold.py     # AlphaFold DB exact 모델 검색
│   ├── esmfold2.py      # ESMFold2 API 예측과 신뢰도 요약
│   ├── http_client.py   # HTTP 세션과 재시도
│   ├── storage.py       # mmCIF 검증과 원자적 저장
│   ├── cache.py         # SHA-256 성공 결과 캐시
│   ├── metadata.py      # 공통 결과 스키마와 검색 이력
│   ├── models.py        # 내부 데이터 모델
│   ├── config.py        # API 주소와 설정
│   └── errors.py        # 단계별 오류
├── tests/               # 모듈별 단위 테스트
├── outputs/             # 실행 결과(버전 관리 제외)
├── sequence_to_structure.py
├── pyproject.toml
└── requirements-esmfold2.txt
```

`sequence_to_structure.py`는 설치 전에도 사용할 수 있는 얇은 호환용
진입점이며 실제 구현은 `src/protein_structure_resolver`에 있습니다.

## 설치

Python 3.11 이상에서:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

ESMFold2까지 사용할 때만 공식 Biohub SDK를 추가로 설치합니다.

```powershell
python -m pip install -r requirements-esmfold2.txt
$env:BIOHUB_API_TOKEN = "발급받은 토큰"
```

API 토큰은 코드나 저장소에 기록하지 않습니다.

실제 RCSB 구조 검색과 ESMFold2 fallback 실행 예시는
[`examples/RUN_EXAMPLE.md`](examples/RUN_EXAMPLE.md)에 입력 FASTA,
명령어, 검색 이력 및 결과 해석과 함께 정리되어 있습니다.

## 실행

서열 직접 입력:

```powershell
python sequence_to_structure.py `
  --sequence "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR" `
  --output outputs\hemoglobin
```

동일한 출력 경로를 다시 사용할 때는 명시적으로 덮어쓰기를 허용해야
합니다.

```powershell
protein-structure-resolver `
  --fasta protein.fasta `
  --output outputs\protein `
  --overwrite
```

설치 후에는 동일한 기능을 명령어로 실행할 수도 있습니다.

```powershell
protein-structure-resolver --fasta protein.fasta --output outputs\protein
```

API 예측 없이 검색 동작만 확인:

```powershell
protein-structure-resolver `
  --fasta protein.fasta `
  --output outputs\search_only `
  --skip-prediction
```

각 실행 폴더에는 다음 파일이 생성됩니다.

```text
outputs/<run-name>/
├── structure.cif
└── metadata.json
```

출력은 같은 상위 디렉터리의 임시 폴더에서 완성된 뒤 최종 경로로
이동합니다. 따라서 실패한 실행의 `metadata.json`과 이전 실행의
`structure.cif`가 섞이지 않습니다.

## 캐시

정규화된 서열의 SHA-256을 캐시 키로 사용합니다. 성공한 구조와
metadata만 OS의 사용자 캐시 디렉터리에 저장하며, 서열 해시,
metadata 스키마, 캐시 형식 및 mmCIF 형식을 검증한 뒤 재사용합니다.

```powershell
# 기존 캐시를 무시하고 다시 검색·예측한 뒤 갱신
protein-structure-resolver `
  --fasta protein.fasta `
  --output outputs\refreshed `
  --refresh-cache

# 캐시를 완전히 사용하지 않음
protein-structure-resolver `
  --fasta protein.fasta `
  --output outputs\uncached `
  --no-cache

# 캐시 위치 직접 지정
protein-structure-resolver `
  --fasta protein.fasta `
  --output outputs\protein `
  --cache-dir cache
```

데이터베이스의 최신 구조를 즉시 다시 평가해야 한다면
`--refresh-cache`를 사용합니다.

## Metadata

`metadata.json`에는 다음 정보가 포함됩니다.

- metadata 및 resolver 버전
- 구조 출처와 `source_class`
- 입력 길이와 SHA-256
- 실험 구조 → AlphaFold DB → ESMFold2 전체 검색 이력
- 실험 후보별 선택 여부와 탈락 단계·이유
- AlphaFold DB 후보별 선택 여부와 탈락 단계·이유
- 선택 구조의 modified monomer와 좌표가 없는 잔기 구간
- 실험 구조 검증 지표 또는 예측 구조 신뢰도
- ESMFold2 반환 서열·pLDDT 개수·pTM 범위·유한 좌표 검증 결과
- 캐시 적중·저장 상태

실험 구조 파일은 annotation보다 먼저 저장합니다. 좌표 결손 또는
modified monomer annotation 조회만 실패하면 구조를 버리지 않고,
metadata의 `retrieval_status`와 `validation.warnings`에 실패를
표시합니다.

## 실험 구조 선택 규칙

1. canonical sequence가 입력과 길이·문자 모두 exact
2. modeled residue count가 가장 큰 후보
3. 동일 실험법 안에서 해상도가 가장 좋은 후보
4. 각 실험법의 우승 후보 중 initial release date가 가장 최신인 후보
5. 완전히 동률이면 PDB ID 오름차순

서로 다른 실험법의 해상도는 직접 비교하지 않습니다.

AlphaFold DB exact 모델이 여러 개면 mean pLDDT 최대 → model version
최신 → model ID 오름차순으로 선택하며, 후보별 선택 근거를 metadata에
기록합니다.

## 테스트

패키지를 editable 모드로 설치했다면:

```powershell
python -m unittest discover -s tests -v
```

설치하지 않은 상태라면:

```powershell
$env:PYTHONPATH = "$PWD\src"
python -m unittest discover -s tests -v
```

## 현재 범위

- 단일 단백질 서열과 20개 표준 아미노산
- 20개 표준 아미노산으로 구성된 비어 있지 않은 서열
- 단일 FASTA 레코드
- 기존 실험 구조는 전체 entry를 mmCIF로 저장

25 residues 미만 서열은 RCSB sequence search와 AlphaFold DB 색인 검색을
건너뛰고 ESMFold2로 전달합니다. `--skip-prediction`도 지정했다면
`not_found` 결과를 반환합니다.

복합체 입력, 유사 서열 구조 대체, CDR/interface 분석, 자동 chain 추출,
biological assembly 선택은 현재 범위에 포함하지 않습니다.
