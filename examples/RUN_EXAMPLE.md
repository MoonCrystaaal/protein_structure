# Protein Structure Resolver 실행 예시

## 목적

실제 FASTA 입력으로 다음 두 경로를 확인한다.

1. RCSB PDB에서 global exact 실험 구조를 찾는 경로
2. RCSB와 AlphaFold DB에 exact 구조가 없을 때 ESMFold2로 예측하는 경로

실행 결과는 저장소 루트의 `outputs/`에 생성된다. `outputs/`는 Git에서
제외되며 예제 입력과 이 문서만 버전 관리한다.

## 공통 준비

Windows PowerShell:

```powershell
cd C:\Users\dayou\immune\protein_structure_resolver

$env:BIOHUB_API_TOKEN = [Environment]::GetEnvironmentVariable(
    "BIOHUB_API_TOKEN",
    "User"
)
```

## 예제 1: RCSB global exact 실험 구조

입력:

```text
examples/inputs/example_rcsb.fasta
```

실행:

```powershell
py -3.12 sequence_to_structure.py `
    --fasta .\examples\inputs\example_rcsb.fasta `
    --output .\outputs\example_rcsb `
    --refresh-cache `
    --verbose `
    2>&1 | Tee-Object -FilePath .\outputs\example_rcsb.log
```

예상 검색 경로:

```text
RCSB 실험 구조 검색 → global exact 후보 선택 → 전체 PDB entry 저장
```

## 예제 2: ESMFold2 fallback

입력:

```text
examples/inputs/example_esmfold2.fasta
```

실행:

```powershell
py -3.12 sequence_to_structure.py `
    --fasta .\examples\inputs\example_esmfold2.fasta `
    --output .\outputs\example_esmfold2 `
    --verbose `
    2>&1 | Tee-Object -FilePath .\outputs\example_esmfold2.log
```

예상 검색 경로:

```text
RCSB not found → AlphaFold DB not found → ESMFold2 신규 예측
```

ESMFold2는 Biohub API 크레딧을 사용할 수 있다. 최초 성공 결과는
서열 SHA-256 캐시에 저장되므로 특별한 이유가 없으면
`--refresh-cache`를 사용하지 않는다.

## 결과 파일

각 실행 폴더:

```text
outputs/<example-name>/
├── structure.cif
└── metadata.json
```

Metadata 핵심 값 확인:

```powershell
$result = Get-Content `
    .\outputs\example_rcsb\metadata.json `
    -Encoding utf8 |
    ConvertFrom-Json

$result.status
$result.source
$result.source_class
$result.pipeline.search_trace
$result.selection
$result.validation
```

## 실제 실행 결과

실행일: 2026-07-29

### RCSB 예제 결과

입력:

| 항목 | 값 |
|---|---|
| FASTA | `examples/inputs/example_rcsb.fasta` |
| 설명 | 사람 hemoglobin alpha 서열 |
| 서열 길이 | 141 residues |
| SHA-256 | `625c708d1e51d4df14c88a2fc08ec59355f96b51c0ea7a0a5eb5a910b92a0675` |

검색 및 선택 결과:

| 항목 | 값 |
|---|---|
| status | `success` |
| source | `RCSB_PDB` |
| source class | `experimental` |
| exact 후보 수 | 248 |
| 선택 PDB | `7PCH` |
| entity / asym ID | `1` / `A` |
| modeled residue count | 141 |
| coordinate coverage | 1.0 |
| 실험법 | `ELECTRON MICROSCOPY` |
| 해상도 | 2.89 Å |
| initial release date | `2022-04-13T00:00:00Z` |
| annotation | `success` |
| 좌표 결손 잔기 | 0 |
| modified monomer | 없음 |
| 캐시 저장 | `stored` |

검색 이력:

```text
rcsb_experimental: found
alphafold_db: not_needed
esmfold2: not_needed
```

생성 구조:

```text
outputs/example_rcsb/structure.cif
크기: 1,251,058 bytes
첫 줄: data_7PCH
```

### ESMFold2 fallback 예제 결과

입력:

| 항목 | 값 |
|---|---|
| FASTA | `examples/inputs/example_esmfold2.fasta` |
| 설명 | fallback 확인용 40 aa 합성 서열 |
| 서열 길이 | 40 residues |
| SHA-256 | `3f58c67c36a585628d33015ae2cd0a227b34f8cbb0e7fcdfaa4c4b3de45c415d` |

검색 및 예측 결과:

| 항목 | 값 |
|---|---|
| status | `success` |
| source | `ESMFold2` |
| source class | `new_prediction` |
| RCSB 검색 | `not_found` |
| AlphaFold DB 검색 | `not_found` |
| model | `esmfold2-fast-2026-05` |
| mean pLDDT | 46.477 |
| min pLDDT | 40.674 |
| pTM | 0.149387 |
| pLDDT < 50 비율 | 0.95 |
| pLDDT ≥ 70 비율 | 0.0 |
| 캐시 저장 | `stored` |

예측 응답 검증:

| 검증 | 결과 |
|---|---|
| 반환 서열 exact | 통과 |
| residue count | 40 |
| pLDDT count | 40 |
| coordinate scalar count | 1,002 |
| 모든 좌표 유한 | 통과 |
| mmCIF | 통과 |

검색 이력:

```text
rcsb_experimental: not_found
alphafold_db: not_found
esmfold2: success
```

생성 구조:

```text
outputs/example_esmfold2/structure.cif
크기: 32,621 bytes
첫 줄: data_complex
```

이 합성 서열은 mean pLDDT와 pTM이 낮으므로 생물학적으로 신뢰할 만한
구조 예제가 아니라 **fallback 및 응답 검증 경로를 확인하기 위한 실행
예제**다. 실제 연구에 사용할 때는 residue별 pLDDT, pTM 및 낮은 신뢰도
구간을 함께 해석해야 한다.
