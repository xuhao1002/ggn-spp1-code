# Analysis plans — index and checksum verification

Each row below is one plan listed in Supplementary Table S1 of the paper. The checksum shown is the first 16 characters of the SHA-256 of the file as released here; every one of them matches the value recorded in Table S1, which was recorded before the corresponding analysis was run.

| Table S1 row | File | SHA-256 (first 16) | Verified |
|---|---|---|---|
| Part 2: precursor-stage target screen | `第二阶段预设方案.md` | 1adf491451868896 | ✓ |
| Part 3: spatial analysis | `第三部分_空间分析预设方案.md` | 5d1a3cc50fd8eb93 | ✓ |
| Part 4: normal→precursor drivers | `第四部分_正常到前体驱动因素预设方案.md` | 63f97c694caa9e90 | ✓ |
| KAC fraction, cross-cohort pooling | `预设方案_KAC比例跨队列合并_2026-10-11.md` | cc63b2d74251927e | ✓ |
| Part 5: human genetics | `第五部分_遗传学证据预设方案.md` | 31da180ab8996952 | ✓ |
| Spatial gene-level exploratory | `预设方案_空间基因层面探索性分析_2026-10-11.md` | 084a60d49f9afd93 | ✓ |
| SPP1 source diagnostic | `预设方案_SPP1表达来源诊断_2026-10-11.md` | af880950dec87d9c | ✓ |
| Sensitivity and protein validation | `预设方案_敏感性与蛋白验证_2026-10-12.md` | b15055f872e0001a | ✓ |
| Xenium independent validation | `预设方案_Xenium独立验证_2026-10-12.md` | 818813d56ab28c6a | ✓ |
| Prognosis robustness and Xenium spillover | `预设方案_预后稳健性与Xenium串扰校准_2026-10-12.md` | ec46b0c72c349cc4 | ✓ |
| Supplementary checks T | `预设方案_补充检查T_2026-10-12.md` | b29834f7de82b644 | ✓ |
| Hypoxia H1/H2 | `预设方案_缺氧H1H2_2026-10-12.md` | d1839e14fc692435 | ✓ |
| Hypoxia H1c | `预设方案_缺氧H1c_2026-10-12.md` | 69038b09d2f7bccf | ✓ |
| TF activity, ligand–receptor, RCTD | `预设方案_调控互作RCTD_2026-10-13.md` | 105809cc8cce637a | ✓ |
| Response to external review (A1–A6), version 1.3 | `预设方案_外审回应A1-A6_v1.3_2026-10-19.md` | 1cde317f3bccac30 | ✓ |

The plans are written in Chinese, the working language of the study.

## Superseded versions of the external-review plan (`superseded/`)

Versions 1.0–1.2 of the external-review plan were replaced before any analysis was run; their checksums are listed in the Table S1 row of version 1.3.

| File | SHA-256 (first 16) |
|---|---|
| `superseded/预设方案_外审回应A1-A6_2026-10-19.md` (v1.0) | 006f67678062defc |
| `superseded/预设方案_外审回应A1-A6_v1.1_2026-10-19.md` | d30f66d2a536471d |
| `superseded/预设方案_外审回应A1-A6_v1.2_2026-10-19.md` | 60ffdf205ce1d3ea |

## Execution notes (`execution_notes/`)

Implementation details not fixed by plan v1.3, written before the corresponding analysis was run; the scripts verify their checksums at start-up.

| File | Analysis | SHA-256 (first 16) |
|---|---|---|
| `execution_notes/执行补充_A1c_A1d_2026-10-19.md` | A1c, A1d | e108a3449e1173b6 |
| `execution_notes/执行补充_A5_2026-10-21.md` | A5 | 74ef4e5cc0e4f9bc |

## Deviation record

`偏离与修订记录_外审回应A1-A6.md` is the append-only record of deviations (D1–D7), process records (P1) and the A3 data-access note for plan v1.3, as it stood when this release was frozen. Entries for post hoc analyses (D2, D5, D6, D7) were written before the analyses were run; the scripts for D5, D6 and D7 refuse to run unless their entry is present and log the checksum of the record at run time. Timestamp corrections are appended as separate lines rather than edited in place. `../ANALYSES.md` maps every analysis ID to its plan or record, script and result directory.

## One checksum in Table S1 has no file here

Table S1 lists two checksums for the TF-activity plan: `87fcde8977805e49` for the version covering the TF-activity analysis (R), and `105809cc8cce637a` for the final version after the registered amendment. Only the final version exists as a file; the earlier checksum documents the state of the plan when the R analysis was run. This is the expected behaviour of an amended plan, not a missing file.

## Why the plans are not annotated

These files are released exactly as they were when their checksums were recorded. Nothing has been appended to them — not even a clarifying note — because any edit would change the SHA-256 and break the correspondence with Table S1. Explanations therefore live in this index, in `../ANALYSES.md` and in `../README.md`, never inside a plan.

## Other project documents

`other_documents/` holds documents from the wider project that are **not** among the Table S1 plans for this paper. They are included so that the release is complete and nothing that was written or run is hidden:

| File | What it is |
|---|---|
| `新研究_靶点优先级筛选_预设方案.md` | Plan for a separate target-prioritisation screen. |
| `第3-5层判定标准草案.md` | Draft decision criteria from an earlier phase. |
| `第一阶段图表清单.md`, `第一阶段结果汇总.md` | Figure list and results summary from phase 1. |
| `STAS分子评分_预设方案_冻结_2026-10-11.md`, `STAS分子评分_预设方案_草稿_2026-10-11.md` | Plan for an unrelated STAS scoring sub-project. The draft and the frozen version are byte-identical (same SHA-256): the plan was frozen without changes, and both names were kept in the project's own records. |

The KAC cross-cohort pooling plan was listed here in versions 1.0.0 and 1.0.1. It was hash-locked before its analysis was run (`scripts/60_kac_pooled.py`) and is therefore a plan used in this study; from version 1.1.0 it is a Table S1 row and is in `plans/`.
