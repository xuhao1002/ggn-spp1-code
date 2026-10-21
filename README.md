# A hypoxia-associated SPP1+ macrophage state distinguishes paired precursor and invasive lung adenocarcinoma lesions — analysis code

This repository contains the analysis code, the pre-specified analysis plans with their checksums, and the source values behind each figure panel for the above study. No new data were generated; every dataset analysed is public and is listed in the Data availability statement of the paper.

## Contents

| Folder | Contents |
|---|---|
| `scripts/` | Analysis scripts in the order they were run (numbered prefixes). Python unless suffixed `.R`. `scripts/review_A/` holds the analyses added in response to external review (plan v1.3, deviation records D1–D7); superseded or failed versions are kept in `scripts/review_A/_archive/`. |
| `plans/` | The 15 hash-locked analysis plans listed in Supplementary Table S1, one file per row, with superseded versions of the external-review plan (`superseded/`), execution notes (`execution_notes/`) and the deviation record. Each plan fixes hypotheses, statistical units, thresholds, tests, multiple-testing families and decision rules, and was hash-locked before the corresponding analysis was run. See `plans/INDEX.md`. |
| `stage1_standards/`, `README_stage1_standards.md` | The first-stage standards, which were written in the project log but not hash-locked, with an account of what the released material can and cannot show. |
| `ANALYSES.md` | Register mapping every analysis to its plan or record, script and result directory. |
| `other_documents/` | Working documents from the wider project that are **not** plans for this paper and are not in Table S1. Listed in `plans/INDEX.md`. |
| `figures/` | `make_figs_npjpo.py`, which produces the six submitted figures (panels labelled a, b, c …). |
| `source_data/` | The values plotted in each figure panel, as flat files. Names are the original result paths with `/` replaced by `__`. |
| `CHECKSUMS.tsv` | SHA-256 of every file in this repository, except `CHECKSUMS.tsv` itself. |
| `LICENSE` | Dual licence; see below. |
| `CITATION.cff` | Machine-readable citation metadata for this repository and for the paper. |
| `.zenodo.json` | Deposit metadata used when the archived release is minted on Zenodo. |

## Analysis plans and checksums

From the second stage onward, each analysis plan was written down before the analysis was run and the SHA-256 checksum of the plan file was recorded. Most analysis scripts recomputed the checksum of their plan at run time and refused to run on a mismatch; the scripts that did not are marked in `ANALYSES.md`, and for them the recorded checksum is the only safeguard. The first-stage standards were written in the project log and not hash-locked; see `README_stage1_standards.md`. Plans were written in stages, so later plans were informed by earlier results; the checksums establish that the content of each plan was fixed before its analysis, not the order in time, which rests on local records. Amendments, deviations and post hoc analyses are recorded at the end of each plan — for the external-review plan v1.3, which was locked and is not edited, in the separate deviation record `plans/偏离与修订记录_外审回应A1-A6.md` — and are listed in Supplementary Table S1.

`plans/INDEX.md` maps every hash-locked Table S1 row to its file and shows that all 15 recorded checksums match the released files byte for byte. The plans are released unmodified: no note has been appended to any of them, because an edit would change the checksum and break that correspondence.

## Reproducing the figures

`figures/make_figs_npjpo.py` takes no hard-coded paths. Four environment variables control it, all optional:

| Variable | Default | Purpose |
|---|---|---|
| `SPP1_SOURCE_DATA` | `../source_data` next to the script | The released source-data files. |
| `SPP1_RESULTS` | unset | The original results tree, if you have it. |
| `SPP1_RAW` | unset | The original raw-input tree, if you have it. |
| `SPP1_OUT` | `../figures_out` | Where the figures are written. |

Every input is looked up first under `SPP1_RESULTS` / `SPP1_RAW` using its original relative path, then in `SPP1_SOURCE_DATA` under its flattened name. A missing file raises an error naming the file and both options.

From a clone, with nothing configured:

```
cd figures && python3 make_figs_npjpo.py
```

This regenerates Figures 1, 2, 3, 5 and 6, and panels a and e of Figure 4, from the released source data alone. Supplementary Figure S1 is drawn by `scripts/review_A/A2b_unit_scatter.py` from the A2 unit tables (released in `source_data/`). The Xenium source table used for Table S3 and Figure 4e is regenerated from the A3 output by `scripts/review_A/xenium_source_fixed_from_A3.py`. Panels b, c and d of Figure 4 read the CPTAC-LUAD proteome tables, which are not redistributed here (see below); download them from the Proteomic Data Commons and point `SPP1_RAW` at the directory holding them, after which all six figures are reproduced:

```
SPP1_RAW=/path/to/raw python3 make_figs_npjpo.py
```

With that one input supplied, the script reproduces the six submitted figures byte-identically.

## What is not included, and why

- **GWAS summary statistics.** The Japanese lung adenocarcinoma GWAS sets were obtained from the J-CGE download site of the National Cancer Center, Japan (https://epi.ncc.go.jp/jcge/) under its conditions of use, which prohibit redistribution. They are available to any researcher who completes the site's application form.
- **The CollecTRI regulatory network.** Retrieved through OmniPath on 13 October 2026. The retrieval step is in the scripts and the version and date are recorded, so the identical snapshot can be regenerated; the network itself is not redistributed here because of its source licence.
- **Primary sequencing, spatial and proteomic data.** All are public under the accessions given in the paper (GEO, SSBD, PDC, GDC, GTEx, HPA, eQTL Catalogue) and are not duplicated here.
- **Two large inputs read by the figure scripts** — the CPTAC proteome tables and the GTEx median-TPM file — are downloaded from their sources rather than redistributed.

## Two pairs of identical files

`CHECKSUMS.tsv` contains two pairs of files with the same SHA-256. Both are intentional, and neither is a packaging error:

- `scripts/13_indirect_path.py` and `scripts/13_indirect_path_v1.py` are byte-identical. The `_v1` suffix marks the version kept when the script was revisited; the revision ended up making no change, and both names were retained in the project's own records.
- `other_documents/STAS分子评分_预设方案_冻结_2026-10-11.md` and `..._草稿_2026-10-11.md` are byte-identical for the same reason: the draft was frozen without edits. Both belong to an unrelated sub-project.

## Software

Python 3.11; scanpy 1.11.5, anndata 0.12.19, harmonypy 2.1.0, pydeseq2 0.5.4, statsmodels 0.15.0, scipy 1.17.1, gseapy, decoupler, cell2location 0.1.5. R 4.5.2 with survival 3.8.3 and coxphf 1.13.4 is used for the Cox diagnostics and Firth models of the external-review analysis A6 (`scripts/review_A/A6_cox_diagnostics.R`) and for the RCTD helper script. The external-review scripts locate the project at the fixed path `~/全部数据库/GGN_target` set in `scripts/review_A/_common.py`, as when they were run; to rerun them elsewhere, edit that path (and the one in `A6_cox_diagnostics.R`).

## Citing

If you use this code, please cite both the archived release and the paper. `CITATION.cff` holds the metadata in machine-readable form; the journal reference and its DOI are added there on acceptance.

Xu H., Zhang H., Zhang L. *Analysis code for: A hypoxia-associated SPP1+ macrophage state distinguishes paired precursor and invasive lung adenocarcinoma lesions*, version 1.2.0 (2026). The Second Affiliated Hospital of Harbin Medical University, Harbin, China. https://github.com/xuhao1002/ggn-spp1-code. Archived at Zenodo: https://doi.org/10.5281/zenodo.23359094 (this DOI always resolves to the latest archived version). Version 1.2.0 is the version of record for the manuscript; later changes, including any made in response to peer review, are released as new versions.

## Licence

Dual-licensed:

- **Code** (`scripts/`, `figures/`) — MIT License.
- **Text and data** (`plans/`, `other_documents/`, `source_data/`, `README.md`, `CHECKSUMS.tsv`) — Creative Commons Attribution 4.0 International (CC BY 4.0).

See `LICENSE` for the full terms. Third-party files keep the licence of their source and are not relicensed here.
