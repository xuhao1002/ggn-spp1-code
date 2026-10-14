# A hypoxia-associated SPP1+ macrophage state marks the transition from precursor lesions to invasive lung adenocarcinoma — analysis code

This repository contains the analysis code, the pre-specified analysis plans with their checksums, and the source values behind each figure panel for the above study. No new data were generated; every dataset analysed is public and is listed in the Data availability statement of the paper.

## Contents

| Folder | Contents |
|---|---|
| `scripts/` | Analysis scripts in the order they were run (numbered prefixes). Python unless suffixed `.R`. |
| `plans/` | The 13 analysis plans listed in Supplementary Table S1, one file per row. Each plan fixes hypotheses, statistical units, thresholds, tests, multiple-testing families and decision rules, and was hash-locked before the corresponding analysis was run. See `plans/INDEX.md`. |
| `other_documents/` | Working documents from the wider project that are **not** plans for this paper and are not in Table S1. Listed in `plans/INDEX.md`. |
| `figures/` | `make_figs_npjpo.py`, which produces the six submitted figures (panels labelled a, b, c …). |
| `source_data/` | The values plotted in each figure panel, as flat files. Names are the original result paths with `/` replaced by `__`. |
| `CHECKSUMS.tsv` | SHA-256 of every file in this repository, except `CHECKSUMS.tsv` itself. |
| `LICENSE` | Dual licence; see below. |
| `CITATION.cff` | Machine-readable citation metadata for this repository and for the paper. |
| `.zenodo.json` | Deposit metadata used when the archived release is minted on Zenodo. |

## Analysis plans and checksums

Each analysis was written down before it was run, and the SHA-256 checksum of the plan file was recorded. At run time the analysis script recomputed the checksum of its plan and refused to run on a mismatch. Plans were written in stages, so later plans were informed by earlier results; the checksums establish that the content of each plan was fixed before its analysis, not the order in time, which rests on local records. Amendments, deviations and post hoc analyses are recorded at the end of each plan and are listed in Supplementary Table S1.

`plans/INDEX.md` maps every Table S1 row to its file and shows that all 13 recorded checksums match the released files byte for byte. The plans are released unmodified: no note has been appended to any of them, because an edit would change the checksum and break that correspondence.

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

This regenerates Figures 1, 2, 3, 5 and 6, and panels a and e of Figure 4, from the released source data alone. Panels b, c and d of Figure 4 read the CPTAC-LUAD proteome tables, which are not redistributed here (see below); download them from the Proteomic Data Commons and point `SPP1_RAW` at the directory holding them, after which all six figures are reproduced:

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

Python 3.11; scanpy 1.11.5, anndata 0.12.19, harmonypy 2.1.0, pydeseq2 0.5.4, statsmodels 0.15.0, scipy 1.17.1, gseapy, decoupler, cell2location 0.1.5. R is used only for the spatial deconvolution helper scripts.

## Citing

If you use this code, please cite both the archived release and the paper. `CITATION.cff` holds the metadata in machine-readable form; the journal reference and its DOI are added there on acceptance.

Xu H., Zhang H., Zhang L. *Analysis code for: A hypoxia-associated SPP1+ macrophage state marks the transition from precursor lesions to invasive lung adenocarcinoma*, version 1.0.1 (2026). The Second Affiliated Hospital of Harbin Medical University, Harbin, China. https://github.com/xuhao1002/ggn-spp1-code. Archived at Zenodo: https://doi.org/10.5281/zenodo.23359094 (this DOI always resolves to the latest archived version).

## Licence

Dual-licensed:

- **Code** (`scripts/`, `figures/`) — MIT License.
- **Text and data** (`plans/`, `other_documents/`, `source_data/`, `README.md`, `CHECKSUMS.tsv`) — Creative Commons Attribution 4.0 International (CC BY 4.0).

See `LICENSE` for the full terms. Third-party files keep the licence of their source and are not relicensed here.
