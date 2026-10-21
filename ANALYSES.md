# Analysis register

Every analysis reported in the paper or its Supplementary Information, mapped to the plan or record that governs it, the script that ran it and the result directory the reported values come from. The register is built from Supplementary Table S1 and the deviation record (`plans/偏离与修订记录_外审回应A1-A6.md`), not from the contents of this repository, so that a missing script would show up as a gap here.

## Stage 1: nomination (standards written in the project log, not hash-locked)

The first stage nominated the candidates that later stages tested. Its standards were recorded in the project log as written and dated (8 October 2026) before the analyses were run, but were not kept as hash-locked files; the order in time rests on the log alone (see `README_stage1_standards.md` and `stage1_standards/`). SPP1 was nominated here, as the myeloid gene with the largest paired increase, and its replication in the two single-cell cohorts was governed by the same standards (one-sided BH q < 0.05 in the same direction in at least one cohort).

| Analysis | Scripts | Governing text |
|---|---|---|
| Nuclei QC, annotation, epithelial subtyping | `01_snRNA_qc_annotate.py`, `01b_refine_annotation.py`, `02_epithelial_infercnv.py` (`_v1` kept) | research design document; project log |
| Paired pseudobulk comparison, epithelial (primary) | `03a_pairing_table.py`, `03_pseudobulk_paired.py`, `04_source_filter.py` (`_v1` kept) | log: "03 伪批量配对分析（预设）" |
| External validation cohorts: preparation and epithelial validation | `05_prepare_validation.py`, `05b_apply_overrides.py`, `06_external_validation.py` | research design document (layer-1 rule) |
| Epithelial candidate layers 3–5 and integration | `07_candidates.py`, `08_layer3_depmap.py`, `09_layer4_clinical.py`, `10_layer5_safety_druggable.py` (`_v1` kept), `11_integrate.py`, `12_layer3_lincs.py`, `13_indirect_path.py` (`_v1` identical) | log: layer 3–5 and supplementary-analysis standards |
| Microenvironment comparisons incl. myeloid SPP1, with external replication | `14_microenv_layer1.py`, `15_microenv_layers45.py` (`_v1` kept), `17_ligand_receptor.py` | log: supplementary-analysis standards |
| Proliferation-adjusted Cox (exploratory; labelled post hoc in Table S1) | `16_prolif_adjusted_cox.py` | log: supplementary-analysis standards |
| Stage sensitivity (17 pairs) | `18_sensitivity_stageI.py` | log: sensitivity analysis entry |

## Stage 2 onward: hash-locked plans

The plan file named in each script's header, as found in the released scripts. Scripts marked † do not recompute the plan checksum at run time; for them the checksum recorded in Table S1 before the analysis was run is the only safeguard (stated in the paper's Methods).

| Table S1 plan | Scripts |
|---|---|
| Part 2: precursor-stage target screen | `21_stage2_discovery.py`†, `22_stage2_validation_layers.py`† |
| Part 3: spatial analysis | `70_build_spatial_reference.py`, `71_fit_reference_model.py`†, `90_merged_reference.py`†, `91_identifiability_sim.py`†, `92_c2l_map.py`†, `93_epoch_check.py`†, `96_spp1_niche_prognosis.py`† (state-score prognosis); helper `c2l_common.py` |
| Part 4: normal→precursor drivers | `40_step0_KAC_counts.py`†, `40_stepwise_tests.py`†, `41_state_matched_KAC.py`†, `43_hypothesis_TIM3.py`† |
| KAC fraction, cross-cohort pooling | `60_kac_pooled.py` |
| Part 5: human genetics | `95_mr_coloc.py`, `97_mr_spp1_macrophage.py` |
| Spatial gene-level exploratory | `80_spatial_gene_exploratory.py` |
| SPP1 source diagnostic | `81_spp1_source_diagnostic.py` |
| Sensitivity and protein validation | `98_niche_sensitivity.py`, `99_il1b_power.py`, `100_cptac_spp1.py` |
| Xenium independent validation | `101_xenium_validation.py` |
| Prognosis robustness and Xenium spillover | `102_prognosis_robustness.py`, `103_xenium_spillover.py` |
| Supplementary checks T | `104_checkT_xenium.py`, `105_checkT_hypoxia.py` |
| Hypoxia H1/H2 | `106_hypoxia_H1.py`, `107_hypoxia_H2.py` |
| Hypoxia H1c | `108_hypoxia_H1c.py` |
| TF activity, ligand–receptor, RCTD | `110_tf_activity.py`, `112a_rctd_prep.py`, `112b_rctd.R`†, `113_rctd_stats.py` (RCTD not executed; see Table S1) |

`102_prognosis_robustness.py` and `105_checkT_hypoxia.py` are kept as run; the defect in their bulk-cohort hypoxia gene sets is corrected in `review_A/A6c_corrections.py` (deviation D4). `prolif_genes.py` is a shared gene list. Scripts `50`–`56` and `70_stas_*`–`78_stas_*` belong to separate projects (`other_documents/`) and support no result in this paper.

## Analyses in response to external review (plan v1.3)

Scripts in `scripts/review_A/`; all share `_common.py`, which verifies the plan checksum at start-up and, from 20 October 2026, writes each run to its own `results/<name>/run_<time>/` directory. Results directories are relative to the project root and are not redistributed; the values reported are in the Supplementary Tables named.

| ID | Governing plan or record | Script | Result directory (valid run) | Reported in |
|---|---|---|---|---|
| A1 (a–c) | v1.3 §2-A1, §6.4; execution note A1c/A1d | `A1_macrophage_strict.py` | `review_A1/` (see P1, D3) | Tables S13a–S13c |
| A1d | v1.3 §2-A1; execution note A1c/A1d | `A1d_visium_markers.py` | `review_A1d/` | Table S13d |
| A2 | v1.3 §2-A2 | `A2_external_units.py` | `review_A2/` | Table S11a |
| A2b | v1.3 §2-A2 | `A2b_unit_scatter.py` | `review_A2b/run_20261021_054837/` | Supplementary Figure S1 |
| D2 (post hoc) | deviation record D2 | `A2s_hyy_sensitivity.py` | `review_A2/` | Table S11b |
| A3 | v1.3 §2-A3 | `A3_xenium_sensitivity.py` | `review_A3/` | Tables S3b, S3c |
| A4 | v1.3 §2-A4, §6.2 | `A4_stage_sensitivity.py` | `review_A4/` | Tables S12a, S12b |
| A5 | v1.3 §2-A5, §6.3; execution note A5 | `A5_external_hypoxia.py` | `review_A5/run_20261021_140428/` | Tables S15a, S15b |
| A6 | v1.3 §2-A6, §6.5 | `A6_prepare.py`, `A6_cox_diagnostics.R`, `A6_assemble.py` | `review_A6/run_20261021_100210/` | Tables S6b, S14a–S14e; Figure 6a |
| D1 | deviation record D1 | `D1_hallmark_source_check.py` ‡ | `review_D1/run_20261021_065321/` | Methods |
| D4 (correction) | deviation record D4 | `A6c_corrections.py` (and corrected models in `A6_prepare.py`) | `review_A6c/run_20261021_093506/` | Table S6b |
| D5 (post hoc) | deviation record D5 | `A6o_GSE37745_OS.py` | `review_A6o/run_20261021_140201/` | Table S14f |
| D6 (post hoc) | deviation record D6 | `D6_discovery_depth_matching.py` | `review_D6/run_20261021_141856/` | Table S15c |
| D7 (post hoc) | deviation record D7 | `D7_external_HIF1A_depth.py` | `review_D7/run_20261021_142723/` | Table S15d |
| P1, D3 | deviation record | — (process records for A1) | `review_A1_run1_old_script_2026-10-19_2342/` | Table S1 note |

Notes:

- ‡ `D1_hallmark_source_check.py` was recreated on 21 October 2026 from the original interactive comparison, which was run in a cloud container and whose script was not preserved; the recreated script reads the archived official gene-set file (checksum asserted) and gives an identical result. A first copy without the overlap-gene assertion was also run (`review_D1/run_20261021_065303/`, same result; script in `_archive/`).
- Run directory names written on the Linux virtual machine (`review_A2b`, `review_D1`) use UTC; all others use China Standard Time.
- Failed runs are kept as produced and are not listed above: `review_A5/run_20261021_140205/` (stopped before the summary step by a key-type error), `review_D6/run_20261021_141318/` (stopped before loading data by a NameError). The superseded scripts are in `scripts/review_A/_archive/`, and each failed run's `RUN_INFO.txt` names the cause.
- `review_A3_cloudrun_2026-10-19/` is an earlier execution of the same A3 script in a cloud environment; it agrees with `review_A3/` to floating-point precision.
- The review scripts locate the project with the fixed path `~/全部数据库/GGN_target` defined in `_common.py`, as when they were run; to rerun them elsewhere, edit `P` in `_common.py` (and `P` in `A6_cox_diagnostics.R`).
- Figure 6 (`figures/make_figs_npjpo.py`) reads `review_A6/run_20261021_100210/A6_pooled.tsv`, released as `source_data/review_A6__run_20261021_100210__A6_pooled.tsv`; Supplementary Figure S1 is drawn by `A2b_unit_scatter.py` from the A2 unit tables, released as `source_data/review_A2__A2_units_*.tsv`.
