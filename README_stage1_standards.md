# First-stage standards: what exists and what it can show

The first analysis stage of this study — nuclei quality control and annotation, and the paired pseudobulk comparisons in the discovery cohort with their replication in GSE164789 and GSE189357, including the myeloid comparison that nominated SPP1 — was governed by written standards, not by a hash-locked plan file. This note lists what was released and states plainly what it can and cannot establish. It corresponds to the "Stage 1" row of Supplementary Table S1.

## What the standards required

The comparison that nominated SPP1 was a paired pseudobulk comparison of invasive with precursor lesions within each cell population (pydeseq2, ~ patient + lesion type, FDR < 0.05, at least 50 cells per side), followed by source and ambient-RNA filters, with external replication counted when a gene reached one-sided BH q < 0.05 in the same direction in at least one of the two validation cohorts (GSE164789, GSE189357). SPP1 met this rule in both cohorts and had the largest discovery-cohort increase among myeloid genes that met it (COMP had a larger increase but did not replicate) (`scripts/14_microenv_layer1.py`; `results/14_microenv/`). Every later test of SPP1 reported in the paper was conducted under a hash-locked plan (see `plans/INDEX.md` and `ANALYSES.md`).

## Files in `stage1_standards/`

| File | What it is |
|---|---|
| `stage1_standards_excerpts.md` | Verbatim excerpts of the standards as recorded in the project log (`docs/项目进展与交接.md`), with log line numbers and the SHA-256 of the log at extraction (21 October 2026). |
| `补充分析判定标准草案.md` | The draft of the supplementary-analysis standards (microenvironment comparisons, external validation rule, exploratory proliferation adjustment) as submitted for confirmation on 8 October 2026; the confirmed version, with six changes, is the log excerpt. |
| `research_design_document_rev26_exported_2026-10-21.md` | The study's research design document (a shared online document, revision 26), exported on 21 October 2026. It records the stage-1 objectives and decision rules dated 8 October 2026 and later progress notes. No checksum of any earlier revision exists. |

The earlier draft of the layer 3–5 standards is in `other_documents/第3-5层判定标准草案.md`.

## What these files can and cannot show

- **No checksum was recorded for any stage-1 standard before the analyses were run.** The project directory was not under version control, so there is no earlier commit to point to. A checksum or commit made now would carry today's date and could not serve as evidence of timing; none has been constructed.
- **The order in time — standards written before the analyses — rests on the project log alone.** The log is a working document that is appended to and was also annotated in place (for example, checklist items were later marked as completed), so its current text is not a contemporaneous snapshot.
- **File modification times do not establish the order either.** The draft `补充分析判定标准草案.md` has a modification time of 8 October 2026, 22:30 (China Standard Time), later than the run of `scripts/14_microenv_layer1.py` (21:48 the same day). The draft carries a note added when the first stage was finalised ("注（2026-10-08 定稿时补充）"), which accounts for the later time; the research design document records the same finalisation at 22:30. Neither time says when the standards themselves were first written.
- The research design document is a live online document; its revision history is not exported here and is not used for any claim about timing.

For these reasons the paper describes the first-stage standards as "written and dated in the project log before the analyses were run; not hash-locked", and confines the claim of hash-locked plans to the second stage onward.
