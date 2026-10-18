# cellprograms — analysis reproducibility bundle

Standalone export of the analysis code, data builders, checkpoint tables, and
figures behind the four data-analysis vignettes of the `cellprograms` R
package (cell-type-specific transcriptional programs via empirical-Bayes
matrix factorization). The package source itself lives at
https://github.com/mdmanurung/cellprograms; this bundle contains everything
needed to reproduce the analyses **bit by bit**, with committed checkpoint
tables to verify every stage against.

Live vignettes:

1. [Real-data case study (Stephenson + COMBAT COVID PBMC)](https://mdmanurung.github.io/cellprograms/articles/real_data_case_study.html)
2. [Immune-aging programs across 8 global PBMC cohorts](https://mdmanurung.github.io/cellprograms/articles/immune_aging_cross_cohort.html)
3. [DIVAS COVID-19 multi-omics reanalysis](https://mdmanurung.github.io/cellprograms/articles/divas_case_study.html)
4. [Compared to MOFAcellulaR](https://mdmanurung.github.io/cellprograms/articles/mofacellular_comparison.html)

## Layout

| Path | Contents |
|---|---|
| `vignettes/articles/` | The four article sources (`.Rmd`) + rendered `.md` + all figures (`figures/`, `figure/`) |
| `analyses/immune_aging/` | Immune-aging pipeline: 20 scripts, checkpoint tables (`tables/`), README (vignette 2) |
| `benchmarks/` | All pre-registered benchmark/study folders: `cmp_detection/`, `cmp_bench/`, `divas_comparison/`, `mofacellular_comparison/`, `biomni_replication/`, `model_selection/`, `coordination/` (vignettes 1, 3, 4) |
| `experiments/README.md` | Study index tying each benchmark folder to scripts/figures/reports |
| `reports/` | Analysis narratives (`immune_aging_analysis.md` is the full immune-aging report) |
| `environment/` | `install.R` (package install) + `sessionInfo.txt` (exact as-run versions) |

Large data (raw h5ads, pseudobulk, fitted models) is **not** bundled — every
stage is regenerated from public sources by the scripts below, and the
committed checkpoint tables let you verify each step.

## Environment

```bash
Rscript environment/install.R        # R >= 4.3; installs cellprograms + pinned deps
```

As-run environment: R 4.4.2, `flashier` 1.0.60 (pinned
`willwerscheid/flashier@fd44811`), `flashier.utils` 1.1.0.9000 — full details
in `environment/sessionInfo.txt`. The Python data builder needs
`scanpy`, `anndata`, `pyarrow`, `cellxgene-census` (see
`benchmarks/biomni_replication/env_fix.yml` for the conda recipe).

**Registered operating point** (immune-aging discovery/refit): `max_factors
= 10`, `loading_prior = "point_laplace"`, `var_type = 1`,
`features = "variable"` (top-2000), `feature_scale = "logcounts"`, seed
`20261016`, no covariate residualization, `canonicalize_programs()` after
fit, log1p-CP10K from `.raw/X`. (The CMP-detection and case-study fits use
`max_factors = 30` — selected by the pre-registered benchmark
`benchmarks/cmp_detection/`.)

---

## Vignette 1 — real_data_case_study (Stephenson + COMBAT COVID PBMC)

**Question**: program discovery, cross-cell-type sharing, V(D)J/genotype
artifact masking, and CMP detection in two real COVID PBMC cohorts.

**Bit-by-bit**:

1. **Pseudobulk both cohorts** from CELLxGENE Census (pinned `2025-11-08`):
   ```bash
   python benchmarks/biomni_replication/01_pseudobulk.py --cohort stephenson
   python benchmarks/biomni_replication/01_pseudobulk.py --cohort combat
   ```
   Stephenson = Census collection `ddfad306-714d-4cc0-9985-d9072820c530`
   (also ArrayExpress E-MTAB-8731), 12 cell types, 72–117 donors/type;
   COMBAT = collection `8f126edf-5405-4731-8374-b5ce11f53e82`. Outputs:
   `data/<cohort>/mat/ct_<ct>.csv` (donor x gene), `donor_meta.csv`,
   `composition.csv`.
2. **Fit + share + mask + detect**, following the article's recipe chunks
   (`vignettes/articles/real_data_case_study.Rmd`, all shown with `eval =
   FALSE`): K = 30 factors per cell type, point-Laplace loadings,
   residualize institute + sex, variable-genes panel (2000), then
   `program_sharing()` (mashr), `flag_artifact_genes()` + masked refit, and
   `find_coordinated_programs()`.
3. **As-run campaign code** (hard-coded paths to the original sandbox —
   edit `config.R` / path constants first, as noted in that folder's README):
   `benchmarks/cmp_bench/` (real K30 fits, sharing, TF triage) and
   `benchmarks/cmp_detection/` (pre-registered caller benchmark whose
   outcome fixed K = 30).
4. **Checkpoints**: `vignettes/articles/figures/fig_sharing_celltypes.png`,
   `fig_masked_vs_unmasked.png`, `fig_cmp_network.png`,
   `fig_cmp_context.png`.

## Vignette 2 — immune_aging_cross_cohort (8 cohorts, 2,609 donors)

**Question**: age-associated programs, cross-cohort conservation, and the
composition-level aging axis across 8 healthy PBMC cohorts.

**Bit-by-bit** (all scripts in `analyses/immune_aging/scripts/`; see that
folder's README for the full command table):

1. **Download the 8 cohort h5ads** (~127 GB total) from the Terekhova et al.
   portal, base `https://artyomovlab.wustl.edu/scn-m/var/datasets/curated/maksim/GZMB_GZMK`:
   `ABF300_new_annotation_clean/files/ABF300_reannotated.h5ad`,
   `AIDA_new_annotation_clean/files/AIDA_reannotated.h5ad`,
   `CIMA_new_annotation_clean/files/CIMA_reannotated.h5ad`,
   `KCL_annotation/files/KCL_reannotated.h5ad`,
   `OneK1K_new_annotation_clean/files/OneK1K_reannotated.h5ad`,
   `SoundLife_new_annotation_clean/files/SoundLife_reannotated.h5ad`,
   `SPAC_new_annotation_clean/files/SPAC_filtered_donors_reannotated.h5ad`,
   `UCSF_new_annotation_clean/files/UCSF_reannotated.h5ad`.
   `download_lane.sh` is the as-run wrapper (downloads + pseudobulks +
   composition counts as each file lands; edit its `/workspace/...` paths).
   Donor counts: OneK1K 981, AIDA 625, CIMA 421, ABF300 166, KCL 188,
   UCSF 99, SoundLife 96, SPAC 33.
2. **Pseudobulk**: `pseudobulk_h5ad.py <h5ad> <out_dir> <COHORT> --min-cells 10`
   per cohort (donor x cell-type pseudobulk from `.raw/X`, log1p-CP10K).
3. **Fit** each cohort at the registered operating point (`fit_pilot.R` /
   per-cohort driver), then `canonicalize_programs()`.
4. **Age association** per cohort: `assoc_programs.R` (per cell type,
   `score ~ age10 + sex [+ CMV]`, FDR across programs).
5. **Cross-cohort transfer**: pick a discovery cohort (registered: ABF300),
   `project_cohort.R` for every query cohort (MR-1 projection),
   `assoc_concordance.R`, `assemble_conservation.R`, `gene_overlap.R`,
   `projected_age_replication_gen.R`.
6. **Composition axis**: `composition_counts.py` (per-donor per-type counts)
   -> `composition_aging.R` (proportion trajectories + GZMK+/GZMB+ Tem ratio,
   zero-filled regressions) -> `make_composition_figures.R`.
7. **Discovery sensitivity**: refreeze on OneK1K (981 donors), reproject all
   7 remaining cohorts, `discovery_sensitivity.R`.
8. **Checkpoints** (verify every stage):
   `analyses/immune_aging/tables/` — `age_by_cohort.csv`,
   `age_by_celltype_consistency.csv`, `conservation_panel.csv`,
   `projected_age_replication.csv`, `gene_overlap_top20.csv`,
   `composition/composition_age_trajectories.csv`,
   `composition/tem_ratio_age_cmv.csv`,
   `sensitivity/discovery_sensitivity.csv`,
   `sensitivity/projected_age_replication_onek1k_disc.csv`. Figures:
   `vignettes/articles/figures/fig_aging_*.png` (7 files).

Headline expectations: CD8+/CD4+ RTE decline in 8/8 cohorts (thymic
involution); matched age-profile Spearman null-to-weak under both discovery
freezes; Tem-ratio age trajectory cohort-inconsistent; CMV lowers the ratio
in SoundLife (-1.07 log, p = 2.5e-11).

## Vignette 3 — divas_case_study (DIVAS COVID multi-omics reanalysis)

**Bit-by-bit** (scripts in `benchmarks/divas_comparison/scripts/`, run in
order):

1. `01_download_data.sh` — input blocks + metadata from
   `ByronSyun/DIVAS_COVID19_CaseStudy` (GitHub) and the precomputed DIVAS
   6-block results from Zenodo record `17430294` (419 MB).
2. `02_fit_cellprograms.R` — per-block cellprograms fits at DIVAS parity.
3. `03_cmp_sharing.R` — CMPs + sharing on the multi-omics blocks.
4. `04_extract_divas.R`, `05_compare.R` — head-to-head against the precomputed
   DIVAS components (91: 44 individual, 46 partial, 1 joint).
5. `06_figures_extra.R`, `07_sensitivity_pairing.R` — figures + pairing
   sensitivity.
6. **Checkpoints**: `vignettes/articles/figures/fig_breadth_concordance.png`,
   `fig_concordance_heatmap.png`, `fig_effect_agreement.png`.

## Vignette 4 — mofacellular_comparison (vs MOFAcellulaR, Stephenson)

**Bit-by-bit** (scripts in `benchmarks/mofacellular_comparison/scripts/`;
that folder's README is the authoritative run book):

1. **Pseudobulk**: same Census builder as vignette 1, `--cohort stephenson`;
   assemble the `stephenson_pb/semi_base.rds` object the scripts expect
   (structure documented in `benchmarks/cmp_detection/scripts/real_data.R`:
   per cell type `list(Y = donor x gene counts, ...)` + `donor_meta` with
   `donor_id, site, sex, Status_on_day_collection_summary`).
2. `01_fit_masked.R` — masked K30 fit (274-gene VDJ + genotype-marker
   blacklist), ~25 min.
3. `02_cmp.R` — CMP detection on the fit, ~1 min.
4. `03_mofa.R` — MOFAcellulaR under parity (needs `MOFAcellulaR` + `MOFA2`
   with `mofapy2` via reticulate), ~5 min.
5. `04_fit_excl.R` + `02_cmp.R` + `03_mofa.R` — artifact-exclusion
   sensitivity rerun (`exclude_artifacts = TRUE`, 313-gene blacklist).
6. `05_build_cache.R` — rebuilds the committed caches
   `vignettes/articles/mofacellular_stephenson_cache{,_excl}.rds` that the
   article renders from.
7. **Checkpoints**: the two cache `.rds` files +
   `vignettes/articles/figure/*.png` (8 knitr-output figures referenced by
   the rendered article).

## Scope notes

- Session-internal planning docs (`ROADMAP.md`, `HANDOFF.md`, `plans/`) and
  the package dev checklist stay in the package repository; this bundle
  carries the analysis-facing reports only.
- `benchmarks/model_selection/` result bulk (per-dataset CSVs) is regenerable
  by its tracked scripts and gitignored, matching the source repo.
- The Stephenson `semi_base.rds` assembler predates the committed benchmark
  layout; its exact object structure is documented in
  `benchmarks/cmp_detection/scripts/real_data.R` and the Census builder in
  `benchmarks/biomni_replication/01_pseudobulk.py` reproduces the underlying
  pseudobulk from pinned public data.

## Update 2026-10-19 — spec-v2 validated package state

The bundle's analyses were run against `cellprograms@413e290`. The package
has since completed its spec-v2 milestone queue (calibrated CMP association
suite, program-sharing sensitivity gates, cross-cohort replication
workflow, metacells folded in from the former `cellprogramsMC` sub-package,
projection query-set dependence documented). All changes are
backward-compatible with the analysis code here — legacy function names and
the registered operating point are unchanged.

`environment/install.R` now pins the validated state:
`flashier.utils@738eba1`, `cellprograms@cce06d9` (limitation registers:
`reports/audit/LIMITATIONS.md` in cellprograms, `reports/LIMITATIONS.md` in
flashier-utils). `environment/sessionInfo.txt` remains the historical
as-run record for the published vignettes and is intentionally not
regenerated.

New since the original bundle: metacells (`MetacellsByGroupsSEACells`,
mcRigor QC) in the main package; cross-cohort replication case-study
article (`articles/cross_cohort_replication.html`); confirm-scale
benchmark evidence (`benchmarks/cmp_calibration/`,
`benchmarks/program_sharing_sensitivity/` in the package repo).
