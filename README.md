# Separating tobacco exposure from ancestry in East Asian and Western lung adenocarcinoma

Code, harmonised data and results for a six-cohort whole-exome comparison of
East Asian (EAS) and Western lung adenocarcinoma (LUAD), stratified by smoking.

| Where | What | Access |
|---|---|---|
| This repository | code, harmonised public-cohort data, results, figures | open |
| Zenodo v1, <https://doi.org/10.5281/zenodo.21201441> | open data release | open |
| Zenodo v2, <https://doi.org/10.5281/zenodo.23050408> | complete in-house (Guo2026) data | restricted, on request |

![Graphical abstract](figures/graphical_abstract.jpg)

## The story in three results

The analysis set is 1,453 tumours from six cohorts: 420 EAS never-smokers, 200
EAS ever-smokers, 239 Western never-smokers and 594 Western ever-smokers.

1. **EAS ever-smokers resemble Western never-smokers.** Median TMB was 2.18
   mutations/Mb in Western never-smokers (n = 239) and 2.39 mutations/Mb in EAS
   ever-smokers (n = 190) (median difference = -0.20 mutations/Mb; adjusted
   P = 0.97), none of the 15 canonical LUAD driver genes differed after
   multiple-testing correction, and KRAS was mutated in 17.6% and 17.8% of these
   two groups (OR = 1.02, 95% CI 0.60-1.75; adjusted P = 1.00). Numbers in
   `results/01_four_groups/tables/`: TMB medians and n in `03_tmb_contrasts.csv`;
   driver frequencies, odds ratios and adjusted P in `06_driver_genomewide.csv`.
   The TMB adjusted P (0.97) is from `15_western_never_vs_eas_ever.csv`, which
   adjusts jointly over all 8,524 tests of this cross-check (TMB, 89 signatures,
   six substitution classes, 8,428 genes). `03_tmb_contrasts.csv` adjusts each
   TMB contrast on its own, so its `p_adjusted` for this row is the unadjusted
   P (0.045). The odds-ratio CI is not written by any script (see "Multiple
   testing and odds ratios" below).
2. **Accounting for smoking narrows the ancestry signal.** In the pooled
   comparison 2,063 of 8,428 genes differed between EAS and Western tumours, but
   only 13 differed in both the never-smoker and the ever-smoker comparison, and
   EGFR and KRAS were the only established LUAD drivers among them. Numbers:
   `results/03_pooled_vs_stratified/tables/04_pooled_vs_stratified.csv`.
3. **The smoking-associated changes in SBS4 and SBS1 differ between the
   groups.** Median SBS4 contribution rose from 0% to 6.2% in EAS tumours and
   from 0.4% to 29.6% in Western tumours (difference in the smoking-associated
   change 23.0 percentage points, 95% CI 11.7-30.5; adjusted P = 3.75 × 10<sup>-3</sup>),
   while SBS1 fell from 6.9% to 2.1% in EAS tumours (adjusted P = 9.67 × 10<sup>-8</sup>)
   but changed little in Western tumours (2.3% to 2.0%; adjusted P = 0.30).
   Numbers in `results/01_four_groups/tables/`: `09_signature_contrasts.csv`,
   `10_signature_interaction.csv`.

## Figure 1, panel by panel

| Panel | What it shows | Script | Source file in `results/figure1_panels/` |
|---|---|---|---|
| 1A, left | Never- and ever-smokers in each of the six cohorts | `01_four_groups.Rmd` | `fig1a_left_cohort_composition.pdf` |
| 1A, right | TMB in the four ancestry-by-smoking groups | `01_four_groups.Rmd` | `fig1a_right_tmb.pdf` |
| 1B, top | Mutation frequency (95% CI) of the 15 canonical LUAD drivers in the four groups | `01_four_groups.Rmd` | `fig1b_top_driver_frequency.pdf` |
| 1B, bottom | Genome-wide EAS vs Western scan (8,428 genes), never-smokers | `01_four_groups.Rmd` | `fig1b_bottom_never_smokers.pdf` |
| 1B, bottom | Genome-wide EAS vs Western scan (8,428 genes), ever-smokers | `01_four_groups.Rmd` | `fig1b_bottom_ever_smokers.pdf` |
| 1C | C>A share of all single-base substitutions in the four groups | `01_four_groups.Rmd` | `fig1c_ca_fraction.pdf` |
| 1C | Relative contribution of the labelled signatures in the four groups | `01_four_groups.Rmd` | `fig1c_signatures.pdf` |
| 1D, left | Smoking-associated change estimated within each of four cohorts | `02_within_cohort.Rmd` | `fig1d_left_within_cohort.pdf` |
| 1D, right | Driver-gene frequencies compared between every pair of the four cohorts, by smoking group | `04_cohort_pairs.Rmd` | `fig1d_right_cohort_pairs.pdf` |

The final Figure 1 (`figures/Figure1.pdf`, `figures/Figure1.jpg`) was assembled
from these panels in Adobe Illustrator. Illustrator was used to select facets
(SBS1 and SBS4 from `fig1c_signatures.pdf`; TMB, C>A, SBS1 and SBS4 from
`fig1d_left_within_cohort.pdf`; EGFR, KEAP1, KRAS and TP53 from
`fig1d_right_cohort_pairs.pdf`) and to add annotations such as the comparison
brackets. No code writes to `figures/`.

The numbers quoted for Figure 1D are in
`results/02_within_cohort/tables/01_within_cohort_all_layers.csv` (left) and
`results/04_cohort_pairs/tables/03_pairwise_driver.csv` (right).

## Supplementary Tables S1-S3

`06_supplementary_tables.Rmd` assembles the tables from saved results; it runs
no test and no bootstrap. Each file is a long CSV with one table cell per line
(`part`, `part_title`, `row`, `column`, `value`).

| Table | File in `results/supplementary_tables/` | Built by 06 from |
|---|---|---|
| S1. Cohort composition and endpoint-specific denominators | `TableS1_cohort_composition_and_denominators.csv` | `data/analysis_metadata.csv`; `results/01_four_groups/tables/01_cell_counts.csv`, `07_driver_hit_counts.csv`, `08_driver_frequency_by_cell.csv` |
| S2. Sensitivity analyses for hypermutation | `TableS2_hypermutation_sensitivity.csv` | `results/05_sensitivity/tables/01_hypermutation_balance.csv`, `03_hypermutation_sensitivity.csv` |
| S3. Western cohorts' ethnicity and sensitivity to the Western-group definition | `TableS3_western_group_definition.csv` | `data/analysis_metadata.csv`, `data/TCGA2014_info.csv`, `data/TRACERx2023_info.csv`; `results/05_sensitivity/tables/05_ancestry_restriction.csv`, with the primary TMB and SBS4 rows taken from `03_hypermutation_sensitivity.csv` (scenario `full`) so that they match Table S2 |

## Extended analyses

These results are not shown in the Letter. They support its robustness
statements. Each folder holds `tables/` and `figures/`.

- `results/01_four_groups/`: the rest of the four-group analysis: age and sex,
  the TMB interaction, the gene scan in all five contrasts, every active
  signature and substitution class, the 96-context profile, and the EAS-ever
  versus Western-never cross-check across all layers.
- `results/02_within_cohort/`: ever- versus never-smoker estimates computed
  inside Chen2020, Xu2020, TCGA2014 and TRACERx2023 for every layer (TMB, all
  active signatures, six substitution classes) and every tested gene.
- `results/03_pooled_vs_stratified/`: the pooled EAS versus Western comparison
  that ignores smoking (32.3% versus 71.3% ever-smokers), and the
  pooled-versus-stratified gene table behind result 2.
- `results/04_cohort_pairs/`: differences between all six pairs of the four
  cohorts, for every continuous layer and every tested gene. The two
  within-ancestry pairs (Chen2020/Xu2020, TCGA2014/TRACERx2023) are negative
  controls for technical differences. With six pairs the smallest attainable
  rank-sum P is 0.133, so complete separation of cross- from within-ancestry
  pairs is the informative statistic.
- `results/05_sensitivity/`: hypermutation (four flag definitions, reruns
  without flagged tumours), leave-one-cohort-out, reported-ancestry
  restriction, variant-class audit, and the copy-number comparability audit.

## Data

`data/` holds the harmonised inputs written by `00_harmonise_cohorts.Rmd`:

- `<Cohort>_info.csv` (clinical: sex, age, smoking, reported race or ethnicity),
  `<Cohort>_snv.csv` (variant calls with a common variant class; five public
  cohorts, see below for Guo2026) and `<Cohort>_tmb.csv` (per-sample TMB).
- `signature/`: SBS96 spectra and COSMIC refits for the five cohorts with full
  mutation tables.
- `analysis_metadata.csv`: one row per harmonised tumour with the flags every
  script filters on (`AncestryExcluded`, `HasTMB`, `HasSignature`,
  `HasFullMAF`, `InPerCohortAnalysis`, four hypermutation flags); plus
  `variant_class_audit.csv`, `hypermutation_audit.csv` and `sessionInfo.txt`.

| Cohort | Group | Harmonised | Analysis set, never / ever | Source |
|---|---|---|---|---|
| Guo2026 (in-house) | EAS | 19 | 11 / 8 | this study |
| Chen2019 | EAS | 197 | 143 / 54 | supplementary tables of Chen H, et al. Nat Commun 2019 |
| Chen2020 | EAS | 305 | 189 / 112 | OncoSG, dataset "Lung Adenocarcinoma, GIS 2019", <https://src.gisapps.org/OncoSG/> |
| Xu2020 | EAS | 103 | 77 / 26 | supplementary tables of Xu JY, et al. Cell 2020 |
| TCGA2014 | Western | 595 | 214 / 373 | UCSC Xena, TCGA-LUAD, <https://xenabrowser.net/datapages/> |
| TRACERx2023 | Western | 248 | 25 / 221 | Zenodo record 7822002, <https://zenodo.org/records/7822002> |
| Total | | 1,467 | 1,453 | |

**The in-house cohort.** The complete Guo2026 data folder (`Guo2026_LUAD19/`,
as delivered by the sequencing pipeline, including the 19 per-patient somatic
MAFs in `00.maf/filtered_maf/`) is under restricted access at
<https://doi.org/10.5281/zenodo.23050408>. Request access on that Zenodo record
and state the intended use; access is granted for non-commercial academic
research. This repository does not include `data/Guo2026_snv.csv`, the
harmonised in-house variant table. To reproduce the published numbers exactly,
place the folder at `raw/Guo2026_LUAD19/` and run the Guo2026 chunk of
`00_harmonise_cohorts.Rmd` (section 2.4, after the setup and section 1 chunks),
which writes `data/Guo2026_snv.csv`. Without it the scripts run, warn, and leave
Guo2026 out of the mutation-frequency analyses; Guo2026 still enters the TMB
and signature analyses. The shipped `results/` were produced with the file
present, so a run without it overwrites them with different mutation-frequency
numbers. Raw sequencing data are available from the corresponding authors on
reasonable request.

### Raw data: `raw/` layout

No raw data are included. `00_harmonise_cohorts.Rmd` reads exactly these files,
relative to the repository root. The spreadsheet names are local; 00 expects
the sheets and column names used in the script.

```
raw/Chen2019_LUAD197/01_cohort_info.xlsx
raw/Chen2019_LUAD197/02_somatic_mutation.xlsx
raw/Chen2019_LUAD197/04_tmb.xlsx
raw/Chen2020_LUAD305/01_cohort_info.xlsx
raw/Chen2020_LUAD305/04_tmb.xlsx                            (sheet "final")
raw/Chen2020_LUAD305/gis031_lung_adeno_nccs_2020/snv_indel.maf
raw/Xu2020_LUAD103/01_cohort_info.xlsx
raw/Xu2020_LUAD103/02_somatic_mutation.xlsx                 (sheet "final")
raw/Guo2026_LUAD19/01_cohort.xlsx
raw/Guo2026_LUAD19/00.maf/filtered_maf/*.maf.gz             (19 MAFs, restricted record)
raw/TCGA2014_LUAD/TCGA-LUAD.clinical.tsv.gz
raw/TCGA2014_LUAD/TCGA-LUAD.somaticmutation_wxs.tsv.gz
raw/TRACERx2023_LUAD421/20221109_TRACERx421_all_tumour_df.rds
raw/TRACERx2023_LUAD421/20221109_TRACERx421_mutation_table.fst
raw/reference/COSMIC_Human_SBS-96_GRCh37_v3.6.tsv
raw/reference/COSMIC_Human_SBS-96_GRCh38_v3.6.tsv
```

The COSMIC v3.6 SBS96 reference signatures (GRCh37 and GRCh38) are from
<https://cancer.sanger.ac.uk/signatures/downloads/>.
`raw/Guo2026_LUAD19/01_cohort.xlsx` is part of the restricted record; its
harmonised content is `data/Guo2026_info.csv`. Running 00 writes `data/Guo2026_snv.csv`, which
`.gitignore` keeps out of version control.

### Terms of use for third-party data

The harmonised files in `data/` are derived from the published sources listed
above and remain under those sources' terms.

- `TRACERx2023_*` files in `data/` and `data/signature/` are derived from Zenodo
  record 7822002, which is released under an Academic-Use Only Data Licence:
  academic non-commercial use only. Commercial users must request access as that
  record describes.
- `Chen2019_*` and `Xu2020_*` (journal supplementary tables), `Chen2020_*`
  (OncoSG) and `TCGA2014_*` (UCSC Xena) files are subject to the terms of those
  sources.
- `Guo2026_*` files come from this study; its variant-level calls are
  restricted (see "The in-house cohort").

## How to run

Developed with R 4.3.1 and checked with R 4.6.1. Figures are written with
`cairo_pdf`. `capabilities("cairo")` can return `TRUE` even when the cairo
library cannot be loaded; on macOS, CRAN builds of R need XQuartz
(<https://www.xquartz.org>). Without it, `ggsave()` only warns ("failed to load
cairo DLL"), the scripts finish without error, and every plot goes to
`scripts/Rplots.pdf` instead of `results/`. Check before running; this must
print `TRUE`:

```r
f = tempfile(fileext = ".pdf"); grDevices::cairo_pdf(f); plot.new(); invisible(grDevices::dev.off()); file.exists(f)
```

```r
# scripts 01 to 06 and run.R
install.packages(c("data.table", "ggplot2", "ggrepel", "patchwork", "nlme", "knitr"))
# 00_harmonise_cohorts.Rmd only
install.packages(c("readxl", "fst", "BiocManager"))
BiocManager::install(c("NMF", "BSgenome", "GenomeInfoDb",
                       "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38"))
```

All paths are relative to `scripts/`, so run everything from there. `run.R`
executes one Rmd's code chunks without pandoc:

```sh
cd scripts
Rscript run.R 01_four_groups.Rmd            # one script
sh run_all.sh                               # 01 to 06 in order, stops at the first failure
RSCRIPT=/path/to/Rscript sh run_all.sh      # the same with a specific R
```

Run order and dependencies:

- `00_harmonise_cohorts.Rmd` turns `raw/` into `data/`. It needs `raw/` and the
  packages `fst`, `readxl`, `NMF` and the two BSgenome builds. Its outputs are
  already in `data/`, so `run_all.sh` skips it.
- 01, 02, 04 and 05 read only `data/` and can run in any order.
- 03 reads `results/01_four_groups/tables/06_driver_genomewide.csv`: run it after 01.
- 06 reads tables written by 01 and 05: run it last.

Scripts 00 to 05 set the seed `20260809`. Each script overwrites its own
outputs in `results/`.

## Repository layout

```
scripts/                      analysis code (run from here)
  _functions.R, run.R, run_all.sh   shared functions and runners
  00_harmonise_cohorts.Rmd    raw/ -> data/
  01_four_groups.Rmd          Figure 1A-C
  02_within_cohort.Rmd        Figure 1D, left
  03_pooled_vs_stratified.Rmd pooled versus stratified gene counts
  04_cohort_pairs.Rmd         Figure 1D, right
  05_sensitivity.Rmd          sources of Tables S2 and S3
  06_supplementary_tables.Rmd Tables S1-S3
data/                         harmonised inputs; signature/ holds SBS96 spectra and refits
results/
  figure1_panels/             the nine source panels of Figure 1
  supplementary_tables/       Tables S1-S3
  01_four_groups/ ... 05_sensitivity/   tables/ and figures/ for each script
figures/                      Figure1.pdf, Figure1.jpg, graphical_abstract.jpg
raw/                          not included (see "Raw data: raw/ layout")
```

## Analysis decisions

- **Cohort inclusion.** Chen2019 enters the cohort description and TMB only:
  its deposited mutation file is a curated gene list, not a full exome table
  (flag `HasFullMAF`). Guo2026 enters the pooled analyses only, because eight
  ever-smokers cannot support within-cohort or cohort-pair estimates (flag
  `InPerCohortAnalysis`). The other four cohorts enter everything.
- **Protein-altering restriction.** MAF, VEP and ANNOVAR annotations are mapped
  to one variant class. The non-coding to protein-altering ratio ranges from
  0.031 to 1.110 across cohorts (`data/variant_class_audit.csv`), which reflects
  annotation depth, so TMB and mutation frequencies use protein-altering
  variants only. Signatures use all single-base substitutions, because the
  refit needs the complete 96-context spectrum.
- **TMB and signatures.** TMB is protein-altering variants per 38 Mb, published
  for Chen2019 and Chen2020 and recomputed for the others. SBS96 spectra are
  built on each cohort's reference build and refitted to COSMIC v3.6 by
  non-negative least squares; the 89 signatures above 1% in at least 20 tumours
  are tested.
- **Hypermutation** is flagged and retained. The primary flag is based on
  mechanism (POLE or mismatch-repair signatures above 20%, SBS2 plus SBS13 above
  50%), because a TMB cut-off of 10/Mb would remove 29.0% of Western
  ever-smokers but only 0.5% of EAS never-smokers. 05 reruns the analyses
  without flagged tumours under four definitions (Table S2).
- **Ancestry and smoking.** EAS is assigned at cohort level. Ten Western
  patients reported as Asian or part-Asian (eight TCGA2014, two TRACERx2023) and
  four Chen2020 patients with unknown smoking status are excluded, taking the
  1,467 harmonised tumours to 1,453. Restricting the Western group to reported
  White is a sensitivity analysis (Table S3).
- **TRACERx** uses the tumour-level mutation table, which has one row per tumour
  and locus, so multi-region sampling does not inflate counts.
- **Multiple testing and odds ratios.** Adjusted P values are
  Benjamini-Hochberg, and the family is the one used in the table that reports
  them. In `results/01_four_groups/tables/`, each of the five contrasts is its
  own family (8,428 genes in `06_driver_genomewide.csv`, 89 signatures in
  `09_signature_contrasts.csv`, six classes in `11_point_class_contrasts.csv`,
  one test in `03_tmb_contrasts.csv`); `10_signature_interaction.csv` adjusts
  over the 15 signatures with a non-zero median in at least one group, and
  `15_western_never_vs_eas_ever.csv` over all 8,524 tests of the cross-check.
  Figure 1D, left, is one family of 384 estimates (TMB, 89 signatures and six
  substitution classes in each of four cohorts); Figure 1D, right, is one family
  of 58,135 tests (gene × cohort pair × smoking group). Odds ratios carry a 0.5
  continuity correction. The 95% CIs quoted in the Letter are Wald intervals on
  that log odds ratio, exp(log OR ± 1.96 × SE) with
  SE = sqrt(1/(a + 0.5) + 1/(n<sub>a</sub> − a + 0.5) + 1/(b + 0.5) + 1/(n<sub>b</sub> − b + 0.5)),
  where a = `mutated_a`, n<sub>a</sub> = `n_a`, b = `mutated_b` and
  n<sub>b</sub> = `n_b` in `06_driver_genomewide.csv`; no script writes them.
- **Copy number and fusions are not compared.** Two cohorts supply only
  aggregate copy-number tables, the three sample-level quantifications (GISTIC
  discrete calls, ASCAT2 absolute copy number, GISTIC peak calls) share no
  common scale, and the Guo2026 encoding is undocumented
  (`results/05_sensitivity/tables/08_cnv_comparability_audit.csv`). Fusion calls
  differ too much between pipelines to be compared across cohorts.

## Citation

Please cite the Letter and the data records:

- Liu L, Yu R, et al. Separating tobacco exposure from ancestry in East Asian
  and Western lung adenocarcinoma. MedComm - Oncology. In press.
- Open data release (Zenodo v1): https://doi.org/10.5281/zenodo.21201441
- In-house data, restricted access (Zenodo v2): https://doi.org/10.5281/zenodo.23050408

## Contact

- Liang-Liang Fan, Central South University: swfanliangliang@csu.edu.cn
- Shuai Guo, Emory University: shuai.guo@emory.edu
