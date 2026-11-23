# blood-methylation-comorbidity

Analysis code, run logs and figure source data for the study

**Benchmarked comparison of blood DNA methylation across five common disorders identifies a recurring leukocyte shift but no reproducible shared methylation signature**

The study compares genome-wide case-control effects of peripheral-blood DNA methylation between Parkinson's disease, schizophrenia, Alzheimer's disease, major depressive disorder and inflammatory bowel disease. Each disorder is represented by two independent public cohorts, so that within-disease replication calibrates every cross-disease comparison.

| Disorder | Cohorts (GEO) |
|---|---|
| Parkinson's disease | GSE145361, GSE111629 |
| Schizophrenia | GSE84727, GSE80417 |
| Alzheimer's disease | GSE338167, GSE144858 |
| Major depressive disorder | GSE198904, GSE273970 |
| Inflammatory bowel disease | GSE87648, GSE112611 |

## Contents

| Path | Content |
|---|---|
| `scripts/` | 72 R scripts and 1 shell script |
| `output/logs/` | Run logs: the last run of every script that writes a log |
| `output/53_figures/` | 17 CSV files with the plotted values of Figures 1-6 and S1-S6 |

Raw data are not included. All data sets are public and are downloaded from GEO by the download scripts (`03_download_data.R`, `03b_download_panel2.R`, `03c_fix_gse198904.R`).

## Environment

R 4.4.1 in RStudio on Windows 11. The versions of the key packages are listed at the end of the logs.

## How the scripts are organised

1. Scripts are numbered in the order in which they were written and run. Files starting with `00_` are shared components (logging, matrix reader, plotting).
2. Scripts 30-57 are the analyses of the ten-cohort, five-disorder panel reported in the manuscript. The header of `55_run_after_dedup.R` lists the final order of the downstream steps (54, 39-47, 50-53) and their dependencies. `53_manuscript_figures.R` produces all figures and their source data.
3. Scripts 56 and 57 are two sensitivity analyses: exclusion of the six sex-discordant samples of GSE87648 (56), and the forty cross-disease comparisons under models that include the sample-quality index (57).
4. Scripts 00-17 cover data acquisition, quality control, covariates and the cohort-level models of the first set of cohorts. Four of these cohorts (GSE145361, GSE111629, GSE87648, GSE112611) are part of the final panel, and their models originate from these steps (re-determined in 37, 38 and 54).
5. Scripts 18 and 20-28 belong to an earlier, exploratory stage restricted to Parkinson's disease and inflammatory bowel disease (with schizophrenia as a negative control); scripts 19-19e and 29 document the search for and selection of additional cohorts. Their results are not reported in the manuscript. They are retained so that the repository is the complete record of the project.
6. The log of `21_composition_probe.R` ends with an error; the cause is described and corrected in `22_compositional.R`.

## Language and paths

The analyses up to step 55 were run with comments, messages and labels in Chinese. For this release:

- Comments and character strings containing Chinese text, and one quoted argument name in `50_overlap_screen.R`, were translated into English. All other code tokens are identical to the versions that produced the results (compared token by token with the R parser). The translated scripts were not re-run. Scripts 56 and 57 were written in English.
- The logs are the logs of the original runs. Their text was translated with the message templates of the translated scripts; no number was altered, and lines without Chinese text are unchanged.
- The text labels in three source data files (Fig2, FigS4, FigS4_sex) were translated in the same way.
- Local directory names were replaced by neutral ones. Before re-running, set `PROJ_ROOT` (`"C:/blood_methylation_project"`) and `DATA_ROOT` (`"D:/blood_methylation_data"`) at the top of each script to existing directories. The three download scripts write to their own `DATA_ROOT` (`"C:/blood_methylation_download"`).

## Citation

Please cite the archived release of this repository (see `CITATION.cff`). The reference to the article will be added on publication.

## License

MIT, see `LICENSE`.
