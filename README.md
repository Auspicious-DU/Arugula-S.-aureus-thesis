# Arugula–S. aureus MSc thesis: R analysis code

This repository contains R scripts for the MSc thesis examining variation in *Staphylococcus aureus* growth responses to *Eruca sativa* leaf extracts and their associations with LC-MS features.

**Data are not included.** This is a code-sharing archive, **not a fully validated one-click reproduction pipeline**. The upstream AUC processing also involved GraphPad Prism.

## Scripts

- `scripts/phenotype/`: 2- and 3-replicate plate processing; threshold crossing, merge/filter; `00_DT_setup.R` and `01_PRE_DT_12h_phenotype_harmonisation.R` for harmonizing archived threshold summaries to the 12-h window.
- `scripts/metabolite_association/`: sample matching, LC-MS data preprocessing, phenotype QC, Pearson and DAM, batch-adjusted sensitivity, PLS/VIP, candidate prioritization.
- `scripts/figures/`: shared plot theme and thesis Figure 3.2, 3.3 and 3.4 scripts.

## Working directory and expected inputs

Open the repository folder in RStudio and set it as the working directory (e.g., RStudio menu **Session → Set Working Directory → Choose Directory**). Scripts require `README.md` in that directory. The scripts use relative paths via `file.path(repo_dir, ...)`.

Create **local-only** `input/` folders for original data. See `PATH_AUDIT.md` for the expected subfolders and filenames. Output files are written to local-only `output/`. Both are excluded from version control by `.gitignore`.

## Limitations

The scripts retain original statistical procedures, but paths were reorganized for portability. No full-data rerun has been performed on this packaged version. Some scripts are intended to be executed sequentially and depend on earlier objects or saved intermediate files. Input data are not publicly available and are retained by the research group.
