# BScThesis

Bachelor thesis project for detecting and explaining outliers using Shapley-Mahalanobis values.

The code is organised as a small research project rather than one long script:

- `analysis/main.Rmd`: main notebook for thesis results, figures, and interpretation.
- `R/original_method.R`: wrappers for the Mayrhofer and Filzmoser Shapley-Mahalanobis method via `ShapleyOutlier`.
- `R/extensions.R`: extension methods such as robust covariance choices, grouped contributions, interactions, and bootstrap stability.
- `R/simulation.R`: synthetic data generation with rowwise and cellwise outliers.
- `R/evaluation.R`: rowwise and cellwise performance metrics.
- `R/plotting.R`: reusable plotting functions.
- `replication_code.r`: command-line runner for a quick end-to-end smoke experiment.

## Setup

Install the project dependencies:

```sh
Rscript scripts/install_packages.R
```

Run the command-line replication script:

```sh
Rscript replication_code.r
```

Render the notebook:

```sh
Rscript -e 'rmarkdown::render("analysis/main.Rmd")'
```

Generated tables and figures are written to `outputs/`. Raw data should go in `data/raw/`, and cleaned/intermediate data can go in `data/processed/`.

## Research Workflow

Use `analysis/main.Rmd` as the main thesis-facing file. Keep methodological code in `R/` files so each component can be tested, reused, and swapped independently. A typical experiment should:

1. Create or load data.
2. Run one or more methods from `R/original_method.R` or `R/extensions.R`.
3. Evaluate results with `R/evaluation.R`.
4. Produce plots with `R/plotting.R`.
5. Save reproducible outputs under `outputs/`.

This keeps the notebook readable while preserving reproducibility.
