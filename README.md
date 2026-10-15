
# Do journal policies reduce the use of bar graphs for visualizing continuous data? An observational study

Robert Schulz<sup>1</sup>, Teshome K. Deressa<sup>2,3</sup>, Adrian Barnett<sup>4</sup>, Yazid Zalia<sup>1</sup>, Teresa Cunha-Oliveira<sup>2,3</sup>, Camila Victoria-Quilla Baselly Heibrich<sup>1</sup>, Parya Abbasi<sup>1,5</sup>, Vladislav Nachev<sup>1</sup>, Tracey L. Weissgerber<sup>2,3</sup>

<sup>1</sup> *QUEST Center for Responsible Research, Berlin Institute of Health at Charité – Universitätsmedizin Berlin, Berlin, Germany*  
<sup>2</sup> *Center for Innovative Biomedicine and Biotechnology (CIBB), University of Coimbra, Coimbra, Portugal*  
<sup>3</sup> *Center for Neuroscience and Cell Biology (CNC-UC), University of Coimbra, Coimbra, Portugal*  
<sup>4</sup> *Australian Centre for Health Services Innovation and Centre for Healthcare Transformation, School of Public Health and Social Work, Faculty of Health, Queensland University of Technology, Brisbane, QLD, Australia*  
<sup>5</sup> *SYNLAB MVZ Berlin GmbH, Berlin, Germany*

---

An observational study of the effect of journal editorial policies on the inappropriate use of bar graphs for displaying continuous data, using Bayesian hierarchical Poisson models fitted with the **nimble** package.

### Quick start

- Clone the repository and create the `data/`, `scripts/`, `results/`, and `diagnostics/` folders locally.
- Open the project as an RStudio Project (required for the **here** package to resolve paths correctly and ensure that the working directory is the repository root).
- Install the required R packages (see [0. Installing packages](#0-installing-packages) below). Developed and tested under R 4.5.1.
- Run scripts 2–6 in `scripts/` in order (see [Pipeline overview](#pipeline-overview) below). Script 6 contains the complete secondary-outcome workflow as well as the cumulative rate-ratio plot, so there is no need to loop back through scripts 4 and 5.
- The pipeline starts from the prepared datasets, so `data/` must already contain `df_descriptive.csv` (read by script 2) and `df_sample.csv` (read by scripts 3 and 6 and by `99_functions.R`).
- Scripts 4 and 5 build the model labels from `mcmc_niter`, which is defined in `3__models.R` (and `6__sensitivity.R`) but not in `99_functions.R`. Run them in the same R session as script 3.

### Repository structure

```text
.
├── data/          # Prepared datasets (df_sample.csv, df_descriptive.csv)
├── scripts/       # All analysis R scripts (numbered pipeline order)
├── results/       # Model outputs, summary tables, figures
└── diagnostics/   # MCMC diagnostics (trace/density/ACF PDFs, DIC summaries)
```

### Pipeline overview

Scripts are numbered and must be run in order, since each step reads files written by the previous one. All scripts use the **here** package for path handling, so the project should be opened as an RStudio Project (or `here::here()` should resolve to the repository root).

| # | Script | Purpose |
|---|---|---|
| 2 | `2__descriptives.R` | Produces descriptive summary tables and exploratory figures. |
| 3 | `3__models.R` | Fits the seven candidate Bayesian Poisson models (primary outcome: number of bar graphs). |
| 4 | `4__results.R` | Loads fitted models and builds parameter summary tables (posterior summaries with ESS and R-hat) and the forest plot for the primary outcome. |
| 5 | `5__diagnostics.R` | Produces MCMC diagnostics (trace, density, ACF plots) and DIC for each primary-outcome model, then collates the DIC comparison table. |
| 6 | `6__sensitivity.R` | Sensitivity analyses: refits models 1A–6 for the secondary outcome (number of more informative graphs) and runs their summaries, diagnostics, DIC table and forest plot; also builds the cumulative rate-ratio plot. |
| — | `99_functions.R` | Shared helper functions (model runner, summaries, diagnostics, DIC). Sourced by scripts 3–6. |

### 0. Installing packages

```r
install.packages(c(
  "nimble", "coda", "posterior", "dplyr", "tidyr", "tidyverse",
  "ggplot2", "readr", "purrr", "stringr", "scales", "here",
  "openxlsx", "ggh4x", "patchwork", "fs", "viridis"
))
```

### 1. Input data

- `df_descriptive.csv`
- `df_sample.csv`

### 2. `2__descriptives.R`

Reads `data/df_descriptive.csv`, recodes field names for publication, and produces descriptive summaries and plots, including journal-sample counts, policy adoption, and follow-up summaries.

#### Outputs written to `results/`

- `summary_table.csv` — Overall and by-policy-group counts/percentages of bar-only, informative-only, both, and excluded articles.
- `policy_timeline_plot.png` — Timeline of journal policy adoption patterns.
- `stacked_by_field.png` — Horizontal stacked bar chart of graph-type usage by field.
- `policy_by_year.png` / `nonpolicy_by_year.png` — Temporal trend plots of graph-type proportions by field, faceted for policy vs. non-policy journals, with policy-adoption year markers.

### 3. `3__models.R`

Fits seven Bayesian hierarchical Poisson models to the primary outcome, reading `data/df_sample.csv`. Uses **nimble** with 4 parallel MCMC chains (60,000 iterations, 10,000 burn-in, thinning of 10).

| Model | Description |
|---|---|
| 1A | Base Poisson, step-change policy effect, journal random intercepts and slopes, with observation-level overdispersion |
| 1B | Same as 1A but without overdispersion |
| 2 | Step-change policy effect, fixed time trend, journal random intercepts only |
| 3 | Step-change policy effect, journals nested within fields, random slopes |
| 4 | Gradual policy effect (`year_since_policy`), fixed time trend, random intercepts |
| 5 | Gradual policy effect + journal-specific time trends |
| 6 | Gradual policy effect + field hierarchy + journal-specific time trends |

#### Outputs written to `results/`

- `samples_<label>.rds` — Per-chain posterior samples and WAIC for each model.
- `mcmc_list_<label>.rds` — The same samples as a `coda` `mcmc.list`.

Labels follow the pattern `model1a_pois_60k`, `model1b_pois_60k`, `model2_pois_60k`, …, `model6_pois_60k` (the number is `mcmc_niter` in thousands).

Models are fitted through `build_and_run_pois()` in `99_functions.R`: one chain per core (up to `detectCores() - 2`; any remaining chains queue), with fixed seeds 20260216–20260219.

Run times recorded in the script comments are roughly 10–19 hours for Models 1B, 2, 3 and 6, and under half an hour for Models 1A, 4 and 5.

### 4. `4__results.R`

Loads the fitted model objects from `results/` and produces final summary tables and comparison figures for the primary outcome. Requires `mcmc_niter` to be defined (see [Quick start](#quick-start)).

- Calls `summarise_model()` for each of the 7 models (posterior mean, SD, 95% credible interval, ESS and \(\hat{R}\) for key parameters) → `results/summary_<model_label>.csv`.
- Builds a forest plot of time and policy rate ratios across all models, ordered by DIC, for the primary outcome → `results/forest_plot_all_models_bar.png` (requires `results/dic_comparison_bar.csv` from script 5).
- The cumulative rate-ratio plot and the secondary-outcome forest plot are not produced here; both are built in `6__sensitivity.R` (see [section 6](#6-6__sensitivityr)).

### 5. `5__diagnostics.R`

For each of the 7 primary-outcome models, loads the MCMC/sample objects from `results/` (`mcmc_niter` must be defined — see [Quick start](#quick-start)) and:

- Saves trace and posterior density plots for key parameters (`save_trace_density()`) → `diagnostics/trace_plots_<label>.pdf` and `diagnostics/density_plots_<label>.pdf`.
- Saves autocorrelation (ACF) plots (`save_acf_plot()`) → `diagnostics/acf_<label>.pdf`.
- Computes DIC (`compute_dic()`; deviance information criterion, with model-specific log-likelihood functions for the non-nested, 3-level nested, and field-hierarchy model structures) → `results/dic_<label>.rds` and `diagnostics/dic_<label>.txt`.
- Collates all saved DIC results with `build_dic_table()`, ranked by DIC with ΔDIC → `results/dic_comparison_bar.csv` (the input to the forest plot in script 4).

### 6. `6__sensitivity.R`

#### Secondary-outcome models

Refits the same 7 model structures as `3__models.R`, but with `n_info` (more informative graph counts) as the outcome instead of `n_bar`. Labels are suffixed `_info_pois_60k` (e.g., `model1a_info_pois_60k`) and outputs follow the same pattern as script 3 (`results/samples_<label>.rds`, `results/mcmc_list_<label>.rds`).

The script then runs the rest of the secondary-outcome workflow itself, so there is no need to loop back through scripts 4 and 5.

#### Results, diagnostics and forest plot (secondary outcome)

- **Results:** `summarise_model()` for each model → `results/summary_<label>.csv`.
- **Diagnostics:** Trace, density and ACF plots and DIC for each model (same functions as script 5) → `diagnostics/`, then `build_dic_table()` → `results/dic_comparison_info.csv`.
- **Forest plot:** Forest plot of time and policy rate ratios, ordered by DIC → `results/forest_plot_all_models_info.png`.

Before running the diagnostics section, edit `99_functions.R`: uncomment the secondary-outcome lines (`base_data` with `y = n_info`, and the `dic_comparison_info.csv` `write_csv`/`message` lines in `build_dic_table()`) and comment out their primary-outcome counterparts. Revert these before re-running scripts 3–5.

Note that `build_dic_table()` collects every `dic_*.rds` file in `results/`, so the info table will also list the primary-outcome models; the forest-plot code keeps only the `_info_` models.

#### Cumulative rate-ratio plot

Builds cumulative rate-ratio curves for the 7 primary-outcome models over −3 to +3 years relative to policy adoption, using the joint posterior draws of the time (`beta`) and policy (`gamma`) effects. Step-change models (1A–3) are drawn solid and gradual-change models (4–6) dashed → `results/cumulative_rr_policy_event_time_joint.png`.
