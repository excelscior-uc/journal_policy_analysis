# Do journal policies reduce the use of bar graphs for visualizing continuous data? An observational study

Robert Schulz<sup>1</sup>, Teshome K. Deressa<sup>2,3</sup>, Adrian Barnett<sup>4</sup>, Yazid Zalia<sup>1</sup>, Teresa Cunha-Oliveira<sup>2,3</sup>, Camila Victoria-Quilla Baselly Heibrich<sup>1</sup>, Parya Abbasi<sup>1,5</sup>, Vladislav Nachev<sup>1</sup>, Tracey L. Weissgerber<sup>2,3</sup>

<sup>1</sup> *QUEST Center for Responsible Research, Berlin Institute of Health at Charité – Universitätsmedizin Berlin, Berlin, Germany*
<sup>2</sup> *CIBB, Center for Innovative Biomedicine and Biotechnology, University of Coimbra, Coimbra, Portugal*
<sup>3</sup> *CNC-UC, Center for Neuroscience and Cell Biology, University of Coimbra, Coimbra, Portugal*
<sup>4</sup> *Australian Centre for Health Services Innovation and Centre for Healthcare Transformation, School of Public Health and Social Work, Faculty of Health, Queensland University of Technology, Brisbane, QLD, Australia*
<sup>5</sup> *SYNLAB MVZ Berlin GmbH, Berlin, Germany*

---

Code for an observational study of the effect of journal editorial policies on the inappropriate use of bar graphs for displaying continuous data. The analysis uses Bayesian hierarchical Poisson models fitted with the **nimble** package.

- **Primary outcome:** number of articles with bar graphs for continuous data (`n_bar`)
- **Secondary outcome (sensitivity analysis):** number of articles with more informative graphs (`n_info`)

## Contents

- [Quick start](#quick-start)
- [Repository structure](#repository-structure)
- [Pipeline overview](#pipeline-overview)
- [Input data](#input-data)
- [Script reference](#script-reference)
- [Models](#models)
- [Outputs](#outputs)
- [Known quirks](#known-quirks)

## Quick start

1. Clone the repository and create the `data/`, `scripts/`, `results/` and `diagnostics/` folders locally.
2. Put the prepared datasets in `data/`: `df_descriptive.csv` (read by script 2) and `df_sample.csv` (read by scripts 3 and 6 and by `99_functions.R`).
3. Open the project as an **RStudio Project**. The **here** package needs this to resolve paths from the repository root.
4. Install the required packages (developed and tested under **R 4.5.1**):

```r
   install.packages(c(
     "nimble", "coda", "posterior", "dplyr", "tidyr", "tidyverse", "ggplot2",
     "readr", "purrr", "stringr", "scales", "here", "openxlsx", "ggh4x",
     "patchwork", "fs", "viridis"
   ))
```

5. Run the scripts in `scripts/` in the order shown below, **in a single R session** (see [Known quirks](#known-quirks)).

> **Heads up:** Models 1B, 2, 3 and 6 take roughly 10–19 hours each, so `3__models.R` and `6__sensitivity.R` are long runs. Fitted objects are saved to `results/`, so later scripts load them instead of refitting.

## Repository structure

```text
.
├── data/          # Prepared datasets (df_sample.csv, df_descriptive.csv)
├── scripts/       # All analysis R scripts (numbered pipeline order)
├── results/       # Model outputs, summary tables, figures
└── diagnostics/   # MCMC diagnostics (trace/density/ACF PDFs, DIC summaries)
```

## Pipeline overview

| Script | Stage | Purpose |
|:---|:---|:---|
| `99_functions.R` | Helpers | Shared helper functions; sourced by the other scripts. |
| `2__descriptives.R` | Explore | Descriptive summary tables and exploratory figures. |
| `3__models.R` | Model (primary) | Fits the seven candidate Bayesian Poisson models for `n_bar`. |
| `4__results.R` | Report (primary) | Posterior summary tables (with ESS and R-hat) and the forest plot. |
| `5__diagnostics.R` | Check (primary) | MCMC diagnostics (trace, density, ACF), DIC per model, and the DIC comparison table. |
| `6__sensitivity.R` | Sensitivity | Refits all seven models for `n_info`, runs their summaries, diagnostics, DIC table and forest plot, and builds the cumulative rate-ratio plot. |

**Run order:** 2 → 3 → 5 → 4 → 6. The forest plot in script 4 reads `results/dic_comparison_bar.csv`, which script 5 writes, so run script 5 first (or rerun script 4 afterwards). Script 6 is self-contained for the secondary outcome, so there is no need to loop back through scripts 4 and 5.

## Input data

| File | Description |
|:---|:---|
| `df_descriptive.csv` | Journal-year-level data for descriptive summaries, including counts of articles excluded from the main analysis. |
| `df_sample.csv` | Journal-year-level data restricted to articles containing bar graphs for continuous data, informative visualizations, or both. Used for the primary and secondary analyses. |

### Data dictionary: `df_sample.csv`

3,211 observations and 17 variables. Each observation is a journal-year record within a biomedical field and policy-reference classification.

| Variable | Type | Definition |
|:---|:---|:---|
| `field` | Character | Biomedical subject field of the journal. |
| `journal` | Character | Standardized journal name. |
| `year` | Integer | Publication year (2010–2025). |
| `policy_ref` | Binary (0/1) | Whether the journal is classified as a policy-reference journal (`1`) or not (`0`). |
| `bar_specificity` | Numeric categorical | Policy-strength variable: specificity of the journal's bar-graph guidance. |
| `recommendation` | Numeric categorical | Policy-strength variable: coded recommendation category. |
| `sample_size` | Numeric categorical | Policy-strength variable: coded sample-size category. |
| `n_bar` | Integer count | Articles with at least one bar graph for continuous data (each article counted once). |
| `n_info` | Integer count | Articles with at least one informative visualization: bar-dot, box, dot, histogram or violin plot. |
| `n_bar_or_info` | Integer count | Unique articles with bar graphs, informative visualizations, or both. Denominator for the primary analysis. |
| `n_both` | Integer count | Articles with both bar graphs and informative visualizations. |
| `n_bar_xor_info` | Integer count | Articles with exactly one category (bar graphs or informative visualizations, not both). |
| `policy` | Binary (0/1) | Whether the policy had been implemented by the publication year (`1` = adoption year or later; `0` = before adoption or no recorded adoption year). |
| `year_policy_started` | Integer year | Year the policy was adopted; missing if none recorded. |
| `year_since_policy` | Integer | Time since adoption: `0` before adoption or if no adoption year; `1` in the adoption year; `2` the following year, and so on. |
| `year_centered` | Numeric | `year - mean(year)` across the aggregated dataset. |
| `log_offset` | Numeric | `log(n_bar_or_info / 100)`. Exposure offset in applicable count models. |

**Notes**

- `n_bar` and `n_info` can include the same article, so their sum may exceed `n_bar_or_info`.
- The numeric codes of `bar_specificity`, `recommendation` and `sample_size` should be interpreted with the original policy-strength codebook.
- `year_since_policy = 0` covers both pre-adoption journals and journals with no recorded adoption year. Use `year_policy_started` to tell them apart.

## Script reference

### `2__descriptives.R`

Reads `data/df_descriptive.csv`, recodes field names for publication, and produces descriptive summaries and plots (journal-sample counts, policy adoption, follow-up).

**Writes to `results/`:**

- `summary_table.csv`: overall and by-policy-group counts and percentages of bar-only, informative-only, both, and excluded articles.
- `policy_timeline_plot.png`: timeline of journal policy adoption patterns.
- `stacked_by_field.png`: horizontal stacked bar chart of graph-type usage by field.
- `policy_by_year.png` / `nonpolicy_by_year.png`: temporal trends in graph-type proportions by field, for policy vs. non-policy journals, with adoption-year markers.

### `3__models.R`

Fits the seven Bayesian hierarchical Poisson models (see [Models](#models)) to the primary outcome, reading `data/df_sample.csv`.

- **Sampler:** nimble, 4 parallel chains, 60,000 iterations, 10,000 burn-in, thinning of 10.
- **Implementation:** `build_and_run_pois()` in `99_functions.R`. One chain per core (up to `detectCores() - 2`; remaining chains queue).
- **Seeds:** fixed, 20260216–20260219.
- **Run times:** roughly 10–19 hours for Models 1B, 2, 3 and 6; under half an hour for Models 1A, 4 and 5.

**Writes to `results/`:**

- `samples_<label>.rds`: per-chain posterior samples and WAIC.
- `mcmc_list_<label>.rds`: the same samples as a `coda` `mcmc.list`.

Labels follow `model1a_pois_60k`, `model1b_pois_60k`, `model2_pois_60k`, …, `model6_pois_60k` (the number is `mcmc_niter` in thousands).

### `4__results.R`

Loads the fitted models from `results/` and produces the summary tables and forest plot for the primary outcome. Requires `mcmc_niter`.

- Calls `summarise_model()` for each of the 7 models (posterior mean, SD, 95% credible interval, ESS and R-hat for key parameters) → `results/summary_<model_label>.csv`.
- Draws a forest plot of time and policy rate ratios across all models, ordered by DIC → `results/forest_plot_bar_two_panels_vertical.png`. **Requires `results/dic_comparison_bar.csv` from script 5.**
- The cumulative rate-ratio plot and the secondary-outcome forest plot are built in `6__sensitivity.R`, not here.

### `5__diagnostics.R`

For each of the 7 primary-outcome models, loads the MCMC and sample objects from `results/` (requires `mcmc_niter`) and:

- saves trace and posterior density plots for key parameters (`save_trace_density()`) → `diagnostics/trace_plots_<label>.pdf`, `diagnostics/density_plots_<label>.pdf`;
- saves autocorrelation plots (`save_acf_plot()`) → `diagnostics/acf_<label>.pdf`;
- computes DIC (`compute_dic()`, with model-specific log-likelihood functions for the non-nested, 3-level nested and field-hierarchy structures) → `results/dic_<label>.rds`, `diagnostics/dic_<label>.txt`;
- collates all DIC results with `build_dic_table()`, ranked by DIC with ΔDIC → `results/dic_comparison_bar.csv`.

### `6__sensitivity.R`

**Secondary-outcome models.** Refits the same 7 model structures as script 3 with `n_info` as the outcome. Labels are suffixed `_info_pois_60k` (e.g. `model1a_info_pois_60k`); outputs follow the same pattern as script 3.

**Results, diagnostics and forest plot (secondary outcome).**

- **Results:** `summarise_model()` for each model → `results/summary_<label>.csv`.
- **Diagnostics:** trace, density and ACF plots and DIC per model (same functions as script 5) → `diagnostics/`, then `build_dic_table()` → `results/dic_comparison_info.csv`.
- **Forest plot:** time and policy rate ratios, ordered by DIC → `results/forest_plot_info_two_panels_vertical.png`.

**Cumulative rate-ratio plot.** Builds cumulative rate-ratio curves for the 7 *primary-outcome* models over −3 to +3 years relative to policy adoption, using the joint posterior draws of the time (`beta`) and policy (`gamma`) effects. Step-change models (1A–3) are solid and gradual-change models (4–6) are dashed → `results/cumulative_rr_policy_event_time_joint.png`. This needs the primary models from script 3.

## Models

| Model | Policy effect | Description |
|:---|:---|:---|
| 1A | Step change | Base Poisson, journal random intercepts and slopes, observation-level overdispersion |
| 1B | Step change | Same as 1A without overdispersion |
| 2 | Step change | Fixed time trend, journal random intercepts only |
| 3 | Step change | Journals nested within fields, random slopes |
| 4 | Gradual (`year_since_policy`) | Fixed time trend, random intercepts |
| 5 | Gradual | Journal-specific time trends |
| 6 | Gradual | Field hierarchy + journal-specific time trends |

## Outputs

| Folder | Contents |
|:---|:---|
| `results/` | Fitted samples (`.rds`), `coda` objects, summary and DIC tables (`.csv`), figures (`.png`) |
| `diagnostics/` | Trace, density and ACF plots (`.pdf`), DIC text summaries (`.txt`) |

## Known quirks

- **`mcmc_niter` is not in `99_functions.R`.** It is defined in `3__models.R` and `6__sensitivity.R`, and scripts 4 and 5 use it to build model labels. Run scripts 4 and 5 in the same R session as script 3, or define `mcmc_niter <- 60000` yourself first.
- **Switching `99_functions.R` for the secondary outcome.** Before the diagnostics section of script 6, edit `99_functions.R`: uncomment the secondary-outcome lines (`base_data` with `y = n_info`, and the `dic_comparison_info.csv` `write_csv`/`message` lines in `build_dic_table()`) and comment out their primary-outcome counterparts. **Revert these before re-running scripts 3–5.**
- **The info DIC table includes primary models.** `build_dic_table()` collects every `dic_*.rds` file in `results/`, so `dic_comparison_info.csv` also lists the primary-outcome models. The forest-plot code keeps only the `_info_` models.
