## ===========
## results.R
## ===========

library(coda)
library(here)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(scales)
library(ggh4x)
library(patchwork)
library(posterior)

source("scripts/99_functions.R")

# -------------------------------------
# Model 1A: Poisson with overdispersion
# -------------------------------------
label_model1A_pois <- paste0("model1A_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1A <- readRDS(here("results", paste0("mcmc_list_", label_model1A_pois, ".rds")))
samples_1A <- readRDS(here("results", paste0("samples_", label_model1A_pois, ".rds")))
summarise_model(mcmc_list1A,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs"),
                label = label_model1A_pois)


# ----------------------------------------
# Model 1B: Poisson without overdispersion
# ----------------------------------------
label_model1B_pois <- paste0("model1B_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1B <- readRDS(here("results", paste0("mcmc_list_", label_model1B_pois, ".rds")))
samples_1B <- readRDS(here("results", paste0("samples_", label_model1B_pois, ".rds")))
summarise_model(mcmc_list1B,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b"),
                label = label_model1B_pois)


# -------------------------------------------------
# Model 2: Fixed time trend, random intercepts only
# -------------------------------------------------
label_model2_pois <- paste0("model2_pois_", mcmc_niter %/% 1000, "k")
mcmc_list2 <- readRDS(here("results", paste0("mcmc_list_", label_model2_pois, ".rds")))
samples_2 <- readRDS(here("results", paste0("samples_", label_model2_pois, ".rds")))
summarise_model(mcmc_list2,
                key_params = c("alpha", "beta", "gamma", "sigma_a"),
                label = label_model2_pois)


# ---------------------------------------------------
# Model 3: Nested intercepts (journals within fields)
# ---------------------------------------------------
label_model3_pois <- paste0("model3_pois_", mcmc_niter %/% 1000, "k")
mcmc_list3 <- readRDS(here("results", paste0("mcmc_list_", label_model3_pois, ".rds")))
samples_3 <- readRDS(here("results", paste0("samples_", label_model3_pois, ".rds")))
summarise_model(mcmc_list3,
                key_params = c("A", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f"),
                label = label_model3_pois)

# -------------------------------------------
# Model 4: Gradual policy effect, fixed slope
# -------------------------------------------
label_model4_pois <- paste0("model4_pois_", mcmc_niter %/% 1000, "k")
mcmc_list4 <- readRDS(here("results", paste0("mcmc_list_", label_model4_pois, ".rds")))
samples_4 <- readRDS(here("results", paste0("samples_", label_model4_pois, ".rds")))
summarise_model(mcmc_list4,
                key_params = c("alpha", "beta", "gamma", "sigma_a"),
                label = label_model4_pois)


# --------------------------------------------------------
# Model 5: Gradual policy effect + journal-specific slopes
# --------------------------------------------------------
label_model5_pois <- paste0("model5_pois_", mcmc_niter %/% 1000, "k")
mcmc_list5 <- readRDS(here("results", paste0("mcmc_list_", label_model5_pois, ".rds")))
samples_5 <- readRDS(here("results", paste0("samples_", label_model5_pois, ".rds")))
summarise_model(mcmc_list5,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b"),
                label = label_model5_pois)


# ------------------------------------------------
# Model 6: Gradual policy effect + field hierarchy
# ------------------------------------------------
label_model6_pois <- paste0("model6_pois_", mcmc_niter %/% 1000, "k")
mcmc_list6 <- readRDS(here("results", paste0("mcmc_list_", label_model6_pois, ".rds")))
samples_6 <- readRDS(here("results", paste0("samples_", label_model6_pois, ".rds")))
summarise_model(mcmc_list6,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f"),
                label = label_model6_pois)



# ----------------------------------------------------------
# Forest plot for main effects (time + policy) - per 5 years
# ----------------------------------------------------------
# NOTE: Run this after running the "diagnostics.R".

# ---- Model names -------------------------------------------------------------
model_names <- tibble(
  Model_short = c("model1A", "model1B", "model2", "model3", "model4", "model5", "model6"),
  model_label = c(
    "Model 1A: Base Poisson model",   # with overdispersion
    "Model 1B: Base Poisson model",   # without overdispersion
    "Model 2: Policy leakage",
    "Model 3: Field specific hierarchy",
    "Model 4: Base Poisson model",
    "Model 5: Journal specific trend",
    "Model 6: Field specific hierarchy"),
  change_type = c(rep("Step change", 4), rep("Gradual change", 3)))

# ---- Import DIC --------------------------------------------------------------
dic <- read_csv(
  here("results", "dic_comparison_bar.csv"),
  show_col_types = FALSE) %>%
  arrange(DIC) %>%
  mutate(
    Model_short = str_remove(Model, "_pois_60k$"),
    selected = Delta_DIC == 0) %>%
  left_join(model_names, by = "Model_short") %>%
  mutate(
    model_label = if_else(
      is.na(model_label),
      str_replace(Model_short, "model", "Model "),
      model_label))

parse_cri <- function(x) {
  x <- str_remove_all(x, "[()]")
  parts <- str_split_fixed(x, ",\\s*", 2)
  tibble(low = as.numeric(parts[, 1]), high = as.numeric(parts[, 2]))
}

# ---- Effects table -----------------------------------------------------------
effects <- dic %>%
  select(Model, model_label, selected, DIC, Delta_DIC, change_type) %>%
  rowwise() %>%
  do({
    model_info <- .
    read_csv(here("results", paste0("summary_", model_info$Model, ".csv")), show_col_types = FALSE) %>%
      filter(Parameter %in% c("beta", "gamma")) %>%
      bind_cols(parse_cri(.$`95% CrI`)) %>%
      mutate(
        Model = model_info$Model,
        model_label = model_info$model_label,
        selected = model_info$selected,
        DIC = model_info$DIC,
        Delta_DIC = model_info$Delta_DIC,
        change_type = model_info$change_type,
        effect = case_when(
          Parameter == "beta" ~ "Time effect (per 5 years)",
          Parameter == "gamma" ~ "Policy effect"),
        
        # Scale time effect from per year to per 5 years
        Mean_scaled = if_else(Parameter == "beta", Mean * 5, Mean),
        low_scaled = if_else(Parameter == "beta", low * 5, low),
        high_scaled = if_else(Parameter == "beta", high * 5, high),
        
        rate_ratio = exp(Mean_scaled),
        rr_low = exp(low_scaled),
        rr_high = exp(high_scaled))
  }) %>%
  ungroup()

# ---- Common axis limits ------------------------------------------------------
x_min <- 0.75
x_max <- 1.07
x_breaks <- pretty(c(x_min, x_max), n = 6)
x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
x_labels <- trimws(format(x_breaks, nsmall = 1))

# ---- Vertical layout ---------------------------------------------------------
# Everything is drawn in ONE plot with a single y-scale: the two panels are two
# vertical blocks of that scale, separated by a gap that holds the top panel's
# x-axis labels and the bottom panel's title. The DIC arrows are a second plot
# with the identical y-scale, so they line up exactly with every model row.

top_type    <- "Step change"       # panel (a)
bottom_type <- "Gradual change"    # panel (b)

# Rank models by DIC within each panel (rank 1 = lowest DIC = top of panel)
model_pos <- effects %>%
  distinct(Model, model_label, change_type, DIC) %>%
  group_by(change_type) %>%
  arrange(DIC, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

n_a <- sum(model_pos$change_type == top_type)
n_b <- sum(model_pos$change_type == bottom_type)

panel_gap <- 1.2
b_lo <- 0.7                    # bottom panel: lower limit
b_hi <- n_b + 0.5              #               upper limit
a_lo <- b_hi + panel_gap       # top panel:    lower limit
a0   <- a_lo + 0.3             #               y of its lowest row
a_hi <- a0 + n_a - 1 + 0.5     #               upper limit

y_top <- a_hi + 0.75           # room for title (a) and legend
y_bot <- b_lo - 2.2            # room for x labels, direction arrows, "Rate ratio"
y_limits <- c(y_bot, y_top)

model_pos <- model_pos %>%
  mutate(y_base = if_else(change_type == bottom_type,
                          n_b - rank + 1,
                          a0 + n_a - rank))

plot_data <- effects %>%
  left_join(select(model_pos, Model, y_base), by = "Model") %>%
  mutate(y_pos = y_base + if_else(effect == "Time effect (per 5 years)", 0.16, -0.16))

y_breaks <- model_pos$y_base
y_labels <- model_pos$model_label

# One arrow per panel: from that panel's best model (top row) to its worst (bottom row)
arrow_df <- model_pos %>%
  group_by(change_type) %>%
  summarise(y_best = max(y_base), y_worst = min(y_base), .groups = "drop")

# ---- Styling constants -------------------------------------------------------
grid_col <- "grey92"
grid_lw  <- 14 / 22                    # same as theme_minimal(base_size = 14) major grid
tick_sz  <- 14 * 0.8 / .pt             # same as theme_minimal axis text
dir_gap  <- (x_max - x_min) * 0.015    # small buffer around RR = 1

col_time   <- "#0072B2"
col_policy <- "#E69F00"

# ---- Per-panel legend --------------------------------------------------------
# Right-aligned at x_max. Offsets are in x-axis units; if labels collide with
# keys or with the panel title, adjust the numbers below (or widen the plot).
panel_legend <- function(y) {
  key_w      <- 0.008                  # half-width of the key's error bar
  x_pol_txt  <- x_max                  # right edge of "Policy effect"
  x_pol_key  <- x_pol_txt - 0.060      # key sits left of its label
  x_time_txt <- x_pol_key - 0.020      # right edge of "Time effect ..."
  x_time_key <- x_time_txt - 0.100     # key sits left of its label
  
  list(
    # Time effect key + label
    annotate("segment", x = x_time_key - key_w, xend = x_time_key + key_w,
             y = y, yend = y, colour = col_time, linewidth = 0.7),
    annotate("point", x = x_time_key, y = y, colour = col_time,
             size = 2.4, shape = 16),
    annotate("text", x = x_time_txt, y = y, label = "Time effect (per 5 years)",
             hjust = 1, vjust = 0.5, size = 11 / .pt, colour = "grey20"),
    # Policy effect key + label
    annotate("segment", x = x_pol_key - key_w, xend = x_pol_key + key_w,
             y = y, yend = y, colour = col_policy, linewidth = 0.7),
    annotate("point", x = x_pol_key, y = y, colour = col_policy,
             size = 2.4, shape = 16),
    annotate("text", x = x_pol_txt, y = y, label = "Policy effect",
             hjust = 1, vjust = 0.5, size = 11 / .pt, colour = "grey20")
  )
}

# ---- Forest plot (both panels) -----------------------------------------------
forest <- ggplot(
  plot_data,
  aes(x = rate_ratio, y = y_pos, xmin = rr_low, xmax = rr_high, colour = effect)) +
  
  # vertical grid lines, drawn per panel so they do not cross the gap
  geom_segment(data = tibble(x = x_breaks), inherit.aes = FALSE,
               aes(x = x, xend = x, y = a_lo, yend = a_hi),
               colour = grid_col, linewidth = grid_lw) +
  geom_segment(data = tibble(x = x_breaks), inherit.aes = FALSE,
               aes(x = x, xend = x, y = b_lo, yend = b_hi),
               colour = grid_col, linewidth = grid_lw) +
  
  # null line (rate ratio = 1), one segment per panel
  annotate("segment", x = 1, xend = 1, y = a_lo, yend = a_hi,
           linetype = "dotted", colour = "red") +
  annotate("segment", x = 1, xend = 1, y = b_lo, yend = b_hi,
           linetype = "dotted", colour = "red") +
  
  geom_errorbarh(height = 0.12, linewidth = 0.7) +
  geom_point(size = 2.4, shape = 16) +
  
  # panel titles
  annotate("text", x = x_min, y = a_hi + 0.3, label = "(a) Step-change models",
           hjust = 0, vjust = 0, fontface = "bold", size = 12 / .pt) +
  annotate("text", x = x_min, y = b_hi + 0.3, label = "(b) Gradual change models",
           hjust = 0, vjust = 0, fontface = "bold", size = 12 / .pt) +
  
  # one legend per panel, on the same line as each panel title
  panel_legend(a_hi + 0.3 + 0.17) +
  panel_legend(b_hi + 0.3 + 0.17) +
  
  # x-axis labels under each panel
  annotate("text", x = x_breaks, y = a_lo - 0.15, label = x_labels,
           vjust = 1, colour = "grey30", size = tick_sz) +
  annotate("text", x = x_breaks, y = b_lo - 0.15, label = x_labels,
           vjust = 1, colour = "grey30", size = tick_sz) +
  
  # shared direction arrows + axis title, below the bottom panel
  annotate("segment", x = 1 - dir_gap, xend = x_min, y = b_lo - 1.0, yend = b_lo - 1.0,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("segment", x = 1 + dir_gap, xend = x_max, y = b_lo - 1.0, yend = b_lo - 1.0,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("text", x = (1 + x_min) / 2, y = b_lo - 1.5, label = "Bar graph use decreasing",
           size = 3.8, colour = "grey30") +
  annotate("text", x = x_max - (x_max - x_min) * 0.015, y = b_lo - 1.5,
           label = "Bar graph use increasing",
           size = 3.8, colour = "grey30", hjust = 1) +
  annotate("text", x = (x_min + x_max) / 2, y = b_lo - 2.0, label = "Rate ratio",
           size = 4.2, colour = "black", fontface = "bold") +
  
  scale_colour_manual(
    values = c("Time effect (per 5 years)" = col_time, "Policy effect" = col_policy),
    name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = y_limits, clip = "off", expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",                      # legends drawn manually per panel
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),          # drawn manually per panel
    axis.text.x = element_blank(),                 # drawn manually per panel
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 0, r = 15, b = 0, l = 10))

# ---- One DIC arrow per panel -------------------------------------------------
# Same y-scale and same top/bottom margins as the forest plot, so each arrow
# runs from its panel's best model (top row) to its worst model (bottom row).
arrow_plot <- ggplot(arrow_df) +
  geom_segment(aes(x = 0, xend = 0, y = y_worst, yend = y_best),
               arrow = arrow(length = unit(0.12, "inches"), ends = "both", type = "closed"),
               colour = "grey30", linewidth = 0.3) +
  geom_text(aes(x = 0, y = y_best), label = "Low DIC (best)",
            angle = 90, vjust = -1, size = 3.5, colour = "grey30") +
  geom_text(aes(x = 0, y = y_worst), label = "High DIC (worst)",
            angle = 90, vjust = -1, size = 3.5, colour = "grey30") +
  coord_cartesian(xlim = c(-1, 1), ylim = y_limits, clip = "off", expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(t = 0, r = 2, b = 0, l = 5))

# ---- Combine -----------------------------------------------------------------
final_plot <- arrow_plot + forest + plot_layout(widths = c(0.06, 1))

final_plot

ggsave(filename = here("results", "forest_plot_bar_two_panels_vertical.png"),
       plot = final_plot, width = 10, height = 8, dpi = 300)


# --------------------------------------------------------------
# Model 3: Policy & time effects, rate ratios & absolute changes
# --------------------------------------------------------------

# 1. Inputs
df_sample <- read.csv(here("data", "df_sample.csv"))
df_sample$journal_idx <- as.integer(factor(df_sample$journal))
df_sample$field_idx   <- as.integer(factor(df_sample$field))

stopifnot("log_offset" %in% names(df_sample), is.numeric(df_sample$log_offset))

journal_i  <- df_sample$journal_idx
yc_i       <- df_sample$year_centered
policy_i   <- df_sample$policy
offset_i   <- exp(df_sample$log_offset)      
J          <- max(df_sample$journal_idx)
N          <- nrow(df_sample)
tot_offset <- sum(offset_i)

# 2. Posterior draws as a plain numeric matrix
label_model3_pois <- paste0("model3_pois_", mcmc_niter %/% 1000, "k")
draws3 <- readRDS(here("results", paste0("mcmc_list_", label_model3_pois, ".rds")))
draws_mat <- unclass(posterior::as_draws_matrix(draws3))
attr(draws_mat, "nchains") <- NULL

gamma_draws <- as.numeric(draws_mat[, "gamma"])
beta_draws <- as.numeric(draws_mat[, "beta"])
a_mat <- draws_mat[, paste0("a[", 1:J, "]")]
b_mat <- draws_mat[, paste0("b[", 1:J, "]")]

rr_policy <- exp(gamma_draws)
rr_time_5yr <- exp(5 * beta_draws)    # typical journal (population-level slope)

# 3. Baseline (policy = 0) and +5 year effect, draw by draw
# No 100x needed: sum(mu) / sum(n_bar_or_info / 100) is already in percent
n_draws <- nrow(a_mat)
baseline_pct_draws <- numeric(n_draws)
abs_change_time_pp <- numeric(n_draws)
chunks <- split(seq_len(n_draws), ceiling(seq_len(n_draws) / 500))

for (ix in chunks) {
  a_i <- a_mat[ix, journal_i, drop = FALSE]
  b_i <- b_mat[ix, journal_i, drop = FALSE]
  mu0 <- sweep(exp(a_i + sweep(b_i, 2, yc_i, `*`)), 2, offset_i, `*`)   # policy = 0
  mu5 <- mu0 * exp(5 * b_i)                                             # +5 years
  baseline_pct_draws[ix] <- rowSums(mu0) / tot_offset
  abs_change_time_pp[ix] <- rowSums(mu5 - mu0) / tot_offset
}

# 4. Policy absolute change (exact: no interaction in the model)
abs_change_policy_pp <- baseline_pct_draws * (rr_policy - 1)

# 5. Summary table
summ <- function(x) {
  q <- quantile(x, c(0.025, 0.975))
  data.frame(Mean = mean(x), Lower = q[[1]], Upper = q[[2]], row.names = NULL)
}

results_table <- rbind(
  cbind(Effect = "Policy",             Metric = "Rate ratio",           summ(rr_policy)),
  cbind(Effect = "Policy",             Metric = "Absolute change (pp)", summ(abs_change_policy_pp)),
  cbind(Effect = "Time (per 5 years)", Metric = "Rate ratio",           summ(rr_time_5yr)),
  cbind(Effect = "Time (per 5 years)", Metric = "Absolute change (pp)", summ(abs_change_time_pp)))

results_table
write.csv(results_table, here("results", "model3_policy_time_impact_summary.csv"), row.names = FALSE)
