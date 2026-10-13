
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

source("scripts/99_functions.R")

# -------------------------------------
# Model 1A: Poisson with overdispersion
# -------------------------------------
label_model1a_pois <- paste0("model1a_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1a <- readRDS(here("results", paste0("mcmc_list_", label_model1a_pois, ".rds")))
samples_1a <- readRDS(here("results", paste0("samples_", label_model1a_pois, ".rds")))
summarise_model(mcmc_list1a,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs"),
                label = label_model1a_pois)


# ----------------------------------------
# Model 1B: Poisson without overdispersion
# ----------------------------------------
label_model1b_pois <- paste0("model1b_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1b <- readRDS(here("results", paste0("mcmc_list_", label_model1b_pois, ".rds")))
samples_1b <- readRDS(here("results", paste0("samples_", label_model1b_pois, ".rds")))
summarise_model(mcmc_list1b,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b"),
                label = label_model1b_pois)


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


# --------------------------------------------
# Forest plot for main effects (time + policy) - per 5 years
# --------------------------------------------
# NOTE: Run this after running the the "diagnostics.R".

# Model names
model_names <- tibble(
  Model_short = c("model1A", "model1B", "model2", "model3", "model4", "model5", "model6"),
  model_label = c(
    "Model 1A: Base Poisson model - step change",   # with overdispersion
    "Model 1B: Base Poisson model - step change",   # without overdispersion
    "Model 2: Policy leakage - step change",
    "Model 3: Field specific hierarchy - step change",
    "Model 4: Base Poisson model - gradual change",
    "Model 5: Journal specific trend - gradual change",
    "Model 6: Field specific hierarchy - gradual change"))

# Import DIC 
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

effects <- dic %>%
  select(Model, model_label, selected, DIC, Delta_DIC) %>%
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
  ungroup() %>%
  mutate(
    Model = factor(Model, levels = rev(dic$Model)),  # use Model to keep 1A and 1B separate
    model_label = factor(model_label, levels = rev(dic$model_label)),
    effect = factor(effect, levels = c("Time effect (per 5 years)", "Policy effect")))

# Combine both effects onto one row per model 
n_models <- nlevels(effects$Model)

effects <- effects %>%
  mutate(
    y_base = as.numeric(Model),
    y_pos  = y_base + if_else(effect == "Time effect (per 5 years)", 0.16, -0.16))

y_breaks <- seq_len(n_models)

# Get the corresponding model labels in the correct order
y_labels <- dic %>%
  mutate(Model = factor(Model, levels = rev(dic$Model))) %>%
  arrange(Model) %>%
  pull(model_label)

#x_min <- floor(min(effects$rr_low) * 100) / 100
x_min <- 0.75
x_max <- 1.07
x_breaks <- pretty(c(x_min, x_max), n = 6)

y_limits <- c(0.7, n_models + 0.5)

# Main forest plot 
forest_plot <- ggplot(
  effects,
  aes(x = rate_ratio, y = y_pos, xmin = rr_low, xmax = rr_high, colour = effect)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "red") +
  geom_errorbarh(height = 0.12, linewidth = 0.7) +
  geom_point(size = 2.4, shape = 16) +
  scale_colour_manual(
    values = c("Time effect (per 5 years)" = "#0072B2", "Policy effect" = "#E69F00"), name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = y_limits, clip = "off", expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(1.07, 1.06),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 5, r = 15, b = 5, l = 10))

gap <- (x_max - x_min) * 0.015   # small buffer around RR = 1

direction_plot <- ggplot() +
  annotate("segment", x = 1 - gap, xend = x_min, y = 1, yend = 1,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("segment", x = 1 + gap, xend = x_max, y = 1, yend = 1,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("text", x = (1 + x_min) / 2, y = 0.35, label = "Bar graph use decreasing",
           size = 3.8, colour = "grey30") +
  annotate("text", x = (1 + x_max) / 2, y = 0.35, label = "Bar graph use increasing",
           size = 3.8, colour = "grey30") +
  annotate("text", x = (x_min + x_max) / 2, y = -0.15, label = "Rate ratio",
           size = 4.2, colour = "black", fontface = "bold") +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(-0.5, 1.5), expand = FALSE) +
  theme_void()

# DIC arrow
arrow_plot <- ggplot(tibble(x = 0), aes(x = x)) +
  geom_segment(aes(x = 0, xend = 0, y = 1, yend = n_models),
               arrow = arrow(length = unit(0.12, "inches"), ends = "both", type = "closed"),
               colour = "grey30", linewidth = 0.3) +
  annotate("text", x = 0, y = n_models, label = "Low DIC (best)",
           angle = 90, vjust = -1, size = 3.8, colour = "grey30") +
  annotate("text", x = 0, y = 1, label = "High DIC (worst)",
           angle = 90, vjust = -1, size = 3.8, colour = "grey30") +
  coord_cartesian(xlim = c(-1, 1), ylim = y_limits, clip = "off", expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 5))

# Combine with patchwork
layout <- "
AB
#C
"

final_plot <- arrow_plot + forest_plot + direction_plot +
  plot_layout(design = layout, widths = c(0.06, 1), heights = c(1, 0.15))

final_plot

ggsave(filename = here("results", "forest_plot_all_models_bar.png"), plot = final_plot, width = 10, height = 6, dpi = 300)


# --------------------------------------------
# Forest plot for main effects (time + policy) - per 5 years
# --------------------------------------------
# NOTE: Run this after running the the "diagnostics.R".

# Model names
model_names <- tibble(
  Model_short = c("model1A", "model1B", "model2", "model3", "model4", "model5", "model6"),
  model_label = c(
    "Model 1A: Base Poisson model - step change",   # with overdispersion
    "Model 1B: Base Poisson model - step change",   # without overdispersion
    "Model 2: Policy leakage - step change",
    "Model 3: Field specific hierarchy - step change",
    "Model 4: Base Poisson model - gradual change",
    "Model 5: Journal specific trend - gradual change",
    "Model 6: Field specific hierarchy - gradual change"))

# Import DIC 
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

effects <- dic %>%
  select(Model, model_label, selected, DIC, Delta_DIC) %>%
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
  ungroup() %>%
  mutate(
    Model = factor(Model, levels = rev(dic$Model)),  # use Model to keep 1A and 1B separate
    model_label = factor(model_label, levels = rev(dic$model_label)),
    effect = factor(effect, levels = c("Time effect (per 5 years)", "Policy effect")))

# Combine both effects onto one row per model 
n_models <- nlevels(effects$Model)

effects <- effects %>%
  mutate(
    y_base = as.numeric(Model),
    y_pos  = y_base + if_else(effect == "Time effect (per 5 years)", 0.16, -0.16))

y_breaks <- seq_len(n_models)

# Get the corresponding model labels in the correct order
y_labels <- dic %>%
  mutate(Model = factor(Model, levels = rev(dic$Model))) %>%
  arrange(Model) %>%
  pull(model_label)

#x_min <- floor(min(effects$rr_low) * 100) / 100
x_min <- 0.75
x_max <- 1.07
x_breaks <- pretty(c(x_min, x_max), n = 6)

y_limits <- c(0.7, n_models + 0.5)

# Main forest plot 
forest_plot <- ggplot(
  effects,
  aes(x = rate_ratio, y = y_pos, xmin = rr_low, xmax = rr_high, colour = effect)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "red") +
  geom_errorbarh(height = 0.12, linewidth = 0.7) +
  geom_point(size = 2.4, shape = 16) +
  scale_colour_manual(
    values = c("Time effect (per 5 years)" = "#0072B2", "Policy effect" = "#E69F00"), name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = y_limits, clip = "off", expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(1.07, 1.06),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 5, r = 15, b = 5, l = 10))

gap <- (x_max - x_min) * 0.015   # small buffer around RR = 1

direction_plot <- ggplot() +
  annotate("segment", x = 1 - gap, xend = x_min, y = 1, yend = 1,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("segment", x = 1 + gap, xend = x_max, y = 1, yend = 1,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("text", x = (1 + x_min) / 2, y = 0.35, label = "Bar graph use decreasing",
           size = 3.8, colour = "grey30") +
  annotate("text", x = x_max - (x_max - x_min) * 0.015, y = 0.35, label = "Bar graph use increasing",
           size = 3.8, colour = "grey30", hjust = 1) +
  annotate("text", x = (x_min + x_max) / 2, y = -0.15, label = "Rate ratio",
           size = 4.2, colour = "black", fontface = "bold") +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(-0.5, 1.5), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 5))

# DIC arrow
arrow_plot <- ggplot(tibble(x = 0), aes(x = x)) +
  geom_segment(aes(x = 0, xend = 0, y = 1, yend = n_models),
               arrow = arrow(length = unit(0.12, "inches"), ends = "both", type = "closed"),
               colour = "grey30", linewidth = 0.3) +
  annotate("text", x = 0, y = n_models, label = "Low DIC (best)",
           angle = 90, vjust = -1, size = 3.8, colour = "grey30") +
  annotate("text", x = 0, y = 1, label = "High DIC (worst)",
           angle = 90, vjust = -1, size = 3.8, colour = "grey30") +
  coord_cartesian(xlim = c(-1, 1), ylim = y_limits, clip = "off", expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 5))

# Combine with patchwork
layout <- "
AB
#C
"

final_plot <- arrow_plot + forest_plot + direction_plot +
  plot_layout(design = layout, widths = c(0.06, 1), heights = c(1, 0.15))

final_plot

ggsave(filename = here("results", "forest_plot_all_models_bar.png"), plot = final_plot, width = 10, height = 6, dpi = 300)
