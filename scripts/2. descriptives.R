

## ================
## descriptive.R
## ================

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(openxlsx)
library(here)

df_descriptive <- read.csv(here("data", "df_descriptive.csv"))
dim(df_descriptive)

# -----------------------------
# Rename fields for publication
# -----------------------------
df_descriptive <- df_descriptive %>%
  mutate(field = recode(
    field,
    "pharmacology & pharmacy"          = "Pharmacology & Pharmacy",
    "physiology scie"                  = "Physiology",
    "endocrinology & metabolism"       = "Endocrinology & Metabolism",
    "genetics & heredity"              = "Genetics & Heredity",
    "oncology"                         = "Oncology",
    "urology & nephrology"             = "Urology & Nephrology",
    "orthopedics"                      = "Orthopedics",
    "cardiac & cardiovascular systems" = "Cardiac & Cardiovascular Systems",
    "neurosciences"                    = "Neurosciences",
    "immunology"                       = "Immunology",
    "clinical neurology"               = "Clinical Neurology",
    "rheumatology"                     = "Rheumatology"))


# -------------------
# Prepare proportions
# -------------------
df_descriptive <- df_descriptive %>%
  mutate(n_total      = n_bar_or_info + n_others,
         prop_bar     = (n_bar - n_both) / n_total,
         prop_info    = (n_info - n_both) / n_total,
         prop_both    = n_both / n_total,
         prop_others  = n_others / n_total,
         policy_group = case_when(
           policy_ref == 1 ~ "Policy",
           policy_ref == 0 ~ "Non-policy"))


# -------------------------
# 1. Journal sample summary
# -------------------------
journal_sample <- df_descriptive %>%
  distinct(journal, field, policy_ref) %>%
  summarise(
    total_journals     = n_distinct(journal),
    policy_journals    = n_distinct(journal[policy_ref == 1]),
    nonpolicy_journals = n_distinct(journal[policy_ref == 0]),
    fields             = n_distinct(field))
journal_sample


# --------------------
# 2. Overall summaries
# --------------------
# Create mutually exclusive categories
df_summary <- df_descriptive %>%
  mutate(
    n_bar_only  = coalesce(n_bar, 0L)  - coalesce(n_both, 0L),
    n_info_only = coalesce(n_info, 0L) - coalesce(n_both, 0L),
    n_both      = coalesce(n_both, 0L),
    n_others    = coalesce(n_others, 0L),
    n_included  = n_bar_only + n_info_only + n_both,
    n_total     = n_included + n_others)

# Create combined summary table
summary_table <- bind_rows(
  
  # By policy group
  df_summary %>%
    group_by(policy_group) %>%
    summarise(
      total_articles = sum(n_total, na.rm = TRUE),
      
      included_n = sum(n_included, na.rm = TRUE),
      included_pct = round(included_n / total_articles * 100, 1),
      
      bar_only_n = sum(n_bar_only, na.rm = TRUE),
      bar_only_pct = round(bar_only_n / included_n * 100, 1),
      
      info_only_n = sum(n_info_only, na.rm = TRUE),
      info_only_pct = round(info_only_n / included_n * 100, 1),
      
      both_n = sum(n_both, na.rm = TRUE),
      both_pct = round(both_n / included_n * 100, 1),
      
      others_n = sum(n_others, na.rm = TRUE),
      others_pct = round(others_n / total_articles * 100, 1),
      
      .groups = "drop") %>%
    rename(Group = policy_group),
  
  # Overall total row
  df_summary %>%
    summarise(
      Group = "Total",
      
      total_articles = sum(n_total, na.rm = TRUE),
      
      included_n = sum(n_included, na.rm = TRUE),
      included_pct = round(included_n / total_articles * 100, 1),
      
      bar_only_n = sum(n_bar_only, na.rm = TRUE),
      bar_only_pct = round(bar_only_n / included_n * 100, 1),
      
      info_only_n = sum(n_info_only, na.rm = TRUE),
      info_only_pct = round(info_only_n / included_n * 100, 1),
      
      both_n = sum(n_both, na.rm = TRUE),
      both_pct = round(both_n / included_n * 100, 1),
      
      others_n = sum(n_others, na.rm = TRUE),
      others_pct = round(others_n / total_articles * 100, 1))) %>%
  mutate(
    
    # Percent contribution to all articles
    total_pct = round(
      total_articles / total_articles[Group == "Total"] * 100, 1),
    
    `Total articles` = ifelse(
      Group == "Total",
      comma(total_articles),
      paste0(
        comma(total_articles),
        " (", total_pct, "%)")),
    
    `Included: bar or informative` = paste0(
      comma(included_n),
      " (", included_pct, "%)"),
    
    `Bar only` = paste0(
      comma(bar_only_n),
      " (", bar_only_pct, "%)"),
    
    `Informative only` = paste0(
      comma(info_only_n),
      " (", info_only_pct, "%)"),
    
    `Both bar & informative` = paste0(
      comma(both_n),
      " (", both_pct, "%)"),
    
    `Excluded: other graph types` = paste0(
      comma(others_n),
      " (", others_pct, "%)")) %>%
  select(Group, `Total articles`, `Included: bar or informative`, `Bar only`,
         `Informative only`, `Both bar & informative`, `Excluded: other graph types`)

summary_table
write.csv(summary_table, here("results", "summary_table.csv"), row.names = FALSE)

# ------------------------------------------------------
# 3. Journal-level policy adoption & follow-up summary
# ------------------------------------------------------
policy_adoption <- df_descriptive %>%
  filter(policy_ref == 1) %>%
  group_by(journal, field) %>%
  summarise(
    adoption_year = if (any(policy == 1, na.rm = TRUE)) {
      max(year_policy_started[policy == 1], na.rm = TRUE)
    } else {
      NA_real_
    },
    year_policy_started = if (any(!is.na(year_policy_started))) {
      max(year_policy_started, na.rm = TRUE)
    } else {
      NA_real_
    },
    first_observed_year = min(year, na.rm = TRUE),
    last_observed_year  = max(year, na.rm = TRUE),
    pre_policy_years    = n_distinct(year[policy == 0 & !is.na(policy)]),
    post_policy_years   = n_distinct(year[policy == 1 & !is.na(policy)]),
    max_followup_years  = if (any(!is.na(year_since_policy))) {
      max(year_since_policy, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop")

# Descriptive summary
policy_adoption_summary <- policy_adoption %>%
  filter(!is.na(adoption_year)) %>%
  summarise(
    policy_journals              = n_distinct(journal),
    earliest_adoption_year       = min(adoption_year, na.rm = TRUE),
    latest_adoption_year         = max(adoption_year, na.rm = TRUE),
    median_adoption_year         = median(adoption_year, na.rm = TRUE),
    median_post_policy_followup  = median(post_policy_years, na.rm = TRUE),
    min_post_policy_followup     = min(post_policy_years, na.rm = TRUE),
    max_post_policy_followup     = max(post_policy_years, na.rm = TRUE))

policy_adoption_summary

# Number of journals adopting policies by year
adoption_by_year <- policy_adoption %>%
  filter(!is.na(adoption_year)) %>%
  count(adoption_year, name = "n_journals") %>%
  arrange(adoption_year)

adoption_by_year

# Number of post-policy years observed per journal
followup_distribution <- policy_adoption %>%
  count(post_policy_years, name = "n_journals") %>%
  arrange(post_policy_years)

followup_distribution

# Journals with no pre-policy data
no_prepolicy <- df_descriptive %>%
  filter(policy_ref == 1) %>%
  group_by(journal, field) %>%
  summarise(
    first_observed_year = min(year, na.rm = TRUE),
    last_observed_year  = max(year, na.rm = TRUE),
    pre_policy_years    = sum(policy == 0, na.rm = TRUE),
    post_policy_years   = sum(policy == 1, na.rm = TRUE),
    .groups = "drop") %>%
  filter(pre_policy_years == 0 & post_policy_years > 0)

no_prepolicy

# -------------------------------------------
# Policy timeline grouped by adoption pattern
# -------------------------------------------
policy_timeline <- df_descriptive %>%
  filter(policy_ref == 1) %>%
  group_by(journal) %>%
  summarise(
    first_observed_year = min(year, na.rm = TRUE),
    last_observed_year = max(year, na.rm = TRUE),
    n_articles = sum(n_bar_or_info, na.rm = TRUE),
    
    adoption_year = {
      adoption_candidates <- year[
        policy == 1 &
          !is.na(policy) &
          !is.na(year_since_policy)] -
        year_since_policy[
          policy == 1 &
            !is.na(policy) &
            !is.na(year_since_policy)] + 1
      
      if (length(adoption_candidates) == 0 ||
          all(is.na(adoption_candidates))) {
        NA_real_
      } else {
        min(adoption_candidates, na.rm = TRUE)
      }
    },
    .groups = "drop") %>%
  filter(!is.na(adoption_year), is.finite(adoption_year)) %>%
  mutate(
    n_pre_policy_years = pmax(0, adoption_year - first_observed_year))

policy_patterns <- policy_timeline %>%
  group_by(adoption_year, n_pre_policy_years) %>%
  summarise(
    first_observed_year = min(first_observed_year),
    last_observed_year = max(last_observed_year),
    n_journals = n_distinct(journal),
    n_articles = sum(n_articles),
    .groups = "drop") %>%
  arrange(adoption_year, n_pre_policy_years) %>%
  mutate(
    pattern_id = paste0("pattern_", row_number()),
    pattern_id = factor(pattern_id, levels = rev(pattern_id)))

y_labels <- setNames(
  paste0(
    "n\u2c7c = ", policy_patterns$n_journals,
    ", n\u2090 = ", format(policy_patterns$n_articles, big.mark = ",")),
  policy_patterns$pattern_id)

p_policy_timeline <- ggplot(policy_patterns, aes(y = pattern_id)) +
  
  geom_segment(
    data = subset(policy_patterns, n_pre_policy_years > 0),
    aes(
      x = first_observed_year,
      xend = adoption_year,
      yend = pattern_id,
      color = "Pre-policy",
      linewidth = n_articles),
    lineend = "round") +
  
  geom_segment(
    aes(
      x = adoption_year,
      xend = last_observed_year,
      yend = pattern_id,
      color = "Post-policy",
      linewidth = n_articles),
    lineend = "round") +
  
  geom_point(
    aes(x = adoption_year),
    color = "#D46A4C",  size = 2.5) +
  
  scale_color_manual(
    values = c("Post-policy" = "#1F7A7A", "Pre-policy" = "#B0B0B0")) +
  
  scale_linewidth_continuous(trans = "sqrt", range = c(0.5, 3),  guide = "none") +
  
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(
      min(policy_patterns$first_observed_year, na.rm = TRUE),
      max(policy_patterns$last_observed_year, na.rm = TRUE),
      by = 1)) +
  
  scale_y_discrete(
    labels = y_labels,
    expand = c(0.02, 0.02)) +
  
  coord_cartesian(clip = "off") +
  
  labs(x = NULL, y = "Journals with the same adoption pattern", color = NULL) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.box.spacing = unit(0, "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1, margin = margin(r = 4)),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 0, l = 5))

p_policy_timeline

# Save plot
ggsave(here("results", "policy_timeline_plot.png"), p_policy_timeline, width = 9, height = 6, dpi = 300)


# ------------------------------------------------------
# Prepare data for horizontal stacked bar chart by field
# ------------------------------------------------------
df_stacked_field <- df_descriptive %>%
  mutate(n_bar_only = n_bar - n_both,
         n_info_only = n_info - n_both) %>%
  group_by(field) %>%
  summarise(
    n_bar_only = sum(n_bar_only),
    n_info_only = sum(n_info_only),
    n_both = sum(n_both),
    n_others = sum(n_others),
    .groups = "drop") %>%
  mutate(
    n_total = n_bar_only + n_info_only + n_both + n_others,
    pct_bar_only = n_bar_only / n_total * 100,
    pct_info_only = n_info_only / n_total * 100,
    pct_both = n_both / n_total * 100,
    pct_others = n_others / n_total * 100) %>%
  arrange(desc(pct_bar_only))

df_labels_field <- df_stacked_field %>%
  mutate(
    field = factor(field, levels = field),
    field_id = as.numeric(field)) %>%
  select(field, field_id, n_total)

df_stacked_field_long <- df_stacked_field %>%
  mutate(field = factor(field, levels = field)) %>%
  pivot_longer(cols = starts_with("pct_"), names_to = "category", values_to = "percentage") %>%
  mutate(
    category = factor(
      category,
      levels = c("pct_bar_only", "pct_both", "pct_info_only", "pct_others"),
      labels = c("Bar only", "Both bar & informative", "Informative only", "Others")))

bar_height <- 0.65

df_stacked_field_rect <- df_stacked_field_long %>%
  arrange(
    field, factor(category,
      levels = c("Bar only", "Both bar & informative", "Informative only", "Others"))) %>%
  group_by(field) %>%
  mutate(
    xmin = lag(cumsum(percentage), default = 0),
    xmax = cumsum(percentage),
    field_id = as.numeric(field),
    ymin = field_id - bar_height / 2,
    ymax = field_id + bar_height / 2) %>%
  ungroup()

df_both_checker <- df_stacked_field_rect %>%
  filter(category == "Both bar & informative", xmax > xmin) %>%
  group_by(field) %>%
  group_modify(~ {
    x_min <- .x$xmin[1]
    x_max <- .x$xmax[1]
    y_min <- .x$ymin[1]
    y_max <- .x$ymax[1]
    
    cell_width <- 0.8
    n_rows <- 4
    cell_height <- (y_max - y_min) / n_rows
    edge_margin <- 0.25
    
    if ((x_max - x_min) <= 2 * edge_margin) {
      return(data.frame(
        xmin = numeric(0),
        xmax = numeric(0),
        ymin = numeric(0),
        ymax = numeric(0)))
    }
    
    x_starts <- seq(from = x_min + edge_margin, to = x_max - edge_margin, by = cell_width)
    
    y_starts <- seq(from = y_min, to = y_max - cell_height, length.out = n_rows)
    
    expand.grid(
      x_id = seq_along(x_starts),
      y_id = seq_along(y_starts)) %>%
      mutate(
        xmin = x_starts[x_id],
        xmax = pmin(xmin + cell_width, x_max - edge_margin),
        ymin = y_starts[y_id],
        ymax = ymin + cell_height) %>%
      filter((x_id + y_id) %% 2 == 0) %>%
      select(xmin, xmax, ymin, ymax)
  }) %>%
  ungroup()

legend_y <- max(df_labels_field$field_id) + 0.7

df_legend <- data.frame(
  label = c("Bar only", "Both bar & informative", "Informative only", "Others"),
  x = c(0, 20, 55, 86), y = legend_y,
  fill = c("#D9897B","#D9897B", "#A8C5E8","#DCE0E4"))

df_legend_both_checker <- expand.grid(
  x_id = 1:6, y_id = 1:3) %>%
  mutate(xmin = 20 + 0.25 + (x_id - 1) * 0.6, xmax = xmin + 0.6,
         ymin = legend_y - 0.12 + (y_id - 1) * 0.08, ymax = ymin + 0.08) %>%
  filter((x_id + y_id) %% 2 == 0)

stacked_by_field <- ggplot() +
  geom_rect(
    data = df_stacked_field_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category)) +
  geom_rect(
    data = df_both_checker,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#A8C5E8") +
  geom_text(
    data = df_labels_field, aes(x = 102, y = field_id, label = paste0("n = ", scales::comma(n_total))),
    size = 4, hjust = 0) +
  geom_rect(data = df_legend, aes(xmin = x, xmax = x + 4, ymin = y - 0.12, ymax = y + 0.12),
    fill = df_legend$fill, inherit.aes = FALSE) +
  geom_rect(
    data = df_legend_both_checker,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#A8C5E8", inherit.aes = FALSE) +
  geom_text(data = df_legend, aes(x = x + 5, y = y, label = label),
    hjust = 0, size = 3.5, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c("Bar only" = "#D9897B", "Both bar & informative" = "#D9897B",
               "Informative only" = "#A8C5E8", "Others" = "#DCE0E4"), guide = "none") +
  scale_x_continuous(
    limits = c(0, 114),
    breaks = seq(0, 100, by = 25),
    labels = scales::label_percent(scale = 1),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    breaks = df_labels_field$field_id,
    labels = df_labels_field$field,
    limits = c(0.5, legend_y + 0.35),
    expand = expansion(mult = c(0, 0))) +
  labs(title = NULL, x = "Percentage of articles",  y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 2, r = 15, b = 5, l = 5))

stacked_by_field

# Save plot
ggsave(here("results", "stacked_by_field.png"), stacked_by_field, width = 10, height = 5, dpi = 300)


# ---------------
# Temporal trends
# ---------------
# Function to prepare plotting data
# Proportions are relative to: n_total = n_bar_only + n_info_only + n_both + n_others
prepare_plot_data <- function(data) {
  data %>%
    mutate(
      n_bar_only  = coalesce(n_bar, 0L)  - coalesce(n_both, 0L),
      n_info_only = coalesce(n_info, 0L) - coalesce(n_both, 0L),
      n_both      = coalesce(n_both, 0L),
      n_others    = coalesce(n_others, 0L)) %>%
    group_by(field, year) %>%
    summarise(
      n_bar_only  = sum(n_bar_only, na.rm = TRUE),
      n_info_only = sum(n_info_only, na.rm = TRUE),
      n_both      = sum(n_both, na.rm = TRUE),
      n_others    = sum(n_others, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(
      n_total = n_bar_only + n_info_only + n_both + n_others) %>%
    pivot_longer(
      cols = c(n_bar_only, n_info_only, n_both, n_others),
      names_to = "graph_type", values_to = "count") %>%
    mutate(
      graph_type = factor(
        graph_type,
        levels = c("n_bar_only", "n_info_only", "n_both", "n_others")),
      proportion = ifelse(n_total > 0, count / n_total, NA_real_))
}

# Split policy vs non-policy journals
policy_df <- df_descriptive %>%
  filter(policy_ref == 1)

nonpolicy_df <- df_descriptive %>%
  filter(policy_ref == 0)

# Total articles by field: policy journals
field_totals_policy <- policy_df %>%
  mutate(
    n_bar_only  = coalesce(n_bar, 0L)  - coalesce(n_both, 0L),
    n_info_only = coalesce(n_info, 0L) - coalesce(n_both, 0L),
    n_both      = coalesce(n_both, 0L),
    n_others    = coalesce(n_others, 0L),
    n_total     = n_bar_only + n_info_only + n_both + n_others) %>%
  group_by(field) %>%
  summarise(
    total_articles = sum(n_total, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    field_label = paste0(field, "\n(n = ", comma(total_articles), ")"))

# Total articles by field: non-policy journals
field_totals_nonpolicy <- nonpolicy_df %>%
  mutate(
    n_bar_only  = coalesce(n_bar, 0L)  - coalesce(n_both, 0L),
    n_info_only = coalesce(n_info, 0L) - coalesce(n_both, 0L),
    n_both      = coalesce(n_both, 0L),
    n_others    = coalesce(n_others, 0L),
    n_total     = n_bar_only + n_info_only + n_both + n_others) %>%
  group_by(field) %>%
  summarise(
    total_articles = sum(n_total, na.rm = TRUE),  .groups = "drop") %>%
  mutate(
    field_label = paste0(field, "\n(n = ", comma(total_articles), ")"))

# Prepare plotting datasets
policy_long <- prepare_plot_data(policy_df) %>%
  left_join(
    field_totals_policy %>%
      select(field, field_label),
    by = "field")

nonpolicy_long <- prepare_plot_data(nonpolicy_df) %>%
  left_join(
    field_totals_nonpolicy %>%
      select(field, field_label),
    by = "field")

# Policy adoption years
# Each vertical line represents an adoption year.
# Line thickness is proportional to the number of journals adopting in that year.
policy_adoption_years <- policy_df %>%
  filter(
    policy == 1, !is.na(year_since_policy)) %>%
  mutate(
    adoption_year = year - year_since_policy + 1) %>%
  distinct(field, journal, adoption_year) %>%
  count(field, adoption_year, name = "n_journals_started") %>%
  left_join(
    field_totals_policy %>%
      select(field, field_label), by = "field")
trend_colors <- c(n_bar_only  = "#B23A2B", n_info_only = "#2F5DA8", n_both = "#1F7A4C", n_others = "#404040")

# Plot function
make_trend_plot <- function(
    plot_data, adoption_data = NULL,
    year_breaks, x_start, x_end, trend_colors, 
    policy_line_color = "#9370DB" , ncol = 3, base_size = 16) {
  
  p <- ggplot(plot_data, aes(x = year, y = proportion, color = graph_type, group = graph_type))
  
  if (!is.null(adoption_data)) {
    
    p <- p +
      geom_vline(
        data = adoption_data,
        aes(xintercept = adoption_year,
          linewidth = n_journals_started,
          linetype = "Policy adoption"), inherit.aes = FALSE,
        color = policy_line_color,  alpha = 0.8, show.legend = FALSE)
  }
  
  p +
    geom_line(aes(alpha = graph_type), linewidth = 0.8) +
    geom_point(aes(alpha = graph_type), size = 1.3) +
    facet_wrap( ~ field_label, ncol = ncol) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 0.90), breaks = seq(0, 1, by = 0.20)) +
    scale_x_continuous(
      breaks = year_breaks,
      limits = c(x_start, x_end)) +
    scale_color_manual(
      values = trend_colors,
      labels = c(
        n_bar_only  = "Bar only", n_info_only = "Informative only",
        n_both      = "Both", n_others = "Others")) +
    scale_alpha_manual(
      values = c(n_bar_only  = 0.80, n_info_only = 0.80, n_both  = 0.80, 
                 n_others = 0.80), guide = "none") +
    guides(color = guide_legend(order = 1,
        override.aes = list(
          alpha = c(0.80, 0.80, 0.80, 0.80),
          linewidth = 1.2, size = 2)),
      linetype = guide_legend(order = 2),
      linewidth = "none") +
    labs(x = NULL, y = "Percentage of articles", color = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position    = "top",
      legend.box.spacing = unit(0, "mm"),
      legend.margin      = margin(t = 0, r = 0, b = 8, l = 0),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1))
}

# Axis settings
year_breaks <- seq(min(df_descriptive$year), max(df_descriptive$year), by = 3)

x_start <- min(df_descriptive$year)
x_end <- max(df_descriptive$year)

# Create policy plot
policy_plot <- make_trend_plot(
  plot_data     = policy_long,
  adoption_data = policy_adoption_years,
  year_breaks   = year_breaks,
  x_start       = x_start,
  x_end         = x_end,
  trend_colors  = trend_colors)

policy_plot

# Create non-policy plot
nonpolicy_plot <- make_trend_plot(
  plot_data     = nonpolicy_long,
  adoption_data = NULL,
  year_breaks   = year_breaks,
  x_start       = x_start,
  x_end         = x_end,
  trend_colors  = trend_colors)

nonpolicy_plot

# Save plots
ggsave(here("results", "policy_by_year.png"), policy_plot, width = 14, height = 12, dpi = 300)
ggsave(here("results", "nonpolicy_by_year.png"), nonpolicy_plot, width = 14, height = 10, dpi = 300)


