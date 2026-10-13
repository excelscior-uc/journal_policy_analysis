## ===============
## sensitivity.R
## ===============

library(coda)
library(nimble)
library(parallel)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(fs)
library(stringr)
library(patchwork)

source("scripts/99_functions.R")

# Spike-and-slab sensitvity analyses
# We used a spike-and-slab approach to test whether the policy had a zero effect in some journals and stronger effect in others.
# However, the model did not converge giving a clear indication that this was not a good approach.

# Policy strength sensitivity analyses (Models 8-10)
# The protocol says "only run the best model if there are at least 20 journals in each policy-strength category."

# Check descriptive stat first

#1. Language used to "ban vs. discourage"
# Tests whether journals that explicitly ban bar graphs show larger reductions than those that merely discourage them.
df_sample %>%
  filter(policy_ref == 1) %>%
  distinct(journal, recommendation) %>%
  count(recommendation, name = "unique_journal_count") #recommendation (0) = 36, requirement (1) = 15, mixed(2) = 20

#2. Language used on sample size
# Tests whether policies that require bar graph replacement regardless of sample size are more 
# effective than those applying only below a specified threshold (e.g. n < 10 per group).
df_sample %>%
  filter(policy_ref == 1) %>%
  distinct(journal, sample_size) %>%
  count(sample_size, name = "unique_journal_count")  #recommendation (0) = 27, requirement (1) = 14, mixed(2) = 30

#3. Bar specificity
# Tests whether policies that specifically name bar graphs as inappropriate for continuous data are 
# more impactful than general data-display guidelines.
df_sample %>%
  filter(policy_ref == 1) %>%
  distinct(journal, bar_specificity) %>%
  count(bar_specificity, name = "unique_journal_count") #0=20, 1=51


# --------------------------------------------------------
# Models 1a - 6 secondary outcome (# of informative graphs)
# --------------------------------------------------------
# Load data
df_sample <- read.csv(here("data", "df_sample.csv"))
df_sample$journal_idx <- as.integer(factor(df_sample$journal)) 
df_sample$field_idx <- as.integer(factor(df_sample$field))    

# Shared data / constants
base_data <- list(y = df_sample$n_info)

base_constants <- list(N = nrow(df_sample),
                       J = max(df_sample$journal_idx),
                       log_offset = df_sample$log_offset, 
                       policy = df_sample$policy,
                       year_centered = df_sample$year_centered,
                       journal = df_sample$journal_idx)
# MCMC settings
mcmc_niter   <- 60000
mcmc_nburnin <- 10000
mcmc_thin    <- 10
mcmc_nchains <- 4
mcmc_seeds <- c(20260216, 20260217, 20260218, 20260219)

## 1. MODELS

# -------------------------------------------------------------
# Model 1A: Base Poisson with observation-level overdispersion
# -------------------------------------------------------------
model1a_overdisp_code_info <- nimbleCode({
  for (i in 1:N) {
    epsilon[i] ~ dnorm(0, sd = sigma_obs)
    
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log_offset[i] +  a[journal[i]] + b[journal[i]] * year_centered[i] +  gamma * policy[i] +  epsilon[i]
    
    loglik[i] <- dpois(y[i], mu[i], log = TRUE)    # monitored for DIC
  }
  
  for (j in 1:J) {
    a[j] ~ dnorm(alpha, sd = sigma_a)
    b[j] ~ dnorm(beta,  sd = sigma_b)
  }
  
  # Vague priors
  alpha ~ dnorm(0, sd = 100)
  beta  ~ dnorm(0, sd = 100)
  gamma ~ dnorm(0, sd = 100)
  
  sigma_a   ~ T(dnorm(0, sd = 10), 0, )
  sigma_b   ~ T(dnorm(0, sd = 10), 0, )
  sigma_obs ~ T(dnorm(0, sd = 0.30), 0, ) 
})

inits_model1a <- function() {
  list(
    alpha     = rnorm(1, 3.5, 0.1),
    beta      = rnorm(1, -0.04, 0.01),
    gamma     = rnorm(1, 0, 0.05),
    sigma_a   = runif(1, 0.4, 1.0),
    sigma_b   = runif(1, 0.01, 0.08),
    sigma_obs = runif(1, 0.05, 0.25),
    a         = rnorm(base_constants$J, 3.5, 0.5),
    b         = rnorm(base_constants$J, -0.04, 0.03),
    epsilon   = rnorm(base_constants$N, 0, 0.1))
}

label_model1a <- paste0("model1a_info_pois_", mcmc_niter %/% 1000, "k")

res1a <- build_and_run_pois(
  code      = model1a_overdisp_code_info,
  data      = base_data,
  constants = base_constants,
  inits_fn  = inits_model1a,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs", "a", "b", "loglik"),
  label     = label_model1a) # Duration: 00:17:01 (hh:mm:ss)


# ----------------------------------------------
# Model 1B: Base Poisson without overdispersion
# ----------------------------------------------
model1b_code_info <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dpois(mu[i])
    
    log(mu[i]) <- log_offset[i] +  a[journal[i]] +  b[journal[i]] * year_centered[i] + gamma * policy[i]
    
    loglik[i] <- dpois(y[i], mu[i], log = TRUE)
  }
  
  for (j in 1:J) {
    a[j] ~ dnorm(alpha, sd = sigma_a)
    b[j] ~ dnorm(beta,  sd = sigma_b)
  }
  
  # Vague priors 
  alpha ~ dnorm(0, sd = 100)
  beta  ~ dnorm(0, sd = 100)
  gamma ~ dnorm(0, sd = 100)
  
  sigma_a ~ T(dnorm(0, sd = 10), 0, )
  sigma_b ~ T(dnorm(0, sd = 10), 0, )
})

inits_model1b <- function() {
  list(
    alpha   = rnorm(1, 3.5, 0.1),
    beta    = rnorm(1, -0.04, 0.01),
    gamma   = rnorm(1, 0, 0.05),
    sigma_a = runif(1, 0.4, 1.0),
    sigma_b = runif(1, 0.01, 0.08),
    a       = rnorm(base_constants$J, 3.5, 0.5),
    b       = rnorm(base_constants$J, -0.04, 0.03))
}

label_model1b <- paste0("model1b_info_pois_", mcmc_niter %/% 1000, "k")

res1b <- build_and_run_pois(
  code      = model1b_code_info,
  data      = base_data,
  constants = base_constants,
  inits_fn  = inits_model1b,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "a", "b", "loglik"),
  label     = label_model1b) # Duration: 00:12:13 (hh:mm:ss)


# ------------------------------------------------------
# Model 2: Poisson, fixed time trend, random intercepts
# ------------------------------------------------------
model2_code_info <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dpois(mu[i])
    
    log(mu[i]) <- log_offset[i] + a[journal[i]] + beta * year_centered[i] +  gamma * policy[i]
    
    loglik[i] <- dpois(y[i], mu[i], log = TRUE)
  }
  
  for (j in 1:J) {
    a[j] ~ dnorm(alpha, sd = sigma_a)
  }
  
  # Vague priors
  alpha ~ dnorm(0, sd = 100)
  beta  ~ dnorm(0, sd = 100)
  gamma ~ dnorm(0, sd = 100)
  
  sigma_a ~ T(dnorm(0, sd = 10), 0, )
})

inits_model2 <- function() {
  list(
    alpha   = rnorm(1, 3.5, 0.1),
    beta    = rnorm(1, -0.04, 0.01),
    gamma   = rnorm(1, 0, 0.05),
    sigma_a = runif(1, 0.4, 1.0),
    a       = rnorm(base_constants$J, 3.5, 0.5))
}

label_model2_pois <- paste0("model2_info_pois_", mcmc_niter %/% 1000, "k")

res2 <- build_and_run_pois(
  code      = model2_code_info,
  data      = base_data,
  constants = base_constants,
  inits_fn  = inits_model2,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "a", "loglik"),
  label     = label_model2_pois) # Duration: 00:11:42 (hh:mm:ss)


# -------------------------------------------------------
# Model 3: Nested: journals within fields, random slopes
# -------------------------------------------------------
model3_code_info <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dpois(mu[i])
    
    log(mu[i]) <- log_offset[i] +  a[journal[i]] + b[journal[i]] * year_centered[i] + gamma * policy[i]
    
    loglik[i] <- dpois(y[i], mu[i], log = TRUE)
  }
  
  for (j in 1:J) {
    a[j] ~ dnorm(alpha_f[field[j]], sd = sigma_a)
    b[j] ~ dnorm(beta, sd = sigma_b)
  }
  
  for (f in 1:NF) {
    alpha_f[f] ~ dnorm(A, sd = sigma_f)
  }
  
  # Vague priors
  A     ~ dnorm(0, sd = 100)
  beta  ~ dnorm(0, sd = 100)
  gamma ~ dnorm(0, sd = 100)
  
  sigma_a ~ T(dnorm(0, sd = 10), 0, )
  sigma_b ~ T(dnorm(0, sd = 10), 0, )
  sigma_f ~ T(dnorm(0, sd = 10), 0, )
})
journal_field <- tapply(df_sample$field_idx, df_sample$journal_idx, function(x) unique(x))
stopifnot(all(lengths(journal_field) == 1))
journal_field <- as.integer(unlist(journal_field))

data3 <- base_data
constants3 <- c(base_constants, list(NF = max(journal_field), field = journal_field))

inits_model3 <- function() {
  NF <- max(constants$NF)   # 'constants' is exported to workers
  J  <- constants$J
  list(
    A       = rnorm(1, 3.5, 0.1),
    beta    = rnorm(1, -0.04, 0.01),
    gamma   = rnorm(1, 0, 0.05),
    sigma_a = runif(1, 0.4, 1.0),
    sigma_b = runif(1, 0.01, 0.08),
    sigma_f = runif(1, 0.4, 1.0),
    alpha_f = rnorm(NF, 3.5, 0.5),
    a       = rnorm(J, 3.5, 0.5),
    b       = rnorm(J, -0.04, 0.03))
}

label_model3 <- paste0("model3_info_pois_", mcmc_niter %/% 1000, "k")

res3 <- build_and_run_pois(
  code      = model3_code_info,
  data      = data3,
  constants = constants3,
  inits_fn  = inits_model3,
  monitors  = c("A", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f", "alpha_f", "a", "b", "loglik"),
  label     = label_model3) # Duration: 00:13:42 (hh:mm:ss)


# ----------------------------------------------------------
# Model 4: Gradual policy effect (policy_year), fixed slope
# ----------------------------------------------------------
base_constants4 <- list(N = nrow(df_sample),
                        J = max(df_sample$journal_idx),
                        log_offset = df_sample$log_offset,
                        policy = df_sample$policy,
                        year_centered = df_sample$year_centered,
                        policy_year = df_sample$year_since_policy,
                        journal = df_sample$journal_idx)
model4_code_info <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dpois(mu[i])
    
    log(mu[i]) <- log_offset[i] +  a[journal[i]] +  beta * year_centered[i] +  gamma * policy_year[i]
    
    loglik[i] <- dpois(y[i], mu[i], log = TRUE)
  }
  
  for (j in 1:J) {
    a[j] ~ dnorm(alpha, sd = sigma_a)
  }
  
  # Vague priors
  alpha ~ dnorm(0, sd = 100)
  beta  ~ dnorm(0, sd = 100)
  gamma ~ dnorm(0, sd = 100)
  
  sigma_a ~ T(dnorm(0, sd = 10), 0, )
})

inits_model4 <- function() {
  list(
    alpha   = rnorm(1, 3.5, 0.1),
    beta    = rnorm(1, -0.04, 0.01),
    gamma   = rnorm(1, 0, 0.05),
    sigma_a = runif(1, 0.4, 1.0),
    a       = rnorm(base_constants$J, 3.5, 0.5))
}

label_model4 <- paste0("model4_info_pois_", mcmc_niter %/% 1000, "k")

res4 <- build_and_run_pois(
  code      = model4_code_info,
  data      = base_data,
  constants = base_constants4,
  inits_fn  = inits_model4,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "a", "loglik"),
  label     = label_model4) # Duration: 00:12:31 (hh:mm:ss)


# --------------------------------------------------------
# Model 5: Gradual policy effect + journal-specific slope
# --------------------------------------------------------
base_constants5 <- list(N = nrow(df_sample),
                        J = max(df_sample$journal_idx),
                        log_offset = df_sample$log_offset,
                        policy = df_sample$policy,
                        year_centered = df_sample$year_centered,
                        policy_year = df_sample$year_since_policy,
                        journal = df_sample$journal_idx)

model5_code_info <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dpois(mu[i])
    
    log(mu[i]) <- log_offset[i] + a[journal[i]] + b[journal[i]] * year_centered[i] + gamma * policy_year[i]
    
    loglik[i] <- dpois(y[i], mu[i], log = TRUE)
  }
  
  for (j in 1:J) {
    a[j] ~ dnorm(alpha, sd = sigma_a)
    b[j] ~ dnorm(beta,  sd = sigma_b)
  }
  
  # Vague priors
  alpha ~ dnorm(0, sd = 100)
  beta  ~ dnorm(0, sd = 100)
  gamma ~ dnorm(0, sd = 100)
  
  sigma_a ~ T(dnorm(0, sd = 10), 0, )
  sigma_b ~ T(dnorm(0, sd = 10), 0, )
})

inits_model5 <- function() {
  list(
    alpha   = rnorm(1, 3.5, 0.1),
    beta    = rnorm(1, -0.04, 0.01),
    gamma   = rnorm(1, 0, 0.05),
    sigma_a = runif(1, 0.4, 1.0),
    sigma_b = runif(1, 0.01, 0.08),
    a       = rnorm(base_constants$J, 3.5, 0.5),
    b       = rnorm(base_constants$J, -0.04, 0.03))
}

label_model5 <- paste0("model5_info_pois_", mcmc_niter %/% 1000, "k")

res5 <- build_and_run_pois(
  code      = model5_code_info,
  data      = base_data,
  constants = base_constants5,
  inits_fn  = inits_model5,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "a", "b", "loglik"),
  label     = label_model5) # Duration: 00:12:08 (hh:mm:ss)


# -------------------------------------------------
# Model 6: Gradual policy effect + field hierarchy
# -------------------------------------------------
model6_code_info <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dpois(mu[i])
    
    log(mu[i]) <- log_offset[i] + a[journal[i]] + b[journal[i]] * year_centered[i] + gamma * policy_year[i]
    
    loglik[i] <- dpois(y[i], mu[i], log = TRUE)
  }
  
  for (j in 1:J) {
    a[j] ~ dnorm(alpha_f[field[j]], sd = sigma_a)
    b[j] ~ dnorm(beta, sd = sigma_b)
  }
  
  for (f in 1:NF) {
    alpha_f[f] ~ dnorm(alpha, sd = sigma_f)
  }
  
  # Vague priors
  alpha ~ dnorm(0, sd = 100)
  beta  ~ dnorm(0, sd = 100)
  gamma ~ dnorm(0, sd = 100)
  
  sigma_a ~ T(dnorm(0, sd = 10), 0, )
  sigma_b ~ T(dnorm(0, sd = 10), 0, )
  sigma_f ~ T(dnorm(0, sd = 10), 0, )
})

# Journal-to-field mapping
journal_field <- tapply(df_sample$field_idx, df_sample$journal_idx, function(x) unique(x))
stopifnot(all(lengths(journal_field) == 1))
journal_field <- as.integer(unlist(journal_field))

data6 <- base_data
constants6 <- c(base_constants, list(policy_year = df_sample$year_since_policy, NF = max(journal_field), field = journal_field))

inits_model6 <- function() {
  NF <- constants$NF
  J  <- constants$J
  
  list(
    alpha   = rnorm(1, 3.5, 0.1),
    beta    = rnorm(1, -0.04, 0.01),
    gamma   = rnorm(1, 0, 0.05),
    sigma_a = runif(1, 0.4, 1.0),
    sigma_b = runif(1, 0.01, 0.08),
    sigma_f = runif(1, 0.4, 1.0),
    alpha_f = rnorm(NF, 3.5, 0.5),
    a       = rnorm(J, 3.5, 0.5),
    b       = rnorm(J, -0.04, 0.03))
}

label_model6 <- paste0("model6_info_pois_", mcmc_niter %/% 1000, "k")

res6 <- build_and_run_pois(
  code      = model6_code_info,
  data      = data6,
  constants = constants6,
  inits_fn  = inits_model6,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f", "a", "b", "alpha_f", "loglik"),
  label     = label_model6) # Duration: 00:13:44 (hh:mm:ss)


## 2. RESULTS

# -------------------------------------
# Model 1A: Poisson with overdispersion
# -------------------------------------
label_model1a_pois <- paste0("model1a_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1a <- readRDS(here("results", paste0("mcmc_list_", label_model1a_pois, ".rds")))
samples_1a <- readRDS(here("results", paste0("samples_", label_model1a_pois, ".rds")))
summarise_model(mcmc_list1a,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs"),
                label = label_model1a_pois)


# ----------------------------------------
# Model 1B: Poisson without overdispersion
# ----------------------------------------
label_model1b_pois <- paste0("model1b_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1b <- readRDS(here("results", paste0("mcmc_list_", label_model1b_pois, ".rds")))
samples_1b <- readRDS(here("results", paste0("samples_", label_model1b_pois, ".rds")))
summarise_model(mcmc_list1b,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b"),
                label = label_model1b_pois)


# -------------------------------------------------
# Model 2: Fixed time trend, random intercepts only
# -------------------------------------------------
label_model2_pois <- paste0("model2_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list2 <- readRDS(here("results", paste0("mcmc_list_", label_model2_pois, ".rds")))
samples_2 <- readRDS(here("results", paste0("samples_", label_model2_pois, ".rds")))
summarise_model(mcmc_list2,
                key_params = c("alpha", "beta", "gamma", "sigma_a"),
                label = label_model2_pois)


# ---------------------------------------------------
# Model 3: Nested intercepts (journals within fields)
# ---------------------------------------------------
label_model3_pois <- paste0("model3_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list3 <- readRDS(here("results", paste0("mcmc_list_", label_model3_pois, ".rds")))
samples_3 <- readRDS(here("results", paste0("samples_", label_model3_pois, ".rds")))
summarise_model(mcmc_list3,
                key_params = c("A", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f"),
                label = label_model3_pois)

# -------------------------------------------
# Model 4: Gradual policy effect, fixed slope
# -------------------------------------------
label_model4_pois <- paste0("model4_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list4 <- readRDS(here("results", paste0("mcmc_list_", label_model4_pois, ".rds")))
samples_4 <- readRDS(here("results", paste0("samples_", label_model4_pois, ".rds")))
summarise_model(mcmc_list4,
                key_params = c("alpha", "beta", "gamma", "sigma_a"),
                label = label_model4_pois)


# --------------------------------------------------------
# Model 5: Gradual policy effect + journal-specific slopes
# --------------------------------------------------------
label_model5_pois <- paste0("model5_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list5 <- readRDS(here("results", paste0("mcmc_list_", label_model5_pois, ".rds")))
samples_5 <- readRDS(here("results", paste0("samples_", label_model5_pois, ".rds")))
summarise_model(mcmc_list5,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b"),
                label = label_model5_pois)


# ------------------------------------------------
# Model 6: Gradual policy effect + field hierarchy
# ------------------------------------------------
label_model6_pois <- paste0("model6_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list6 <- readRDS(here("results", paste0("mcmc_list_", label_model6_pois, ".rds")))
samples_6 <- readRDS(here("results", paste0("samples_", label_model6_pois, ".rds")))
summarise_model(mcmc_list6,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f"),
                label = label_model6_pois)


## 3. DIAGNOSTICS (don't forget to uncomment few lines in "99_functions.R" before you run)

# ---------
# Model 1A
# ---------
label_model1a_pois <- paste0("model1a_info_pois_", mcmc_niter %/% 1000, "k")  # secondary outcome
mcmc_list1a <- readRDS(here("results", paste0("mcmc_list_", label_model1a_pois, ".rds")))
samples_1a <- readRDS(here("results", paste0("samples_", label_model1a_pois, ".rds")))
all_samples1a <- do.call(rbind, samples_1a$samples)

key_params1a <- c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs")
labels1a <- c(alpha     = "Intercepts",
              beta      = "Time effects",
              gamma     = "Policy effects",
              sigma_a   = "SD of journal intercepts",
              sigma_b   = "SD of journal slopes",
              sigma_obs = "SD of observation-level errors")

save_trace_density(mcmc_list1a, key_params1a, labels1a, label_model1a_pois, nrow_panels = 3)
save_acf_plot(mcmc_list1a, key_params1a, labels1a, label_model1a_pois)

# DIC calculation
compute_dic(samples_1a, label_model1a_pois)


# ---------
# Model 1B
# ---------
label_model1b_pois <- paste0("model1b_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1b <- readRDS(here("results", paste0("mcmc_list_", label_model1b_pois, ".rds")))
samples_1b <- readRDS(here("results", paste0("samples_", label_model1b_pois, ".rds")))
all_samples1b <- do.call(rbind, samples_1b$samples)

key_params1b <- c("alpha", "beta", "gamma", "sigma_a", "sigma_b")
labels1b <- c(alpha   = "Intercepts",
              beta    = "Time effects",
              gamma   = "Policy effects",
              sigma_a = "SD of journal intercepts",
              sigma_b = "SD of journal time slopes")

save_trace_density(mcmc_list1b, key_params1b, labels1b, label_model1b_pois, nrow_panels = 3)
save_acf_plot(mcmc_list1b, key_params1b, labels1b, label_model1b_pois)

# DIC calculation
compute_dic(samples_1b, label_model1b_pois)


# --------
# Model 2
# --------
label_model2_pois <- paste0("model2_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list2 <- readRDS(here("results", paste0("mcmc_list_", label_model2_pois, ".rds")))
samples_2 <- readRDS(here("results", paste0("samples_", label_model2_pois, ".rds")))
all_samples2 <- do.call(rbind, samples_2$samples)

key_params2 <- c("alpha", "beta", "gamma", "sigma_a")
labels2 <- c(alpha   = "Intercepts",
             beta    = "Time effects",
             gamma   = "Policy effects",
             sigma_a = "SD of journal intercepts")

save_trace_density(mcmc_list2, key_params2, labels2, label_model2_pois, nrow_panels = 2)
save_acf_plot(mcmc_list2, key_params2, labels2, label_model2_pois)

# DIC calculation
compute_dic(samples_2,  label_model2_pois)


# --------
# Model 3
# --------
label_model3_pois <- paste0("model3_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list3 <- readRDS(here("results", paste0("mcmc_list_", label_model3_pois, ".rds")))
samples_3 <- readRDS(here("results", paste0("samples_", label_model3_pois, ".rds")))
all_samples3 <- do.call(rbind, samples_3$samples)

key_params3 <- c("A", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f")
labels3 <- c(A       = "Intercepts", 
             beta    = "Time effects",
             gamma   = "Policy effects",
             sigma_a = "SD of journal intercepts",
             sigma_b = "SD of journal slopes",
             sigma_f = "SD of field-level intercepts")

save_trace_density(mcmc_list3, key_params3, labels3, label_model3_pois, nrow_panels = 3)
save_acf_plot(mcmc_list3, key_params3, labels3, label_model3_pois)

# DIC calculation
compute_dic(samples_3, label_model3_pois)


# --------
# Model 4
# --------
label_model4_pois <- paste0("model4_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list4 <- readRDS(here("results", paste0("mcmc_list_", label_model4_pois, ".rds")))
samples_4 <- readRDS(here("results", paste0("samples_", label_model4_pois, ".rds")))
all_samples4 <- do.call(rbind, samples_4$samples)

key_params4 <- c("alpha", "beta", "gamma", "sigma_a")
labels4 <- c(alpha   = "Intercepts",
             beta    = "Time effects",
             gamma   = "Policy effects",
             sigma_a = "SD of journal-specific intercepts")

save_trace_density(mcmc_list4, key_params4, labels4, label_model4_pois, nrow_panels = 2)
save_acf_plot(mcmc_list4, key_params4, labels4, label_model4_pois)

# DIC calculation
compute_dic(samples_4, label_model4_pois)


# --------
# Model 5
# --------
label_model5_pois <- paste0("model5_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list5 <- readRDS(here("results", paste0("mcmc_list_", label_model5_pois, ".rds")))
samples_5 <- readRDS(here("results", paste0("samples_", label_model5_pois, ".rds")))
all_samples5 <- do.call(rbind, samples_5$samples)

key_params5 <- c("alpha", "beta", "gamma", "sigma_a", "sigma_b")
labels5 <- c(alpha   = "Intercepts",
             beta    = "Time effect",
             gamma   = "Policy effect",
             sigma_a = "SD of journal-specific intercepts",
             sigma_b = "SD of journal-specific time slopes")

save_trace_density(mcmc_list5, key_params5, labels5, label_model5_pois, nrow_panels = 3)
save_acf_plot(mcmc_list5, key_params5, labels5, label_model5_pois)

# DIC calculation
compute_dic(samples_5, label_model5_pois)


# --------
# Model 6
# --------
label_model6_pois <- paste0("model6_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list6 <- readRDS(here("results", paste0("mcmc_list_", label_model6_pois, ".rds")))
samples_6 <- readRDS(here("results", paste0("samples_", label_model6_pois, ".rds")))
all_samples6 <- do.call(rbind, samples_6$samples)

key_params6 <- c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f")
labels6 <- c(alpha   = "Intercepts",
             beta    = "Time effects",
             gamma   = "Policy effects",
             sigma_a = "SD of journal-specific intercepts",
             sigma_b = "SD of journal-specific time slopes",
             sigma_f = "SD of field-specific intercepts")

save_trace_density(mcmc_list6, key_params6, labels6, label_model6_pois, nrow_panels = 3)
save_acf_plot(mcmc_list6, key_params6, labels6, label_model6_pois)

# DIC calculation
compute_dic(samples_6, label_model6_pois)


# -------------------------------------------------------
# DIC comparison table — parse all .txt files and collate 
# -------------------------------------------------------
build_dic_table()



## 4. FOREST PLOT
model_names <- tibble(
  Model_short = c("model1A", "model1B", "model2", "model3", "model4", "model5",  "model6"),
  model_label = c(
    "Model 1A: Base Poisson model - step change",  # with overdispersion
    "Model 1B: Base Poisson model - step change",  # without overdispersion
    "Model 2: Policy leakage - step change",
    "Model 3: Field specific hierarchy - step change",
    "Model 4: Base Poisson model - gradual change",
    "Model 5: Journal specific trend - gradual change",
    "Model 6: Field specific hierarchy - gradual change"))

dic <- read_csv(here("results", "dic_comparison_info.csv"), show_col_types = FALSE) %>%
  arrange(DIC) %>%
  mutate(Model_short = str_remove(Model, "_info_pois_60k$"), selected = Delta_DIC == 0) %>%
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
  filter(grepl("_info_", Model)) %>%
  select(Model, model_label, selected, DIC, Delta_DIC) %>%
  rowwise() %>%
  do({
    model_info <- .
    
    summary_df <- read_csv(here("results", paste0("summary_", model_info$Model, ".csv")), show_col_types = FALSE) %>%
      filter(Parameter %in% c("beta", "gamma"))
    
    summary_df %>%
      bind_cols(parse_cri(summary_df$`95% CrI`)) %>%
      mutate(
        Model = model_info$Model,
        model_label = model_info$model_label,
        selected = model_info$selected,
        DIC = model_info$DIC,
        Delta_DIC = model_info$Delta_DIC,
        effect = case_when(
          Parameter == "beta" ~ "Time effect (per 5 years)",
          Parameter == "gamma" ~ "Policy effect",
          TRUE ~ NA_character_),
        
        # Scale time effect from per year to per 5 years
        Mean_scaled = if_else(Parameter == "beta", Mean * 5, Mean),
        low_scaled = if_else(Parameter == "beta", low * 5, low),
        high_scaled = if_else(Parameter == "beta", high * 5, high),
        
        rate_ratio = exp(Mean_scaled),
        rr_low = exp(low_scaled),
        rr_high = exp(high_scaled))}) %>%
  ungroup() %>%
  filter(effect %in% c("Time effect (per 5 years)", "Policy effect"),
         is.finite(rate_ratio), is.finite(rr_low), is.finite(rr_high)) %>%
  mutate(
    model_label = factor(model_label, levels = rev(dic$model_label)),  
    effect = factor(effect, levels = c("Time effect (per 5 years)", "Policy effect")))

# Combine both effects onto one row per model
n_models <- nlevels(effects$model_label)

effects <- effects %>%
  mutate(
    y_base = as.numeric(model_label),
    y_pos  = y_base + if_else(effect == "Time effect (per 5 years)", 0.16, -0.16))

y_breaks <- seq_len(n_models)
y_labels <- levels(effects$model_label)

# Dynamic x-axis limits based on scaled confidence intervals
x_min <- floor(min(effects$rr_low, na.rm = TRUE) * 100) / 100
x_max <- ceiling(max(effects$rr_high, na.rm = TRUE) * 100) / 100

# Add a small margin to both sides
x_range <- x_max - x_min
x_min <- max(0, x_min - x_range * 0.03)
x_max <- x_max + x_range * 0.03

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
    values = c("Time effect (per 5 years)" = "#0072B2", "Policy effect" = "#E69F00"),
    name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = y_limits, clip = "off", expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position      = c(0.60, 1.046),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = "white", colour = NA, linewidth = 0.3),
    legend.title         = element_blank(),
    legend.key.size      = unit(0.85, "lines"),
    panel.grid.minor     = element_blank(),
    axis.text.y          = element_text(size = 12),
    plot.margin          = margin(t = 8, r = 12, b = 5, l = 5))

# Direction-of-effect strip
gap <- (x_max - x_min) * 0.015

direction_plot <- ggplot() +
  annotate("segment", x = 1 - gap, xend = x_min, y = 1, yend = 1,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("segment", x = 1 + gap, xend = x_max, y = 1, yend = 1,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           colour = "grey30", linewidth = 0.3) +
  annotate("text", x = (1 + x_min) / 2, y = 0.35, label = "Informative graph use decreasing",
           size = 0.42, colour = "grey30") +
  annotate("text", x = (1 + x_max) / 2, y = 0.35, label = "Informative graph use increasing",
           size = 3.8, colour = "grey30") +
  annotate("text", x = (x_min + x_max) / 2, y = -0.15, label = "Rate ratio",
           size = 4.2, colour = "black", fontface = "bold") +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(-0.5, 1.5), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 5))

# DIC arrow strip to the left of the model rows
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

ggsave(filename = here("results", "forest_plot_all_models_info.png"), plot = final_plot, width = 10, height = 6, dpi = 300)



## 5. CUMULATIVE RATIOS
# 1. Color palette + line types
pal <- c("Model 1A" = "#E69F00", "Model 1B" = "#CC79A7", "Model 2" = "#56B4E9", "Model 3" = "#009E73",
         "Model 4"  = "#F0E442", "Model 5"  = "#0072B2", "Model 6" = "#D55E00")
ltys <- c("Model 1A" = "solid",    "Model 1B" = "solid",    "Model 2" = "solid", "Model 3" = "solid",
          "Model 4"  = "longdash", "Model 5"  = "longdash", "Model 6" = "longdash")

# 2. Event-time grid: pre-policy piece ends AT 0, post-policy piece starts AT 0
#    (integers / 10 avoid floating-point drift, e.g. -1.4e-16 instead of 0)
grid <- data.frame(
  event_time = c(seq(-30, 0) / 10, seq(0, 30) / 10),
  period     = rep(c("pre", "post"), each = 31))

# 3. Helper: extract beta & gamma from the list-of-chains object
extract_beta_gamma <- function(rds_path) {
  obj <- readRDS(rds_path)
  mat <- do.call(rbind, obj$samples)          # combine chains
  data.frame(beta = mat[, "beta"], gamma = mat[, "gamma"])
}

# 4. Helper: joint-posterior cumulative RR curves
make_curves_joint <- function(draws, model, step = TRUE) {
  n_draws <- nrow(draws)
  n_t     <- nrow(grid)
  
  t_mat    <- matrix(grid$event_time, nrow = n_draws, ncol = n_t, byrow = TRUE)
  post_mat <- matrix(as.numeric(grid$period == "post"), nrow = n_draws, ncol = n_t, byrow = TRUE)
  
  # step: indicator (0 in the pre piece incl. t = 0, 1 in the post piece incl. t = 0)
  # gradual: years since adoption (continuous at 0)
  pol_mat <- if (step) post_mat else pmax(t_mat, 0)
  
  lp <- draws$beta * t_mat + draws$gamma * pol_mat   # beta/gamma recycle down columns
  
  rr_draws <- exp(lp)
  
  data.frame(
    model      = model,
    model_type = ifelse(step, "Step change", "Gradual change"),
    event_time = grid$event_time,
    period     = grid$period,
    rr         = colMeans(rr_draws),
    rr_lo      = apply(rr_draws, 2, quantile, probs = 0.025),
    rr_hi      = apply(rr_draws, 2, quantile, probs = 0.975))
}

# 5. Load draws for every model
draws_1a <- extract_beta_gamma(here("results", "samples_model1a_pois_60k.rds"))
draws_1b <- extract_beta_gamma(here("results", "samples_model1b_pois_60k.rds"))
draws_2  <- extract_beta_gamma(here("results", "samples_model2_pois_60k.rds"))
draws_3  <- extract_beta_gamma(here("results", "samples_model3_pois_60k.rds"))
draws_4  <- extract_beta_gamma(here("results", "samples_model4_pois_60k.rds"))
draws_5  <- extract_beta_gamma(here("results", "samples_model5_pois_60k.rds"))
draws_6  <- extract_beta_gamma(here("results", "samples_model6_pois_60k.rds"))

# 6. Combined data frame
df_event <- bind_rows(
  make_curves_joint(draws_1a, "Model 1A", step = TRUE),
  make_curves_joint(draws_1b, "Model 1B", step = TRUE),
  make_curves_joint(draws_2,  "Model 2",  step = TRUE),
  make_curves_joint(draws_3,  "Model 3",  step = TRUE),
  make_curves_joint(draws_4,  "Model 4",  step = FALSE),
  make_curves_joint(draws_5,  "Model 5",  step = FALSE),
  make_curves_joint(draws_6,  "Model 6",  step = FALSE)) %>%
  mutate(model = factor(model, levels = c("Model 1A", "Model 1B", "Model 2",
                                          "Model 3", "Model 4", "Model 5", "Model 6")))

# 6b. Create jump connectors at event time = 0
df_jump <- df_event %>%
  filter(event_time == 0) %>%
  select(model, period, rr) %>%
  pivot_wider(names_from = period, values_from = rr)

# 7. Plot
theme_rr <- theme_bw(base_size = 14) +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.key.width = unit(2, "cm"),
    legend.spacing.x = unit(0.4, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.3))

p_event <- ggplot(df_event, aes(x = event_time, y = rr,
                                colour = model, fill = model, linetype = model,
                                group = interaction(model, period))) +     # <- breaks the line at t = 0
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.6) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  # jump connector (thin, dotted, same colour) -- remove this block if you prefer a visible gap
  geom_segment(data = df_jump,
               aes(x = 0, xend = 0, y = pre, yend = post, colour = model),
               inherit.aes = FALSE, linewidth = 0.5, linetype = "dotted") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_linetype_manual(values = ltys) +
  scale_x_continuous(
    breaks = -3:3,
    labels = c("\u22123", "\u22122", "\u22121", "0", "+1", "+2", "+3"),
    expand = c(0.01, 0)) +
  scale_y_continuous(
    breaks = seq(0.65, 1.10, by = 0.05),
    labels = function(x) sprintf("%.2f", x)) +
  coord_cartesian(ylim = c(0.65, 1.10)) +       # lower limit widened: Models 5-6 fall below 0.70 at +3
  labs(x = "Years relative to policy adoption",
       y = "Cumulative rate ratio",
       colour = NULL, fill = NULL, linetype = NULL) +
  theme_rr +
  guides(
    colour = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(linewidth = 1,
                          linetype = unname(ltys))),
    fill = "none", linetype = "none")

p_event

# 8. Save
ggsave(here("results", "cumulative_rr_policy_event_time_joint.png"),
       p_event, width = 12, height = 7.5, dpi = 300)
