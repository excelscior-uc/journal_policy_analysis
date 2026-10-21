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
# Models 1A - 6 secondary outcome (# of informative graphs)
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
model1A_overdisp_code_info <- nimbleCode({
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

inits_model1A <- function() {
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

label_model1A <- paste0("model1A_info_pois_", mcmc_niter %/% 1000, "k")

res1A <- build_and_run_pois(
  code      = model1A_overdisp_code_info,
  data      = base_data,
  constants = base_constants,
  inits_fn  = inits_model1A,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs", "a", "b", "loglik"),
  label     = label_model1A) # Duration: 00:17:01 (hh:mm:ss)


# ----------------------------------------------
# Model 1B: Base Poisson without overdispersion
# ----------------------------------------------
model1B_code_info <- nimbleCode({
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

inits_model1B <- function() {
  list(
    alpha   = rnorm(1, 3.5, 0.1),
    beta    = rnorm(1, -0.04, 0.01),
    gamma   = rnorm(1, 0, 0.05),
    sigma_a = runif(1, 0.4, 1.0),
    sigma_b = runif(1, 0.01, 0.08),
    a       = rnorm(base_constants$J, 3.5, 0.5),
    b       = rnorm(base_constants$J, -0.04, 0.03))
}

label_model1B <- paste0("model1B_info_pois_", mcmc_niter %/% 1000, "k")

res1B <- build_and_run_pois(
  code      = model1B_code_info,
  data      = base_data,
  constants = base_constants,
  inits_fn  = inits_model1B,
  monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "a", "b", "loglik"),
  label     = label_model1B) # Duration: 00:12:13 (hh:mm:ss)


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
label_model1A_pois <- paste0("model1A_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1A <- readRDS(here("results", paste0("mcmc_list_", label_model1A_pois, ".rds")))
samples_1A <- readRDS(here("results", paste0("samples_", label_model1A_pois, ".rds")))
summarise_model(mcmc_list1A,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs"),
                label = label_model1A_pois)


# ----------------------------------------
# Model 1B: Poisson without overdispersion
# ----------------------------------------
label_model1B_pois <- paste0("model1B_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1B <- readRDS(here("results", paste0("mcmc_list_", label_model1B_pois, ".rds")))
samples_1B <- readRDS(here("results", paste0("samples_", label_model1B_pois, ".rds")))
summarise_model(mcmc_list1B,
                key_params = c("alpha", "beta", "gamma", "sigma_a", "sigma_b"),
                label = label_model1B_pois)


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


## 3. DIAGNOSTICS

# ---------
# Model 1A
# ---------
label_model1A_pois <- paste0("model1A_info_pois_", mcmc_niter %/% 1000, "k")  # secondary outcome
mcmc_list1A <- readRDS(here("results", paste0("mcmc_list_", label_model1A_pois, ".rds")))
samples_1A <- readRDS(here("results", paste0("samples_", label_model1A_pois, ".rds")))
all_samples1A <- do.call(rbind, samples_1A$samples)

key_params1A <- c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs")
labels1A <- c(alpha     = "Intercepts",
              beta      = "Time effects",
              gamma     = "Policy effects",
              sigma_a   = "SD of journal intercepts",
              sigma_b   = "SD of journal slopes",
              sigma_obs = "SD of observation-level errors")

save_trace_density(mcmc_list1A, key_params1A, labels1A, label_model1A_pois, nrow_panels = 3)
save_acf_plot(mcmc_list1A, key_params1A, labels1A, label_model1A_pois)

# DIC calculation
compute_dic_info(samples_1A, label_model1A_pois)


# ---------
# Model 1B
# ---------
label_model1B_pois <- paste0("model1B_info_pois_", mcmc_niter %/% 1000, "k")
mcmc_list1B <- readRDS(here("results", paste0("mcmc_list_", label_model1B_pois, ".rds")))
samples_1B <- readRDS(here("results", paste0("samples_", label_model1B_pois, ".rds")))
all_samples1B <- do.call(rbind, samples_1B$samples)

key_params1B <- c("alpha", "beta", "gamma", "sigma_a", "sigma_b")
labels1B <- c(alpha   = "Intercepts",
              beta    = "Time effects",
              gamma   = "Policy effects",
              sigma_a = "SD of journal intercepts",
              sigma_b = "SD of journal time slopes")

save_trace_density(mcmc_list1B, key_params1B, labels1B, label_model1B_pois, nrow_panels = 3)
save_acf_plot(mcmc_list1B, key_params1B, labels1B, label_model1B_pois)

# DIC calculation
compute_dic_info(samples_1B, label_model1B_pois)


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
compute_dic_info(samples_2,  label_model2_pois)


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
compute_dic_info(samples_3, label_model3_pois)


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
compute_dic_info(samples_4, label_model4_pois)


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
compute_dic_info(samples_5, label_model5_pois)


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
compute_dic_info(samples_6, label_model6_pois)


# -------------------------------------------------------
# DIC comparison table — parse all .txt files and collate 
# -------------------------------------------------------
build_dic_info_table()



## 4. FOREST PLOT
# ----------------------------------------------------------
# Forest plot for main effects (time + policy) - per 5 years
# ----------------------------------------------------------

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
dic <- read_csv(here("results", "dic_comparison_info.csv"), show_col_types = FALSE) %>%
  arrange(DIC) %>%
  mutate(
    Model_short = str_remove(Model, "_info_pois_60k$"),
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
  filter(grepl("_info_", Model)) %>%
  select(Model, model_label, selected, DIC, Delta_DIC, change_type) %>%
  rowwise() %>%
  do({
    model_info <- .
    
    summary_df <- read_csv(
      here("results", paste0("summary_", model_info$Model, ".csv")),
      show_col_types = FALSE) %>%
      filter(Parameter %in% c("beta", "gamma"))
    
    summary_df %>%
      bind_cols(parse_cri(summary_df$`95% CrI`)) %>%
      mutate(
        Model       = model_info$Model,
        model_label = model_info$model_label,
        selected    = model_info$selected,
        DIC         = model_info$DIC,
        Delta_DIC   = model_info$Delta_DIC,
        change_type = model_info$change_type,
        effect = case_when(
          Parameter == "beta"  ~ "Time effect (per 5 years)",
          Parameter == "gamma" ~ "Policy effect",
          TRUE ~ NA_character_),
        
        # Scale time effect from per year to per 5 years
        Mean_scaled = if_else(Parameter == "beta", Mean * 5, Mean),
        low_scaled  = if_else(Parameter == "beta", low * 5, low),
        high_scaled = if_else(Parameter == "beta", high * 5, high),
        
        rate_ratio = exp(Mean_scaled),
        rr_low     = exp(low_scaled),
        rr_high    = exp(high_scaled))}) %>%
  ungroup() %>%
  filter(effect %in% c("Time effect (per 5 years)", "Policy effect"),
         is.finite(rate_ratio), is.finite(rr_low), is.finite(rr_high)) %>%
  mutate(
    effect = factor(effect, levels = c("Time effect (per 5 years)", "Policy effect")))

# ---- Dynamic x-axis limits ---------------------------------------------------
x_min <- floor(min(effects$rr_low,  na.rm = TRUE) * 100) / 100
x_max <- ceiling(max(effects$rr_high, na.rm = TRUE) * 100) / 100

x_range <- x_max - x_min
x_min <- max(0, x_min - x_range * 0.03)
x_max <- x_max + x_range * 0.03

x_breaks <- pretty(c(x_min, x_max), n = 6)
x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
x_labels <- trimws(format(x_breaks, nsmall = 1))

# ---- Vertical layout ---------------------------------------------------------
# Everything is drawn in ONE plot with a single y-scale: the two panels are two
# vertical blocks of that scale, separated by a gap that holds the top panel's
# x-axis labels and the bottom panel's title. The DIC arrows are a second plot
# with the identical y-scale, so they line up exactly with every model row.

top_type    <- "Gradual change"    # panel (a)
bottom_type <- "Step change"       # panel (b)

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
# Right-aligned at x_max. Offsets are fractions of the x-axis range, so the
# layout adapts to the data-driven limits. If labels collide with keys or with
# the panel title, adjust the fractions below (or widen the plot in ggsave).
panel_legend <- function(y) {
  rng        <- x_max - x_min
  key_w      <- 0.025 * rng            # half-width of the key's error bar
  x_pol_txt  <- x_max                  # right edge of "Policy effect"
  x_pol_key  <- x_pol_txt - 0.19 * rng # key sits left of its label
  x_time_txt <- x_pol_key - 0.06 * rng # right edge of "Time effect ..."
  x_time_key <- x_time_txt - 0.31 * rng # key sits left of its label
  
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
  annotate("text", x = x_min, y = a_hi + 0.3, label = "(a) Gradual change models",
           hjust = 0, vjust = 0, fontface = "bold", size = 12 / .pt) +
  annotate("text", x = x_min, y = b_hi + 0.3, label = "(b) Step-change models",
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
  # left label is right-anchored at the null line: the span from x_min to 1 is
  # narrow, so a centered label would be cut off
  annotate("text", x = 1 - dir_gap, y = b_lo - 1.5,
           label = "Informative graph use decreasing",
           hjust = 1, size = 3.8, colour = "grey30") +
  annotate("text", x = x_max, y = b_lo - 1.5,
           label = "Informative graph use increasing",
           hjust = 1, size = 3.8, colour = "grey30") +
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

ggsave(filename = here("results", "forest_plot_info_two_panels_vertical.png"), plot = final_plot, width = 10, height = 8, dpi = 300)



## 5. CUMULATIVE RATE RATIOS

# 1. Color palette + line types
pal <- c("Model 1A" = "#E69F00", "Model 1B" = "#CC79A7", "Model 2" = "#56B4E9",
         "Model 3"  = "#009E73", "Model 4"  = "#F0E442", "Model 5" = "#0072B2",
         "Model 6"  = "#D55E00")

ltys <- c("Model 1A" = "solid", "Model 1B" = "solid", "Model 2" = "solid",
          "Model 3"  = "solid", "Model 4"  = "solid", "Model 5" = "solid",
          "Model 6"  = "solid")

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
  post_mat <- matrix(as.numeric(grid$period == "post"),
                     nrow = n_draws, ncol = n_t, byrow = TRUE)
  
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
draws_1A <- extract_beta_gamma(here("results", "samples_model1A_pois_60k.rds"))
draws_1B <- extract_beta_gamma(here("results", "samples_model1B_pois_60k.rds"))
draws_2  <- extract_beta_gamma(here("results", "samples_model2_pois_60k.rds"))
draws_3  <- extract_beta_gamma(here("results", "samples_model3_pois_60k.rds"))
draws_4  <- extract_beta_gamma(here("results", "samples_model4_pois_60k.rds"))
draws_5  <- extract_beta_gamma(here("results", "samples_model5_pois_60k.rds"))
draws_6  <- extract_beta_gamma(here("results", "samples_model6_pois_60k.rds"))

# 6. Combined data frame
df_event <- bind_rows(
  make_curves_joint(draws_1A, "Model 1A", step = TRUE),
  make_curves_joint(draws_1B, "Model 1B", step = TRUE),
  make_curves_joint(draws_2,  "Model 2",  step = TRUE),
  make_curves_joint(draws_3,  "Model 3",  step = TRUE),
  make_curves_joint(draws_4,  "Model 4",  step = FALSE),
  make_curves_joint(draws_5,  "Model 5",  step = FALSE),
  make_curves_joint(draws_6,  "Model 6",  step = FALSE)) %>%
  mutate(model = factor(model, levels = c("Model 1A", "Model 1B", "Model 2",
                                          "Model 3", "Model 4", "Model 5", "Model 6")))

# 6b. Create jump connectors at event time = 0 (used in the step panel only)
df_jump <- df_event %>%
  filter(event_time == 0) %>%
  select(model, period, rr) %>%
  pivot_wider(names_from = period, values_from = rr)

# 6c. Data-driven y-limits (include credible bands), rounded outward to 0.05
y_step <- 0.05
y_min  <- floor(min(df_event$rr_lo)   / y_step) * y_step
y_max  <- ceiling(max(df_event$rr_hi) / y_step) * y_step
y_lims   <- c(y_min, y_max)
y_breaks <- seq(y_min, y_max, by = y_step)

# 7. Plot 
theme_rr <- theme_bw(base_size = 14) +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
    plot.title       = element_text(face = "bold", hjust = 0))   # left-aligned

# Panel builder: one panel per model type
make_panel <- function(type, title, show_jump = FALSE) {
  d <- df_event %>%
    filter(model_type == type) %>%
    mutate(model = droplevels(model))
  
  mods <- levels(d$model)
  
  p <- ggplot(d, aes(x = event_time, y = rr,
                     colour = model, fill = model, linetype = model,
                     group = interaction(model, period))) +   # breaks the line at t = 0
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.6) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.15, colour = NA)
  
  # jump connector (thin, dotted, same colour); remove if you prefer a visible gap
  if (show_jump) {
    p <- p +
      geom_segment(data = df_jump %>% filter(model %in% mods),
                   aes(x = 0, xend = 0, y = pre, yend = post, colour = model),
                   inherit.aes = FALSE, linewidth = 0.5, linetype = "dotted")
  }
  
  p +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = pal[mods]) +
    scale_fill_manual(values = pal[mods]) +
    scale_linetype_manual(values = ltys[mods]) +
    scale_x_continuous(
      breaks = -3:3,
      labels = c("\u22123", "\u22122", "\u22121", "0", "+1", "+2", "+3"),
      expand = c(0.01, 0)) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = function(x) sprintf("%.2f", x)) +
    coord_cartesian(ylim = y_lims) +             # data-driven: shows full pre-period lines + bands
    labs(title = title,
         x = "Years relative to policy adoption",
         y = "Cumulative rate ratio",
         colour = NULL, fill = NULL, linetype = NULL) +
    theme_rr +
    guides(
      colour = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(linewidth = 1,
                            linetype = unname(ltys[mods]))),
      fill = "none", linetype = "none")
}

p_step <- make_panel("Step change",    "(a) Step change",    show_jump = TRUE)
p_gradual <- make_panel("Gradual change", "(b) Gradual change", show_jump = FALSE)

p_event <- p_step + p_gradual +
  plot_layout(ncol = 2)

p_event

# 8. Save 
ggsave(here("results", "cumulative_rr_policy_event_time_joint.png"), p_event, width = 14, height = 6.5, dpi = 300)
