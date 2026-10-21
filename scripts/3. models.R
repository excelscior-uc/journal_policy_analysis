## ==========
## models.R
## ==========

library(here)
library(nimble)
library(coda)
library(parallel)
library(tidyverse)
library(openxlsx)

source("scripts/99_functions.R")

# Load data
df_sample <- read.csv(here("data", "df_sample.csv"))
dim(df_sample)

df_sample$journal_idx <- as.integer(factor(df_sample$journal)) 
df_sample$field_idx <- as.integer(factor(df_sample$field))    

# Shared data / constants
base_data <- list(y = df_sample$n_bar)

base_constants <- list(N = nrow(df_sample),
                       J = max(df_sample$journal_idx),
                       log_offset = df_sample$log_offset, 
                       policy = df_sample$policy,
                       year_centered = df_sample$year_centered,
                       journal = df_sample$journal_idx)
# MCMC settings
mcmc_niter <- 60000
mcmc_nburnin <- 10000
mcmc_thin <- 10
mcmc_nchains <- 4
mcmc_seeds <- c(20260216, 20260217, 20260218, 20260219)

# -------------------------------------------------------------
# Model 1A: Base Poisson with observation-level overdispersion
# -------------------------------------------------------------
model1A_overdisp_code <- nimbleCode({
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
  list(alpha     = rnorm(1, 3.5, 0.1),
       beta      = rnorm(1, -0.04, 0.01),
       gamma     = rnorm(1, 0, 0.05),
       sigma_a   = runif(1, 0.4, 1.0),
       sigma_b   = runif(1, 0.01, 0.08),
       sigma_obs = runif(1, 0.05, 0.25),
       a         = rnorm(base_constants$J, 3.5, 0.5),
       b         = rnorm(base_constants$J, -0.04, 0.03),
       epsilon   = rnorm(base_constants$N, 0, 0.1))
}

label_model1A <- paste0("model1A_pois_", mcmc_niter %/% 1000, "k")

res1A <- build_and_run_pois(code      = model1A_overdisp_code,
                            data      = base_data,
                            constants = base_constants,
                            inits_fn  = inits_model1A,
                            monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_obs", "a", "b", "loglik"),
                            label     = label_model1A) 


# ----------------------------------------------
# Model 1B: Base Poisson without overdispersion
# ----------------------------------------------
model1B_code <- nimbleCode({
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
  list(alpha   = rnorm(1, 3.5, 0.1),
       beta    = rnorm(1, -0.04, 0.01),
       gamma   = rnorm(1, 0, 0.05),
       sigma_a = runif(1, 0.4, 1.0),
       sigma_b = runif(1, 0.01, 0.08),
       a       = rnorm(base_constants$J, 3.5, 0.5),
       b       = rnorm(base_constants$J, -0.04, 0.03))
}

label_model1B <- paste0("model1B_pois_", mcmc_niter %/% 1000, "k")

res1B <- build_and_run_pois(code      = model1B_code,
                            data      = base_data,
                            constants = base_constants,
                            inits_fn  = inits_model1B,
                            monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "a", "b", "loglik"),
                            label     = label_model1B) 


# -------------------------------------------------------
# Model 2 — Poisson, fixed time trend, random intercepts
# -------------------------------------------------------
model2_code <- nimbleCode({
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
  list(alpha   = rnorm(1, 3.5, 0.1),
       beta    = rnorm(1, -0.04, 0.01),
       gamma   = rnorm(1, 0, 0.05),
       sigma_a = runif(1, 0.4, 1.0),
       a       = rnorm(base_constants$J, 3.5, 0.5))
}

label_model2_pois <- paste0("model2_pois_", mcmc_niter %/% 1000, "k")

res2 <- build_and_run_pois(code      = model2_code,
                           data      = base_data,
                           constants = base_constants,
                           inits_fn  = inits_model2,
                           monitors  = c("alpha", "beta", "gamma", "sigma_a", "a", "loglik"),
                           label     = label_model2_pois) 


# ------------------------------------------------------
# Model 3: Nested: journals within fields, random slopes
# ------------------------------------------------------
model3_code <- nimbleCode({
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
  list(A       = rnorm(1, 3.5, 0.1),
       beta    = rnorm(1, -0.04, 0.01),
       gamma   = rnorm(1, 0, 0.05),
       sigma_a = runif(1, 0.4, 1.0),
       sigma_b = runif(1, 0.01, 0.08),
       sigma_f = runif(1, 0.4, 1.0),
       alpha_f = rnorm(NF, 3.5, 0.5),
       a       = rnorm(J, 3.5, 0.5),
       b       = rnorm(J, -0.04, 0.03))
}

label_model3 <- paste0("model3_pois_", mcmc_niter %/% 1000, "k")

res3 <- build_and_run_pois(code      = model3_code,
                           data      = data3,
                           constants = constants3,
                           inits_fn  = inits_model3,
                           monitors  = c("A", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f", "alpha_f", "a", "b", "mu", "loglik"),
                           label     = label_model3) 

# ---------------------------------------------------------
# Model 4: Gradual policy effect (policy_year), fixed slope
# ---------------------------------------------------------
base_constants4 <- list(N = nrow(df_sample),
                        J = max(df_sample$journal_idx),
                        log_offset = df_sample$log_offset,
                        policy = df_sample$policy,
                        year_centered = df_sample$year_centered,
                        policy_year = df_sample$year_since_policy,
                        journal = df_sample$journal_idx)

model4_code <- nimbleCode({
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
  list(alpha   = rnorm(1, 3.5, 0.1),
       beta    = rnorm(1, -0.04, 0.01),
       gamma   = rnorm(1, 0, 0.05),
       sigma_a = runif(1, 0.4, 1.0),
       a       = rnorm(base_constants$J, 3.5, 0.5))
}

label_model4 <- paste0("model4_pois_", mcmc_niter %/% 1000, "k")

res4 <- build_and_run_pois(code      = model4_code,
                           data      = base_data,
                           constants = base_constants4,
                           inits_fn  = inits_model4,
                           monitors  = c("alpha", "beta", "gamma", "sigma_a", "a", "loglik"),
                           label     = label_model4) 


# -------------------------------------------------------
# Model 5: Gradual policy effect + journal-specific slope
# -------------------------------------------------------
base_constants5 <- list(N = nrow(df_sample),
                        J = max(df_sample$journal_idx),
                        log_offset = df_sample$log_offset,
                        policy = df_sample$policy,
                        year_centered = df_sample$year_centered,
                        policy_year = df_sample$year_since_policy,
                        journal = df_sample$journal_idx)

model5_code <- nimbleCode({
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
  list(alpha   = rnorm(1, 3.5, 0.1),
       beta    = rnorm(1, -0.04, 0.01),
       gamma   = rnorm(1, 0, 0.05),
       sigma_a = runif(1, 0.4, 1.0),
       sigma_b = runif(1, 0.01, 0.08),
       a       = rnorm(base_constants$J, 3.5, 0.5),
       b       = rnorm(base_constants$J, -0.04, 0.03))
}

label_model5 <- paste0("model5_pois_", mcmc_niter %/% 1000, "k")

res5 <- build_and_run_pois(code      = model5_code,
                           data      = base_data,
                           constants = base_constants5,
                           inits_fn  = inits_model5,
                           monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "a", "b", "loglik"),
                           label     = label_model5) 


# ------------------------------------------------
# Model 6: Gradual policy effect + field hierarchy
# ------------------------------------------------
model6_code <- nimbleCode({
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
  
  list(alpha   = rnorm(1, 3.5, 0.1),
       beta    = rnorm(1, -0.04, 0.01),
       gamma   = rnorm(1, 0, 0.05),
       sigma_a = runif(1, 0.4, 1.0),
       sigma_b = runif(1, 0.01, 0.08),
       sigma_f = runif(1, 0.4, 1.0),
       alpha_f = rnorm(NF, 3.5, 0.5),
       a       = rnorm(J, 3.5, 0.5),
       b       = rnorm(J, -0.04, 0.03))
}

label_model6 <- paste0("model6_pois_", mcmc_niter %/% 1000, "k")

res6 <- build_and_run_pois(code      = model6_code,
                           data      = data6,
                           constants = constants6,
                           inits_fn  = inits_model6,
                           monitors  = c("alpha", "beta", "gamma", "sigma_a", "sigma_b", "sigma_f", "a", "b", "alpha_f", "loglik"),
                           label     = label_model6) 