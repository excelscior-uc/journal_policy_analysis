
## ===============
## diagnostics.R
## ===============

library(coda)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(fs)
library(stringr)
library(nimble)

source("scripts/99_functions.R")
#source("3_models.R")


# ---------
# Model 1A
# ---------
label_model1a_pois <- paste0("model1a_pois_", mcmc_niter %/% 1000, "k")

mcmc_list1a <- readRDS(here("results", paste0("mcmc_list_", label_model1a_pois, ".rds")))
samples_1a  <- readRDS(here("results", paste0("samples_", label_model1a_pois, ".rds")))
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
label_model1b_pois <- paste0("model1b_pois_", mcmc_niter %/% 1000, "k")

mcmc_list1b <- readRDS(here("results", paste0("mcmc_list_", label_model1b_pois, ".rds")))
samples_1b  <- readRDS(here("results", paste0("samples_", label_model1b_pois, ".rds")))
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
label_model2_pois <- paste0("model2_pois_", mcmc_niter %/% 1000, "k")

mcmc_list2 <- readRDS(here("results", paste0("mcmc_list_", label_model2_pois, ".rds")))
samples_2  <- readRDS(here("results", paste0("samples_", label_model2_pois, ".rds")))
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
label_model3_pois <- paste0("model3_pois_", mcmc_niter %/% 1000, "k")

mcmc_list3 <- readRDS(here("results", paste0("mcmc_list_", label_model3_pois, ".rds")))
samples_3  <- readRDS(here("results", paste0("samples_", label_model3_pois, ".rds")))
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
label_model4_pois <- paste0("model4_pois_", mcmc_niter %/% 1000, "k")

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
label_model5_pois <- paste0("model5_pois_", mcmc_niter %/% 1000, "k")

mcmc_list5 <- readRDS(here("results", paste0("mcmc_list_", label_model5_pois, ".rds")))
samples_5  <- readRDS(here("results", paste0("samples_", label_model5_pois, ".rds")))
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
label_model6_pois <- paste0("model6_pois_", mcmc_niter %/% 1000, "k")

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
