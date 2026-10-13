

# ----------------
# Helper functions
# ----------------

build_and_run_pois <- function(code, data, constants, inits_fn, monitors, label) {
  start_time <- Sys.time()
  message("Starting: ", label, " at ", format(start_time, "%H:%M:%S"))
  
  dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)
  
  detected_cores <- parallel::detectCores()
  if (is.na(detected_cores)) {
    detected_cores <- 1L
  }
  
  available_cores <- max(1L, detected_cores - 2L)
  n_cores <- min(mcmc_nchains, available_cores)
  
  if (n_cores < mcmc_nchains) {
    message(sprintf("Note: only %d core(s) available; %d chain(s) will queue.", n_cores, mcmc_nchains - n_cores))
  }
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(
    cl,
    varlist = c("code", "data", "constants", "inits_fn", "monitors",
                "mcmc_niter", "mcmc_nburnin", "mcmc_thin", "mcmc_seeds"),
    envir = environment())
  
  if (exists("base_constants", envir = .GlobalEnv, inherits = FALSE)) {
    parallel::clusterExport(cl, varlist = "base_constants", envir = .GlobalEnv)
  }
  
  parallel::clusterEvalQ(cl, {
    library(nimble)
    library(coda)
    NULL
  })
  
  chain_results <- parallel::parLapply(cl, seq_len(mcmc_nchains), function(chain_id) {
    tryCatch({
      inits <- inits_fn()
      
      model <- nimble::nimbleModel(code = code, data = data, constants = constants, inits = inits, check = TRUE)
      
      mcmc_conf <- nimble::configureMCMC(model, monitors = monitors, useConjugacy = TRUE, enableWAIC = TRUE, autoBlock = FALSE)
      
      mcmc_obj <- nimble::buildMCMC(mcmc_conf)
      
      c_model <- nimble::compileNimble(model, showCompilerOutput = FALSE)
      c_mcmc <- nimble::compileNimble(mcmc_obj, project = model, showCompilerOutput = FALSE)
      
      nimble::runMCMC(mcmc = c_mcmc, niter = mcmc_niter, nburnin = mcmc_nburnin, thin = mcmc_thin,
                      nchains = 1, setSeed = mcmc_seeds[chain_id], WAIC = TRUE, summary = FALSE, progressBar = FALSE)
    }, error = function(e) {
      msg <- conditionMessage(e)
      message("====================================")
      message(sprintf("Chain %d FAILED", chain_id))
      message(msg)
      message("====================================")
      NULL
    })
  })
  
  failed <- vapply(chain_results, is.null, logical(1L))
  
  if (any(failed)) {
    warning(sprintf(
      "%d chain(s) failed: chains %s",
      sum(failed),
      paste(which(failed), collapse = ", ")))
  }
  
  chain_results <- chain_results[!failed]
  
  if (length(chain_results) == 0L) {
    stop("All chains failed for model: ", label)
  }
  
  waic_list <- lapply(chain_results, `[[`, "WAIC")
  samples_list <- lapply(chain_results, `[[`, "samples")
  
  mcmc_list <- coda::as.mcmc.list(lapply(samples_list, coda::as.mcmc))
  
  saveRDS(list(samples = samples_list, WAIC = waic_list), here("results", paste0("samples_", label, ".rds")))
  
  saveRDS(mcmc_list, here("results", paste0("mcmc_list_", label, ".rds")))
  
  end_time <- Sys.time()
  elapsed_sec <- as.integer(difftime(end_time, start_time, units = "secs"))
  elapsed_fmt <- sprintf(
    "%02d:%02d:%02d",
    elapsed_sec %/% 3600L,
    (elapsed_sec %% 3600L) %/% 60L,
    elapsed_sec %% 60L)
  
  message(sprintf("Finished: %s | Duration: %s (hh:mm:ss)", label, elapsed_fmt))
  
  invisible(list(samples = list(samples = samples_list, WAIC = waic_list), mcmc_list = mcmc_list))
}

# -------------------------------------------------------
# Summarize key parameters from an mcmc.list,
# print a table, save to CSV, and return invisibly.
# -------------------------------------------------------
summarise_model <- function(mcmc_list, key_params, label) {
  sum_coda <- summary(mcmc_list[, key_params])
  ess_values <- effectiveSize(mcmc_list[, key_params])
  rhat_values <- gelman.diag(mcmc_list[, key_params])$psrf[, "Point est."]
  
  tbl <- data.frame(Parameter = key_params,
                    Mean      = round(sum_coda$statistics[, "Mean"], 3),
                    SD        = round(sum_coda$statistics[, "SD"],   3),
                    `95% CrI` = paste0("(", round(sum_coda$quantiles[, "2.5%"],  3), ", ", round(sum_coda$quantiles[, "97.5%"], 3), ")"),
                    ESS       = round(ess_values,  0),
                    `Rhat`    = round(rhat_values, 3),
                    check.names = FALSE)
  
  print(tbl, digits = 3, row.names = FALSE)
  out_path <- here("results", paste0("summary_", label, ".csv"))
  write.csv(tbl, out_path, row.names = FALSE)
}


# ---------------------------------------
# Shared helper functions for diagnostics
# ---------------------------------------

#' Compute per-chain ACF averaged across chains
#' @param mcmc_list  coda mcmc.list object
#' @param param_name character, parameter column name
#' @param max_lag    integer, maximum lag to compute
#' @return data.frame with columns lag, acf, parameter
#' 
acf_per_param <- function(mcmc_list, param_name, max_lag = 50) {
  chain_acfs <- lapply(mcmc_list, function(chain)
    acf(chain[, param_name], plot = FALSE, lag.max = max_lag)$acf[-1])
  data.frame(lag = 1:max_lag,
             acf = rowMeans(do.call(cbind, chain_acfs)),
             parameter = param_name)
}

#' Save trace + density PDFs for a set of key parameters
save_trace_density <- function(mcmc_list, key_params, param_labels,
                               model_label, nrow_panels) {
  
  chain_cols <- c("#1E90FF", "#9932CC", "#008B8B", "#E05555")
  chain_ltys  <- c(1, 2, 3, 4)
  n_chains    <- nchain(mcmc_list)
  chain_names <- paste("Chain", seq_len(n_chains))
  
  # ── Trace plots ─────────────────────────────────────────────────────────────
  pdf(here("diagnostics", paste0("trace_plots_", model_label, ".pdf")),  width = 10, height = 8)
  par(mfrow = c(nrow_panels, 2), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2.5, 0))
  
  for (param in key_params) {
    traceplot(mcmc_list[, param],
              main = param_labels[param],
              col  = chain_cols[seq_len(n_chains)],
              lwd  = 1.2,
              lty  = chain_ltys[seq_len(n_chains)])
  }
  mtext(paste("Trace plots"), outer = TRUE, cex = 1, line = 0.8)
  dev.off()   
  
  # ── Density plots ────────────────────────────────────────────────────────────
  pdf(here("diagnostics", paste0("density_plots_", model_label, ".pdf")), width = 10, height = 8)
  par(mfrow = c(nrow_panels, 2), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2.5, 0))
  
  for (param in key_params) {
    densplot(mcmc_list[, param],
             main = param_labels[param],
             show.obs = FALSE)
  }
  mtext(paste("Posterior density plots"), outer = TRUE, cex = 1, line = 0.8)
  dev.off()
  
  message("Saved trace + density PDFs for: ", model_label)
}

#' Save ACF plot as PDF
save_acf_plot <- function(mcmc_list, key_params, param_labels, model_label) {
  acf_data <- bind_rows(
    lapply(key_params, function(p) acf_per_param(mcmc_list, p, max_lag = 50))
  )
  acf_data$parameter <- factor(
    acf_data$parameter,
    levels = names(param_labels),
    labels = unname(param_labels)
  )
  n_per_chain <- nrow(as.matrix(mcmc_list[[1]]))
  acf_ci <- 1.96 / sqrt(n_per_chain)
  
  p <- ggplot(acf_data, aes(x = lag, y = acf)) +
    geom_hline(yintercept =  0,       linetype = "solid",  color = "gray50") +
    geom_hline(yintercept =  acf_ci,  linetype = "dashed", color = "gray70") +
    geom_hline(yintercept = -acf_ci,  linetype = "dashed", color = "gray70") +
    geom_line(color  = "#2c7bb6") +
    geom_point(size  = 1.2, color = "#2c7bb6") +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 13) +
    labs(title = paste("ACF"),
         x = "Lag", y = "Autocorrelation") +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "gray95"),
          strip.text       = element_text(face = "bold"))
  
  ggsave(here("diagnostics", paste0("acf_", model_label, ".pdf")), plot = p, width = 10, height = 8, dpi = 320)
  message("Saved ACF plot for: ", model_label)
}


# ---------------------
# dic_computation.R
# ---------------------

# Load shared data / constants
df_sample <- read.csv(here::here("data", "df_sample.csv"))  
dim(df_sample)

df_sample$journal_idx <- as.integer(factor(df_sample$journal))
df_sample$field_idx <- as.integer(factor(df_sample$field))

# for primary outcome
base_data <- list(y = df_sample$n_bar, log_offset = df_sample$log_offset, 
                  policy = df_sample$policy, policy_year = df_sample$year_since_policy,
                  year_centered = df_sample$year_centered, journal = df_sample$journal_idx)

# for secondary outcome
#base_data <- list(y = df_sample$n_info, log_offset = df_sample$log_offset, 
#                  policy = df_sample$policy, policy_year = df_sample$year_since_policy,
#                  year_centered = df_sample$year_centered, journal = df_sample$journal_idx)

base_constants <- list(N = nrow(df_sample), J = max(df_sample$journal_idx))

data3 <- data6 <- c(base_data, list(field = df_sample$field_idx))
constants3 <- constants6 <- c(base_constants, list(F = max(df_sample$field_idx)))

# ----------------------------------------------------------------
# Log-likelihood helpers
#   Each function takes a named vector `params` (one MCMC draw)
#   and returns the scalar total log-likelihood  sum_i log p(y_i|params).
# ----------------------------------------------------------------

## Poisson log-lik — models 1.0, 1.1, 2-pois, 4, 5  (base data)
ll_poisson_base <- function(params, data) {
  J  <- max(data$journal)
  N  <- length(data$y)
  
  alpha   <- if ("alpha" %in% names(params)) params["alpha"] else 0
  beta    <- params["beta"]
  gamma   <- params["gamma"]
  
  # journal random intercepts  a[j]
  a_cols  <- grep("^a\\[", names(params), value = TRUE)
  a       <- if (length(a_cols) > 0) params[a_cols] else rep(0, J)
  
  # journal random slopes  b[j]  (models 1.0, 1.1, 5)
  b_cols  <- grep("^b\\[", names(params), value = TRUE)
  b       <- if (length(b_cols) > 0) params[b_cols] else rep(0, J)
  
  # observation-level overdispersion  eps[i]  (model 1.0)
  eps_cols <- grep("^eps\\[", names(params), value = TRUE)
  eps      <- if (length(eps_cols) > 0) params[eps_cols] else rep(0, N)
  
  j <- data$journal
  log_lambda <- data$log_offset +
    alpha + a[j] +
    (beta  + b[j]) * data$year_centered +
    gamma * data$policy_year +
    eps
  
  sum(dpois(data$y, lambda = exp(log_lambda), log = TRUE))
}


## Poisson log-lik — model 3 (nested: journals within fields)
ll_poisson_nested3 <- function(params, data) {
  J <- max(data$journal)
  A      <- params["A"]
  beta   <- params["beta"]
  gamma  <- params["gamma"]
  
  # field random intercepts  f[field]
  f_cols <- grep("^f\\[", names(params), value = TRUE)
  f      <- if (length(f_cols) > 0) params[f_cols] else rep(0, max(data$field))
  
  # journal-within-field random intercepts  a[j]
  a_cols <- grep("^a\\[", names(params), value = TRUE)
  a      <- if (length(a_cols) > 0) params[a_cols] else rep(0, J)
  
  j  <- data$journal
  fi <- data$field
  log_lambda <- data$log_offset +
    A + f[fi] + a[j] +
    beta  * data$year_centered +
    gamma * data$policy_year
  
  sum(dpois(data$y, lambda = exp(log_lambda), log = TRUE))
}


## Poisson log-lik — model 6 (random intercepts + slopes + field hierarchy)
ll_poisson_nested6 <- function(params, data) {
  J <- max(data$journal)
  alpha  <- params["alpha"]
  beta   <- params["beta"]
  gamma  <- params["gamma"]
  
  f_cols <- grep("^f\\[", names(params), value = TRUE)
  f      <- if (length(f_cols) > 0) params[f_cols] else rep(0, max(data$field))
  
  a_cols <- grep("^a\\[", names(params), value = TRUE)
  a      <- if (length(a_cols) > 0) params[a_cols] else rep(0, J)
  
  b_cols <- grep("^b\\[", names(params), value = TRUE)
  b      <- if (length(b_cols) > 0) params[b_cols] else rep(0, J)
  
  j  <- data$journal
  fi <- data$field
  log_lambda <- data$log_offset +
    alpha + f[fi] + a[j] +
    (beta + b[j]) * data$year_centered +
    gamma * data$policy_year
  
  sum(dpois(data$y, lambda = exp(log_lambda), log = TRUE))
}


# ----------------------------------------------------------------
# DIC = D_bar + pD
#   D_bar  = mean deviance over all posterior samples
#   D_hat  = deviance at posterior mean of parameters
#   pD     = D_bar - D_hat   (effective number of parameters)
#   DIC    = D_bar + pD  =  2*D_bar - D_hat
#
# Arguments:
#   samples_obj   nimble samples object  (list with $samples element,
#                 each element being a matrix [iter x params])
#   label         character label used for the saved .rds / .txt file
#   ll_fun        log-likelihood function  (one of the helpers above)
#   data          named list of data passed to ll_fun
# ----------------------------------------------------------------
compute_dic <- function(samples_obj, label, ll_fun = NULL,  data   = NULL) {
  
  # ---- auto-select ll_fun and data from label ----
  if (is.null(ll_fun)) {
    if (grepl("model1a|model1b|model1[._]0|model1[._]1", label, ignore.case = TRUE)) {
      ll_fun <- ll_poisson_base
      data   <- base_data
    } else if (grepl("model2", label, ignore.case = TRUE)) {
      ll_fun <- ll_poisson_base
      data   <- base_data
    } else if (grepl("model3", label, ignore.case = TRUE)) {
      ll_fun <- ll_poisson_nested3
      data   <- data3
    } else if (grepl("model4|model5", label, ignore.case = TRUE)) {
      ll_fun <- ll_poisson_base
      data   <- base_data
    } else if (grepl("model6", label, ignore.case = TRUE)) {
      ll_fun <- ll_poisson_nested6
      data   <- data6
    } else {
      stop("compute_dic: could not auto-detect ll_fun for label '", label,
           "'. Supply ll_fun and data explicitly.")
    }
  }
  
  # ---- pool all chains into one matrix ----
  all_samples <- do.call(rbind, samples_obj$samples)   # [S x params]
  
  S <- nrow(all_samples)
  
  # ---- per-sample deviance  D(theta_s) = -2 * ll(theta_s) ----
  message("Computing DIC for: ", label, "  (", S, " posterior draws) ...")
  deviances <- apply(all_samples, 1, function(row) {
    ll <- tryCatch(ll_fun(row, data), error = function(e) NA_real_)
    -2 * ll
  })
  
  n_failed <- sum(is.na(deviances))
  if (n_failed > 0)
    warning(n_failed, " draws returned NA log-likelihood for ", label)
  
  D_bar <- mean(deviances, na.rm = TRUE)
  
  # ---- deviance at posterior mean ----
  theta_bar <- colMeans(all_samples)
  D_hat     <- -2 * ll_fun(theta_bar, data)
  
  pD  <- D_bar - D_hat
  DIC <- D_bar + pD          # equivalently: 2*D_bar - D_hat
  
  result <- list(label  = label,
                 D_bar  = round(D_bar, 0),
                 D_hat  = round(D_hat, 0),
                 pD     = round(pD,    0),
                 DIC    = round(DIC,   0))
  
  # ---- save ----
  out_rds <- here("results", paste0("dic_", label, ".rds"))
  saveRDS(result, out_rds)
  
  out_txt <- here("diagnostics", paste0("dic_", label, ".txt"))
  sink(out_txt)
  cat("DIC summary:", label, "\n")
  cat("  D_bar (mean deviance)        :", result$D_bar, "\n")
  cat("  D_hat (deviance at post.mean):", result$D_hat, "\n")
  cat("  pD    (eff. # parameters)    :", result$pD,    "\n")
  cat("  DIC                          :", result$DIC,   "\n")
  sink()
  
  message("  D_bar=", result$D_bar, "  D_hat=", result$D_hat,
          "  pD=", result$pD, "  DIC=", result$DIC)
  invisible(result)
}


# ----------------------------------------------------------------
# Collate all saved DIC results into one comparison table
#   Call this AFTER all compute_dic() calls have run.
# ----------------------------------------------------------------
build_dic_table <- function(labels = NULL) {
  
  if (is.null(labels)) {
    # auto-discover from results/ folder
    rds_files <- list.files(here("results"), pattern = "^dic_.*\\.rds$", full.names = TRUE)
    labels    <- sub("^dic_", "", tools::file_path_sans_ext(basename(rds_files)))
  }
  
  dic_list <- lapply(labels, function(lbl) {
    path <- here("results", paste0("dic_", lbl, ".rds"))
    if (!file.exists(path)) {
      warning("DIC file not found: ", path)
      return(NULL)
    }
    r <- readRDS(path)
    data.frame(Model  = r$label, D_bar  = r$D_bar, D_hat  = r$D_hat,
               pD     = r$pD, DIC    = r$DIC, stringsAsFactors = FALSE)
  })
  
  tbl <- bind_rows(Filter(Negate(is.null), dic_list)) %>%
    arrange(DIC) %>%
    mutate(Delta_DIC = round(DIC - min(DIC, na.rm = TRUE), 3),
           across(where(is.numeric), ~ round(.x, 3)))
  
  print(tbl, row.names = FALSE)
  write_csv(tbl, here("results", "dic_comparison_bar.csv"))
  message("DIC comparison table saved to results/DIC_comparison_bar.csv")
  #write_csv(tbl, here("results", "dic_comparison_info.csv"))
  #message("DIC comparison table saved to results/DIC_comparison_info.csv")
  invisible(tbl)
}
