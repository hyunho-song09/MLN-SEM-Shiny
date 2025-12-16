################################################################################
# MLN-SEM Framework v2.1 - STRONGLY BALANCED
#
# Critical fixes:
# 1. Distinct module names (MX_ vs MZ_)
# 2. Proper multi-IV formula
# 3. Correct n_instruments reporting
#
# IMPROVED: Strongly balanced module construction
# - Aggressive splitting of large modules (max 30%)
# - Initial k = 1.5x larger
# - Iterative split until all modules < 30%
# - Eigengene-based merging for small modules
################################################################################

rm(list = ls())

# Required packages
required_packages <- c("tidyverse", "dplyr", "tidyr", "purrr", "tibble",
                      "AER", "stringr", "Matrix", "glmnet", "ivreg",
                      "sandwich", "lmtest", "lavaan", "igraph", "ggplot2")

new_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if(length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(AER)
  library(glmnet)
  library(ivreg)
  library(sandwich)
  library(lmtest)
  library(lavaan)
  library(igraph)
})

set.seed(42)

################################################################################
# UTILITY FUNCTIONS
################################################################################

clr_transform <- function(A, pseudo = 1e-6) {
  A2 <- as.matrix(A)
  A2[A2 <= 0 | is.na(A2)] <- pseudo
  gm <- exp(rowMeans(log(A2)))
  out <- log(A2 / gm)
  out <- scale(out)
  colnames(out) <- colnames(A)
  return(out)
}

as_num_df <- function(df) {
  df2 <- df
  for(col in names(df2)) {
    if(is.factor(df2[[col]])) {
      df2[[col]] <- as.numeric(df2[[col]])
    } else if(!is.numeric(df2[[col]])) {
      df2[[col]] <- as.numeric(as.character(df2[[col]]))
    }
  }
  return(df2)
}

merge_by_index <- function(pheno, metabolome, microbiome, key = "Index") {
  colnames(pheno) <- make.names(colnames(pheno))
  colnames(metabolome) <- make.names(colnames(metabolome))
  colnames(microbiome) <- make.names(colnames(microbiome))
  
  stopifnot(all(key %in% colnames(pheno), 
                key %in% colnames(metabolome), 
                key %in% colnames(microbiome)))
  
  d <- pheno %>% 
    select(all_of(key), everything()) %>%
    left_join(metabolome, by = key) %>%
    left_join(microbiome, by = key)
  
  return(d)
}

scale_matrix <- function(df, cols) {
  X <- as.matrix(df[, cols, drop = FALSE])
  col_names <- colnames(df)[cols]
  
  for(j in seq_len(ncol(X))) {
    if(anyNA(X[, j])) {
      med <- median(X[, j], na.rm = TRUE)
      X[is.na(X[, j]), j] <- med
    }
  }
  
  X <- scale(X)
  colnames(X) <- col_names
  
  return(X)
}

################################################################################
# MODULE CONSTRUCTION - STRONGLY BALANCED
################################################################################

make_modules <- function(mat, k = NULL, min_size = 10, prefix = "M_") {
  p <- ncol(mat)
  
  if(p < 3) {
    warning("Too few features (<3) for clustering")
    return(list(
      membership = setNames(rep(1, p), colnames(mat)),
      MEs = as.data.frame(rowMeans(mat)),
      dend = NULL
    ))
  }
  
  sds <- apply(mat, 2, sd, na.rm = TRUE)
  keep <- which(sds > 1e-9)
  if(length(keep) < 2) {
    stop("No features with variance after filtering")
  }
  mat <- mat[, keep, drop = FALSE]
  p <- ncol(mat)
  
  cr <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
  cr[is.na(cr)] <- 0
  dist_mat <- as.dist(1 - abs(cr))
  
  hc <- hclust(dist_mat, method = "average")
  
  # Determine k with LARGER initial value to prevent huge modules
  if(is.null(k)) {
    k <- max(3, min(floor(p/min_size * 1.5), 25))  # 1.5x more modules initially
  }
  
  # Maximum size allowed for a single module (30% of total)
  max_size <- floor(p * 0.30)
  
  # Initial clustering
  clust <- cutree(hc, k = k)
  
  # AGGRESSIVE module splitting loop
  max_split_iterations <- 20
  for(split_iter in 1:max_split_iterations) {
    tab <- table(clust)
    large_modules <- names(tab)[tab > max_size]
    
    if(length(large_modules) == 0) break
    
    # Split each large module
    for(lm in large_modules) {
      lm_idx <- which(clust == lm)
      n_features <- length(lm_idx)
      
      if(n_features <= max_size) next
      
      # Calculate how many sub-modules needed
      n_sub <- ceiling(n_features / (max_size * 0.8))
      
      # Sub-clustering
      sub_mat <- mat[, lm_idx, drop = FALSE]
      sub_cr <- cor(sub_mat, use = "pairwise.complete.obs", method = "spearman")
      sub_cr[is.na(sub_cr)] <- 0
      sub_dist <- as.dist(1 - abs(sub_cr))
      sub_hc <- hclust(sub_dist, method = "average")
      sub_clust <- cutree(sub_hc, k = n_sub)
      
      # Assign new module IDs
      max_id <- max(clust)
      new_ids <- max_id + seq_along(unique(sub_clust))
      clust[lm_idx] <- new_ids[sub_clust]
    }
  }
  
  # Now handle small modules - merge based on eigengenes
  max_merge_iterations <- 15
  for(merge_iter in 1:max_merge_iterations) {
    tab <- table(clust)
    small_modules <- names(tab)[tab < min_size]
    
    if(length(small_modules) == 0) break
    
    # Calculate module eigengenes for all modules
    unique_mods <- unique(clust)
    MEs_temp <- matrix(NA, nrow = nrow(mat), ncol = length(unique_mods))
    
    for(i in seq_along(unique_mods)) {
      mod_idx <- which(clust == unique_mods[i])
      if(length(mod_idx) == 1) {
        MEs_temp[, i] <- mat[, mod_idx]
      } else {
        svd_res <- svd(mat[, mod_idx, drop = FALSE])
        MEs_temp[, i] <- svd_res$u[, 1] * svd_res$d[1]
      }
    }
    
    # Merge smallest module with most correlated large module
    large_modules <- setdiff(unique_mods, small_modules)
    
    if(length(large_modules) == 0) {
      # No large modules, relax min_size
      if(min_size > 5) {
        min_size <- floor(min_size * 0.8)
      } else {
        break
      }
      next
    }
    
    # Find smallest module
    smallest <- small_modules[which.min(tab[small_modules])]
    s_idx_in_unique <- which(unique_mods == smallest)
    s_me <- MEs_temp[, s_idx_in_unique]
    
    # Calculate correlations with large modules
    cors <- sapply(large_modules, function(m) {
      m_idx_in_unique <- which(unique_mods == m)
      cor(s_me, MEs_temp[, m_idx_in_unique], use = "pairwise.complete.obs")
    })
    
    # Merge with most correlated
    if(any(!is.na(cors))) {
      target <- large_modules[which.max(abs(cors))]
      clust[clust == smallest] <- target
    } else {
      break
    }
  }
  
  # Renumber modules consecutively
  clust <- match(clust, sort(unique(clust)))
  names(clust) <- colnames(mat)
  
  k_final <- max(clust)
  MEs <- matrix(NA, nrow = nrow(mat), ncol = k_final)
  colnames(MEs) <- paste0(prefix, seq_len(k_final))
  
  for(i in seq_len(k_final)) {
    mod_features <- which(clust == i)
    if(length(mod_features) == 1) {
      MEs[, i] <- mat[, mod_features]
    } else {
      svd_res <- svd(mat[, mod_features, drop = FALSE])
      MEs[, i] <- svd_res$u[, 1] * svd_res$d[1]
    }
  }
  
  return(list(
    membership = clust,
    MEs = as.data.frame(MEs),
    dend = hc
  ))
}

################################################################################
# IV SELECTION
################################################################################

select_instruments <- function(x, Z, covars_df = NULL,
                               alpha = 1.0, max_inst = 5, seed = 1) {
  if(!is.numeric(x)) stop("x must be numeric vector")
  
  Zm <- as_num_df(Z)
  Xmat <- as.matrix(Zm)
  
  x_use <- x
  if(!is.null(covars_df) && ncol(covars_df) > 0) {
    covars_clean <- covars_df
    for(col in names(covars_clean)) {
      if(is.factor(covars_clean[[col]])) {
        covars_clean[[col]] <- as.numeric(covars_clean[[col]])
      }
    }
    
    var_check <- sapply(covars_clean, function(v) sd(v, na.rm = TRUE) > 0)
    if(any(var_check)) {
      covars_clean <- covars_clean[, var_check, drop = FALSE]
      
      ok <- complete.cases(x, covars_clean)
      if(sum(ok) >= 10) {
        tryCatch({
          X_cov <- as.matrix(covars_clean[ok, , drop = FALSE])
          fit_cov <- lm(x[ok] ~ X_cov)
          x_use <- rep(NA, length(x))
          x_use[ok] <- residuals(fit_cov)
        }, error = function(e) {
          warning("Covariate residualization failed")
          x_use <<- x
        })
      }
    }
  }
  
  keep <- which(apply(Xmat, 2, function(v) sd(v, na.rm = TRUE)) > 0)
  if(length(keep) == 0) {
    warning("All microbiome modules have zero variance")
    return(character(0))
  }
  Xmat <- Xmat[, keep, drop = FALSE]
  
  ok <- complete.cases(x_use, Xmat)
  if(sum(ok) < 10) {
    warning(sprintf("Too few complete cases (n=%d)", sum(ok)))
    cr <- as.numeric(cor(Xmat, x, use = "pairwise.complete.obs"))
    names(cr) <- colnames(Xmat)
    ord <- order(abs(cr), decreasing = TRUE)
    sel <- names(cr)[ord][seq_len(min(max_inst, length(cr)))]
    return(sel[!is.na(sel)])
  }
  
  x_ok <- x_use[ok]
  X_ok <- Xmat[ok, , drop = FALSE]
  
  set.seed(seed)
  tryCatch({
    cvfit <- glmnet::cv.glmnet(x = X_ok, y = x_ok, alpha = alpha,
                               nfolds = max(3, min(5, sum(ok) - 1)))
    lam <- cvfit$lambda.1se
    fit <- glmnet::glmnet(x = X_ok, y = x_ok, alpha = alpha, lambda = lam)
    
    b <- as.numeric(fit$beta)
    nm <- rownames(fit$beta)
    sel <- nm[abs(b) > 0]
    
    if(length(sel) == 0) {
      cr <- as.numeric(cor(X_ok, x_ok, use = "pairwise.complete.obs"))
      names(cr) <- colnames(X_ok)
      sel <- names(sort(abs(cr), decreasing = TRUE))
      sel <- sel[seq_len(min(max_inst, length(sel)))]
    } else {
      ord <- order(abs(b[match(sel, nm)]), decreasing = TRUE)
      sel <- sel[ord]
      if(length(sel) > max_inst) sel <- sel[seq_len(max_inst)]
    }
    
    if(length(sel) < 2) {
      cr <- as.numeric(cor(X_ok, x_ok, use = "pairwise.complete.obs"))
      names(cr) <- colnames(X_ok)
      sel <- names(sort(abs(cr), decreasing = TRUE))
      sel <- sel[seq_len(min(max(2, max_inst), length(sel)))]
    }
    
    return(sel[!is.na(sel)])
    
  }, error = function(e) {
    warning(sprintf("Elastic Net failed: %s", e$message))
    cr <- as.numeric(cor(X_ok, x_ok, use = "pairwise.complete.obs"))
    names(cr) <- colnames(X_ok)
    sel <- names(sort(abs(cr), decreasing = TRUE))
    sel <- sel[seq_len(min(max_inst, length(sel)))]
    return(sel[!is.na(sel)])
  })
}

################################################################################
# 2SLS - FIXED TO USE ALL IVs TOGETHER
################################################################################

fit_iv <- function(df, y, x, z_vec, covars = NULL, verbose = FALSE) {
  covs_str <- if(!is.null(covars) && length(covars) > 0) {
    paste("+", paste(covars, collapse = " + "))
  } else {
    ""
  }
  
  ivs_str <- paste(z_vec, collapse = " + ")
  
  # CRITICAL: ALL IVs in one formula
  fml_str <- sprintf("%s ~ %s %s | %s %s", y, x, covs_str, ivs_str, covs_str)
  
  if(verbose) {
    cat(sprintf("\nFormula: %s\n", fml_str))
  }
  
  fml <- as.formula(fml_str)
  
  fit <- tryCatch({
    ivreg::ivreg(fml, data = df)
  }, error = function(e) {
    if(verbose) cat(sprintf("ERROR: %s\n", e$message))
    return(NULL)
  })
  
  if(is.null(fit)) return(NULL)
  
  sm <- tryCatch({
    summary(fit, diagnostics = TRUE)
  }, error = function(e) {
    return(summary(fit, diagnostics = FALSE))
  })
  
  weakF <- tryCatch({
    if(!is.null(sm$diagnostics) && "Weak instruments" %in% rownames(sm$diagnostics)) {
      diag_df <- as.data.frame(sm$diagnostics)
      diag_df["Weak instruments", "statistic"]
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  
  sargP <- tryCatch({
    if(!is.null(sm$diagnostics) && "Sargan" %in% rownames(sm$diagnostics)) {
      diag_df <- as.data.frame(sm$diagnostics)
      diag_df["Sargan", "p-value"]
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  
  hausP <- tryCatch({
    if(!is.null(sm$diagnostics) && "Wu-Hausman" %in% rownames(sm$diagnostics)) {
      diag_df <- as.data.frame(sm$diagnostics)
      diag_df["Wu-Hausman", "p-value"]
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  
  if(verbose) {
    cat(sprintf("  Weak-F: %.2f\n", weakF))
    cat(sprintf("  Sargan p: %.4f\n", sargP))
  }
  
  return(list(
    fit = fit,
    coeftest = coeftest(fit, vcov = vcovHC(fit, type = "HC1")),
    weakF = weakF,
    sargP = sargP,
    hausP = hausP,
    formula = fml_str
  ))
}

################################################################################
# SEM VALIDATION
################################################################################

fit_sem_model <- function(df, mediator, outcome, instruments, covariates = NULL) {
  tryCatch({
    med_on_iv <- paste(mediator, "~", paste(instruments, collapse = " + "))
    out_on_med <- paste(outcome, "~", mediator)
    
    if(!is.null(covariates) && length(covariates) > 0) {
      cov_str <- paste(covariates, collapse = " + ")
      med_on_iv <- paste(med_on_iv, "+", cov_str)
      out_on_med <- paste(out_on_med, "+", cov_str)
    }
    
    model_str <- paste(med_on_iv, out_on_med, sep = "\n")
    
    fit <- lavaan::sem(model_str, data = df, missing = "ml")
    
    list(
      fit = fit,
      summary = summary(fit, fit.measures = TRUE, standardized = TRUE),
      fit_measures = fitMeasures(fit),
      success = TRUE
    )
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

suggest_parameters <- function(n, p_metab, p_micro) {
  if(n >= 100) {
    metab_k <- max(5, min(10, floor(p_metab / 25)))
    micro_k <- max(8, min(15, floor(p_micro / 40)))
    min_size <- 10
    max_inst <- 5
    weak_f_threshold <- 10
  } else if(n >= 50) {
    metab_k <- max(4, min(8, floor(p_metab / 30)))
    micro_k <- max(6, min(12, floor(p_micro / 50)))
    min_size <- 8
    max_inst <- 4
    weak_f_threshold <- 5
  } else {
    metab_k <- max(3, min(5, floor(p_metab / 40)))
    micro_k <- max(5, min(8, floor(p_micro / 60)))
    min_size <- 6
    max_inst <- 3
    weak_f_threshold <- 3
  }
  
  min_size <- min(min_size, floor(p_metab / (metab_k * 1.5)))
  
  return(list(
    metab_k = metab_k,
    micro_k = micro_k,
    min_size = min_size,
    max_inst = max_inst,
    weak_f_threshold = weak_f_threshold,
    n_top_modules = max(3, min(6, floor(metab_k * 0.6)))
  ))
}

################################################################################
# MAIN PIPELINE - COMPLETELY FIXED
################################################################################

run_mln_sem <- function(pheno, metabolome, microbiome, 
                        outcomes, covariates = c(),
                        n_top_modules = 10,
                        metab_k = 20, micro_k = 25, min_size = 10,
                        weak_f_threshold = 10, p_threshold = 0.05, sargan_threshold = 0.05,
                        use_clr = TRUE, alpha = 1.0, max_inst = 5) {
  
  cat("\n====================================\n")
  cat("MLN-SEM v2.1 (COMPLETELY FIXED)\n")
  cat("====================================\n\n")
  
  cat("Step 1: Merging datasets...\n")
  d <- merge_by_index(pheno, metabolome, microbiome)
  cat(sprintf("  Merged: %d samples × %d variables\n", nrow(d), ncol(d)))
  
  cat("\nStep 2: Extracting features...\n")
  metab_cols <- grep("^M_", colnames(d))
  micro_cols <- grep("^s_", colnames(d))
  
  X_metab <- scale_matrix(d, metab_cols)
  
  if(use_clr) {
    cat("  CLR transformation (microbiome)...\n")
    Z_micro <- clr_transform(d[, micro_cols, drop = FALSE])
  } else {
    Z_micro <- scale_matrix(d, micro_cols)
  }
  
  cat(sprintf("  Metabolome: %d × %d\n", nrow(X_metab), ncol(X_metab)))
  cat(sprintf("  Microbiome: %d × %d%s\n", 
              nrow(Z_micro), ncol(Z_micro),
              if(use_clr) " (CLR)" else ""))
  
  cat("\nStep 3: Constructing modules...\n")
  
  # CRITICAL: Use distinct prefixes
  metab_mod <- make_modules(X_metab, k = metab_k, min_size = min_size, prefix = "MX_")
  micro_mod <- make_modules(Z_micro, k = micro_k, min_size = min_size, prefix = "MZ_")
  
  ME_metab <- metab_mod$MEs
  ME_micro <- micro_mod$MEs
  
  cat(sprintf("  Metabolome: %d modules (%s)\n", 
              ncol(ME_metab), paste(head(colnames(ME_metab), 3), collapse = ", ")))
  cat(sprintf("  Microbiome: %d modules (%s)\n", 
              ncol(ME_micro), paste(head(colnames(ME_micro), 3), collapse = ", ")))
  
  cat("\nStep 4: Preparing covariates...\n")
  if(length(covariates) > 0) {
    covars <- d %>%
      select(all_of(covariates)) %>%
      mutate(across(everything(), ~ {
        if(is.factor(.)) {
          if(nlevels(.) == 1) as.numeric(as.character(.)) else as.numeric(.)
        } else {
          as.numeric(.)
        }
      }))
    
    var_check <- sapply(covars, function(v) sd(v, na.rm = TRUE) > 0)
    if(!all(var_check)) {
      removed <- names(covars)[!var_check]
      cat(sprintf("  ⚠ Removing zero-variance: %s\n", paste(removed, collapse = ", ")))
      covars <- covars[, var_check, drop = FALSE]
      covariates <- intersect(covariates, names(covars))
    }
    
    cat(sprintf("  Covariates: %s\n", paste(covariates, collapse = ", ")))
  } else {
    covars <- NULL
    cat("  No covariates\n")
  }
  
  cat("\nStep 5: Screening metabolome modules...\n")
  Y_targets <- intersect(outcomes, colnames(d))
  
  cand_scores <- sapply(colnames(ME_metab), function(me) {
    cors <- sapply(Y_targets, function(y) {
      abs(cor(ME_metab[[me]], d[[y]], use = "pairwise.complete.obs"))
    })
    max(cors, na.rm = TRUE)
  })
  
  cand_modules <- names(sort(cand_scores, decreasing = TRUE))[seq_len(min(n_top_modules, length(cand_scores)))]
  cat(sprintf("  Top %d modules selected\n", length(cand_modules)))
  
  cat("\nStep 6: Running 2SLS IV analysis...\n")
  cat(sprintf("  alpha=%.1f, max_inst=%d\n", alpha, max_inst))
  
  iv_results <- list()
  iv_count <- 0
  
  for(me in cand_modules) {
    cat(sprintf("  Analyzing %s...\n", me))
    
    sel_iv <- try(
      select_instruments(ME_metab[[me]], ME_micro, 
                        covars_df = covars,
                        alpha = alpha,
                        max_inst = max_inst),
      silent = TRUE
    )
    
    if(inherits(sel_iv, "try-error") || length(sel_iv) < 2) {
      cat(sprintf("    ✗ IV selection failed\n"))
      next
    }
    
    cat(sprintf("    ✓ IVs: %s\n", paste(sel_iv, collapse = ", ")))
    
    # CRITICAL: Use ALL IVs together in ONE analysis per outcome
    for(y in Y_targets) {
      df_iv <- cbind(
        d[, y, drop = FALSE],
        ME_metab[, me, drop = FALSE],
        ME_micro[, sel_iv, drop = FALSE],
        if(!is.null(covars)) covars else NULL
      )
      
      out <- try(
        fit_iv(df = df_iv, y = y, x = me,
               z_vec = sel_iv, covars = covariates),
        silent = TRUE
      )
      
      if(inherits(out, "try-error") || is.null(out)) next
      
      beta <- tryCatch(out$coeftest[me, "Estimate"], error = function(e) NA_real_)
      se <- tryCatch(out$coeftest[me, "Std. Error"], error = function(e) NA_real_)
      pval <- tryCatch(out$coeftest[me, "Pr(>|t|)"], error = function(e) NA_real_)
      
      # CRITICAL: Store as single result with all IVs
      iv_results[[paste(me, y, sep = "__")]] <- list(
        module = me,
        y = y,
        IVs = sel_iv,  # Vector, not string
        n_IVs = length(sel_iv),  # CRITICAL
        beta = beta,
        se = se,
        p = pval,
        weakF = out$weakF,
        sargP = out$sargP,
        hausP = out$hausP
      )
      
      iv_count <- iv_count + 1
    }
  }
  
  cat(sprintf("  Completed %d IV analyses\n", iv_count))
  
  cat("\nStep 7: Summarizing...\n")
  
  # Initialize sem_results
  sem_results <- list()
  
  if(length(iv_results) == 0) {
    iv_table <- tibble(
      module = character(0),
      outcome = character(0),
      n_instruments = integer(0),
      instruments = character(0),
      beta = numeric(0),
      se = numeric(0),
      p_value = numeric(0),
      weak_F = numeric(0),
      sargan_p = numeric(0),
      hausman_p = numeric(0)
    )
  } else {
    iv_table <- purrr::map_dfr(iv_results, function(res) {
      tibble(
        module = res$module,
        outcome = res$y,
        n_instruments = res$n_IVs,  # CORRECT
        instruments = paste(res$IVs, collapse = ", "),
        beta = res$beta,
        se = res$se,
        p_value = res$p,
        weak_F = res$weakF,
        sargan_p = res$sargP,
        hausman_p = res$hausP
      )
    })
  }
  
  if(nrow(iv_table) > 0) {
    cat(sprintf("  Total: %d pathways\n", nrow(iv_table)))
    
    has_weakF <- iv_table %>% filter(!is.na(weak_F))
    cat(sprintf("  With weak_F: %d\n", nrow(has_weakF)))
    
    if(nrow(has_weakF) > 0) {
      cat(sprintf("    Range: %.2f - %.2f\n", 
                  min(has_weakF$weak_F), max(has_weakF$weak_F)))
    }
    
    passed <- iv_table %>%
      filter(!is.na(weak_F), weak_F > weak_f_threshold,
             p_value < p_threshold,
             !is.na(sargan_p), sargan_p > sargan_threshold)
    
    cat(sprintf("  Passed all criteria: %d\n", nrow(passed)))
  
  # Step 8: SEM validation
  cat("\nStep 8: SEM validation...\n")
  
  if(nrow(passed) > 0) {
    for(i in seq_len(nrow(passed))) {
      me <- passed$module[i]
      y <- passed$outcome[i]
      ivs <- strsplit(passed$instruments[i], ",\\s*")[[1]]
      
      cat(sprintf("  Validating %s → %s...\n", me, y))
      
      # Build SEM dataframe
      df_sem <- cbind(
        d[, y, drop = FALSE],
        ME_metab[, me, drop = FALSE],
        ME_micro[, ivs, drop = FALSE],
        if(!is.null(covars)) covars else NULL
      )
      
      # Fit SEM
      sem_out <- fit_sem_model(
        df = df_sem,
        mediator = me,
        outcome = y,
        instruments = ivs,
        covariates = covariates
      )
      
      if(sem_out$success) {
        cat(sprintf("    ✓ SEM fit successful\n"))
      } else {
        cat(sprintf("    ✗ SEM fit failed: %s\n", sem_out$error))
      }
      
      sem_results[[paste(me, y, sep = "__")]] <- sem_out
    }
    
    cat(sprintf("  Completed %d SEM validations\n", length(sem_results)))
  } else {
    cat("  No pathways to validate\n")
  }
  } else {
    passed <- iv_table
  }
  
  cat("\n====================================\n")
  cat("Analysis complete!\n")
  cat("====================================\n\n")
  
  return(list(
    merged_data = d,
    metabolome_modules = metab_mod,
    microbiome_modules = micro_mod,
    iv_results = iv_results,
    iv_table = iv_table,
    passed_pathways = passed,
    sem_results = sem_results,  # FIXED: return actual results
    parameters = list(
      use_clr = use_clr,
      alpha = alpha,
      max_inst = max_inst,
      weak_f_threshold = weak_f_threshold
    )
  ))
}

cat("MLN-SEM v2.1 COMPLETELY FIXED loaded!\n")
cat("All bugs fixed:\n")
cat("  ✓ Distinct module names (MX_ vs MZ_)\n")
cat("  ✓ Proper multi-IV formula\n")
cat("  ✓ Correct n_instruments\n\n")
