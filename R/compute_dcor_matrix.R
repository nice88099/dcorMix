# =============================================================================
# Distance Correlation Matrix with Permutation Tests
#
# Computes pairwise distance correlation coefficients and permutation test
# p-values for variables of mixed types (continuous, binary/ternary nominal,
# binary/ternary ordinal). Handles missing values via pairwise deletion.
#
# Dependencies: energy package (install.packages("energy"))
# =============================================================================

#' Compute Pairwise Distance Correlation Matrix with Permutation Tests
#'
#' @description
#' Calculates pairwise distance correlation coefficients and their permutation
#' test p-values for all pairs of variables in a data frame. Handles mixed
#' variable types (continuous, binary/ternary nominal, binary/ternary ordinal)
#' and missing values.
#'
#' @param data A data frame containing the variables to analyze.
#' @param var_types Optional named character vector specifying variable types.
#'   Supported types: "continuous", "binary_nominal", "binary_ordinal",
#'   "ternary_nominal", "ternary_ordinal", "nominal", "ordinal".
#'   If NULL (default), types are auto-detected.
#' @param R Integer. Number of permutations for the test (default: 1000).
#'   Set \code{R = 0} to use the fast asymptotic t-test
#'   (\code{energy::dcorT.test}, based on the bias-corrected dcorT) instead
#'   of the permutation test — recommended for large samples (n > 200),
#'   where results are nearly identical but computation is orders of
#'   magnitude faster.
#' @param p_adjust Character. P-value adjustment method for multiple
#'   comparisons (default: "none", no adjustment). Any method supported by
#'   \code{stats::p.adjust} may be used, e.g. "BH" (Benjamini-Hochberg FDR),
#'   "bonferroni", "holm", "hochberg", "BY". When not "none", two extra
#'   columns are appended: \code{Adjusted P Value} and
#'   \code{Adjusted Significance}. Adjustment is applied over the unique
#'   variable pairs (upper triangle) only.
#' @param verbose Logical. If TRUE, prints progress messages (default: TRUE).
#' @param include_diag Logical. If TRUE, includes self-pairs (diagonal,
#'   dcor = 1) in the output (default: TRUE).
#'
#' @return A data frame in long format with 6 columns (8 when
#'   \code{p_adjust != "none"}):
#'   \itemize{
#'     \item \code{Column Variable} - Name of the column variable
#'     \item \code{Row Variable} - Name of the row variable
#'     \item \code{Distance Correlation} - Distance correlation coefficient
#'     \item \code{P Value} - Permutation test p-value (raw)
#'     \item \code{Significance} - Significance marker (based on raw p)
#'     \item \code{N Pair} - Number of complete observations used for this
#'       pair (after pairwise deletion)
#'     \item \code{Adjusted P Value} - (only if \code{p_adjust != "none"})
#'     \item \code{Adjusted Significance} - (only if \code{p_adjust != "none"})
#'   }
#'
#' @details
#' Variable type encoding strategies:
#' \itemize{
#'   \item Continuous: used as-is (Euclidean distance on raw values)
#'   \item Binary (nominal/ordinal): 0/1 encoding (matching distance)
#'   \item Ternary+ nominal: one-hot/dummy encoding (matching distance)
#'   \item Ternary+ ordinal: integer encoding (ordered distance, preserving
#'     ordinal structure: adjacent categories are closer than non-adjacent)
#' }
#'
#' Missing values are handled by pairwise deletion: for each pair of variables,
#' only observations with complete (non-missing) data for both variables are
#' used in the computation. The actual number of observations used is reported
#' in the \code{N Pair} column for full transparency.
#'
#' Hypothesis testing:
#' \itemize{
#'   \item \code{R > 0} (default): permutation test via
#'     \code{energy::dcor.test} with R permutations
#'   \item \code{R = 0}: asymptotic t-test via \code{energy::dcorT.test}
#'     (bias-corrected dcorT), which is much faster for large samples
#' }
#' If \code{dcor.test} is not available in the installed version, the function
#' falls back to \code{energy::dcov.test}.
#'
#' Significance markers:
#' \tabular{ll}{
#'   \code{***} \tab p < 0.001
#'   \code{**}  \tab p < 0.01
#'   \code{*}   \tab p < 0.05
#'   \code{ns}  \tab p >= 0.05
#'   \code{-}   \tab not applicable (self-pair or not computable)
#' }
#'
#' @examples
#' \dontrun{
#' library(energy)
#'
#' # Example data with mixed variable types
#' set.seed(42)
#' n <- 200
#' dat <- data.frame(
#'   age = rnorm(n, 55, 12),
#'   bmi = rnorm(n, 24, 4),
#'   sex = factor(sample(c("Male", "Female"), n, TRUE)),
#'   hypertension = factor(sample(c("No", "Yes"), n, TRUE),
#'                         levels = c("No", "Yes"), ordered = TRUE),
#'   blood_type = factor(sample(c("A", "B", "O"), n, TRUE)),
#'   tumor_stage = factor(sample(c("I", "II", "III"), n, TRUE),
#'                        levels = c("I", "II", "III"), ordered = TRUE)
#' )
#' # Add some missing values
#' dat$age[sample(1:n, 5)] <- NA
#' dat$bmi[sample(1:n, 3)] <- NA
#'
#' result <- compute_dcor_matrix(dat, R = 1000)
#' print(head(result, 20))
#' }
#'
#' @export
compute_dcor_matrix <- function(data, var_types = NULL, R = 1000,
                                p_adjust = "none",
                                verbose = TRUE, include_diag = TRUE) {

  # ---- 1. Input validation ----
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("The 'energy' package is required. Install it with: install.packages('energy')")
  }

  data <- as.data.frame(data)
  var_names <- names(data)
  n_vars <- length(var_names)

  if (n_vars < 2) {
    stop("At least two variables are required.")
  }

  # Validate p_adjust early (fail fast)
  if (!is.character(p_adjust) || length(p_adjust) != 1) {
    stop("'p_adjust' must be a single character string.")
  }
  if (!p_adjust %in% c("none", stats::p.adjust.methods)) {
    stop(sprintf("Invalid p_adjust '%s'. Use 'none' or one of: %s",
                 p_adjust, paste(stats::p.adjust.methods, collapse = ", ")))
  }

  test_label <- if (R > 0) {
    sprintf("%d permutations", R)
  } else {
    "asymptotic t-test"
  }

  # ---- 2. Detect or use provided variable types ----
  if (is.null(var_types)) {
    var_types <- detect_var_types(data)
    if (verbose) {
      message("\n========== Auto-detected variable types ==========")
      for (vn in var_names) {
        message(sprintf("  %-25s : %s", vn, var_types[vn]))
      }
      message("===================================================\n")
    }
  } else {
    # Validate provided types
    missing_vars <- setdiff(var_names, names(var_types))
    if (length(missing_vars) > 0) {
      stop(paste("Variable types not specified for:",
                 paste(missing_vars, collapse = ", ")))
    }
    valid_types <- c("continuous", "binary_nominal", "binary_ordinal",
                     "ternary_nominal", "ternary_ordinal",
                     "nominal", "ordinal", "constant")
    invalid <- setdiff(var_types, valid_types)
    if (length(invalid) > 0) {
      stop(paste("Invalid variable types:",
                 paste(invalid, collapse = ", ")))
    }
    if (verbose) {
      message("\n========== User-specified variable types ==========")
      for (vn in var_names) {
        message(sprintf("  %-25s : %s", vn, var_types[vn]))
      }
      message("=====================================================\n")
    }
  }

  # ---- 3. Encode each variable according to its type ----
  encoded_list <- lapply(var_names, function(vn) {
    encode_variable(data[[vn]], var_types[vn], vn)
  })
  names(encoded_list) <- var_names

  # ---- 4. Identify constant (zero-variance) variables ----
  is_constant <- sapply(var_names, function(vn) {
    enc <- encoded_list[[vn]]
    vals <- as.vector(as.matrix(enc))
    vals <- vals[!is.na(vals)]
    length(unique(vals)) <= 1
  })

  # ---- 5. Compute pairwise distance correlations ----
  # Use symmetry: dcor(X,Y) = dcor(Y,X), and the permutation test
  # p-value is also symmetric. Only compute upper triangle (i < j).

  n_pairs_upper <- n_vars * (n_vars - 1) / 2

  if (verbose) {
    message(sprintf("Computing %d unique pairs (%s each)...\n",
                    n_pairs_upper, test_label))
  }

  dcor_matrix <- matrix(NA_real_, nrow = n_vars, ncol = n_vars,
                        dimnames = list(var_names, var_names))
  pval_matrix <- matrix(NA_real_, nrow = n_vars, ncol = n_vars,
                        dimnames = list(var_names, var_names))
  npair_matrix <- matrix(NA_integer_, nrow = n_vars, ncol = n_vars,
                         dimnames = list(var_names, var_names))

  pair_count <- 0

  for (i in seq_len(n_vars)) {
    for (j in seq(i, n_vars)) {

      if (i == j) {
        # Self-pair: distance correlation of a variable with itself is 1
        dcor_matrix[i, j] <- 1
        pval_matrix[i, j] <- NA_real_
        next
      }

      pair_count <- pair_count + 1

      if (verbose && (pair_count == 1 || pair_count %% 3 == 0 ||
                      pair_count == n_pairs_upper)) {
        message(sprintf("  [%d/%d] %s  vs  %s ...",
                        pair_count, n_pairs_upper,
                        var_names[i], var_names[j]))
      }

      x <- encoded_list[[i]]
      y <- encoded_list[[j]]

      # ---- Pairwise deletion: keep only complete cases ----
      complete_idx <- complete.cases(x, y)
      x_complete <- x[complete_idx, , drop = FALSE]
      y_complete <- y[complete_idx, , drop = FALSE]
      n_complete <- sum(complete_idx)

      # ---- Edge case: constant variable ----
      if (is_constant[i] || is_constant[j]) {
        if (verbose) {
          message(sprintf("      -> Skipped (constant variable): %s or %s",
                          var_names[i], var_names[j]))
        }
        next
      }

      # ---- Edge case: too few complete observations ----
      if (n_complete < 5) {
        if (verbose) {
          message(sprintf("      -> Skipped (only %d complete obs): %s vs %s",
                          n_complete, var_names[i], var_names[j]))
        }
        next
      }

      # ---- Convert to matrix ----
      x_mat <- as.matrix(x_complete)
      y_mat <- as.matrix(y_complete)

      # ---- Check for zero variance after deletion ----
      x_has_var <- any(apply(x_mat, 2, function(col) {
        v <- stats::var(col, na.rm = TRUE)
        !is.na(v) && v > 0
      }))
      y_has_var <- any(apply(y_mat, 2, function(col) {
        v <- stats::var(col, na.rm = TRUE)
        !is.na(v) && v > 0
      }))

      if (!x_has_var || !y_has_var) {
        if (verbose) {
          message(sprintf("      -> Skipped (zero variance after deletion): %s vs %s",
                          var_names[i], var_names[j]))
        }
        next
      }

      # ---- Compute distance correlation (energy::dcor) ----
      dcor_val <- tryCatch({
        result <- energy::dcor(x_mat, y_mat)
        if (is.na(result) || is.nan(result)) NA_real_ else result
      }, error = function(e) {
        if (verbose) {
          message(sprintf("      -> Error computing dcor for %s vs %s: %s",
                          var_names[i], var_names[j], e$message))
        }
        NA_real_
      }, warning = function(w) {
        result <- suppressWarnings(energy::dcor(x_mat, y_mat))
        if (is.na(result) || is.nan(result)) NA_real_ else result
      })

      # ---- Compute hypothesis test ----
      # R > 0: permutation test (dcor.test)
      # R = 0: asymptotic t-test (dcorT.test on bias-corrected dcor)
      #        — much faster for large samples
      p_val <- NA_real_
      if (!is.na(dcor_val)) {
        test_result <- tryCatch({
          if (R > 0) {
            energy::dcor.test(x_mat, y_mat, R = R)
          } else {
            energy::dcorT.test(x_mat, y_mat)
          }
        }, error = function(e1) {
          # Fallback for R > 0: dcov.test (older versions of energy)
          if (R > 0) {
            tryCatch({
              energy::dcov.test(x_mat, y_mat, R = R)
            }, error = function(e2) {
              if (verbose) {
                message(sprintf("      -> Error in hypothesis test for %s vs %s: %s",
                                var_names[i], var_names[j], e2$message))
              }
              NULL
            })
          } else {
            if (verbose) {
              message(sprintf("      -> Error in asymptotic test for %s vs %s: %s",
                              var_names[i], var_names[j], e1$message))
            }
            NULL
          }
        }, warning = function(w) {
          # Retry suppressing warnings
          tryCatch(
            {
              if (R > 0) {
                suppressWarnings(energy::dcor.test(x_mat, y_mat, R = R))
              } else {
                suppressWarnings(energy::dcorT.test(x_mat, y_mat))
              }
            },
            error = function(e) {
              tryCatch(
                {
                  if (R > 0) {
                    suppressWarnings(energy::dcov.test(x_mat, y_mat, R = R))
                  } else {
                    suppressWarnings(energy::dcorT.test(x_mat, y_mat))
                  }
                },
                error = function(e2) NULL
              )
            }
          )
        })

        if (!is.null(test_result) && !is.null(test_result$p.value)) {
          p_val <- test_result$p.value
        }
      }

      # ---- Store symmetric results (dcor, p, and pair sample size) ----
      dcor_matrix[i, j] <- dcor_val
      dcor_matrix[j, i] <- dcor_val
      pval_matrix[i, j] <- p_val
      pval_matrix[j, i] <- p_val
      npair_matrix[i, j] <- n_complete
      npair_matrix[j, i] <- n_complete
    }
  }

  # ---- 6. Adjust p-values for multiple comparisons (optional) ----
  # Raw p-values are preserved in 'P Value'; adjusted values go to a
  # separate matrix so both can be reported side by side.
  adj_pval_matrix <- NULL
  if (p_adjust != "none") {
    # Adjust over unique pairs only (upper triangle, non-NA), then mirror.
    # This avoids counting each pair twice (i != j and j != i rows).
    adj_pval_matrix <- pval_matrix
    ut <- upper.tri(adj_pval_matrix)
    idx_pairs <- which(ut & !is.na(adj_pval_matrix))
    if (length(idx_pairs) > 0) {
      adj_pval_matrix[idx_pairs] <- stats::p.adjust(adj_pval_matrix[idx_pairs],
                                                    method = p_adjust)
      # Mirror adjusted values to the lower triangle
      adj_pval_matrix[lower.tri(adj_pval_matrix)] <-
        t(adj_pval_matrix)[lower.tri(adj_pval_matrix)]
    }
  }

  # ---- 7. Assemble long-format result ----
  result_list <- list()
  idx <- 1

  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {

      if (!include_diag && i == j) next

      dcor_val <- dcor_matrix[i, j]
      p_val <- pval_matrix[i, j]
      sig <- get_significance(p_val, is_self = (i == j))
      n_pair_val <- if (i == j) NA_integer_ else npair_matrix[i, j]

      if (!is.null(adj_pval_matrix)) {
        adj_p_val <- adj_pval_matrix[i, j]
        result_list[[idx]] <- data.frame(
          col_variable  = var_names[j],
          row_variable  = var_names[i],
          dcor          = dcor_val,
          p_value       = p_val,
          significance  = sig,
          n_pair        = n_pair_val,
          adj_p_value   = adj_p_val,
          adj_signif    = get_significance(adj_p_val, is_self = (i == j)),
          stringsAsFactors = FALSE
        )
      } else {
        result_list[[idx]] <- data.frame(
          col_variable  = var_names[j],
          row_variable  = var_names[i],
          dcor          = dcor_val,
          p_value       = p_val,
          significance  = sig,
          n_pair        = n_pair_val,
          stringsAsFactors = FALSE
        )
      }
      idx <- idx + 1
    }
  }

  result_df <- do.call(rbind, result_list)

  # ---- 8. Set professional English column names ----
  if (!is.null(adj_pval_matrix)) {
    colnames(result_df) <- c(
      "Column Variable",
      "Row Variable",
      "Distance Correlation",
      "P Value",
      "Significance",
      "N Pair",
      "Adjusted P Value",
      "Adjusted Significance"
    )
  } else {
    colnames(result_df) <- c(
      "Column Variable",
      "Row Variable",
      "Distance Correlation",
      "P Value",
      "Significance",
      "N Pair"
    )
  }

  if (verbose) {
    message("\n========== Computation Complete ==========")
    n_total <- nrow(result_df)
    n_sig   <- sum(result_df$Significance %in% c("***", "**", "*"),
                   na.rm = TRUE)
    n_na    <- sum(is.na(result_df$`Distance Correlation`))
    message(sprintf("  Total pairs in output : %d", n_total))
    message(sprintf("  Significant pairs      : %d", n_sig))
    message(sprintf("  Non-computable pairs   : %d", n_na))
    if (!is.null(adj_pval_matrix)) {
      n_sig_adj <- sum(result_df$`Adjusted Significance` %in%
                         c("***", "**", "*"), na.rm = TRUE)
      message(sprintf("  Significant (adjusted) : %d  [%s]",
                      n_sig_adj, p_adjust))
    }
    message(sprintf("  Test type              : %s", test_label))
    message("==========================================\n")
  }

  return(result_df)
}


# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Detect Variable Types from Data
#'
#' Auto-detects variable types based on R object class and number of unique
#' values. Ordered factors -> ordinal; unordered factors -> nominal; numeric
#' with <= 2 unique -> binary; 3 unique -> ternary ordinal (default);
#' > 3 unique -> continuous.
#'
#' @param data A data frame.
#' @return A named character vector of variable types.
#' @noRd
detect_var_types <- function(data) {
  var_names <- names(data)
  types <- sapply(var_names, function(vn) {
    v <- data[[vn]]
    v_non_na <- stats::na.omit(v)

    if (length(v_non_na) == 0) return("constant")

    if (is.factor(v)) {
      n_levels <- nlevels(droplevels(v_non_na))
      if (n_levels <= 1) return("constant")
      if (n_levels == 2) {
        if (is.ordered(v)) return("binary_ordinal")
        else               return("binary_nominal")
      }
      if (n_levels == 3) {
        if (is.ordered(v)) return("ternary_ordinal")
        else               return("ternary_nominal")
      }
      if (is.ordered(v)) return("ordinal")
      else               return("nominal")

    } else if (is.numeric(v)) {
      n_unique <- length(unique(v_non_na))
      if (n_unique <= 1) return("constant")
      if (n_unique == 2) return("binary_nominal")   # binary numeric, nominal by default
      if (n_unique == 3) return("ternary_ordinal")   # 3-level numeric, ordinal by default
      return("continuous")

    } else if (is.character(v)) {
      n_unique <- length(unique(v_non_na))
      if (n_unique <= 1) return("constant")
      if (n_unique == 2) return("binary_nominal")
      if (n_unique == 3) return("ternary_nominal")
      return("nominal")

    } else {
      return("nominal")  # fallback
    }
  })
  return(types)
}


#' Encode Variable Based on Type
#'
#' Applies the appropriate encoding for each variable type so that Euclidean
#' distances computed on the encoded data yield meaningful distances:
#'
#' - Continuous: as-is (raw numeric values)
#' - Binary nominal/ordinal: 0/1 encoding (same distance regardless of type)
#' - Ternary+ nominal: one-hot encoding (all category pairs equidistant)
#' - Ternary+ ordinal: integer encoding (adjacent categories closer)
#'
#' @param x A vector (the variable to encode).
#' @param var_type Character string specifying the variable type.
#' @param var_name Character string, the variable name (for one-hot naming).
#' @return A data frame with the encoded variable(s).
#' @noRd
encode_variable <- function(x, var_type, var_name) {

  switch(var_type,

    "continuous" = {
      data.frame(v = as.numeric(x), stringsAsFactors = FALSE)
    },

    "binary_nominal" = {
      # 0/1 encoding: distance = 0 (same) or 1 (different)
      if (is.factor(x)) {
        v <- as.numeric(x) - 1
      } else {
        uniq <- sort(unique(stats::na.omit(x)))
        v <- match(x, uniq) - 1
      }
      data.frame(v = as.numeric(v), stringsAsFactors = FALSE)
    },

    "binary_ordinal" = {
      # For binary, 0/1 encoding preserves order and gives matching distance
      if (is.factor(x)) {
        v <- as.numeric(x) - 1
      } else {
        uniq <- sort(unique(stats::na.omit(x)))
        v <- match(x, uniq) - 1
      }
      data.frame(v = as.numeric(v), stringsAsFactors = FALSE)
    },

    "ternary_nominal" = {
      # One-hot encoding: all category pairs equidistant (sqrt(2) for different)
      if (is.factor(x)) {
        levels_vec <- levels(droplevels(x))
      } else {
        levels_vec <- sort(unique(stats::na.omit(x)))
      }
      mat <- sapply(levels_vec, function(lv) {
        as.numeric(x == lv)  # NA stays NA, which is correct
      })
      if (is.null(dim(mat))) {
        mat <- matrix(mat, ncol = 1)
      }
      colnames(mat) <- paste0(var_name, "_", levels_vec)
      as.data.frame(mat)
    },

    "ternary_ordinal" = {
      # Integer encoding: |level_i - level_j| preserves ordinal structure
      if (is.ordered(x)) {
        v <- as.numeric(x)
      } else if (is.factor(x)) {
        v <- as.numeric(x)
      } else {
        uniq <- sort(unique(stats::na.omit(x)))
        v <- match(x, uniq)
      }
      data.frame(v = as.numeric(v), stringsAsFactors = FALSE)
    },

    "nominal" = {
      # One-hot encoding for general (4+ category) nominal
      if (is.factor(x)) {
        levels_vec <- levels(droplevels(x))
      } else {
        levels_vec <- sort(unique(stats::na.omit(x)))
      }
      mat <- sapply(levels_vec, function(lv) {
        as.numeric(x == lv)
      })
      if (is.null(dim(mat))) {
        mat <- matrix(mat, ncol = 1)
      }
      colnames(mat) <- paste0(var_name, "_", levels_vec)
      as.data.frame(mat)
    },

    "ordinal" = {
      # Integer encoding for general (4+ category) ordinal
      if (is.ordered(x)) {
        v <- as.numeric(x)
      } else if (is.factor(x)) {
        v <- as.numeric(x)
      } else {
        uniq <- sort(unique(stats::na.omit(x)))
        v <- match(x, uniq)
      }
      data.frame(v = as.numeric(v), stringsAsFactors = FALSE)
    },

    "constant" = {
      data.frame(v = as.numeric(x), stringsAsFactors = FALSE)
    },

    # Default fallback: treat as continuous
    data.frame(v = as.numeric(x), stringsAsFactors = FALSE)
  )
}


#' Get Significance Marker from P-value
#'
#' @param p P-value (numeric).
#' @param is_self Logical, whether this is a self-pair (diagonal).
#' @return Character string: "***", "**", "*", "ns", or "-".
#' @noRd
get_significance <- function(p, is_self = FALSE) {
  if (is_self) {
    return("-")
  }
  if (is.na(p)) {
    return("-")
  }
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}
