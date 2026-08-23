# =============================================================================
# Partial Distance Correlation for Mixed-Type Variables
#
# Computes the partial distance correlation (Székely & Rizzo, 2014) between
# two variables while controlling for one or more covariates, using the same
# mixed-type encoding machinery as compute_dcor_matrix(). The permutation
# test is delegated to energy::pdcor.test().
#
# Reference:
#   Székely, G. J., & Rizzo, M. L. (2014). Partial distance correlation with
#   methods for dissimilarities. The Annals of Statistics, 42(6), 2382-2412.
# =============================================================================

#' Partial Distance Correlation with Mixed-Type Variables
#'
#' @description
#' Computes the partial distance correlation between two variables
#' \code{x} and \code{y} after removing the linear and nonlinear dependence
#' shared with a set of control variables \code{z} (e.g. age, sex —
#' confounders). Unlike partial Pearson correlation, partial distance
#' correlation captures any form of dependence, not just linear.
#'
#' The method follows Székely & Rizzo (2014): the double-centered (U-centered)
#' distance matrices of \code{x} and \code{y} are projected away from the
#' U-centered distance matrix of \code{z}, and the distance correlation of the
#' projected matrices is computed. The test is performed by
#' \code{energy::pdcor.test()} (permutation of the partial distance
#' covariance).
#'
#' @param data A data frame containing \code{x}, \code{y}, and all of
#'   \code{z}.
#' @param x,y Single character strings naming the two variables of interest.
#' @param z Character vector naming the control (confounding) variables.
#'   At least one variable is required.
#' @param var_types Optional named character vector specifying variable types
#'   (same codes as \code{\link{compute_dcor_matrix}}):
#'   "continuous", "binary_nominal", "binary_ordinal", "ternary_nominal",
#'   "ternary_ordinal", "nominal", "ordinal". If NULL (default), types are
#'   auto-detected for the involved variables.
#' @param R Integer. Number of permutations for the test (default: 1000).
#'   Set \code{R = 0} to skip the permutation test (the partial distance
#'   correlation estimate is still returned; p-value will be NA).
#' @param verbose Logical. If TRUE, prints progress messages (default TRUE).
#'
#' @return A data frame with one row and columns:
#'   \itemize{
#'     \item \code{X Variable} - name of x
#'     \item \code{Y Variable} - name of y
#'     \item \code{Controlled Variables} - comma-separated names of z
#'     \item \code{Partial Distance Correlation} - the partial dcor estimate
#'     \item \code{P Value} - permutation test p-value (NA if R = 0)
#'     \item \code{Significance} - significance marker
#'     \item \code{N Complete} - number of complete observations used
#'       (listwise deletion across x, y, and all z)
#'   }
#'
#' @details
#' Variable encoding is identical to \code{\link{compute_dcor_matrix}}:
#' continuous variables enter as-is; nominal multi-category variables are
#' one-hot encoded (equidistant categories); ordinal variables are integer
#' encoded (adjacent categories closer). This makes partial distance
#' correlation directly applicable to mixed-type medical data.
#'
#' Missing values are handled by listwise deletion (rows with any missing
#' value among x, y, z are dropped), because the partial correlation is a
#' single joint computation over all variables involved. At least 10 complete
#' observations are required; more are strongly recommended (the U-centering
#' estimators are unstable for small n).
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 300
#' crp <- rnorm(n, 5, 2)
#' dat <- data.frame(
#'   age = rnorm(n, 55, 12),
#'   crp = crp,
#'   esr = 0.6 * crp + rnorm(n),   # correlated with crp
#'   sex = factor(sample(c("Male", "Female"), n, TRUE))
#' )
#'
#' # Partial dcor of crp vs esr, controlling for age and sex
#' compute_pdcor(dat, x = "crp", y = "esr", z = c("age", "sex"), R = 999)
#' }
#'
#' @references
#' Székely, G. J., & Rizzo, M. L. (2014). Partial distance correlation with
#' methods for dissimilarities. \emph{The Annals of Statistics}, 42(6),
#' 2382-2412.
#'
#' @export
compute_pdcor <- function(data, x, y, z, var_types = NULL, R = 1000,
                          verbose = TRUE) {

  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("The 'energy' package is required. Install it with: install.packages('energy')")
  }

  # ---- 1. Input validation ----
  data <- as.data.frame(data)

  if (!is.character(x) || length(x) != 1 ||
      !is.character(y) || length(y) != 1) {
    stop("'x' and 'y' must each be a single variable name (character string).")
  }
  if (!is.character(z) || length(z) < 1) {
    stop("'z' must be a character vector with at least one control variable.")
  }

  all_vars <- c(x, y, z)
  if (anyDuplicated(all_vars) > 0) {
    stop("'x', 'y' and 'z' must not overlap: ",
         paste(all_vars[duplicated(all_vars)], collapse = ", "))
  }
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in 'data': ", paste(missing_vars, collapse = ", "))
  }

  # ---- 2. Detect or use provided variable types (for involved vars only) ----
  if (is.null(var_types)) {
    var_types <- detect_var_types(data[all_vars])
    if (verbose) {
      message("\n========== Auto-detected variable types ==========")
      for (vn in all_vars) {
        message(sprintf("  %-25s : %s", vn, var_types[vn]))
      }
      message("===================================================\n")
    }
  } else {
    missing_types <- setdiff(all_vars, names(var_types))
    if (length(missing_types) > 0) {
      stop("Variable types not specified for: ",
           paste(missing_types, collapse = ", "))
    }
  }

  # ---- 3. Encode variables ----
  encoded <- lapply(all_vars, function(vn) {
    encode_variable(data[[vn]], var_types[vn], vn)
  })
  names(encoded) <- all_vars

  # ---- 4. Listwise deletion across all involved variables ----
  complete_idx <- Reduce(`&`, lapply(encoded, function(e) stats::complete.cases(e)))
  n_complete <- sum(complete_idx)

  if (n_complete < 10) {
    stop(sprintf(
      "Only %d complete observations across x, y and z; at least 10 are required.",
      n_complete))
  }
  if (verbose) {
    message(sprintf("Complete observations used (listwise deletion): %d of %d\n",
                    n_complete, nrow(data)))
  }

  # ---- 5. Build encoded matrices ----
  x_mat <- as.matrix(encoded[[x]][complete_idx, , drop = FALSE])
  y_mat <- as.matrix(encoded[[y]][complete_idx, , drop = FALSE])
  z_mat <- as.matrix(do.call(cbind, lapply(z, function(zn) {
    encoded[[zn]][complete_idx, , drop = FALSE]
  })))

  # ---- 6. Compute partial distance correlation ----
  if (R > 0) {
    test_result <- tryCatch({
      energy::pdcor.test(x_mat, y_mat, z_mat, R = R)
    }, warning = function(w) {
      tryCatch(suppressWarnings(energy::pdcor.test(x_mat, y_mat, z_mat, R = R)),
               error = function(e) NULL)
    }, error = function(e) {
      if (verbose) {
        message("Error in pdcor.test: ", e$message)
      }
      NULL
    })

    if (is.null(test_result)) {
      pdcor_val <- tryCatch(energy::pdcor(x_mat, y_mat, z_mat),
                            error = function(e) NA_real_)
      p_val <- NA_real_
    } else {
      pdcor_val <- as.numeric(test_result$estimate)
      p_val <- as.numeric(test_result$p.value)
    }
  } else {
    pdcor_val <- tryCatch(energy::pdcor(x_mat, y_mat, z_mat),
                          error = function(e) {
                            if (verbose) message("Error in pdcor: ", e$message)
                            NA_real_
                          })
    p_val <- NA_real_
    if (verbose) {
      message("R = 0: permutation test skipped (p-value reported as NA).")
    }
  }

  if (verbose) {
    message(sprintf("Partial dcor(%s, %s | %s) = %.4f",
                    x, y, paste(z, collapse = ", "),
                    ifelse(is.na(pdcor_val), NA, pdcor_val)))
  }

  # ---- 7. Assemble output ----
  out <- data.frame(
    x_variable      = x,
    y_variable      = y,
    controlled      = paste(z, collapse = ", "),
    pdcor           = pdcor_val,
    p_value         = p_val,
    significance    = get_significance(p_val),
    n_complete      = n_complete,
    stringsAsFactors = FALSE
  )
  colnames(out) <- c(
    "X Variable",
    "Y Variable",
    "Controlled Variables",
    "Partial Distance Correlation",
    "P Value",
    "Significance",
    "N Complete"
  )

  out
}
