# =============================================================================
# Convert Long-Format dcor Results to Matrix Form
#
# Converts the long-format data frame returned by compute_dcor_matrix() into
# a symmetric matrix (dcor / p-value / significance / n), suitable for
# corrplot, heatmap(), or direct publication tables.
# =============================================================================

#' Convert Long-Format Distance Correlation Results to a Matrix
#'
#' @description
#' Converts the long-format data frame produced by
#' \code{\link{compute_dcor_matrix}} into a symmetric matrix. This is useful
#' for feeding results into \code{corrplot::corrplot()}, \code{heatmap()},
#' or for pasting a compact matrix into a manuscript table.
#'
#' @param data A data frame in the long format produced by
#'   \code{compute_dcor_matrix()}. Column names are auto-detected
#'   (\code{Column Variable}, \code{Row Variable}, \code{Distance
#'   Correlation}, \code{P Value}, \code{Significance}, \code{N Pair},
#'   \code{Adjusted P Value}, \code{Adjusted Significance}); alternatively
#'   supply \code{col_names}.
#' @param value Character. Which quantity to return as the matrix:
#'   \describe{
#'     \item{\code{"dcor"}}{distance correlation coefficients (default)}
#'     \item{\code{"p"}}{raw p-values}
#'     \item{\code{"adj_p"}}{adjusted p-values (only if computed with
#'       \code{p_adjust != "none"})}
#'     \item{\code{"sig"}}{significance markers (\code{***}, \code{**}, ...)}
#'     \item{\code{"adj_sig"}}{adjusted significance markers}
#'     \item{\code{"n"}}{number of complete observations per pair}
#'   }
#' @param col_names Optional named character vector mapping standard roles to
#'   actual column names in \code{data}:
#'   \code{c(col = "...", row = "...", dcor = "...", p = "...", sig = "...",
#'   n = "...", adj_p = "...", adj_sig = "...")}. Only \code{col}, \code{row}
#'   and the column selected via \code{value} are required.
#' @param fill_lower Logical. If TRUE (default), the matrix is made
#'   symmetrical by mirroring upper-triangle entries to the lower triangle.
#'   Set FALSE to keep only the entries explicitly present in \code{data}
#'   (e.g. when \code{include_diag = FALSE} and only unique pairs were
#'   returned... note the full output contains both directions by default).
#' @param diag_value Value to place on the diagonal for \code{value = "dcor"}
#'   (default 1). Use \code{NA} to leave the diagonal empty.
#'
#' @return A symmetric matrix with variables as row and column names. The
#'   diagonal contains 1 (for \code{value = "dcor"}) or \code{NA} (other
#'   values), unless \code{include_diag = FALSE} rows are absent.
#'
#' @details
#' If the input contains both (i, j) and (j, i) rows — as in the default
#' \code{compute_dcor_matrix()} output — duplicated entries are verified to
#' be identical and only one copy is used. If only one direction is present,
#' it is mirrored to the other.
#'
#' @examples
#' \dontrun{
#' result <- compute_dcor_matrix(dat, R = 0)
#' m <- dcor_as_matrix(result)                  # dcor matrix
#' pm <- dcor_as_matrix(result, value = "p")    # p-value matrix
#' nm <- dcor_as_matrix(result, value = "n")    # sample size matrix
#'
#' # Use with corrplot
#' # corrplot::corrplot(m, method = "color", p.mat = pm, sig.level = 0.05)
#' }
#'
#' @export
dcor_as_matrix <- function(data,
                           value = c("dcor", "p", "adj_p", "sig", "adj_sig", "n"),
                           col_names = NULL,
                           fill_lower = TRUE,
                           diag_value = 1) {
  data <- as.data.frame(data)
  if (nrow(data) == 0) stop("'data' contains no rows.")

  value <- match.arg(value)

  # ---- 1. Identify columns ----
  if (is.null(col_names)) {
    find_col <- function(patterns, what, optional = FALSE) {
      for (p in patterns) {
        hit <- grep(p, names(data), ignore.case = TRUE, value = TRUE)
        if (length(hit) > 0) return(hit[1])
      }
      if (optional) return(NA_character_)
      stop("Could not identify the '", what, "' column. ",
           "Provide 'col_names' explicitly.")
    }
    col_names <- c(
      col    = find_col(c("^Column Variable$", "Column[ ._-]?Var", "^col$"), "Column Variable"),
      row    = find_col(c("^Row Variable$", "Row[ ._-]?Var", "^row$"), "Row Variable"),
      dcor   = find_col(c("^Distance Correlation$", "Distance[ ._-]?Cor"), "Distance Correlation"),
      p      = find_col(c("^P Value$", "P[ ._-]?Value"), "P Value", optional = TRUE),
      adj_p  = find_col(c("^Adjusted P Value$", "Adjusted[ ._-]?P"), "Adjusted P Value", optional = TRUE),
      sig    = find_col(c("^Significance$", "^Sig"), "Significance", optional = TRUE),
      adj_sig = find_col(c("^Adjusted Significance$", "Adjusted[ ._-]?Sig"), "Adjusted Significance", optional = TRUE),
      n      = find_col(c("^N Pair$", "^N[ ._-]?Pair", "^n$"), "N Pair", optional = TRUE)
    )
  } else {
    valid_roles <- c("col", "row", "dcor", "p", "adj_p", "sig", "adj_sig", "n")
    bad_roles <- setdiff(names(col_names), valid_roles)
    if (length(bad_roles) > 0) {
      stop("Invalid roles in 'col_names': ", paste(bad_roles, collapse = ", "),
           ". Valid roles: ", paste(valid_roles, collapse = ", "), ".")
    }
    for (r in names(col_names)) {
      if (!is.na(col_names[[r]]) && !col_names[[r]] %in% names(data)) {
        stop("Column '", col_names[[r]], "' (role '", r, "') not found in 'data'.")
      }
    }
  }

  # ---- 2. Extract the requested values ----
  role <- switch(value,
                 dcor = "dcor", p = "p", adj_p = "adj_p",
                 sig = "sig", adj_sig = "adj_sig", n = "n")

  vcol <- col_names[[role]]
  if (is.null(vcol) || is.na(vcol)) {
    stop(sprintf(
      "The '%s' column is not present in 'data' (value = \"%s\"). %s",
      role, value,
      if (value %in% c("adj_p", "adj_sig"))
        "Re-run compute_dcor_matrix() with p_adjust set to a method other than 'none'."
      else "Provide it via 'col_names'."
    ))
  }

  x_ <- as.character(data[[col_names[["col"]]]])
  y_ <- as.character(data[[col_names[["row"]]]])
  v_ <- data[[vcol]]
  if (value %in% c("dcor", "p", "adj_p", "n")) {
    v_ <- suppressWarnings(as.numeric(v_))
  } else {
    v_ <- as.character(v_)
  }

  vars <- unique(c(x_, y_))
  n <- length(vars)
  if (n < 1) stop("No variables found in 'data'.")

  m <- matrix(if (value %in% c("dcor", "p", "adj_p", "n")) NA_real_ else NA_character_,
              nrow = n, ncol = n, dimnames = list(vars, vars))

  # ---- 3. Fill the matrix ----
  for (k in seq_along(x_)) {
    xi <- x_[k]; yi <- y_[k]
    if (xi == yi) {
      # Self-pair rows carry dcor = 1 but p = NA; keep what the data says
      if (!is.na(v_[k])) m[yi, xi] <- v_[k]
      next
    }
    # Check consistency if both directions are present
    existing <- m[yi, xi]
    if (!is.null(existing) && !is.na(existing) && !is.na(v_[k])) {
      if (!isTRUE(all.equal(existing, v_[k], tolerance = 1e-12))) {
        warning(sprintf(
          "Conflicting values for pair (%s, %s): '%s' vs '%s'. Keeping the first.",
          yi, xi, format(existing), format(v_[k])))
        next
      }
    }
    if (is.na(m[yi, xi])) m[yi, xi] <- v_[k]
    if (fill_lower && is.na(m[xi, yi])) m[xi, yi] <- v_[k]
  }

  # ---- 4. Diagonal for dcor matrices ----
  if (value == "dcor" && !is.na(diag_value)) {
    diag(m) <- diag_value
  }

  m
}
