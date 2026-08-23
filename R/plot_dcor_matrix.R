# =============================================================================
# Distance Correlation Matrix Plot (ggplot2)
#
# Visualizes the long-format output of compute_dcor_matrix() as a color-coded
# matrix heatmap with optional correlation values, significance markers, and
# p-values.
#
# Dependencies: ggplot2 package (install.packages("ggplot2"))
# =============================================================================

# Avoid R CMD check "no visible binding" notes for ggplot2 aesthetic columns
utils::globalVariables(c("x_", "y_", "fill_", "label_text", "dark_text"))

#' Plot a Distance Correlation Matrix Heatmap
#'
#' @description
#' Creates a ggplot2 heatmap of pairwise distance correlation coefficients from
#' the long-format data frame returned by \code{\link{compute_dcor_matrix}}.
#' Tiles are colored by the distance correlation, with optional numeric
#' annotations, significance stars, and p-values.
#'
#' @param data A data frame in the long format produced by
#'   \code{compute_dcor_matrix()}, i.e. with columns \code{Column Variable},
#'   \code{Row Variable}, \code{Distance Correlation}, \code{P Value} and
#'   \code{Significance}. Column names are auto-detected; alternatively
#'   supply \code{col_names}.
#' @param col_names Optional named character vector mapping the standard roles
#'   to actual column names in \code{data}:
#'   \code{c(col = "...", row = "...", dcor = "...", p = "...", sig = "...")}.
#'   Use this when your data frame uses renamed columns. \code{p} and \code{sig}
#'   are optional; they are simply not displayed when absent.
#' @param labels Optional named character vector for renaming variables on the
#'   axes, e.g. \code{c(age = "Age (years)", bmi = "BMI")}. Unmentioned
#'   variables keep their original names.
#' @param var_order Character vector of variable names defining the display
#'   order of rows/columns. Alternatively \code{"asis"} (default; order of
#'   first appearance in the input), \code{"alphabetical"}, or
#'   \code{"cluster"} — hierarchical clustering on the distance
#'   \code{1 - |dcor|}, so that strongly correlated variables are grouped
#'   together and correlation blocks become visually apparent.
#' @param show_values Logical. If TRUE (default), prints the distance
#'   correlation coefficient in each tile.
#' @param show_significance Logical. If TRUE (default), appends significance
#'   markers (\code{***}, \code{**}, \code{*}) to the printed values.
#' @param show_ns Logical. If TRUE, also appends \code{"ns"} for non-significant
#'   pairs (default FALSE).
#' @param show_p Logical. If TRUE, prints the permutation-test p-value in
#'   parentheses below the coefficient (default FALSE). For p < 0.001 it is
#'   shown as \code{"p<0.001"}.
#' @param show_diag Logical. If TRUE (default), the diagonal (self-pairs,
#'   dcor = 1) is drawn; set FALSE for a cleaner look.
#' @param digits Integer. Number of decimal places for printed coefficients
#'   (default 2).
#' @param low,high Colors for the ends of the fill gradient
#'   (defaults: \code{"#2166AC"} blue and \code{"#B2182B"} red).
#' @param fill_limits Numeric vector of length 2. Fill scale limits
#'   (default \code{c(0, 1)}, the theoretical range of distance correlation).
#' @param na_fill Color used for missing (non-computable) tiles
#'   (default \code{"grey90"}).
#' @param text_colors Vector of two colors: text color for light tiles and for
#'   dark tiles (default \code{c("black", "white")}).
#' @param text_threshold Numeric. Distance correlations above this value are
#'   considered "dark tiles" and receive the second text color
#'   (default 0.6).
#' @param text_size Numeric. Size of the tile annotations (default 3.5).
#' @param legend_name Character. Legend title (default \code{"Distance\nCorrelation"}).
#' @param title Optional plot title.
#' @param base_size Numeric. Base font size passed to
#'   \code{ggplot2::theme_minimal} (default 11).
#'
#' @return A ggplot object.
#'
#' @details
#' The x axis shows \code{Column Variable} and the y axis shows
#' \code{Row Variable}, with the first variable placed at the top-left so the
#' figure reads like the printed matrix. Text color automatically switches to
#' white on strongly-correlated (dark) tiles. Missing pairs (e.g. excluded when
#' \code{include_diag = FALSE} was used in \code{compute_dcor_matrix()}) are
#' shown in \code{na_fill} grey.
#'
#' When \code{var_order = "cluster"}, variables are ordered by hierarchical
#' clustering (average linkage) on the distance \code{1 - |dcor|}. Missing
#' pairwise dcors are imputed with the median observed dcor for the clustering
#' step only (they remain grey in the plot itself).
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 200
#' dat <- data.frame(
#'   age  = rnorm(n, 55, 12),
#'   bmi  = rnorm(n, 24, 4),
#'   sex  = factor(sample(c("Male", "Female"), n, TRUE)),
#'   hypertension = factor(sample(c("No", "Yes"), n, TRUE),
#'                          levels = c("No", "Yes"), ordered = TRUE)
#' )
#'
#' result <- compute_dcor_matrix(dat, R = 1000)
#' p <- plot_dcor_matrix(result, title = "Distance Correlation Matrix")
#' print(p)
#' # Show p-values instead of stars
#' plot_dcor_matrix(result, show_significance = FALSE, show_p = TRUE)
#' # Rename axis variables
#' plot_dcor_matrix(result, labels = c(age = "Age (years)", bmi = "BMI"))
#' }
#'
#' @export
plot_dcor_matrix <- function(data,
                             col_names = NULL,
                             labels = NULL,
                             var_order = c("asis", "alphabetical", "cluster"),
                             show_values = TRUE,
                             show_significance = TRUE,
                             show_ns = FALSE,
                             show_p = FALSE,
                             show_diag = TRUE,
                             digits = 2,
                             low = "#2166AC",
                             high = "#B2182B",
                             fill_limits = c(0, 1),
                             na_fill = "grey90",
                             text_colors = c("black", "white"),
                             text_threshold = 0.6,
                             text_size = 3.5,
                             legend_name = "Distance\nCorrelation",
                             title = NULL,
                             base_size = 11) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required. Install it with: install.packages('ggplot2')")
  }

  data <- as.data.frame(data)
  if (nrow(data) == 0) {
    stop("'data' contains no rows.")
  }

  # ---- 1. Identify the required columns ----
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
      col  = find_col(c("^Column Variable$", "Column[ ._-]?Var", "^col$"), "Column Variable"),
      row  = find_col(c("^Row Variable$", "Row[ ._-]?Var", "^row$"), "Row Variable"),
      dcor = find_col(c("^Distance Correlation$", "Distance[ ._-]?Cor", "dcor", "value"), "Distance Correlation"),
      p    = find_col(c("^P Value$", "P[ ._-]?Value", "^p$", "pval"), "P Value", optional = TRUE),
      sig  = find_col(c("^Significance$", "Sig"), "Significance", optional = TRUE)
    )
  } else {
    valid_roles <- c("col", "row", "dcor", "p", "sig")
    bad_roles <- setdiff(names(col_names), valid_roles)
    if (length(bad_roles) > 0) {
      stop("Invalid roles in 'col_names': ", paste(bad_roles, collapse = ", "),
           ". Valid roles: ", paste(valid_roles, collapse = ", "), ".")
    }
    missing_cols <- setdiff(col_names[intersect(names(col_names), c("col", "row", "dcor"))],
                            names(data))
    if (length(missing_cols) > 0) {
      stop("Columns not found in 'data': ", paste(missing_cols, collapse = ", "))
    }
    if (!"p" %in% names(col_names)) col_names <- c(col_names, p = NA_character_)
    if (!"sig" %in% names(col_names)) col_names <- c(col_names, sig = NA_character_)
  }

  # ---- 2. Extract into standard columns ----
  df <- data.frame(
    x_ = as.character(data[[col_names[["col"]]]]),
    y_ = as.character(data[[col_names[["row"]]]]),
    fill_ = suppressWarnings(as.numeric(data[[col_names[["dcor"]]]])),
    stringsAsFactors = FALSE
  )
  p_col <- if (!is.na(col_names[["p"]])) col_names[["p"]] else NA_character_
  s_col <- if (!is.na(col_names[["sig"]])) col_names[["sig"]] else NA_character_
  df$p_ <- if (!is.na(p_col)) suppressWarnings(as.numeric(data[[p_col]])) else NA_real_
  df$s_ <- if (!is.na(s_col)) as.character(data[[s_col]]) else NA_character_

  # ---- 3. Rename variables if requested ----
  if (!is.null(labels)) {
    unknown <- setdiff(names(labels), unique(c(df$x_, df$y_)))
    if (length(unknown) > 0) {
      warning("labels for unknown variables ignored: ",
              paste(unknown, collapse = ", "))
    }
    map <- stats::setNames(names(labels), unname(labels))  # new -> old
    for (v in names(map)) {
      df$x_[df$x_ == map[[v]]] <- v
      df$y_[df$y_ == map[[v]]] <- v
    }
  }

  # ---- 4. Variable order ----
  if (is.character(var_order) &&
      all(var_order %in% c("asis", "alphabetical", "cluster"))) {
    vars <- unique(c(df$x_, df$y_))
    if (var_order[1] == "alphabetical") {
      vars <- sort(vars)
    } else if (var_order[1] == "cluster") {
      vars <- cluster_variable_order(df, vars)
    }
  } else {
    vars <- var_order
    unknown <- setdiff(unique(c(df$x_, df$y_)), vars)
    if (length(unknown) > 0) {
      stop("var_order is missing variables: ", paste(unknown, collapse = ", "))
    }
  }
  df$x_ <- factor(df$x_, levels = vars)
  df$y_ <- factor(df$y_, levels = rev(vars))  # first variable at the top

  # ---- 5. Filter diagonal ----
  if (!show_diag) {
    keep <- as.character(df$x_) != as.character(df$y_)
    df <- df[keep, , drop = FALSE]
  }

  # ---- 6. Build tile annotation text ----
  sig_marks <- c("***", "**", "*")
  txt <- character(nrow(df))

  for (i in seq_len(nrow(df))) {
    parts <- character(0)
    if (show_values && !is.na(df$fill_[i])) {
      parts <- c(parts, formatC(df$fill_[i], format = "f", digits = digits))
    }
    if (show_significance && !is.na(df$s_[i])) {
      if (df$s_[i] %in% sig_marks) {
        parts <- c(parts, df$s_[i])
      } else if (show_ns && df$s_[i] == "ns") {
        parts <- c(parts, "ns")
      }
    }
    line1 <- paste(parts, collapse = "")
    if (show_p && !is.na(df$p_[i])) {
      p_txt <- if (df$p_[i] < 0.001) "p<0.001" else
        paste0("p=", formatC(df$p_[i], format = "f", digits = 3))
      line1 <- if (line1 == "") p_txt else paste0(line1, "\n(", p_txt, ")")
    }
    txt[i] <- line1
  }
  df$label_text <- txt
  df$dark_text <- !is.na(df$fill_) & df$fill_ > text_threshold

  # ---- 7. Assemble the ggplot ----
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_, y = y_, fill = fill_)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradient(
      low = low, high = high, limits = fill_limits, na.value = na_fill,
      name = legend_name
    )

  if (any(nzchar(df$label_text))) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = label_text, color = dark_text),
        size = text_size, lineheight = 0.9, na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = stats::setNames(text_colors, c(FALSE, TRUE)), guide = "none"
      )
  }

  p <- p +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1,
        margin = ggplot2::margin(t = 2, unit = "pt")
      ),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.title.position = "plot"
    )

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  p
}


# =============================================================================
# Internal: cluster-based variable ordering
# =============================================================================

#' Order Variables by Hierarchical Clustering of dcors
#'
#' Builds a symmetric dcor matrix from the (possibly half-filled) long-format
#' pairs, converts to the distance 1 - |dcor|, and returns variables ordered
#' by average-linkage hierarchical clustering. Missing dcors are imputed with
#' the median observed dcor for the clustering step only.
#'
#' @param df Data frame with columns x_, y_, fill_ (standardized by the
#'   caller).
#' @param vars Character vector of all variable names.
#' @return Character vector of variable names in clustered order.
#' @noRd
cluster_variable_order <- function(df, vars) {
  n <- length(vars)
  if (n < 3) return(vars)  # nothing to cluster

  m <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(vars, vars))

  for (k in seq_len(nrow(df))) {
    xi <- df$x_[k]; yi <- df$y_[k]; vi <- df$fill_[k]
    if (is.na(vi)) next
    m[yi, xi] <- vi
    m[xi, yi] <- vi
  }

  # Impute missing pairwise dcors with the median observed value
  offdiag <- m[upper.tri(m)]
  fill_val <- if (all(is.na(offdiag))) 0 else stats::median(offdiag, na.rm = TRUE)
  m[is.na(m)] <- fill_val
  diag(m) <- 1

  # Distance between variables: 1 - |dcor|
  d <- stats::as.dist(1 - abs(m))

  hc <- tryCatch(
    stats::hclust(d, method = "average"),
    error = function(e) NULL
  )
  if (is.null(hc)) return(vars)

  vars[hc$order]
}
