\name{plot_dcor_matrix}
\alias{plot_dcor_matrix}
\title{Plot a Distance Correlation Matrix Heatmap}
\description{
  Creates a ggplot2 heatmap of pairwise distance correlation coefficients
  from the long-format data frame returned by
  \code{\link{compute_dcor_matrix}}. Tiles are colored by the distance
  correlation, with optional numeric annotations, significance stars, and
  p-values.
}
\usage{
plot_dcor_matrix(data, col_names = NULL, labels = NULL,
  var_order = c("asis", "alphabetical", "cluster"),
  triangle = c("full", "lower", "upper"), show_values = TRUE,
  show_significance = TRUE, show_ns = FALSE, show_p = FALSE,
  show_diag = TRUE, digits = 2, low = "#2166AC", high = "#B2182B",
  fill_limits = c(0, 1), na_fill = "grey90",
  text_colors = c("black", "white"), text_threshold = 0.6,
  text_size = 3.5, legend_name = "Distance\nCorrelation",
  title = NULL, base_size = 11)
}
\arguments{
  \item{data}{A data frame in the long format produced by
    \code{compute_dcor_matrix()}, i.e. with columns \code{Column Variable},
    \code{Row Variable}, \code{Distance Correlation}, \code{P Value} and
    \code{Significance}. Column names are auto-detected; alternatively
    supply \code{col_names}.}
  \item{col_names}{Optional named character vector mapping the standard roles
    to actual column names in \code{data}:
    \code{c(col = "...", row = "...", dcor = "...", p = "...", sig = "...")}.
    Use this when your data frame uses renamed columns. \code{p} and
    \code{sig} are optional; they are simply not displayed when absent.}
  \item{labels}{Optional named character vector for renaming variables on the
    axes, e.g. \code{c(age = "Age (years)", bmi = "BMI")}. Unmentioned
    variables keep their original names.}
  \item{var_order}{Character vector of variable names defining the display
    order of rows/columns. Alternatively \code{"asis"} (default; order of
    first appearance in the input), \code{"alphabetical"}, or
    \code{"cluster"} — hierarchical clustering on the distance
    \code{1 - |dcor|}, so that strongly correlated variables are grouped
    together and correlation blocks become visually apparent.}
  \item{triangle}{Character. Which part of the matrix to display:
    \code{"full"} (default; the complete matrix with both triangles),
    \code{"lower"} (lower-left triangle only), or \code{"upper"}
    (upper-right triangle only). In \code{"lower"}/\code{"upper"} modes the
    diagonal is included when \code{show_diag = TRUE}. Half-filled inputs
    (only one orientation per pair) are handled automatically by mirroring
    each pair onto the displayed side.}
  \item{show_values}{Logical. If TRUE (default), prints the distance
    correlation coefficient in each tile.}
  \item{show_significance}{Logical. If TRUE (default), appends significance
    markers (\code{***}, \code{**}, \code{*}) to the printed values.}
  \item{show_ns}{Logical. If TRUE, also appends \code{"ns"} for
    non-significant pairs (default FALSE).}
  \item{show_p}{Logical. If TRUE, prints the permutation-test p-value in
    parentheses below the coefficient (default FALSE). For p < 0.001 it is
    shown as \code{"p<0.001"}.}
  \item{show_diag}{Logical. If TRUE (default), the diagonal (self-pairs,
    dcor = 1) is drawn; set FALSE for a cleaner look.}
  \item{digits}{Integer. Number of decimal places for printed coefficients
    (default 2).}
  \item{low, high}{Colors for the ends of the fill gradient (defaults:
    \code{"#2166AC"} blue and \code{"#B2182B"} red).}
  \item{fill_limits}{Numeric vector of length 2. Fill scale limits
    (default \code{c(0, 1)}, the theoretical range of distance correlation).}
  \item{na_fill}{Color used for missing (non-computable) tiles
    (default \code{"grey90"}).}
  \item{text_colors}{Vector of two colors: text color for light tiles and for
    dark tiles (default \code{c("black", "white")}).}
  \item{text_threshold}{Numeric. Distance correlations above this value are
    considered "dark tiles" and receive the second text color
    (default 0.6).}
  \item{text_size}{Numeric. Size of the tile annotations (default 3.5).}
  \item{legend_name}{Character. Legend title
    (default \code{"Distance\nCorrelation"}).}
  \item{title}{Optional plot title.}
  \item{base_size}{Numeric. Base font size passed to
    \code{ggplot2::theme_minimal} (default 11).}
}
\value{
  A ggplot object.
}
\details{
  The x axis shows \code{Column Variable} and the y axis shows
  \code{Row Variable}, with the first variable placed at the top-left so the
  figure reads like the printed matrix. Text color automatically switches to
  white on strongly-correlated (dark) tiles. Missing pairs (e.g. excluded
  when \code{include_diag = FALSE} was used in
  \code{compute_dcor_matrix()}) are shown in \code{na_fill} grey.

  When \code{triangle = "lower"} or \code{"upper"}, only the corresponding
  half of the matrix is drawn; each pair is placed on the displayed side
  regardless of its original (row, column) orientation in \code{data}. Use
  this for compact publication figures, e.g.
  \code{plot_dcor_matrix(result, triangle = "lower")}.
}
\examples{
\dontrun{
set.seed(42)
n <- 200
dat <- data.frame(
  age  = rnorm(n, 55, 12),
  bmi  = rnorm(n, 24, 4),
  sex  = factor(sample(c("Male", "Female"), n, TRUE)),
  hypertension = factor(sample(c("No", "Yes"), n, TRUE),
                         levels = c("No", "Yes"), ordered = TRUE)
)

result <- compute_dcor_matrix(dat, R = 1000)
p <- plot_dcor_matrix(result, title = "Distance Correlation Matrix")
print(p)
plot_dcor_matrix(result, show_significance = FALSE, show_p = TRUE)
plot_dcor_matrix(result, labels = c(age = "Age (years)", bmi = "BMI"))
plot_dcor_matrix(result, triangle = "lower")
plot_dcor_matrix(result, triangle = "upper", show_diag = FALSE)
}
}
\seealso{
  \code{\link{compute_dcor_matrix}}
}
