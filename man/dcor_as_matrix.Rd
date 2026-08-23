\name{dcor_as_matrix}
\alias{dcor_as_matrix}
\title{Convert Long-Format Distance Correlation Results to a Matrix}
\usage{
dcor_as_matrix(data, value = c("dcor", "p", "adj_p", "sig", "adj_sig", "n"),
               col_names = NULL, fill_lower = TRUE, diag_value = 1)
}
\arguments{
\item{data}{A data frame in the long format produced by
\code{\link{compute_dcor_matrix}}. Column names are auto-detected
(\code{Column Variable}, \code{Row Variable}, \code{Distance Correlation},
\code{P Value}, \code{Significance}, \code{N Pair}, \code{Adjusted P Value},
\code{Adjusted Significance}); alternatively supply \code{col_names}.}
\item{value}{Character. Which quantity to return as the matrix:
\code{"dcor"} (default), \code{"p"} (raw p-values), \code{"adj_p"}
(adjusted p-values), \code{"sig"} (significance markers), \code{"adj_sig"}
(adjusted significance markers), or \code{"n"} (observations per pair).}
\item{col_names}{Optional named character vector mapping standard roles to
actual column names in \code{data}:
\code{c(col = "...", row = "...", dcor = "...", p = "...", sig = "...",
n = "...", adj_p = "...", adj_sig = "...")}. Only \code{col}, \code{row}
and the column selected via \code{value} are required.}
\item{fill_lower}{Logical. If TRUE (default), the matrix is made
symmetrical by mirroring upper-triangle entries to the lower triangle.}
\item{diag_value}{Value to place on the diagonal for
\code{value = "dcor"} (default 1). Use \code{NA} to leave the diagonal
empty.}
}
\value{
A symmetric matrix with variables as row and column names. The diagonal
contains 1 (for \code{value = "dcor"}) or \code{NA} (other values).
}
\description{
Converts the long-format data frame produced by
\code{\link{compute_dcor_matrix}} into a symmetric matrix. This is useful
for feeding results into \code{corrplot::corrplot()}, \code{heatmap()},
or for pasting a compact matrix into a manuscript table.
}
\details{
If the input contains both (i, j) and (j, i) rows — as in the default
\code{compute_dcor_matrix()} output — duplicated entries are verified to
be identical and only one copy is used. If only one direction is present,
it is mirrored to the other.
}
\examples{
\dontrun{
result <- compute_dcor_matrix(dat, R = 0)
m  <- dcor_as_matrix(result)                 # dcor matrix
pm <- dcor_as_matrix(result, value = "p")    # p-value matrix
nm <- dcor_as_matrix(result, value = "n")    # sample size matrix

# Use with corrplot
# corrplot::corrplot(m, method = "color", p.mat = pm, sig.level = 0.05)
}
}
\author{
Maintainer: Surface User
}
\keyword{multivariate}
\keyword{utilities}
