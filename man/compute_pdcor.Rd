\name{compute_pdcor}
\alias{compute_pdcor}
\title{Partial Distance Correlation with Mixed-Type Variables}
\usage{
compute_pdcor(data, x, y, z, var_types = NULL, R = 1000, verbose = TRUE)
}
\arguments{
\item{data}{A data frame containing \code{x}, \code{y}, and all of
\code{z}.}
\item{x, y}{Single character strings naming the two variables of
interest.}
\item{z}{Character vector naming the control (confounding) variables.
At least one variable is required.}
\item{var_types}{Optional named character vector specifying variable
types (same codes as \code{\link{compute_dcor_matrix}}):
\code{"continuous"}, \code{"binary_nominal"}, \code{"binary_ordinal"},
\code{"ternary_nominal"}, \code{"ternary_ordinal"}, \code{"nominal"},
\code{"ordinal"}. If \code{NULL} (default), types are auto-detected for
the involved variables.}
\item{R}{Integer. Number of permutations for the test (default: 1000).
Set \code{R = 0} to skip the permutation test (the partial distance
correlation estimate is still returned; p-value will be \code{NA}).}
\item{verbose}{Logical. If \code{TRUE}, prints progress messages
(default: \code{TRUE}).}
}
\value{
A data frame with one row and columns:
\item{X Variable}{name of x}
\item{Y Variable}{name of y}
\item{Controlled Variables}{comma-separated names of z}
\item{Partial Distance Correlation}{the partial dcor estimate}
\item{P Value}{permutation test p-value (\code{NA} if \code{R = 0})}
\item{Significance}{significance marker}
\item{N Complete}{number of complete observations used (listwise
deletion across x, y, and all z)}
}
\description{
Computes the partial distance correlation between two variables
\code{x} and \code{y} after removing the linear and nonlinear dependence
shared with a set of control variables \code{z} (e.g. age, sex —
confounders). Unlike partial Pearson correlation, partial distance
correlation captures any form of dependence, not just linear.
}
\details{
The method follows Szekely & Rizzo (2014): the U-centered distance
matrices of \code{x} and \code{y} are projected away from the U-centered
distance matrix of \code{z}, and the distance correlation of the
projected matrices is computed. The test is performed by
\code{energy::pdcor.test()} (permutation of the partial distance
covariance).

Variable encoding is identical to
\code{\link{compute_dcor_matrix}}: continuous variables enter as-is;
nominal multi-category variables are one-hot encoded (equidistant
categories); ordinal variables are integer encoded (adjacent categories
closer). This makes partial distance correlation directly applicable to
mixed-type medical data.

Missing values are handled by listwise deletion (rows with any missing
value among x, y, z are dropped), because the partial correlation is a
single joint computation over all variables involved. At least 10
complete observations are required; more are strongly recommended (the
U-centering estimators are unstable for small n).
}
\examples{
\dontrun{
set.seed(42)
n <- 300
crp <- rnorm(n, 5, 2)
dat <- data.frame(
  age = rnorm(n, 55, 12),
  crp = crp,
  esr = 0.6 * crp + rnorm(n),   # correlated with crp
  sex = factor(sample(c("Male", "Female"), n, TRUE))
)

# Partial dcor of crp vs esr, controlling for age and sex
compute_pdcor(dat, x = "crp", y = "esr", z = c("age", "sex"), R = 999)
}
}
\references{
Szekely, G. J., & Rizzo, M. L. (2014). Partial distance correlation with
methods for dissimilarities. \emph{The Annals of Statistics}, 42(6),
2382-2412.
}
\author{
Maintainer: Surface User
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{htest}
