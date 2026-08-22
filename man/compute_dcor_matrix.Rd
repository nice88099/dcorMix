\name{compute_dcor_matrix}
\alias{compute_dcor_matrix}
\title{Compute Pairwise Distance Correlation Matrix with Permutation Tests}
\usage{
compute_dcor_matrix(data, var_types = NULL, R = 1000,
                    verbose = TRUE, include_diag = TRUE)
}
\arguments{
\item{data}{A data frame containing the variables to analyze.}

\item{var_types}{Optional named character vector specifying variable
types. Supported types: \code{"continuous"}, \code{"binary_nominal"},
\code{"binary_ordinal"}, \code{"ternary_nominal"},
\code{"ternary_ordinal"}, \code{"nominal"}, \code{"ordinal"}.
If \code{NULL} (default), types are auto-detected.}

\item{R}{Integer. Number of permutations for the test (default: 1000).}

\item{verbose}{Logical. If \code{TRUE}, prints progress messages
(default: \code{TRUE}).}

\item{include_diag}{Logical. If \code{TRUE}, includes self-pairs
(diagonal, dcor = 1) in the output (default: \code{TRUE}).}
}
\value{
A data frame in long format with 5 columns:
\item{Column Variable}{Name of the column variable}
\item{Row Variable}{Name of the row variable}
\item{Distance Correlation}{Distance correlation coefficient (0 to 1)}
\item{P Value}{Permutation test p-value}
\item{Significance}{Significance marker: \code{"***"} (p < 0.001),
\code{"**"} (p < 0.01), \code{"*"} (p < 0.05), \code{"ns"}
(not significant), or \code{"-"} (self-pair or not computable)}
}
\description{
Calculates pairwise distance correlation coefficients and their
permutation test p-values for all pairs of variables in a data frame.
Handles mixed variable types (continuous, binary/ternary nominal,
binary/ternary ordinal) and missing values.
}
\details{
\strong{Variable type encoding strategies:}

\itemize{
\item \strong{Continuous}: used as-is (Euclidean distance on raw values)
\item \strong{Binary (nominal/ordinal)}: 0/1 encoding (matching distance)
\item \strong{Ternary+ nominal}: one-hot/dummy encoding (all category
pairs equidistant, eliminating spurious ordinal structure)
\item \strong{Ternary+ ordinal}: integer encoding (ordered distance,
preserving ordinal structure: adjacent categories are closer than
non-adjacent)
}

\strong{Missing values} are handled by pairwise deletion: for each pair
of variables, only observations with complete (non-missing) data for
both variables are used in the computation.

\strong{Statistical inference} uses \code{energy::dcor.test} with R
permutations. If \code{dcor.test} is not available in the installed
version, the function falls back to \code{energy::dcov.test}.

\strong{Auto-detection rules} (when \code{var_types = NULL}):
ordered factors are classified as ordinal; unordered factors as nominal;
numeric variables with 2 unique values as binary, 3 unique values as
ternary ordinal (default), and more than 3 as continuous. It is
recommended to manually specify \code{var_types} to ensure correct
classification, especially for 3-level numeric variables where the
nominal/ordinal distinction cannot be inferred automatically.

\strong{Technical documentation}: A detailed methodology document
suitable for writing the Methods section of a scientific manuscript is
included with the package and can be located via
\code{system.file("doc", "technical_documentation.md", package = "dcorMix")}.
}
\examples{
\dontrun{
library(dcorMix)

# Example data with mixed variable types
set.seed(42)
n <- 200
dat <- data.frame(
  age = rnorm(n, 55, 12),
  bmi = rnorm(n, 24, 4),
  sex = factor(sample(c("Male", "Female"), n, TRUE)),
  hypertension = factor(sample(c("No", "Yes"), n, TRUE),
                        levels = c("No", "Yes"), ordered = TRUE),
  blood_type = factor(sample(c("A", "B", "O"), n, TRUE)),
  tumor_stage = factor(sample(c("I", "II", "III"), n, TRUE),
                       levels = c("I", "II", "III"), ordered = TRUE)
)
# Add some missing values
dat$age[sample(1:n, 5)] <- NA
dat$bmi[sample(1:n, 3)] <- NA

# Auto-detect variable types
result <- compute_dcor_matrix(dat, R = 1000)
print(head(result, 20))

# Or manually specify variable types (recommended)
result <- compute_dcor_matrix(dat,
  var_types = c(
    age        = "continuous",
    bmi        = "continuous",
    sex        = "binary_nominal",
    hypertension = "binary_ordinal",
    blood_type = "ternary_nominal",
    tumor_stage = "ternary_ordinal"
  ),
  R = 1000
)

# Access the technical documentation path
tech_doc <- system.file("doc", "technical_documentation.md",
                        package = "dcorMix")
cat("Technical documentation:", tech_doc, "\n")
}
}
\references{
Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and
testing dependence by correlation of distances. \emph{The Annals of
Statistics}, 35(6), 2769-2794.

Szekely, G. J., & Rizzo, M. L. (2009). Brownian distance covariance.
\emph{The Annals of Applied Statistics}, 3(4), 1236-1265.

Rizzo, M. L., & Szekely, G. J. (2024). energy: E-Statistics:
Multivariate Inference via the Energy of Data. R package.
\url{https://CRAN.R-project.org/package=energy}
}
\author{
Maintainer: Surface User
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{htest}
