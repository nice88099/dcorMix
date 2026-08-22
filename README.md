# dcorMix

Distance Correlation Analysis for Mixed-Type Variables

`dcorMix` provides `compute_dcor_matrix()`, a ready-to-use tool for computing
pairwise **distance correlation coefficients** (Szekely et al., 2007) between
variables of **mixed types** — continuous, binary/ternary nominal, and
binary/ternary ordinal — with permutation-test p-values, designed for medical
research datasets.

## Features

- Automatic variable-type detection (continuous / binary / ternary,
  nominal / ordinal), with optional manual override
- Statistically principled encoding strategies:
  - nominal (>=3 levels): one-hot encoding (all categories equidistant)
  - ordinal (>=3 levels): integer encoding (adjacent categories closer)
  - binary: 0/1 encoding
- Missing values handled by pairwise deletion
- Permutation test via `energy::dcor.test` (default R = 1000 replicates)
- Significance markers: `***` (p < 0.001), `**` (p < 0.01), `*` (p < 0.05)
- Returns a tidy long-format data frame (column var, row var, dcor, p, sig)

## Installation

```r
# install.packages("remotes")   # if not already installed
remotes::install_github("<your-github-username>/dcorMix")
```

Requires R >= 4.0.0. The `energy` package is installed automatically as a
dependency.

## Usage

```r
library(dcorMix)

# Example: mixed-type medical data
dat <- data.frame(
  age    = rnorm(200, 55, 12),                      # continuous
  sex    = factor(sample(c("M", "F"), 200, TRUE)),  # binary nominal
  stage  = factor(sample(1:3, 200, TRUE), levels = 1:3, ordered = TRUE),
                                                    # ternary ordinal
  blood  = factor(sample(c("A", "B", "O"), 200, TRUE))  # ternary nominal
)

result <- compute_dcor_matrix(dat, R = 1000)
head(result)

# Manual type specification (recommended for numeric-coded nominal vars)
result2 <- compute_dcor_matrix(
  dat,
  var_types = c(age = "continuous", sex = "binary_nominal",
                stage = "ternary_ordinal", blood = "ternary_nominal"),
  R = 1000
)
```

Output columns: `Column Variable`, `Row Variable`, `Distance Correlation`,
`P Value`, `Significance`.

## Documentation

A full technical document (in Chinese) covering the mathematical background,
encoding rationale, and permutation-test methodology is included:

```r
cat(system.file("doc", "technical_documentation.md", package = "dcorMix"))
```

Function reference: `?compute_dcor_matrix` after installation.

## References

- Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
  dependence by correlation of distances. *Annals of Statistics*, 35(6), 2769-2814.
- Szekely, G. J., & Rizzo, M. L. (2009). Brownian distance covariance.
  *Annals of Applied Statistics*, 3(4), 1236-1265.
- Rizzo, M. L., & Szekely, G. J. (2024). _energy: E-Statistics: Multivariate
  Inference via the Energy of Data_. R package.

## License

MIT
