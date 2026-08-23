# dcorMix

Distance Correlation Analysis for Mixed-Type Variables

`dcorMix` provides `compute_dcor_matrix()`, a ready-to-use tool for computing
pairwise **distance correlation coefficients** (Szekely et al., 2007) between
variables of **mixed types** — continuous, binary/ternary nominal, and
binary/ternary ordinal — with permutation or asymptotic tests, designed for
medical research datasets.

The package now also provides:
- `plot_dcor_matrix()` — ggplot2 heatmap with optional clustered ordering
- `dcor_as_matrix()` — convert results to a matrix for `corrplot` / tables
- `compute_pdcor()` — **partial** distance correlation controlling for
  confounders (e.g. age, sex)

## Features

- Automatic variable-type detection (continuous / binary / ternary,
  nominal / ordinal), with optional manual override
- Statistically principled encoding strategies:
  - nominal (>=3 levels): one-hot encoding (all categories equidistant)
  - ordinal (>=3 levels): integer encoding (adjacent categories closer)
  - binary: 0/1 encoding
- Missing values handled by pairwise deletion, with per-pair **N reported**
  for full transparency
- Two test options:
  - permutation test via `energy::dcor.test` (R = 1000 by default)
  - **asymptotic t-test** via `energy::dcorT.test` (set `R = 0`,
    orders of magnitude faster for large samples)
- **Multiple-comparison correction**: `p_adjust = "BH"`, `"bonferroni"`,
  `"holm"`, `"BY"`, etc. (FDR is the typical choice for medical research)
- **Partial distance correlation** (`compute_pdcor()`) controls for
  confounders via the Szekely-Rizzo 2014 U-centered projection
- `plot_dcor_matrix()`: ggplot2 heatmap with values, significance stars,
  optional p-values, variable renaming, **hierarchical-clustered ordering**
  (`var_order = "cluster"`), automatic white/black text on dark/light tiles
- `dcor_as_matrix()`: convert long-format results to a symmetric matrix
  (dcor / p / sig / adj_p / n) for `corrplot`, `heatmap()` or manuscript
  tables

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
`P Value`, `Significance`, `N Pair`. With `p_adjust` set, two extra columns
are appended: `Adjusted P Value` and `Adjusted Significance`.

### Fast asymptotic test (large samples)

```r
# For large samples, the asymptotic t-test is orders of magnitude faster
# and produces nearly identical p-values to the permutation test
result_fast <- compute_dcor_matrix(dat, R = 0)   # uses dcorT.test
```

### Multiple-comparison correction

```r
# Benjamini-Hochberg FDR — typical default in medical research
result_fdr <- compute_dcor_matrix(dat, R = 0, p_adjust = "BH")

# Bonferroni (conservative, controls family-wise error)
result_bonf <- compute_dcor_matrix(dat, R = 0, p_adjust = "bonferroni")
```

### Partial distance correlation (controlling for confounders)

```r
# Distance correlation of crp vs esr after removing the effect of age and sex
result_partial <- compute_pdcor(dat,
                                x = "crp", y = "esr",
                                z = c("age", "sex"),
                                R = 999)
```

### Converting to a matrix (for `corrplot`, `heatmap()`, or a manuscript table)

```r
m_dcor <- dcor_as_matrix(result_fdr)              # distance correlation matrix
m_p    <- dcor_as_matrix(result_fdr, value = "p") # p-value matrix
m_n    <- dcor_as_matrix(result_fdr, value = "n") # sample-size matrix

# corrplot::corrplot(m_dcor, method = "color", p.mat = m_p, sig.level = 0.05)
```

### Plotting the matrix

```r
# Basic heatmap: tiles colored by dcor, values + significance stars annotated
p <- plot_dcor_matrix(result, title = "Distance Correlation Matrix")
print(p)

# Show p-values instead of significance stars
plot_dcor_matrix(result, show_significance = FALSE, show_p = TRUE)

# Rename variables on the axes and hide the diagonal
plot_dcor_matrix(result,
                 labels = c(age = "Age (years)", bmi = "BMI (kg/m2)"),
                 show_diag = FALSE)

# Custom gradient colors and variable order
plot_dcor_matrix(result,
                 low = "white", high = "#D73027",
                 var_order = c("age", "bmi", "sex", "stage", "blood"))

# Hierarchical clustering (groups strongly correlated variables together)
plot_dcor_matrix(result, var_order = "cluster",
                 title = "Clustered Distance Correlation Matrix")
```

The function returns a standard ggplot object, so you can further customize it
(e.g. `+ ggplot2::theme_minimal(14)`) or save it with `ggplot2::ggsave()`.

## Documentation

A full technical document (in Chinese) covering the mathematical background,
encoding rationale, and permutation-test methodology is included:

```r
cat(system.file("doc", "technical_documentation.md", package = "dcorMix"))
```

Function reference: `?compute_dcor_matrix`, `?plot_dcor_matrix`,
`?dcor_as_matrix` and `?compute_pdcor` after installation.

## References

- Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
  dependence by correlation of distances. *Annals of Statistics*, 35(6), 2769-2814.
- Szekely, G. J., & Rizzo, M. L. (2009). Brownian distance covariance.
  *Annals of Applied Statistics*, 3(4), 1236-1265.
- Szekely, G. J., & Rizzo, M. L. (2014). Partial distance correlation with
  methods for dissimilarities. *Annals of Statistics*, 42(6), 2382-2412.
- Rizzo, M. L., & Szekely, G. J. (2024). _energy: E-Statistics: Multivariate
  Inference via the Energy of Data_. R package.

## License

MIT
