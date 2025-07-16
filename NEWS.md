# oneinfl 1.0.2 (2025-02-12)

Made the `makeXZy()` function external.

# oneinfl 1.0.1 (2025-02-12)

## Bug Fixes

- Resolved an error that necessitated changes to `margins.R`, `dEdq_nb.R`, and `dEdq_pois.R`.
- Changed the spacing in the printed table from the `margins.R` function.

# oneinfl 1.0.0 (2025-02-04)

## Initial CRAN Release

This is the first release of the **oneinfl** package, providing functions for **one-inflated zero-truncated truncated count regression models**. It supports **Poisson** and **Negative Binomial** distributions with estimation and ancillary tools.

### **Features**
- Implements `oneinfl()` for fitting **one-inflated zero-truncated count regression models**.
- Implements `truncreg()` for fitting **zero-truncated count regression models**.
- Provides `margins()` for estimating marginal effects.
- Provides `oneLRT()`, `oneWald()`, and `signifWald()` for various tests.
- Provides `predict()` to compute **expected responses**.
- Provides `summary()` to display model estimates and various diagnostics.
- Provides `oneplot()` to visually compare actual counts and any fitted models.
- Provides random variate generation using a matrix of covariates in `rpp()`, `roipp()`, and `roiztnb()`.
- Compatible with standard R modeling functions (`formula` interface).
- Provides documentation.
- Provides a detailed README.md with an application.

### **Testing & Compatibility**
- Successfully tested on:
  - ✅ **Windows** (`windows-latest`)
  - ✅ **Linux** (`ubuntu-latest`)
  - ✅ **macOS (Intel & ARM64)**
- Passed **`R CMD check`** and **GitHub Actions CI checks**.
- No known issues.