# hrfdecode Vignettes

This directory contains comprehensive documentation for the `hrfdecode` package, following the Hadley-grade house style with the Albers Minimalist visual system.

## Vignette Overview

### Tutorials (Getting Started)

1. **01-getting-started.Rmd** — Quick tour of the core workflow
   - Basic decoder fitting from fMRI time series
   - Making predictions on new data
   - Understanding trial vs. TR-level predictions
   - Read time: ~8-10 minutes

### How-tos (Task-Oriented Recipes)

2. **02-ar-prewhitening.Rmd** — Handle temporal autocorrelation
   - AR(1), ARMA, and automatic order selection
   - Multi-run data with run-specific AR
   - Impact on prediction accuracy
   - Read time: ~12-15 minutes

3. **03-rmvpa-integration.Rmd** — Cross-validation with rMVPA
   - Creating continuous MVPA designs
   - Searchlight and ROI-based analysis
   - Cross-validation structure for time series
   - Read time: ~10-12 minutes

### Explanations (Concepts & Design)

4. **04-hrf-estimation.Rmd** — Understanding joint HRF learning
   - Why estimate HRF alongside decoder?
   - HRF basis functions (SPM canonical + derivative)
   - Trial aggregation with HRF weighting
   - Read time: ~8-10 minutes

5. **05-weakly-supervised.Rmd** — Algorithm internals
   - Alternating least squares (ALS) optimization
   - Regularization parameters and tuning
   - Convergence diagnostics
   - Soft labels interpretation
   - Read time: ~12-15 minutes

## Visual System

All vignettes use the **Albers Minimalist** theme:

- **albers.css** — Clean, accessible styling with dark mode support
- **theme_albers()** — Matching ggplot2 theme (see `R/theme-albers.R`)
- Colorblind-safe Okabe-Ito palette
- System fonts, generous white space, single accent color

## Building Vignettes

To build all vignettes:

```r
# Install package first
devtools::install()

# Build vignettes
devtools::build_vignettes()

# Or build pkgdown site
pkgdown::build_site()
```

To render a single vignette:

```r
rmarkdown::render("vignettes/01-getting-started.Rmd")
```

## Navigation Structure

The `_pkgdown.yml` file organizes vignettes into logical sections:

- **Tutorials** → Start here for guided learning
- **How-tos** → Task-specific recipes
- **Explanations** → Conceptual deep dives

Each vignette includes:
- Table of contents (depth: 2)
- Stable anchors on all H2 headings
- 3-5 cross-links to related articles
- Session info footer for reproducibility

## Dependencies

Vignettes require these suggested packages:

- **knitr** — Vignette engine
- **rmarkdown** — Document rendering
- **ggplot2** — Plotting
- **viridisLite** — Color scales (fallback)
- **sessioninfo** — Session information

All are listed in `DESCRIPTION` under `Suggests`.

## Chunk Defaults

All vignettes use consistent chunk options:

```r
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.align = "center",
  fig.retina = 2,
  out.width = "100%",
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.asp = 0.618  # Golden ratio
)
set.seed(123)
options(pillar.sigfig = 7, width = 80)
```

## Coverage Matrix

| Feature | Vignettes |
|---------|-----------|
| Basic fitting | 01 |
| Prediction (TR & trial) | 01, 03 |
| AR prewhitening | 02 |
| Multi-run data | 02, 03 |
| rMVPA integration | 03 |
| Searchlight | 03 |
| HRF estimation | 04 |
| HRF basis functions | 04 |
| ALS algorithm | 05 |
| Regularization | 05 |
| Soft labels | 05 |
| Parameter tuning | 05 |

## File Organization

```
vignettes/
├── 01-getting-started.Rmd       [Tutorial]
├── 02-ar-prewhitening.Rmd       [How-to]
├── 03-rmvpa-integration.Rmd     [How-to]
├── 04-hrf-estimation.Rmd        [Explanation]
├── 05-weakly-supervised.Rmd     [Explanation]
├── albers.css                   [Visual theme]
└── README.md                    [This file]
```

## Next Steps

After package installation, vignettes will be available via:

```r
browseVignettes("hrfdecode")
vignette("getting-started", package = "hrfdecode")
```

Or view the pkgdown site:

```r
pkgdown::build_site()
```
