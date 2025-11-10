# Estimate nuisance rank from whitened residuals

Uses a Marčenko–Pastur bulk-edge threshold together with a
parallel-analysis null (voxel-wise circular shifts) to select the
low-rank nuisance dimension after whitening and optional baseline
removal.

## Usage

``` r
estimate_rank_auto(
  Y,
  X_list,
  X_base = NULL,
  runs = NULL,
  control = list(),
  alpha = 0.05,
  B = 50L,
  r_max = 50L,
  block_shifts = TRUE
)
```

## Arguments

- Y:

  Numeric matrix (T x V) of voxel time series.

- X_list:

  List of design blocks (each T x K). If no design exists, pass a single
  all-zero column to keep the code path identical to the fitter.

- X_base:

  Optional baseline design (T x B) subtracted after whitening.

- runs:

  Optional integer vector of length T with run identifiers.

- control:

  Optional list with entries \`lambda_A_init_ridge\` and \`lambda_base\`
  for the initial ridge fits.

- alpha:

  Family-wise level for the parallel-analysis quantile.

- B:

  Number of circular-shift null draws (set to 0 to disable).

- r_max:

  Maximum number of singular values to probe.

- block_shifts:

  Logical; if TRUE, circular shifts are done voxel-wise.

## Value

List with fields \`r\`, \`r_mp\`, \`r_pa\`, \`sigma2\`, \`lambda_plus\`,
\`ev\`.
