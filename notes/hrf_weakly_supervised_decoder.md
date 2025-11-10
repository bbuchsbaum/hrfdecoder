# HRF-Aware Weakly Supervised Decoder — Notes

## Section 1: Joint estimation framework

- Goal: directly model TR-level class mixtures by jointly estimating (i) decoder weights `W`, (ii) latent soft labels `z(t)`, and (iii) a compact HRF `h`, removing single-trial beta estimation while keeping standard MVPA train/test workflow (train on runs, evaluate on held-out runs).
- Treat each TR as a mixture of classes with smooth, nonnegative weights encouraged to match HRF-convolved designs.
- Variables per ROI/searchlight:
  - TR index `t=1..T`, voxel vector `y_t ∈ R^V` (or PCA/ICA reduced).
  - Classes `c=1..C` plus optional background `c=0`.
  - Event trains `s_c(t)=Σ_i δ(t-τ_{c,i})` on TR grid.
  - Latent labels `z_t ∈ Δ^C`, decoder probabilities `q_t = softmax(W y_t + b)`.
  - HRF `h ∈ R^K`, nonnegative, normalized; design prior per class `g_c = s_c * h`, per-TR prior `π_t = softmax([g_{t,1},…,g_{t,C}, β₀])`.

### Objective

Minimize over `W, z, h`:

```
L = Σ_t CE(z_t, q_t)
    + λ ||D² z||_F²         # smoothness on z
    + ρ Σ_t ||z_t - π_t||₂² # stay close to HRF-convolved design
    + α ||W||_F²            # decoder regularization
    + β ||L h||₂²           # HRF roughness
```

Subject to `z_t ∈ Δ^C`, `h ≥ 0`, `Σ_k h_k = 1`. `D²` is second-difference operator on time; `L` penalizes HRF roughness. Hyperparameters (`λ, ρ, α, β`) control smoothness, prior adherence, and regularization.

### Intuition

- Cross-entropy aligns decoder outputs with latent soft labels.
- `ρ` term ties labels to HRF-convolved design prior; `λ` enforces temporal smoothness.
- HRF penalty maintains physiological plausibility.
- Extreme cases:
  - `ρ → ∞`, large `λ`: `z` ≈ `π`, equivalent to canonical GLM-style soft labels.
  - Small `ρ`: data-driven corrections to overlap/timing, with smoothness still enforced.
  - Fixing `z = π`, skipping Z/H steps yields a fast baseline (logistic regression on fixed soft labels).

### Alternating optimization

1. **W-step (convex)**: weighted multinomial logistic regression using fractional targets `z_t`. Use LBFGS/coordinate descent; PCA keeps it fast.
2. **Z-step**: solve banded QP/ADMM enforcing simplex constraints:
   - Minimize `Σ_t -Σ_c z_{t,c} log q_{t,c} + ρ ||z_t - π_t||₂² + λ ||D² z||_F²`.
   - Complexity near `O(TC)`; banded due to `D²`.
3. **H-step**: nonnegative ridge regression to fit HRF:
   - `min_{h ≥ 0, 1ᵀh = 1} Σ_c ||z_{:,c} - s_c * h||₂² + β ||L h||₂²`.
   - Tiny NNLS problem per ROI/searchlight; optional low-dimensional HRF basis (canonical + derivatives or FLOBS) for simpler QP.

Repeat 5–10 iterations (often 2–3 suffice) until loss stabilizes.

### Training/testing workflow

Training:
- Preprocess runs (motion, high-pass, nuisance). For each ROI/searchlight, reduce dimensionality (~30–100 PCs).
- Initialize `h` with canonical HRF, normalize; compute `π` via convolution; set `z ← π`.
- Alternate W → Z → H steps. Use validation fold to tune `λ, ρ, α, β` (defaults: `λ` giving ~1–2 s smoothness; `ρ` so `z` follows design except where data evidence; `β` small but >0).

Testing:
- Freeze `W` and `h`. For TR-level decoding, output `q_t = softmax(W y_t)`.
- For trial-wise accuracy, aggregate log-probabilities over HRF window after onset `τ`: `score_c(τ) = Σ_{k=0}^{K-1} h_k log q_{τ+k, c}`; predict `argmax_c score_c`.

### Practical considerations

- Learning `h`: one HRF per ROI/searchlight (share across classes unless expecting class-specific timing).
- Background class `c=0` absorbs ongoing activity; prior weight handled via `β₀` in `π_t`.
- Decoder regularization: default L2; consider group-lasso for high-dimensional searchlights.
- Hyperparameter search via nested CV (`λ ∈ {1e-3…1e1}`, `ρ ∈ {1e-2…1e2}`, `β ∈ {1e-3…1e-1}`).
- Speed tips: precompute `-log q` for Z-step; small K (≤20) keeps H-step trivial; W-step dominates, so use dimensionality reduction and warm starts.

### Pseudocode sketch

```python
h = canonical_hrf(K); h /= h.sum()
pi = softmax(stack([conv(s_c, h) for c in classes] + [beta0_background]))
z = pi.copy()

for iter in range(Iters):
    W, b = fit_logistic(Y_train, targets=z, reg=alpha)
    q = softmax(Y_train @ W + b)
    z = solve_banded_qp(
        linear_cost=-np.log(q),
        quadratic_terms=rho * I + lambda_ * (D2.T @ D2),
        prior=rho * pi,
        constraints="z_t in simplex"
    )
    h = nnls_ridge_fit(z[:, 1:C], event_trains, roughness=beta, K=K)
    h = project_nonneg_normalize(h)

# Testing
q_test = softmax(Y_test @ W + b)
score_c(tau) = sum_{k=0}^{K-1} h[k] * log q_test[tau + k, c]
```

### Outcomes & variants

- Outputs: TR-wise decoder probabilities `q_t`, interpretable soft labels `z_t`, estimated HRF `h`, and trial-wise predictions via HRF-weighted evidence.
- Variants:
  - Basis HRF (`h = B θ`) for small QP.
  - TV penalty on `z` (`||Dz||₁`) for edge-preserving smoothness (ADMM-friendly).
  - Fix `z = π`, augment decoder with temporal features (lagged `y_t`) for fast baseline.
  - Multi-task penalties on `W` to share structure across ROIs.

### Advantages

- Elegant single objective coupling data, design, physiology; eliminates beta extraction.
- Innovative: learns smooth, HRF-consistent soft labels with data-driven HRF adjustment.
- Fast: each subproblem convex/small; banded structure enables near-linear-time updates, practical at searchlight scale.

## Section 2: R package scaffold (hrfdecode)

- Goal: deliver HRF-aware weakly supervised MVPA as an R package backed by Rcpp/RcppArmadillo for core computations (W-, Z-, H-steps). Package name `hrfdecode`; maintains MVPA-style training/testing but bypasses trial betas.
- Skeleton setup: run usethis commands to create package, git, MIT license, add core dependencies (`Rcpp`, `RcppArmadillo`, `Matrix`, `stats`), enable Rcpp, testthat.
- Directory layout:
  - `DESCRIPTION`: metadata, linking to Rcpp/RcppArmadillo, imports (`Rcpp`, `Matrix`, `stats`), SystemRequirements C++17.
  - `NAMESPACE` exporting fit/predict/utility functions, registering native routines via `useDynLib`.
  - `R/`: `fit.R`, `predict.R`, `utils.R`, `simulate.R`, `zzz.R`.
  - `src/`: `softmax_lbfgs.cpp` (W-step), `admm_z.cpp` (Z-step), `hrf_fit.cpp` (H-step), `Makevars`, plus auto-generated `RcppExports.cpp`.
  - `tests/testthat/test-basic.R` for smoke test; optional `vignettes/`.

### R functions

- `fit_hrfdecoder()` (R/fit.R):
  - Inputs: matrix `Y` (T×V), list of onset times per class, TR, HRF support `K`, hyperparams (`lambda`, `rho`, `alpha`, `beta`, `mu`, `admm_iter`, `outer_iter`), flags (`background`, `standardize`), verbosity.
  - Workflow: optional z-score Y; build event matrix (`make_event_matrix`); initialize HRF via `canonical_hrf` and design prior `pi`; iterate W/Z/H:
    1. W-step via C++ `softmax_lbfgs` (L2-regularized multiclass logistic with soft targets).
    2. Compute decoder probabilities `Q`, log-probs `logQ`; Z-step via C++ `admm_z` using ADMM to enforce smoothness (second-diff), adherence to design prior, and simplex constraints.
    3. H-step via C++ `hrf_fit`, using projected-gradient NNLS on simplex with roughness; recompute `pi`.
  - Returns list with `W`, `b`, HRF `h`, settings, and training `z`/`pi`.

- `predict_hrfdecoder()` (R/predict.R):
  - Computes logits `Y_test %*% W + b`, returns TR-wise probabilities by default.
  - For `mode="trial"`, requires `onsets_test`; accumulates HRF-weighted log probabilities over window of length `K` to produce per-trial scores and predicted class.

- `utils.R` helpers:
  - `row_softmax`, `toeplitz_conv`, exported `make_event_matrix`, exported `canonical_hrf` (double-gamma, clipped nonnegative, normalized).
- `simulate.R`: `simulate_rapid_er()` to generate synthetic ROI data (makes event trains, HRF, latent soft labels, voxel weights, noise).
- `zzz.R`: registers native routines (`@useDynLib`, `@importFrom Rcpp sourceCpp`).

### C++ cores (src/)

- `softmax_lbfgs.cpp`:
  - Implements row-wise softmax, objective/gradient, and L-BFGS optimizer with Armijo line search and limited memory `m`.
  - Parameters flattened as `[vec(W); b]`; gradient includes L2 penalty on W only.
  - Returns fitted `W`, `b`.

- `admm_z.cpp`:
  - Constructs sparse second-difference `D2`, uses `D2ᵀD2` (pentadiagonal) in smoothing term.
  - ADMM variables: primal `Z`, smoothed `U`, simplex-projected `V`, duals `A`, `B`.
  - Z-update closed-form; U-update solves `(μI + 2λ D2ᵀD2)U = μ(Z+A)` per class via sparse solve (SuperLU) with CG fallback; V-update projects each TR row onto simplex.
  - Returns `V` (feasible `z`).

- `hrf_fit.cpp`:
  - Builds Toeplitz design per event class, sums to `XtX` and `Xᵀz`.
  - Adds roughness penalty via second-difference matrix.
  - Solves small K-dimensional NNLS with simplex constraint using projected gradient descent with step size from spectral norm.

- `src/Makevars`: sets `CXX_STD = CXX17`.

### Testing & build steps

- Minimal test (`tests/testthat/test-basic.R`): simulate data, run `fit_hrfdecoder`, ensure HRF length correct and predictions have right dimensions.
- Build pipeline: `Rcpp::compileAttributes()`, `devtools::document()`, `devtools::build()`, `devtools::install()`.

### Usage example

```r
sim <- simulate_rapid_er(T=300, V=100, TR=1, onsets=list(seq(10,290,20), seq(14,294,20)))
fit <- fit_hrfdecoder(sim$Y, sim$onsets, TR=1,
                      K=16, lambda=0.5, rho=1.0, alpha=0.1, beta=0.01,
                      mu=1.0, admm_iter=40, outer_iter=4,
                      background=TRUE, verbose=1)
Q_tr <- predict_hrfdecoder(fit, sim$Y, mode="tr")
pred <- predict_hrfdecoder(fit, sim$Y, onsets_test=sim$onsets, mode="trial")
```

- `pred` returns per-trial TR index, true/predicted classes, and log-prob scores.

### Extensions & notes

- Potential additions: searchlight/ROI runner with shared HRFs, HRF basis parameterization (`h = Bθ`), total-variation penalty on `z`, group penalties on `W`, banded Cholesky for Z-step smoothing.
- Design decisions: SPD smoothing matrix solved via SuperLU; intercept (`b`) unregularized; background class optional; recommend PCA to ~50–100 dims for high-V searchlights.
- Additional guidance: use `future.apply`/`RcppParallel` for ROI-level parallelism; deterministic pentadiagonal solver can replace SuperLU if needed.

## Section 3: Integration with fmridesign & fmrihrf

- Objective: refactor `hrfdecode` package to consume design/HRF infrastructure from `fmridesign` and `fmrihrf`, keeping W/Z steps in Rcpp but delegating onset handling, nuisance regression, and HRF basis/penalties to those packages.

### Package metadata updates

- `DESCRIPTION`: bump to version 0.2.0; add Imports `Rcpp`, `RcppArmadillo`, `Matrix`, `stats`, `fmridesign (>=0.1.0)`, `fmrihrf (>=0.1.0)`; add `Remotes: bbuchsbaum/fmridesign, bbuchsbaum/fmrihrf`; note LinkingTo (Rcpp/RcppArmadillo) and `SystemRequirements: C++17`.
- Description text updated to emphasize fmrihrf basis and fmridesign models.

### Interop helpers (`R/interop_fmri.R`)

- `.get_sframe(ev_model)`: derives sampling frame (block lengths, TR) from event model attributes; falls back to median sample spacing via `fmrihrf::samples`.
- `build_condition_basis(ev_model, hrf = NULL)`:
  - Ensures `event_model` input; obtains sampling frame if missing.
  - Selects HRF basis (explicit arg or default via `fmrihrf::getHRF("spmg1")`).
  - Extracts first event term, splits into per-condition T×P basis matrices using `fmridesign::condition_basis_list()`.
  - Returns list with `X_list` (named by condition), HRF object, sampling frame, condition names.
- `residualize_baseline(Y, base_model)`: calls `fmridesign::residualize` on `baseline_model` design to project out nuisance regressors; identity if `NULL`.

### HRF basis solver (`src/hrf_fit_basis.cpp`)

- New C++ exported function `hrf_fit_basis_cpp` solving `(XtX + beta * R) θ = Xᵀz` where `XtX` sums over per-condition basis cross-products, `R` is fmrihrf penalty matrix; uses `solve` with `solve_opts::likely_sympd`, general solve fallback.
- Replaces previous sample-wise HRF `hrf_fit.cpp` (now obsolete); smoothness controlled via `fmrihrf::penalty_matrix`.
- After fitting `θ`, convert to HRF object via `fmrihrf::hrf_from_coefficients`.

- Existing `softmax_lbfgs.cpp` and `admm_z.cpp` remain unchanged for W/Z.

### Training refactor (`R/fit.R`)

- Signature now expects `ev_model` (required `event_model`) and optional `base_model` (for nuisance); optional `hrf` object.
- Steps:
  1. Validate matrix input, optional z-scoring; residualize via `residualize_baseline`.
  2. Call `build_condition_basis` to get condition-specific regressors `X_list`, HRF object, sampling frame, condition names.
  3. Initialize HRF coefficients `θ` (e.g., first basis column) and design prior `π` by projecting basis with θ; optionally append background column.
  4. Iterate W/Z/H: W-step via `softmax_lbfgs`; Z-step via `admm_z`; H-step builds `XtX = Σ_c X_cᵀX_c`, `Xᵀz = Σ_c X_cᵀ z_c`, penalty `R = fmrihrf::penalty_matrix(hrf)`, solves θ with `hrf_fit_basis_cpp`, rebuilds `π`.
  5. Store final HRF as `fmrihrf::hrf_from_coefficients(hrf_obj, θ)` plus metadata (conditions, TR from sampling frame median spacing).
- Model object now contains `theta`, `hrf` (fmrihrf), `conditions`, and training priors.

### Prediction refactor (`R/predict.R`)

- TR mode unchanged (matrix multiply, softmax).
- Trial mode now requires `ev_model_test` (event_model):
  - Extract event table via `fmridesign::events`, mapping conditions to columns learned during training.
  - Evaluate fitted HRF kernel via `fmrihrf::evaluate(hrf, tgrid)` where `tgrid` spans `span` seconds at TR spacing (default 24s).
  - For each event, accumulate HRF-weighted log-probs from `Q` (dropping background column) to produce per-trial scores and predicted class; output data frame with TR index, onset, true/pred class, and score columns.

### Utilities cleanup

- Removed obsolete functions `canonical_hrf`, `make_event_matrix`, `toeplitz_conv`; design/HRF ops now handled entirely by fmridesign/fmrihrf.

### Example workflow (interop)

1. Build sampling frame (`fmrihrf::sampling_frame`) with multiple runs, TR.
2. Create event table with conditions/onsets/runs; feed into `fmridesign::event_model` using `hrf(condition, basis="spmg3")`; create baseline model for nuisance.
3. Simulate or load ROI data `Y`.
4. Fit via `fit_hrfdecoder(Y, ev_model=evmod, base_model=basem, hrf=fmrihrf::getHRF("spmg3"), ...)`.
5. Obtain TR-wise probabilities (`predict(..., mode="tr")`) or trial-wise predictions using same/held-out event model and HRF span.

- Demonstrates seamless use of sampling_frame, event_model, baseline_model, fmrihrf basis evaluation, and `hrfdecode` training.

### Notes & design rationale

- `event_model` + `sampling_frame` govern timing, block structure, TR mapping; `condition_basis_list` yields authoritative regressors per condition (handles multi-run onsets, high-res timing).
- `baseline_model` + `residualize` ensure nuisance regressors removed before decoding, aligning with fmridesign pipeline.
- HRF estimation done in basis space with penalty matrices specialized to chosen basis (FIR, SPMG, splines, etc.), enabling easy swapping and maintaining compatibility with fmrihrf plotting/evaluation tools.
- Model stores fmrihrf HRF object so downstream code can `evaluate`, `plot`, or reuse it in design construction.
- Column metadata accessible via `design_colmap` if needed; `conditions` tracked within model for matching trial predictions.
- Core computation remains fast: W-step (L-BFGS softmax), Z-step (ADMM with banded SPD solve), H-step (tiny SPD solve in P dimensions).

### Tests

- New test ensures fmridesign/fmrihrf interop: build sampling frame + event/baseline models, simulate data, run `fit_hrfdecoder` with few iterations, confirm θ length >0 and predict returns matrix of proper size.

### Future extensions

- ROI/searchlight runner using shared `ev_model`/`base_model`, per-ROI θ estimation.
- Class-specific HRFs by fitting θ per condition (avoid summing XtX across conditions).
- Alternative penalties (e.g., TV on `z`) by swapping Z-step; fmridesign/fmrihrf integration unchanged.

## Section 4: rMVPA integration plan

- Goal: add a continuous-time MVPA path that plugs into rMVPA’s existing searchlight and cross-validation machinery using the HRF-aware soft-label decoder implemented in `hrfdecoder`. Leverages fmridesign/fmrihrf for event modeling, HRF bases, and nuisance handling while keeping rMVPA iterators untouched.

### hrfdecoder package roadmap

1. **Data prep glue**
   - `prepare_event_design()`: wrap `fmridesign::event_model()` to map onsets/durations to TR grid per run; pull HRF basis via `fmrihrf` (canonical/FLOBS); optionally build `baseline_model()` and return nuisance matrix for projection.
   - `build_softlabel_prior()`: construct initial TR×K soft-label grid from HRF basis + event sticks (no trial betas), producing `DBβ` prior for optimization.

2. **Core solver (Rcpp/Eigen)**
   - Alternating least squares minimizing `||XW - P||_F^2 + λ_HRF||P - DBβ||_F^2 + λ_smooth ∑_k ||L p_k||_2^2 + λ_W||W||_F^2`:
     - `X`: T×V ROI data (baseline-projected).
     - `P`: TR×K soft labels; constrained nonnegative (optionally row-sum≤1 via “rest” class).
     - `W`: V×K decoder weights (multiresponse ridge).
     - `DBβ`: HRF-convolved prior per class (from fmridesign + fmrihrf basis).
     - `L`: block-diagonal second-difference Laplacian per run (no cross-run bleed).
   - Closed-form updates:
     - `W ← (XᵀX + λ_W I)⁻¹ Xᵀ P` (ridge per class or joint multiresponse).
     - `p_k ← (I + λ_HRF I + λ_smooth L)⁻¹ (X w_k + λ_HRF (DBβ)_k)` solved with SPD solver (Eigen SimplicialLDLT); optional β update for HRF basis coefficients.
     - Project `P` to ≥0 each iteration; optional row-sum constraint.
   - Rcpp API:
     ```cpp
     List fit_softlabels_als(const Eigen::Map<MatrixXd>& X,
                             const Eigen::Map<MatrixXd>& P0,
                             const Eigen::Map<SparseMatrix<double>>& L,
                             const Eigen::Map<MatrixXd>& DBbeta,
                             double lambda_W, double lambda_HRF, double lambda_smooth,
                             int max_iter, double tol, bool nonneg, bool threads);
     MatrixXd predict_softlabels(const Eigen::Map<MatrixXd>& Xtest,
                                 const Eigen::Map<MatrixXd>& W);
     ```
   - Return `hrfdecoder_fit` S3 object storing `W`, final `P`, HRF basis coefficients β, metadata.

3. **R wrappers**
   - `hrfdecoder_fit(X, event_model, basis, baseline=NULL, ...)`: orchestrates data prep, builds Laplacian/prior, calls C++ solver; returns `hrfdecoder_fit`.
   - `predict.hrfdecoder_fit(fit, Xnew, ...)`: outputs TR×K soft labels (decoder scores).
   - `aggregate_events(P, event_table, window=c(4,8), weights="hrf"|"flat")`: converts TR-level soft labels to event-level probabilities using HRF-weighted windows post-onset.

### Rationale

- Single alternating ridge + SPD solves per class/run (fast in Eigen), no per-trial GLMs.
- HRF/baseline handling delegated to existing packages; overlap resolved via `P` updates; smoothness via Laplacian prevents “wiggles”.

### rMVPA integration

1. **Adapter (`as_mvpa_dataset.fmri_dataset`)**
   - Convert `fmridataset::fmri_dataset` to `mvpa_dataset` by extracting `NeuroVec` and mask so rMVPA’s ROI/searchlight extraction works without change.

2. **Continuous design (`continuous_mvpa_design`)**
   - Wrap `mvpa_design` with dummy `y_train` but store `event_model`, per-TR block/run IDs, and trial table.
   - CV remains blocked by run via `block_var`; `y_train` used only for fold construction.

3. **New model class (`hrfdecoder_model`)**
   - Use `create_model_spec()` to define spec with hyperparameters (`lambda_W`, `lambda_HRF`, `lambda_smooth`, HRF basis, window, nonneg, max_iter, tol), cross-validation strategy (default blocked via runs), and optional performance function.
   - S3 hooks:
     - `y_train.hrfdecoder_model`: returns dummy vector (length T).
     - `train_model.hrfdecoder_model`: for each ROI/fold, extract TR×V matrix, prepare decoder inputs (event design, Laplacian, priors), call `hrfdecoder_fit`, return wrap containing fit + context.
     - `format_result.hrfdecoder_model`: predict TR-level soft labels on test TRs, aggregate to event-level predictions, construct classification-style tibble (`class`, `probs`, `y_true`, `test_ind`, optional fit) so rMVPA’s combiner/performance machinery applies.
     - `merge_results.hrfdecoder_model`: combine fold outputs per ROI, wrap as classification_result, compute performance (either custom or default multiclass).
   - Optional `run_searchlight.hrfdecoder_model` to choose specific combiners (`combine_randomized`, `pool_randomized`), though default runner works.

4. **Searchlight/CV usage**
   - Core iterators (mvpa_iterate, run_searchlight) remain untouched; they call the new model’s hooks.
   - Performance aggregation uses existing combiners; cross-validation uses blocked splits (or custom).

5. **Example workflow**

```r
dset <- as_mvpa_dataset(fmridataset::fmri_dataset(...))
evm  <- fmridesign::event_model(onsets ~ hrf(condition, basis="spmg3"), data=..., block=~run,
                                sampling_frame=fmrihrf::sampling_frame(...))
des  <- continuous_mvpa_design(event_model=evm,
                               block_var=evm$tr_grid$run,
                               design_df_events=evm$events)
mspec <- hrfdecoder_model(dataset=dset, design=des,
                          basis=fmrihrf::spmg1(),
                          lambda_W=10, lambda_HRF=1, lambda_smooth=5,
                          window=c(4,8))
sl_res <- run_searchlight(mspec, radius=8, method="randomized", niter=4)
```

- Searchlight outputs event-level metrics per ROI using existing rMVPA infrastructure.

### Implementation details to enforce stability

- Laplacian `L` block-diagonal per run to avoid inter-run smoothing.
- Nonnegativity enforced via projection after `P` updates; optional row-sum≤1 by adding rest class.
- Nuisance removal: either residualize `X` with baseline regressors or include baseline columns in `P` prior update.
- Factorizations reused within fold since `(I + λ_HRF I + λ_smooth L)` constant; parallelize class updates and folds (OpenMP).

### File layout guidance

- `hrfdecoder`:
  - `R/fit.R`, `R/prep.R`, `R/adapters.R`, etc.
  - `src/softlabels_als.cpp`, `src/laplacian.cpp/h`, `src/init.cpp`.
- `rMVPA` (or `hrfdecoder` depending on ownership):
  - `R/hrfdecoder_model.R` defining model spec and S3 hooks.

### Optional enhancements

- AR(1) prewhitening per run prior to fitting.
- Basis learning via β updates.
- Diagnostics exposing smoothness/HRF fidelity per fold.
- Unit tests: convergence on synthetic ROIs, overlap recovery, smoothness energy reduction.

- Ready path to extend to searchlight-scale analyses without modifying rMVPA core; relies on documented “new analysis” pattern.

## Section 5: Implementation plan

1. **Package scaffolding & deps**
   - Bump `hrfdecode` metadata (DESCRIPTION/NAMESPACE) to include `fmridesign`, `fmrihrf`, Remotes, and remove legacy HRF utilities.
   - Add interop helpers for event/baseline handling and ensure baseline residualization utilities are in place.

2. **Eigen solver integration**
   - Implement `fit_softlabels_als`/`predict_softlabels` plus Laplacian builders in C++ (Rcpp/Eigen), exposing them via R wrappers.
   - Update `fit_hrfdecoder`, `predict_hrfdecoder`, and `aggregate_events` to use fmrihrf/fmridesign bases, penalties, and event models.

3. **Dataset/design adapters**
   - Provide adapters (`as_mvpa_dataset.fmri_dataset`, `continuous_mvpa_design`, preparation helpers) that map fmri_dataset inputs and event models into the solver-ready Laplacian/prior structures.

4. **rMVPA model hooks**
   - Define the `hrfdecoder_model` spec and required S3 methods (`y_train`, `train_model`, `format_result`, `merge_results`) so the decoder slots into existing rMVPA searchlight and cross-validation flows.

5. **Tests & examples**
   - Add synthetic tests for convergence/interop, fmrihrf/fmridesign compatibility checks, and an rMVPA searchlight example to verify end-to-end wiring.
