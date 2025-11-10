# Engine Test Plan

We will implement the following tests to raise confidence in the hrfdecoder ALS engine. Each item will be checked off once implemented.

1. [x] **Objective vs. helper parity**  
   Compare `als_objective()` (R helper) against the C++ `engine_objective()` for random inputs; assert equality within tolerance.

2. [x] **Finite-difference gradient check**  
   Tiny synthetic problem (e.g., T=6, V=3): perturb entries of `W`/`P` and confirm objective changes match gradient predictions.

3. [x] **Ground-truth recovery**  
   Simulate data with known `W*`, `P*`, `θ*` and verify fitted parameters correlate strongly with the truth (≥0.95).

4. [x] **HRF basis stress test**  
   Use a multi-basis HRF (e.g., SPMG3/FLOBS) and ensure `estimate_hrf_theta()` recovers random θ with low relative error.

5. [x] **Run-boundary isolation**  
   Events in run 1 only; confirm smoothing + HRF penalties do not leak mass into run 2 (`max|P_run2| < 1e-6`).

6. [x] **Laplacian limit cases**  
   - λ_smooth → 0 ⇒ `P ≈ XW`.  
   - λ_smooth → large ⇒ each run’s `P` is nearly constant (second differences ≈ 0).

7. [x] **HRF limit cases**  
   - λ_HRF → 0 ⇒ `P` follows `XW`.  
   - λ_HRF → ∞ ⇒ `P` follows `DBβ`.

8. [x] **Edge-case aggregation**  
   Events at TR boundaries; spans longer than available TRs; ensure aggregation truncates correctly and stays run-local.

9. [x] **Monte Carlo stability sweep**  
   Random λ_W/H/smooth combinations on small problems; assert convergence (final `rel_dP` < 1e-2) and monotone objective.

10. [x] **Warm-start continuity**  
    Run ALS for k iterations, resume with `P0 = P_k`, and verify objectives/traces continue without jumps.
