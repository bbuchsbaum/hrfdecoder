# hrfdecode

<div class="callout">
  <strong>Primary metric</strong>: In rMVPA integrations, results flag a primary metric and value (default <em>AUC (macro one‑vs‑one)</em>). Configure via <code>hrfdecoder_model(metrics = c('auc_ovo','accuracy'), primary_metric = 'auc_ovo')</code>. See the <a href="articles/06-performance-metrics.html">Performance metrics</a> article.
</div>

HRF‑aware weakly supervised decoding directly from fMRI time series and event timing. Joint estimation of soft labels, HRF parameters, and decoder weights via alternating least squares with temporal smoothing and optional AR prewhitening.

- Get started: articles/01-getting-started.html
- Metrics and confidence: articles/06-performance-metrics.html
- Reference: reference/index.html
