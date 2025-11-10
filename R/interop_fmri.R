#' @keywords internal
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' @keywords internal
.get_sframe <- function(ev_model) {
  stopifnot(inherits(ev_model, "event_model"))
  sf <- attr(ev_model, "sampling_frame", exact = TRUE)
  if (!is.null(sf)) return(sf)
  blocklens <- attr(ev_model, "blocklens", exact = TRUE)
  if (is.null(blocklens)) {
    dm <- fmridesign::design_matrix(ev_model)
    blocklens <- nrow(dm)
  }
  tr_attr <- attr(ev_model, "TR", exact = TRUE)
  tr <- tr_attr %||% 2
  fmrihrf::sampling_frame(blocklens = blocklens, TR = tr)
}

#' Build per-condition basis matrices for decoding
#' @param ev_model fmridesign::event_model result
#' @param hrf optional fmrihrf HRF object
#' @keywords internal
build_condition_basis <- function(ev_model, hrf = NULL) {
  stopifnot(inherits(ev_model, "event_model"))
  sframe <- attr(ev_model, "sampling_frame", exact = TRUE)
  if (is.null(sframe)) sframe <- .get_sframe(ev_model)
  if (is.null(hrf)) {
    hrf <- fmrihrf::getHRF("spmg1")
  }

  # Rebuild a consistent event_model from its events + sframe to ensure
  # block/run metadata and timing align with the sampling frame.
  events <- get_event_table(ev_model)
  if (is.null(events$run)) {
    events$run <- 1L
  }
  # Basic validation of runs and onsets
  # Prefer fmrihrf accessor for robust block lengths
  blocklens <- tryCatch(as.integer(fmrihrf::blocklens(sframe)), error = function(e) NULL)
  if (is.null(blocklens)) {
    blocklens <- attr(sframe, "blocklens", exact = TRUE) %||% sframe$blocklens
  }
  if (is.null(blocklens)) {
    # Fallback: single block of full length (best-effort)
    blocklens <- c(nrow(fmridesign::design_matrix(ev_model)))
  }
  nblocks <- length(blocklens)
  if (any(events$run < 1L | events$run > nblocks)) {
    stop("prepare_decoder_inputs: events$run contains values outside 1..", nblocks)
  }
  TR <- tryCatch(as.numeric(attr(sframe, "TR", exact = TRUE)), error = function(e) NULL)
  if (is.null(TR)) TR <- sframe$TR %||% 2
  total_time <- sum(blocklens) * TR
  if (any(events$onset > total_time)) {
    warning("Some event onsets exceed total acquisition time; they will be ignored downstream.")
  }

  basis_name <- attr(hrf, "name") %||% "spmg1"
  # Canonicalize basis name if coming from fmrihrf::hrf_from_coefficients
  available_hrfs <- try(fmrihrf::list_available_hrfs(), silent = TRUE)
  if (!inherits(available_hrfs, "try-error")) {
    if (!(basis_name %in% available_hrfs)) {
      basis_name <- sub("_from_coef$", "", basis_name)
      if (!(basis_name %in% available_hrfs)) basis_name <- "spmg1"
    }
  } else {
    basis_name <- sub("_from_coef$", "", basis_name)
  }
  ev_model2 <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = basis_name),
    data = events,
    block = ~ run,
    sampling_frame = sframe
  )

  terms <- fmridesign::event_terms(ev_model2)
  if (length(terms) != 1L) {
    warning("Multiple event terms detected; using the first term for decoding.")
  }
  term <- terms[[1L]]
  X_list <- fmridesign::condition_basis_list(term, hrf, sframe, output = "condition_list")
  conds <- names(X_list)
  list(X_list = X_list, hrf = hrf, sframe = sframe, conditions = conds)
}

#' Residualize Y against a baseline model
#' @param Y numeric matrix (T x V)
#' @param base_model optional fmridesign::baseline_model
#' @keywords internal
residualize_baseline <- function(Y, base_model = NULL) {
  if (is.null(base_model)) return(Y)
  fmridesign::residualize(base_model, Y)
}

#' @keywords internal
get_event_table <- function(ev_model) {
  out <- tryCatch(
    fmridesign::events(ev_model),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    terms <- tryCatch(fmridesign::event_terms(ev_model), error = function(e) NULL)
    if (!is.null(terms) && length(terms) > 0) {
      tables <- lapply(terms, function(term) {
        tbl <- term$event_table
        if (is.null(tbl)) return(NULL)
        df <- as.data.frame(tbl)
        if (!is.null(term$onsets)) df$onset <- term$onsets
        if (!is.null(term$blockids)) df$run <- term$blockids
        if (!is.null(term$durations)) df$duration <- term$durations
        df
      })
      tables <- Filter(Negate(is.null), tables)
      if (length(tables) > 0) {
        tab <- tables[[1]]
        if (!is.null(tab)) return(tab)
      }
    }
    tab <- attr(ev_model, "events", exact = TRUE)
    if (!is.null(tab)) return(tab)
    stop(out)
  }
  out
}
