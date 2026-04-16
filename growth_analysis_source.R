# Core growth-analysis functions consolidated from project source files:
# - LT1_funcs.R: get_growth, get_growth_dist
# - scripts_growth.R: processGrowth
# - special_helpers.R: spline/normalization/stat helpers
# - funcs.R: rowCVs
# - PhenoFunctions.R: adaptiveGaussianSNN

rowCVs <- function(df) {
  apply(df, 1, function(x) stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)
}

cv <- function(x, ..., aszero = FALSE, na.rm = FALSE) {
  x <- c(x, ...)
  z <- x[!is.na(x)]
  if (length(z) == 0) {
    return(NA_real_)
  } else if (!na.rm && (length(z) < length(x))) {
    return(NA_real_)
  } else if (length(z) == 1 && aszero) {
    return(0)
  } else {
    x <- mean(abs(z))
    if (x == 0) {
      return(0)
    } else {
      return(100 * stats::sd(z) / x)
    }
  }
}

get_colFracts <- function(data) {
  t(t(data) / colSums(data, na.rm = TRUE))
}

min_max_normalization <- function(x, new_min = 0, new_max = 1) {
  if (length(x) == 1) return(new_min)
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) *
    (new_max - new_min) + new_min
}

min_max_normalization_cols_matrix <- function(mat, new_min, new_max) {
  col_mins <- matrixStats::colMins(mat, na.rm = TRUE)
  col_maxs <- matrixStats::colMaxs(mat, na.rm = TRUE)
  ((mat - matrix(col_mins, nrow(mat), ncol(mat), byrow = TRUE)) /
     matrix(col_maxs - col_mins, nrow(mat), ncol(mat), byrow = TRUE)) *
    (new_max - new_min) + new_min
}

scale2 <- function(data, center = TRUE, scale = TRUE) {
  data <- scale(data, center = center, scale = scale)
  attr(data, "scaled:center") <- NULL
  attr(data, "scaled:scale") <- NULL
  data
}

geometric_mean <- function(x, na.rm = TRUE, epsilon = 1e-12) {
  x <- as.numeric(x)
  if (na.rm) {
    x <- x[is.finite(x)]
  } else if (any(!is.finite(x))) {
    return(NA_real_)
  }
  if (!length(x)) return(NA_real_)
  x[x <= 0] <- epsilon
  exp(mean(log(x)))
}

resolve_curve_times <- function(curves, times = NULL) {
  if (!is.null(times)) return(as.numeric(times))
  if (is.null(colnames(curves))) {
    stop("times is NULL and curves have no colnames to parse.")
  }
  parsed <- sub("^.*_", "", colnames(curves))
  times_num <- suppressWarnings(as.numeric(parsed))
  if (any(!is.finite(times_num))) {
    stop("Could not parse numeric times from curve colnames.")
  }
  times_num
}

earliness_score <- function(curves, gr_centmass, times = NULL) {
  inp <- as.matrix(curves)
  if (!is.numeric(inp)) storage.mode(inp) <- "double"
  if (is.null(rownames(inp))) stop("curves must have rownames.")
  if (!nrow(inp)) stop("curves has zero rows.")
  if (!ncol(inp)) stop("curves has zero columns.")

  times <- resolve_curve_times(inp, times = times)
  if (length(times) != ncol(inp)) stop("Length of times must match ncol(curves).")

  inp <- inp + find_min_except_zero(inp)
  if (!is.null(names(gr_centmass))) {
    gr_centmass <- gr_centmass[rownames(inp)]
  }
  gr_centmass <- as.numeric(gr_centmass)
  if (length(gr_centmass) != nrow(inp)) stop("Length of gr_centmass must match nrow(curves).")

  r2_values <- apply(inp, 1, function(pop) summary(stats::lm(log(pop) ~ times))$r.squared)
  intercept_values <- apply(inp, 1, function(pop) stats::coef(stats::lm(log(pop) ~ times))[1])
  plotdata <- cbind(gr_centmass = gr_centmass, r2_values = r2_values, inv_intercept = 1 / intercept_values)
  rownames(plotdata) <- rownames(inp)
  plotdata <- plotdata[stats::complete.cases(plotdata), , drop = FALSE]

  out_scores <- setNames(rep(NA_real_, nrow(inp)), rownames(inp))
  if (nrow(plotdata)) {
    plotdata_scaled <- scale2(plotdata)
    plotdata_scaled <- min_max_normalization_cols_matrix(plotdata_scaled, 0.01, 1)
    gm <- apply(plotdata_scaled, 1, geometric_mean)
    out_scores[names(gm)] <- gm
  }

  list(scores = out_scores, plotdata = plotdata)
}

lateness_score <- function(curves, gr_centmass, times = NULL) {
  inp <- as.matrix(curves)
  if (!is.numeric(inp)) storage.mode(inp) <- "double"
  if (is.null(rownames(inp))) stop("curves must have rownames.")
  if (!nrow(inp)) stop("curves has zero rows.")
  if (!ncol(inp)) stop("curves has zero columns.")

  times <- resolve_curve_times(inp, times = times)
  if (length(times) != ncol(inp)) stop("Length of times must match ncol(curves).")

  inp <- inp + find_min_except_zero(inp)
  if (!is.null(names(gr_centmass))) {
    gr_centmass <- gr_centmass[rownames(inp)]
  }
  gr_centmass <- as.numeric(gr_centmass)
  if (length(gr_centmass) != nrow(inp)) stop("Length of gr_centmass must match nrow(curves).")

  plotdata <- cbind(inv_gr_centmass = 1 / gr_centmass, inv_rowsum = 1 / rowSums(inp))
  rownames(plotdata) <- rownames(inp)
  plotdata <- plotdata[stats::complete.cases(plotdata), , drop = FALSE]

  out_scores <- setNames(rep(NA_real_, nrow(inp)), rownames(inp))
  if (nrow(plotdata)) {
    plotdata_scaled <- scale2(plotdata)
    plotdata_scaled <- min_max_normalization_cols_matrix(plotdata_scaled, 0.01, 1)
    gm <- apply(plotdata_scaled, 1, geometric_mean)
    out_scores[names(gm)] <- gm
  }

  list(scores = out_scores, plotdata = plotdata)
}

MinMaxQuant <- function(data, q = 0.1, low_q = NULL, high_q = NULL) {
  if (!is.null(q)) {
    if (sign(0.5 - q) == -1) {
      low_q <- 1 - q
      high_q <- q
    } else {
      low_q <- q
      high_q <- 1 - q
    }
    if (is.null(low_q)) low_q <- q
    if (is.null(high_q)) high_q <- q
  } else if (!is.null(low_q) && !is.null(high_q)) {
    if (low_q > high_q) stop("Error low_q > high_q")
  } else {
    stop("Error: q or low_q and high_q should be provided")
  }
  min_v <- stats::quantile(data, low_q, na.rm = TRUE)
  max_v <- stats::quantile(data, high_q, na.rm = TRUE)
  data[data > max_v] <- max_v
  data[data < min_v] <- min_v
  data
}

summarize_over_rows_given_colgroup <- function(mat, groups, func1, retain_colnames = FALSE) {
  if (!is.matrix(mat)) stop("The first argument 'mat' must be a matrix.")
  if (length(groups) != ncol(mat)) stop("Length of 'groups' must match the number of columns in 'mat'.")
  unique_groups <- unique(groups)
  result <- matrix(NA_real_, nrow = nrow(mat), ncol = length(unique_groups))
  colnames(result) <- unique_groups
  for (i in seq_len(nrow(mat))) {
    row_values <- mat[i, ]
    summarized_values <- sapply(unique_groups, function(g) func1(row_values[groups == g]))
    result[i, ] <- summarized_values
  }
  if (retain_colnames && !is.null(rownames(mat))) rownames(result) <- rownames(mat)
  result
}

smooth_timeseries <- function(data, new_time, time, suffix = "T_", deg.fr = 5, spars = NULL,
                              prefit_linear = FALSE, log_transform = FALSE) {
  rownames_data <- rownames(data)
  data <- summarize_over_rows_given_colgroup(data, time, mean)
  time <- unique(time)
  norm_data_fit <- matrix(0, ncol = length(new_time), nrow = nrow(data))
  for (i in seq_len(nrow(norm_data_fit))) {
    if (prefit_linear) {
      fit_data <- approx(time, data[i, ], xout = new_time, method = "linear", rule = 2)$y
      fit_time <- new_time
    } else {
      fit_data <- data[i, ]
      fit_time <- time
    }
    if (log_transform) fit_data <- log10(fit_data)
    if (is.null(spars)) {
      fit <- stats::smooth.spline(fit_time, fit_data, df = deg.fr, spar = NULL)
    } else {
      fit <- stats::smooth.spline(fit_time, fit_data, df = deg.fr, spar = spars[i])
    }
    pred <- stats:::predict.smooth.spline(fit, new_time)$y
    if (log_transform) {
      pred <- 10^pred
      pred[pred < 0] <- 0
    }
    norm_data_fit[i, ] <- pred
  }
  rownames(norm_data_fit) <- rownames_data
  colnames(norm_data_fit) <- paste0(suffix, new_time)
  norm_data_fit
}

find_min_except_zero <- function(mat) {
  non_zero_elements <- mat[mat != 0]
  min(non_zero_elements, na.rm = TRUE)
}

centermass_for_rows2 <- function(inp) {
  inp_gr <- data.matrix(inp)
  inp_gr_rev <- t(apply(inp_gr, 1, rev))
  if (requireNamespace("riskRegression", quietly = TRUE)) {
    left <- riskRegression::rowCumSum(inp_gr)
    right <- t(apply(riskRegression::rowCumSum(inp_gr_rev), 1, rev))
  } else {
    row_cumsum <- function(x) t(apply(x, 1, cumsum))
    left <- row_cumsum(inp_gr)
    right <- t(apply(row_cumsum(inp_gr_rev), 1, rev))
  }
  inp_gr_centmass <- -abs(right - left)
  inp_gr_centmass <- which.colMaxs(t(inp_gr_centmass))
  names(inp_gr_centmass) <- rownames(inp_gr)
  inp_gr_centmass
}

pseudostage_for_rows <- function(input, coeff = 1.1) {
  input_norm <- t(min_max_normalization_cols_matrix(t(input), 0, 1))
  psdstg_mx <- which.colMaxs(t(input_norm))
  psdstg_mx[psdstg_mx == max(psdstg_mx)] <- round(max(psdstg_mx[psdstg_mx != max(psdstg_mx)]) * coeff)
  psdstg_mx
}

spline_smooth <- function(counts, times, fit_interval = 1, deg.fr = 4, prefit_linear = FALSE,
                          log_transform = TRUE, q_rat = 0.01, spar_range = c(0.3, 0.65)) {
  new_times <- seq(min(times) + 1, max(times) + 1, fit_interval)
  counts <- counts + find_min_except_zero(counts) / 2
  cvs <- summarize_over_rows_given_colgroup(counts, as.factor(times), cv)
  means <- summarize_over_rows_given_colgroup(counts, as.factor(times), mean)
  cvs_means <- rowCVs(means)
  cvs_repl <- rowMeans(cvs, na.rm = TRUE)
  rat <- log(cvs_means) / log(cvs_repl)
  rat <- MinMaxQuant(rat, q_rat)
  if (!is.null(spar_range)) {
    spars <- min_max_normalization(rat, spar_range[1], spar_range[2])
  } else {
    spars <- NULL
  }
  cnts_gr_sm <- smooth_timeseries(
    counts, new_times, times, suffix = "T_", deg.fr = deg.fr, spars = spars,
    prefit_linear = prefit_linear, log_transform = log_transform
  )
  list(times = new_times, counts = cnts_gr_sm)
}

get_growth_rates <- function(counts, times) {
  counts <- as.matrix(counts)
  if (any(duplicated(times))) {
    mean_cnts <- summarize_over_rows_given_colgroup(counts, as.factor(times), mean)
    mean_cnts <- mean_cnts + 0.25
    times <- unique(times)
  } else {
    mean_cnts <- counts
  }
  if (is.unsorted(times)) {
    ord <- order(times)
    mean_cnts <- mean_cnts[, ord, drop = FALSE]
    times <- times[ord]
  }
  delta_time <- diff(times)
  tms_inp2d <- t(replicate(nrow(mean_cnts), delta_time))
  gr_rate_log <- log(mean_cnts[, 2:ncol(mean_cnts), drop = FALSE] / mean_cnts[, 1:(ncol(mean_cnts) - 1), drop = FALSE])
  gr_rate_raw <- gr_rate_log / tms_inp2d
  list(growth_rates = gr_rate_raw, delta_time = delta_time)
}

processGrowth <- function(counts, times, fit_interval = 2, deg.fr = 4, prefit_linear = FALSE,
                          log_transform = TRUE, q_rat = 0.01, spar_range = c(0.25, 0.65),
                          metadata = NULL, simulate_no_response = FALSE) {
  counts <- counts[, order(times), drop = FALSE]
  times <- times[order(times)]
  counts <- counts[!(rowSums(counts) == 0), , drop = FALSE]
  cnts_smooth_obj <- spline_smooth(
    counts, times = times, fit_interval = fit_interval, deg.fr = deg.fr,
    prefit_linear = prefit_linear, log_transform = log_transform,
    q_rat = q_rat, spar_range = spar_range
  )
  times_smooth <- cnts_smooth_obj$times
  cnts_smooth <- cnts_smooth_obj$counts
  cnts_smooth_mm <- t(min_max_normalization_cols_matrix(t(cnts_smooth), 0, 1))
  cnts_smooth_maxN <- cnts_smooth / matrixStats::rowMaxs(cnts_smooth)
  gr_centmass <- centermass_for_rows2(cnts_smooth_maxN)
  psds <- pseudostage_for_rows(cnts_smooth, coeff = 1.1)
  cnts_smooth_sort <- cnts_smooth[names(sort(psds, decreasing = FALSE)), , drop = FALSE]
  cnts_smooth_sort <- cnts_smooth_sort[names(sort(gr_centmass, decreasing = FALSE)), , drop = FALSE]
  cnts_smooth_sort_mm <- t(min_max_normalization_cols_matrix(t(cnts_smooth_sort), 0, 1))
  raw_growth_obj <- get_growth_rates(counts, times)
  raw_growth <- raw_growth_obj$growth_rates
  rownames(raw_growth) <- rownames(cnts_smooth)
  delta_time_raw <- raw_growth_obj$delta_time
  growth_smooth_obj <- get_growth_rates(cnts_smooth, times_smooth)
  growth_smooth <- growth_smooth_obj$growth_rates
  delta_time_smooth <- growth_smooth_obj$delta_time
  list(
    counts_raw = counts,
    times_raw = times,
    counts_smooth = cnts_smooth,
    cnts_smooth_maxN = cnts_smooth_maxN,
    cnts_smooth_sort = cnts_smooth_sort,
    cnts_smooth_sort_mm = cnts_smooth_sort_mm,
    gr_raw = raw_growth,
    gr_smooth = growth_smooth,
    gr_deltatime_raw = delta_time_raw,
    gr_deltatime_smooth = delta_time_smooth,
    times_smooth = times_smooth,
    gr_centmass = gr_centmass,
    psds = psds
  )
}

crossdist_to_matrix <- function(crossdist) {
  if (inherits(crossdist, "dist")) return(as.matrix(crossdist))
  mat <- `attr<-`(as.matrix(as.data.frame.matrix(crossdist)), "dimnames", NULL)
  dimnames(mat) <- dimnames(crossdist)
  mat
}

get_growth_dist <- function(gr_obj, overwrite = FALSE, type = "counts_smooth_mm") {
  if (!requireNamespace("dtwclust", quietly = TRUE)) {
    stop("Please install dtwclust package")
  }
  if (is.null(gr_obj$growth_sdtw) || overwrite) {
    if (overwrite) print("Overwriting growth_sdtw")
    growth_sdtw <- proxy::dist(data.matrix(gr_obj[[type]]), method = "sdtw")
    growth_sdtw <- crossdist_to_matrix(growth_sdtw)
    growth_sdtw <- growth_sdtw - find_min_except_zero(growth_sdtw)
    diag(growth_sdtw) <- 0
    gr_obj$growth_sdtw <- growth_sdtw
    gr_obj$dist_type <- type
  }
  gr_obj
}

get_growth <- function(counts, brcs, smpls, metadata, fit_interval = 2, spar_range = c(0.25, 0.65)) {
  if (!all(smpls %in% rownames(metadata))) stop("Not all samples are in rownames metadata")
  if (!all(smpls %in% colnames(counts))) stop("Not all samples are in colnames counts")
  brcs <- intersect(brcs, rownames(counts))
  print("n barcodes overlap: ")
  print(length(brcs))
  metadata_time <- metadata[smpls, , drop = FALSE]
  metadata_time <- metadata_time[order(metadata_time$time), , drop = FALSE]
  smpls <- rownames(metadata_time)
  fractions <- get_colFracts(counts)
  fractions <- fractions[brcs, smpls, drop = FALSE]
  gr_obj <- processGrowth(
    fractions, metadata_time$time, fit_interval = fit_interval, prefit_linear = FALSE,
    log_transform = TRUE, q_rat = 0.01, spar_range = spar_range
  )
  gr_obj$counts_smooth_mm <- t(min_max_normalization_cols_matrix(t(gr_obj$counts_smooth), 0, 1))
  gr_obj
}

keep_inds_in_square <- function(index_table, square_mat, replace = 0) {
  index_table <- as.matrix(index_table)
  if (!is.matrix(square_mat) || nrow(square_mat) != ncol(square_mat)) {
    stop("`square_mat` must be a square matrix.")
  }
  if (nrow(index_table) != nrow(square_mat)) {
    stop("Number of rows in `index_table` must match the number of rows/columns in `square_mat`.")
  }
  n <- nrow(square_mat)
  row_indices <- rep(seq_len(n), times = ncol(index_table))
  col_indices <- as.vector(index_table)
  valid <- !is.na(col_indices) & col_indices >= 1 & col_indices <= n
  rows_keep <- row_indices[valid]
  cols_keep <- col_indices[valid]
  mask <- matrix(FALSE, nrow = n, ncol = n)
  unique_pairs <- unique(cbind(rows_keep, cols_keep))
  mask[cbind(unique_pairs[, 1], unique_pairs[, 2])] <- TRUE
  square_mat_filtered <- square_mat
  square_mat_filtered[!mask] <- replace
  square_mat_filtered
}

set_utf8_locale <- function() {
  candidates <- c("C.UTF-8", "en_US.UTF-8", "UTF-8")
  for (loc in candidates) {
    out <- suppressWarnings(try(Sys.setlocale("LC_CTYPE", loc), silent = TRUE))
    if (!inherits(out, "try-error") && !is.na(out)) return(invisible(out))
  }
  invisible(NULL)
}

read_dslt_rds <- function(path) {
  withCallingHandlers(
    readRDS(path),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("cannot be translated from 'US-ASCII' to UTF-8, but is valid UTF-8", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
      if (grepl("strings not representable in native encoding will be translated to UTF-8", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

parse_times_from_sample_names <- function(sample_names) {
  m <- regexec("^T([0-9]+)(?:_([0-9]+))?$", sample_names, perl = TRUE)
  parts <- regmatches(sample_names, m)
  bad <- lengths(parts) < 2
  if (any(bad)) {
    stop("Invalid sample names for time parsing: ", paste(sample_names[bad], collapse = ", "))
  }
  as.numeric(vapply(parts, `[`, character(1), 2))
}

resolve_existing_file <- function(paths, label = "file") {
  paths <- unique(paths)
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop(
      "Missing ", label, ". Checked: ",
      paste(paths, collapse = " | ")
    )
  }
  existing[[1]]
}

growth_similarity_to_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "dist")) return(as.matrix(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x) && is.numeric(x)) return(x)
  if (is.list(x)) {
    candidates <- c("similaritySNN", "similarity", "distance_matrix", "raw_data")
    for (nm in candidates) {
      val <- x[[nm]]
      if (is.null(val)) next
      if (inherits(val, "dist")) val <- as.matrix(val)
      if (is.data.frame(val)) val <- as.matrix(val)
      if (is.matrix(val) && is.numeric(val)) return(val)
    }
  }
  stop(
    "growth_similarity must be a numeric matrix/dist or a list containing ",
    "a numeric matrix in one of: similaritySNN, similarity, distance_matrix, raw_data."
  )
}

parse_growth_metadata <- function(sample_names) {
  m <- regexec("^T([0-9]+)(?:_([0-9]+))?$", sample_names, perl = TRUE)
  parts <- regmatches(sample_names, m)
  bad <- lengths(parts) < 2
  if (any(bad)) {
    stop("Invalid sample names. Expected T<time>_<rep>: ", paste(sample_names[bad], collapse = ", "))
  }
  time <- as.numeric(vapply(parts, `[`, character(1), 2))
  replica <- as.numeric(vapply(
    parts,
    function(x) if (length(x) >= 3 && nzchar(x[3])) x[3] else "1",
    character(1)
  ))
  metadata <- data.frame(
    sample = sample_names,
    time = time,
    replica = replica,
    treatment = "control_prepared",
    experiment_id = "prepared",
    stringsAsFactors = FALSE,
    row.names = sample_names
  )
  metadata[order(metadata$time, metadata$replica), , drop = FALSE]
}

prepare_numeric_matrix <- function(x, layer_name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("Layer ", layer_name, " must be a matrix or data.frame.")
  if (!is.numeric(x)) stop("Layer ", layer_name, " must be numeric.")
  if (is.null(rownames(x)) || anyDuplicated(rownames(x))) {
    stop("Layer ", layer_name, " must have unique rownames.")
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0(layer_name, "_", seq_len(ncol(x)))
  }
  if (anyDuplicated(colnames(x))) {
    stop("Layer ", layer_name, " has duplicated colnames.")
  }
  storage.mode(x) <- "double"
  x
}

vector_to_single_col_matrix <- function(vec, layer_name, col_name) {
  if (!is.atomic(vec) || is.null(vec)) {
    stop("Layer ", layer_name, " must be an atomic vector.")
  }
  row_ids <- names(vec)
  if (is.null(row_ids) || any(!nzchar(row_ids)) || anyDuplicated(row_ids)) {
    stop("Vector layer ", layer_name, " must have unique non-empty names.")
  }
  mat <- matrix(as.numeric(vec), ncol = 1, dimnames = list(row_ids, col_name))
  prepare_numeric_matrix(mat, layer_name)
}

pad_to_reference_rows <- function(mat, reference_rows, layer_name) {
  extra_rows <- setdiff(rownames(mat), reference_rows)
  if (length(extra_rows)) {
    stop(
      "Layer ", layer_name, " has rows absent from cGR: ",
      paste(utils::head(extra_rows, 10), collapse = ", ")
    )
  }
  padded <- matrix(
    NA_real_,
    nrow = length(reference_rows),
    ncol = ncol(mat),
    dimnames = list(reference_rows, colnames(mat))
  )
  padded[rownames(mat), ] <- mat
  padded
}

add_prepared_growth_to_dslt <- function(
  ds,
  gr_obj,
  layer_prefix = "prepared_",
  dist_type_fallback = "counts_smooth_mm"
) {
  dist_type_fallback <- as.character(dist_type_fallback)[1]
  if (is.na(dist_type_fallback) || !nzchar(dist_type_fallback)) {
    stop("dist_type_fallback must be a non-empty character scalar.")
  }
  report <- data.frame(
    layer_type = character(),
    layer_name = character(),
    status = character(),
    details = character(),
    n_rows = integer(),
    n_cols = integer(),
    stringsAsFactors = FALSE
  )
  append_report <- function(
    layer_type,
    layer_name,
    status,
    details = "",
    n_rows = NA_integer_,
    n_cols = NA_integer_
  ) {
    report <<- rbind(
      report,
      data.frame(
        layer_type = layer_type,
        layer_name = layer_name,
        status = status,
        details = details,
        n_rows = as.integer(n_rows),
        n_cols = as.integer(n_cols),
        stringsAsFactors = FALSE
      )
    )
  }
  safe_add_assay <- function(layer_name, mat) {
    tryCatch(
      {
        ds$addAssay(layer_name, mat)
        append_report("assay", layer_name, "added", n_rows = nrow(mat), n_cols = ncol(mat))
      },
      error = function(e) {
        append_report(
          "assay", layer_name, "skipped", conditionMessage(e),
          n_rows = if (is.null(dim(mat))) NA_integer_ else nrow(mat),
          n_cols = if (is.null(dim(mat))) NA_integer_ else ncol(mat)
        )
      }
    )
  }
  safe_add_column_metadata <- function(assay_name, metadata_name, values, value_col) {
    layer_name <- paste0(assay_name, "::", metadata_name)
    assay_mat <- ds$assays[[assay_name]]
    if (is.null(assay_mat)) {
      append_report("column_metadata", layer_name, "skipped", "Target assay not present in dslt.")
      return(invisible(NULL))
    }
    if (is.null(values)) {
      append_report("column_metadata", layer_name, "skipped", "Source vector is missing.")
      return(invisible(NULL))
    }
    assay_cols <- colnames(assay_mat)
    if (length(values) != length(assay_cols)) {
      append_report(
        "column_metadata", layer_name, "skipped",
        paste0("Length mismatch: vector=", length(values), ", assay cols=", length(assay_cols), ".")
      )
      return(invisible(NULL))
    }
    md <- data.frame(value = as.numeric(values), row.names = assay_cols, stringsAsFactors = FALSE)
    colnames(md) <- value_col
    tryCatch(
      {
        ds$addColumnMetadata(assay_name, metadata_name, md)
        append_report("column_metadata", layer_name, "added", n_rows = nrow(md), n_cols = ncol(md))
      },
      error = function(e) {
        append_report("column_metadata", layer_name, "skipped", conditionMessage(e), n_rows = nrow(md), n_cols = ncol(md))
      }
    )
  }
  safe_add_graph <- function(assay_name, graph_name, value) {
    layer_name <- paste0(assay_name, "::", graph_name)
    assay_mat <- ds$assays[[assay_name]]
    if (is.null(assay_mat)) {
      append_report("graph", layer_name, "skipped", "Target assay not present in dslt.")
      return(invisible(NULL))
    }
    if (is.null(value)) {
      append_report("graph", layer_name, "skipped", "Source graph matrix is missing.")
      return(invisible(NULL))
    }
    graph_mat <- tryCatch(
      prepare_numeric_matrix(value, graph_name),
      error = function(e) {
        append_report("graph", layer_name, "skipped", conditionMessage(e))
        NULL
      }
    )
    if (is.null(graph_mat)) return(invisible(NULL))
    if (nrow(graph_mat) != ncol(graph_mat)) {
      append_report(
        "graph", layer_name, "skipped",
        paste0("Graph must be square: ", nrow(graph_mat), "x", ncol(graph_mat), "."),
        n_rows = nrow(graph_mat), n_cols = ncol(graph_mat)
      )
      return(invisible(NULL))
    }
    if (!identical(rownames(graph_mat), colnames(graph_mat))) {
      append_report(
        "graph", layer_name, "skipped",
        "Graph rownames/colnames must match in the same order.",
        n_rows = nrow(graph_mat), n_cols = ncol(graph_mat)
      )
      return(invisible(NULL))
    }
    assay_rows <- rownames(assay_mat)
    if (!is.null(assay_rows)) {
      extra_rows <- setdiff(rownames(graph_mat), assay_rows)
      if (length(extra_rows) > 0) {
        append_report(
          "graph", layer_name, "skipped",
          paste0("Graph rows absent from assay: ", paste(utils::head(extra_rows, 10), collapse = ", ")),
          n_rows = nrow(graph_mat), n_cols = ncol(graph_mat)
        )
        return(invisible(NULL))
      }
    }
    tryCatch(
      {
        ds$addGraph(assay = assay_name, input = graph_mat, names = graph_name, copy_assay_on_mismatch = TRUE)
        append_report("graph", layer_name, "added", n_rows = nrow(graph_mat), n_cols = ncol(graph_mat))
      },
      error = function(e) {
        append_report("graph", layer_name, "skipped", conditionMessage(e), n_rows = nrow(graph_mat), n_cols = ncol(graph_mat))
      }
    )
  }

  reference_rows <- rownames(ds$getAssay("cGR", force = TRUE))
  if (is.null(reference_rows)) stop("cGR assay has no rownames; cannot align prepared growth layers.")

  matrix_layer_values <- list(
    gr_obj$counts_raw,
    gr_obj$counts_smooth,
    gr_obj$cnts_smooth_maxN,
    gr_obj$cnts_smooth_sort,
    gr_obj$cnts_smooth_sort_mm,
    gr_obj$counts_smooth_mm,
    gr_obj$gr_raw,
    gr_obj$gr_smooth
  )
  matrix_layer_names <- paste0(
    layer_prefix,
    c(
      "counts_raw",
      "counts_smooth",
      "counts_smooth_maxN",
      "counts_smooth_sort",
      "counts_smooth_sort_mm",
      "counts_smooth_mm",
      "gr_raw",
      "gr_smooth"
    )
  )
  matrix_layers <- stats::setNames(matrix_layer_values, matrix_layer_names)

  for (layer_name in names(matrix_layers)) {
    value <- matrix_layers[[layer_name]]
    if (is.null(value)) {
      append_report("assay", layer_name, "skipped", "Layer missing from gr_obj.")
      next
    }
    mat <- tryCatch(
      prepare_numeric_matrix(value, layer_name),
      error = function(e) {
        append_report("assay", layer_name, "skipped", conditionMessage(e))
        NULL
      }
    )
    if (is.null(mat)) next
    padded <- tryCatch(
      pad_to_reference_rows(mat, reference_rows, layer_name),
      error = function(e) {
        append_report("assay", layer_name, "skipped", conditionMessage(e), n_rows = nrow(mat), n_cols = ncol(mat))
        NULL
      }
    )
    if (is.null(padded)) next
    safe_add_assay(layer_name, padded)
  }

  vector_layer_values <- list(gr_obj$gr_centmass, gr_obj$psds)
  vector_layer_names <- paste0(layer_prefix, c("gr_centmass", "psds"))
  for (idx in seq_along(vector_layer_names)) {
    layer_name <- vector_layer_names[[idx]]
    value <- vector_layer_values[[idx]]
    if (is.null(value)) {
      append_report("assay", layer_name, "skipped", "Vector layer missing from gr_obj.")
      next
    }
    mat <- tryCatch(
      vector_to_single_col_matrix(value, layer_name, col_name = layer_name),
      error = function(e) {
        append_report("assay", layer_name, "skipped", conditionMessage(e))
        NULL
      }
    )
    if (is.null(mat)) next
    padded <- pad_to_reference_rows(mat, reference_rows, layer_name)
    safe_add_assay(layer_name, padded)
  }

  safe_add_column_metadata(paste0(layer_prefix, "counts_raw"), "time_raw", gr_obj$times_raw, "time_raw")
  safe_add_column_metadata(paste0(layer_prefix, "counts_smooth"), "time_smooth", gr_obj$times_smooth, "time_smooth")
  safe_add_column_metadata(paste0(layer_prefix, "counts_smooth_maxN"), "time_smooth", gr_obj$times_smooth, "time_smooth")
  safe_add_column_metadata(paste0(layer_prefix, "counts_smooth_sort"), "time_smooth", gr_obj$times_smooth, "time_smooth")
  safe_add_column_metadata(paste0(layer_prefix, "counts_smooth_sort_mm"), "time_smooth", gr_obj$times_smooth, "time_smooth")
  safe_add_column_metadata(paste0(layer_prefix, "counts_smooth_mm"), "time_smooth", gr_obj$times_smooth, "time_smooth")
  safe_add_column_metadata(paste0(layer_prefix, "gr_raw"), "delta_time_raw", gr_obj$gr_deltatime_raw, "delta_time_raw")
  safe_add_column_metadata(paste0(layer_prefix, "gr_smooth"), "delta_time_smooth", gr_obj$gr_deltatime_smooth, "delta_time_smooth")

  graph_dist_type <- gr_obj$dist_type
  graph_dist_type <- as.character(graph_dist_type)[1]
  if (is.na(graph_dist_type) || !nzchar(graph_dist_type)) graph_dist_type <- dist_type_fallback
  graph_assay_name <- paste0(layer_prefix, graph_dist_type)
  graph_name <- paste0(layer_prefix, "growth_sdtw")
  safe_add_graph(graph_assay_name, graph_name, gr_obj$growth_sdtw)
  sim_graph_name <- paste0(layer_prefix, "growth_similarity")
  safe_add_graph(
    graph_assay_name,
    sim_graph_name,
    value = tryCatch(
      growth_similarity_to_matrix(gr_obj$growth_similarity),
      error = function(e) {
        append_report(
          "graph",
          paste0(graph_assay_name, "::", sim_graph_name),
          "skipped",
          conditionMessage(e)
        )
        NULL
      }
    )
  )

  append_report(
    layer_type = "scalar",
    layer_name = paste0(layer_prefix, "dist_type"),
    status = "skipped",
    details = "Not added: scalar field has no native DatasetLT layer type."
  )

  list(ds = ds, report = report)
}

derive_dslt_output_path <- function(
  input_path,
  overwrite = FALSE,
  suffix = "_preparedGrowth",
  file_suffix = ""
) {
  if (overwrite && !nzchar(file_suffix)) {
    return(input_path)
  }
  input_name <- basename(input_path)
  input_stem <- sub("[.]rds$", "", input_name, ignore.case = TRUE)
  file.path(dirname(input_path), paste0(input_stem, suffix, file_suffix, ".rds"))
}

build_gr_obj <- function(
  cell_line,
  counts_all,
  cgr_file,
  fit_interval = 3,
  test_run = FALSE,
  test_barcodes_per_line = NULL,
  growth_dist_type = "counts_smooth_mm",
  calculate_growth_similarity = TRUE,
  growth_sim_k_per_1k = 75,
  growth_sim_kt_fraction = 0.25,
  growth_sim_gamma = 2
) {
  growth_dist_type <- as.character(growth_dist_type)[1]
  if (is.na(growth_dist_type) || !nzchar(growth_dist_type)) {
    stop("growth_dist_type must be a non-empty character scalar.")
  }
  growth_sim_k_per_1k <- as.numeric(growth_sim_k_per_1k)[1]
  growth_sim_kt_fraction <- as.numeric(growth_sim_kt_fraction)[1]
  growth_sim_gamma <- as.numeric(growth_sim_gamma)[1]
  if (!is.finite(growth_sim_k_per_1k) || growth_sim_k_per_1k <= 0) stop("growth_sim_k_per_1k must be a positive number.")
  if (!is.finite(growth_sim_kt_fraction) || growth_sim_kt_fraction <= 0) stop("growth_sim_kt_fraction must be a positive number.")
  if (!is.finite(growth_sim_gamma) || growth_sim_gamma <= 0) stop("growth_sim_gamma must be a positive number.")

  if (!cell_line %in% names(counts_all)) {
    stop("Cell line ", cell_line, " missing in prepared counts list.")
  }

  ds <- read_dslt_rds(cgr_file)
  cgr <- ds$getAssay("cGR", force = TRUE)
  if (is.null(rownames(cgr))) stop("cGR assay in ", basename(cgr_file), " has no rownames.")

  counts_line <- as.matrix(counts_all[[cell_line]])
  barcodes <- intersect(rownames(cgr), rownames(counts_line))
  if (length(barcodes) == 0) stop("No overlapping barcodes for ", cell_line)

  if (test_run && !is.null(test_barcodes_per_line)) {
    test_barcodes_per_line <- as.integer(test_barcodes_per_line)
    if (!is.na(test_barcodes_per_line) && test_barcodes_per_line > 0) {
      keep_n <- min(length(barcodes), test_barcodes_per_line)
      barcode_strength <- rowSums(counts_line[barcodes, , drop = FALSE], na.rm = TRUE)
      barcodes <- names(sort(barcode_strength, decreasing = TRUE))[seq_len(keep_n)]
      cat("  test_run: keeping", keep_n, "barcodes\n")
    }
  }

  counts_sel <- counts_line[barcodes, , drop = FALSE]
  counts_sel[is.na(counts_sel)] <- 0
  if (any(!is.finite(counts_sel))) stop("Non-finite values in counts for ", cell_line)
  if (any(counts_sel < 0)) stop("Negative values in counts for ", cell_line)
  if (any(colSums(counts_sel) <= 0)) stop("At least one empty sample in ", cell_line)

  counts_sel <- counts_sel[rowSums(counts_sel) > 0, , drop = FALSE]
  barcodes <- rownames(counts_sel)
  metadata <- parse_growth_metadata(colnames(counts_sel))
  if (length(unique(metadata$time)) < 3) stop("Need at least 3 distinct time points for ", cell_line)

  smpls <- rownames(metadata)
  counts_sel <- counts_sel[, smpls, drop = FALSE]

  gr_obj <- get_growth(
    counts = counts_sel,
    brcs = barcodes,
    smpls = smpls,
    metadata = metadata,
    fit_interval = fit_interval
  )
  if (is.null(gr_obj[[growth_dist_type]])) {
    stop("get_growth_dist type '", growth_dist_type, "' is missing from gr_obj.")
  }
  gr_obj <- get_growth_dist(gr_obj, type = growth_dist_type)

  if (isTRUE(calculate_growth_similarity)) {
    if (!exists("adaptiveGaussianSNN", mode = "function")) {
      stop("adaptiveGaussianSNN function is not available in the current R session.")
    }
    k <- as.integer(round(length(barcodes) * growth_sim_k_per_1k / 1000))
    k <- max(2L, min(k, max(2L, length(barcodes) - 1L)))
    kt <- k * growth_sim_kt_fraction
    gr_obj$growth_similarity <- adaptiveGaussianSNN(
      as.dist(gr_obj$growth_sdtw),
      k = k,
      kt = kt,
      gamma = growth_sim_gamma
    )
  }

  summary_row <- data.frame(
    cell_line = cell_line,
    test_run = test_run,
    test_barcodes_per_line = if (is.null(test_barcodes_per_line)) NA_integer_ else as.integer(test_barcodes_per_line),
    n_barcodes_counts = nrow(counts_line),
    n_barcodes_cgr = nrow(cgr),
    n_active_barcodes = length(barcodes),
    n_samples = ncol(counts_sel),
    n_timepoints = length(unique(metadata$time)),
    dist_type = gr_obj$dist_type,
    growth_similarity = isTRUE(calculate_growth_similarity),
    counts_smooth_rows = nrow(gr_obj$counts_smooth),
    counts_smooth_cols = ncol(gr_obj$counts_smooth),
    stringsAsFactors = FALSE
  )

  list(gr_obj = gr_obj, summary = summary_row)
}

to_long_df <- function(mat, x_name = "value", group_name = "condition") {
  out <- data.frame(
    rep(colnames(mat), each = nrow(mat)),
    as.numeric(mat),
    stringsAsFactors = FALSE
  )
  colnames(out) <- c(group_name, x_name)
  out
}

plot_ridge_distribution <- function(
  df,
  title,
  subtitle,
  xlab,
  value_col = "value",
  group_col = "condition",
  ridge_bins = 120L,
  ridge_scale = 1.2,
  ridge_alpha = 0.55
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")
  if (!requireNamespace("ggridges", quietly = TRUE)) stop("Package ggridges is required.")
  if (!all(c(value_col, group_col) %in% colnames(df))) {
    stop("Missing required columns in df: ", value_col, " and/or ", group_col)
  }
  d <- df[is.finite(df[[value_col]]), , drop = FALSE]
  if (!nrow(d)) stop("No finite transformed values to plot.")
  group_levels <- names(table(d[[group_col]]))
  d[[group_col]] <- factor(d[[group_col]], levels = group_levels, ordered = TRUE)

  xr <- range(d[[value_col]])
  baseline_df <- data.frame(
    group = factor(group_levels, levels = group_levels, ordered = TRUE),
    xmin = xr[[1]],
    xmax = xr[[2]],
    stringsAsFactors = FALSE
  )
  cols <- grDevices::hcl.colors(length(group_levels), palette = "Spectral")
  names(cols) <- group_levels

  ggplot2::ggplot(d, ggplot2::aes(x = .data[[value_col]], y = .data[[group_col]])) +
    ggridges::geom_density_ridges(
      ggplot2::aes(fill = .data[[group_col]]),
      stat = "binline",
      bins = as.integer(ridge_bins),
      color = NA,
      alpha = ridge_alpha,
      scale = ridge_scale
    ) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::geom_segment(
      data = baseline_df,
      ggplot2::aes(x = xmin, xend = xmax, y = group, yend = group),
      inherit.aes = FALSE,
      color = "#2F2F2F",
      linewidth = 0.16,
      lineend = "butt"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = "") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "#303030", linewidth = 0.35),
      axis.ticks.x = ggplot2::element_line(color = "#303030", linewidth = 0.35),
      axis.text.x = ggplot2::element_text(color = "#303030"),
      axis.title.x = ggplot2::element_text(color = "#303030"),
      axis.text.y = ggplot2::element_text(size = 6.5),
      plot.title = ggplot2::element_text(size = 9.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8)
    )
}
NN_similarity_feature_enrichment= function(nn_similarity,feature_similarity){
  
  ### cols values - col for cell A values indicate enrichment of growth pattern of cell A in the sc neighborhood of
  # of cells A-X
  ### rows values - row for cell A values indicate enrichment of growth patterns of of cells A-X in the sc neighborhood of cell A
  
  nn_similarity = nn_similarity/rowSums(nn_similarity)
  feature_similarity= get_colFracts(feature_similarity)
  
  nn_feature_enrichment= nn_similarity %*% feature_similarity
  return(nn_feature_enrichment)
  
}
