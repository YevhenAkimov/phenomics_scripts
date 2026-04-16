
# ------------------------------------------------------------------ #
#                         history helper                              #
# ------------------------------------------------------------------ #

record_history <- function(dslt, step, params = list(), note = NULL) {
    hist <- attr(dslt, "pipeline_history")
    if (is.null(hist)) hist <- list()
    hist[[length(hist) + 1L]] <- list(
        step      = step,
        timestamp = Sys.time(),
        params    = params,
        note      = note
    )
    attr(dslt, "pipeline_history") <- hist
    invisible(dslt)
}

print_history <- function(dslt) {
    hist <- attr(dslt, "pipeline_history")
    if (is.null(hist) || !length(hist)) {
        cat("No pipeline history recorded.\n")
        return(invisible(NULL))
    }
    for (i in seq_along(hist)) {
        h <- hist[[i]]
        cat(sprintf(
            "[%d] %s  %s\n",
            i, format(h$timestamp, "%Y-%m-%d %H:%M:%S"), h$step
        ))
        if (length(h$params)) {
            for (nm in names(h$params)) {
                v <- h$params[[nm]]
                if (length(v) > 4) v <- paste0(paste(head(v, 4), collapse = ","), ",...")
                cat(sprintf("      %s = %s\n", nm, paste(v, collapse = ",")))
            }
        }
        if (!is.null(h$note)) cat("      note: ", h$note, "\n", sep = "")
    }
    invisible(hist)
}


# ================================================================== #
#                   (B) WORKFLOW-SPECIFIC FUNCTIONS                   #
# ================================================================== #

# --- collect barcodes from Seurat meta.data -------------------------
read_seurat_barcodes <- function(rds_paths,
                                 barcode_column = "barcodes_in_lineage") {
    rds_paths <- unique(as.character(rds_paths))
    rds_paths <- rds_paths[nzchar(rds_paths)]
    missing <- rds_paths[!file.exists(rds_paths)]
    if (length(missing))
        stop("Seurat RDS files not found: ", paste(missing, collapse = ", "))
    if (!length(rds_paths)) return(character())
    per <- lapply(rds_paths, function(p) {
        obj <- readRDS(p)
        md  <- slot(obj, "meta.data")
        b   <- as.character(md[[barcode_column]])
        b   <- unique(trimws(unlist(strsplit(b[!is.na(b) & nzchar(b)], ",\\s*"))))
        b[!is.na(b) & nzchar(b)]
    })
    Reduce(union, per)
}


# Sets the active barcode set on a DSLT by combining:
#   - barcodes passing an average-fraction threshold (fraction + fraction_control)
#   - barcodes coming from one or more Seurat RDS files
# Defaults reflect lapl_smooth.R / lapl_smoothH160.R / lapl_smoothH180.R.
set_active_barcodes <- function(dslt,
                                seurat_barcodes   = character(),
                                fraction_assays   = c("fraction", "fraction_control"),
                                reference_assay   = "cGR",
                                fraction_threshold = 3e-5,
                                fraction_op       = "mean",   # "mean" or "max"
                                seurat_mode       = "add") {  # "add" or "only"
    ref <- dslt$assays[[reference_assay]]
    if (is.null(ref)) stop("Missing reference assay: ", reference_assay)

    fraction_brcs <- character()
    if (length(fraction_assays) && all(fraction_assays %in% names(dslt$assays))) {
        mats <- lapply(fraction_assays, function(a) dslt$getAssay(a))
        agg <- Reduce("+", mats) / length(mats)
        fraction_brcs <- rownames(agg)[matrixStats::rowMaxs(agg) > fraction_threshold]
    }

    seurat_brcs <- intersect(seurat_barcodes, rownames(ref))

    brcs <- if (identical(seurat_mode, "only") && length(seurat_brcs)) {
        seurat_brcs
    } else {
        union(fraction_brcs, seurat_brcs)
    }
    brcs <- intersect(brcs, rownames(ref))
    dslt$setActiveSamples(brcs)

    dslt <- record_history(dslt, "set_active_barcodes",
        params = list(
            reference_assay = reference_assay,
            fraction_threshold = fraction_threshold,
            seurat_mode = seurat_mode,
            n_fraction = length(fraction_brcs),
            n_seurat = length(seurat_brcs),
            n_active = length(brcs)
        ))
    invisible(dslt)
}


# Add a time-ordered counts matrix (barcodes x samples) as a new assay with
# a numeric `times` vector stored as column metadata.
# Generic: caller is responsible for deriving `times` from sample names.
integrate_longterm_growth_curves <- function(dslt,
                                             counts,
                                             times,
                                             assay_name    = "growth_counts_longterm",
                                             time_metadata = "time") {
    counts <- as.matrix(counts)
    if (length(times) != ncol(counts))
        stop("length(times) must equal ncol(counts)")

    active <- dslt$activeSamples
    missing_brcs <- setdiff(active, rownames(counts))
    if (length(missing_brcs))
        warning(sprintf(
            "integrate_longterm_growth_curves: %d / %d active barcodes not in counts",
            length(missing_brcs), length(active)))

    dslt$addAssay(assay_name, counts)
    time_df <- data.frame(time = times, row.names = colnames(counts),
                          stringsAsFactors = FALSE)
    names(time_df) <- time_metadata
    dslt$addColumnMetadata(assay_name, time_metadata, time_df)

    dslt <- record_history(dslt, "integrate_longterm_growth_curves",
        params = list(assay = assay_name, n_rows = nrow(counts), n_cols = ncol(counts),
                      n_active_missing = length(missing_brcs)))
    invisible(dslt)
}


# ---- utility: named-vector → single-column assay matrix ----
scores_to_assay <- function(scores, row_ids, col_name) {
    out <- matrix(NA_real_, nrow = length(row_ids), ncol = 1,
                  dimnames = list(row_ids, col_name))
    keep <- intersect(row_ids, names(scores))
    out[keep, 1] <- as.numeric(scores[keep])
    out
}


# Fit + smooth growth curves.  Times are read from dslt column metadata
# (written by integrate_longterm_growth_curves), so no sample-name parsing here.
smooth_growth_curves <- function(dslt,
                                 counts_assay      = "growth_counts_longterm",
                                 active_assay_name = "growth_counts_active",
                                 time_metadata     = "time",
                                 fit_interval      = 3,
                                 drop_zero_rows    = FALSE,
                                 assay_names = list(
                                     counts_smooth    = "growth_counts_smooth",
                                     counts_smooth_mm = "growth_counts_smooth_mm",
                                     gr_smooth        = "growth_gr_smooth",
                                     gr_centmass      = "growth_gr_centmass",
                                     psds             = "growth_psds"
                                 )) {
    counts   <- as.matrix(dslt$getAssay(counts_assay, force = TRUE))
    barcodes <- intersect(dslt$activeSamples, rownames(counts))
    counts   <- counts[barcodes, , drop = FALSE]
    counts[is.na(counts)] <- 0
    if (drop_zero_rows) {
        counts <- counts[rowSums(counts) > 0, , drop = FALSE]
    }
    barcodes <- rownames(counts)

    # times from column metadata (set by integrate_longterm_growth_curves)
    time_meta <- dslt$getColumnMetadata(counts_assay, time_metadata)
    smpls     <- colnames(counts)
    metadata  <- data.frame(time = time_meta[smpls, 1], row.names = smpls,
                            stringsAsFactors = FALSE)
    metadata  <- metadata[order(metadata$time), , drop = FALSE]
    smpls     <- rownames(metadata)
    counts    <- counts[, smpls, drop = FALSE]

    gr_obj <- get_growth(counts = counts, brcs = barcodes, smpls = smpls,
                         metadata = metadata, fit_interval = fit_interval)

    # --- push assays ---
    dslt$addAssay(active_assay_name, counts)
    dslt$addColumnMetadata(active_assay_name, "time",
        data.frame(time = metadata[smpls, "time"], row.names = smpls,
                   stringsAsFactors = FALSE))
    dslt$addAssay(assay_names$counts_smooth,    gr_obj$counts_smooth)
    dslt$addAssay(assay_names$counts_smooth_mm, gr_obj$counts_smooth_mm)
    dslt$addAssay(assay_names$gr_smooth,        gr_obj$gr_smooth)

    row_ids <- rownames(gr_obj$counts_smooth_mm)
    dslt$addAssay(assay_names$gr_centmass,
                  scores_to_assay(gr_obj$gr_centmass, row_ids, assay_names$gr_centmass))
    dslt$addAssay(assay_names$psds,
                  scores_to_assay(gr_obj$psds, row_ids, assay_names$psds))

    dslt$addColumnMetadata(assay_names$counts_smooth, "time_smooth",
        data.frame(time_smooth = gr_obj$times_smooth,
                   row.names = colnames(gr_obj$counts_smooth), stringsAsFactors = FALSE))
    dslt$addColumnMetadata(assay_names$counts_smooth_mm, "time_smooth",
        data.frame(time_smooth = gr_obj$times_smooth,
                   row.names = colnames(gr_obj$counts_smooth_mm), stringsAsFactors = FALSE))
    dslt$addColumnMetadata(assay_names$gr_smooth, "delta_time_smooth",
        data.frame(delta_time_smooth = gr_obj$gr_deltatime_smooth,
                   row.names = colnames(gr_obj$gr_smooth), stringsAsFactors = FALSE))

    dslt <- record_history(dslt, "smooth_growth_curves",
        params = list(counts_assay = counts_assay, n_barcodes = length(barcodes),
                      fit_interval = fit_interval))
    invisible(dslt)
}


# Earliness / lateness scores from smoothed curves + centre-of-mass.
compute_growth_scores <- function(dslt,
                                  curves_assay   = "growth_counts_smooth_mm",
                                  centmass_assay = "growth_gr_centmass",
                                  assay_names = list(
                                      earliness = "growth_earliness_score",
                                      lateness  = "growth_lateness_score"
                                  )) {
    curves      <- as.matrix(dslt$getAssay(curves_assay))
    centmass_m  <- as.matrix(dslt$getAssay(centmass_assay))
    gr_centmass <- setNames(centmass_m[, 1], rownames(centmass_m))
    time_meta   <- dslt$getColumnMetadata(curves_assay, "time_smooth")
    times       <- time_meta[colnames(curves), 1]

    earl <- earliness_score(curves = curves, gr_centmass = gr_centmass, times = times)
    late <- lateness_score(curves  = curves, gr_centmass = gr_centmass, times = times)

    row_ids <- rownames(curves)
    dslt$addAssay(assay_names$earliness,
                  scores_to_assay(earl$scores, row_ids, assay_names$earliness))
    dslt$addAssay(assay_names$lateness,
                  scores_to_assay(late$scores, row_ids, assay_names$lateness))

    dslt <- record_history(dslt, "compute_growth_scores",
        params = list(curves_assay = curves_assay, n_barcodes = length(row_ids)))
    invisible(dslt)
}


# Soft-DTW distance matrix + optional SNN graph + MDS embedding.
compute_growth_distances <- function(dslt,
                                     curves_assay         = "growth_counts_smooth_mm",
                                     growth_dist_type     = "counts_smooth_mm",
                                     calculate_similarity = TRUE,
                                     sim_k_per_1k         = 20,
                                     sim_kt_fraction      = 0.15,
                                     sim_gamma            = 1,
                                     graph_names = list(
                                         sdtw       = "growth_sdtw",
                                         similarity = "growth_similarity_snn",
                                         embedding  = "growth_mds2"
                                     )) {
    curves   <- as.matrix(dslt$getAssay(curves_assay))
    barcodes <- rownames(curves)

    gr_obj <- get_growth_dist(list(counts_smooth_mm = curves), type = growth_dist_type)

    dslt$addGraph(assay = curves_assay,
                  input = gr_obj$growth_sdtw, names = graph_names$sdtw)

    if (calculate_similarity && length(barcodes) > 2) {
        k  <- max(2L, min(as.integer(round(length(barcodes) * sim_k_per_1k / 1000)),
                          length(barcodes) - 1L))
        kt <- k * sim_kt_fraction
        cat("k_neighbors =", k, "snn_threshold =", kt, "\n")
        sim <- adaptiveGaussianSNN(as.dist(gr_obj$growth_sdtw),
                                   k_neighbors = k, snn_threshold = kt, gamma = sim_gamma)
        dslt$addGraph(assay = curves_assay,
                      input = growth_similarity_to_matrix(sim),
                      names = graph_names$similarity)
    }

    mds <- tryCatch(stats::cmdscale(as.dist(gr_obj$growth_sdtw), k = 2),
                    error = function(e) NULL)
    if (!is.null(mds)) {
        if (ncol(mds) < 2) mds <- cbind(MDS1 = mds[, 1], MDS2 = 0)
        colnames(mds) <- c("MDS1", "MDS2")
        dslt$addEmbedding(assay = curves_assay, input = mds, names = graph_names$embedding)
    }

    dslt <- record_history(dslt, "compute_growth_distances",
        params = list(curves_assay = curves_assay, n_barcodes = length(barcodes),
                      calculate_similarity = calculate_similarity))
    invisible(dslt)
}


# ================================================================== #
#                   (A) GENERIC FUNCTIONS                             #
# ================================================================== #

# Laplacian adaptive denoising of any assay. Result is stored as a new assay.
# Condition: `assay` exists in dslt and is a numeric matrix.
laplace_denoize <- function(dslt,
                            assay           = "cGR",
                            new_assay       = paste0(assay, "_smoothed"),
                            center          = FALSE,
                            scale           = TRUE,
                            k_neighbors_prop = NULL,   # default: 80/nrow
                            snn_threshold_prop = 0.15,
                            gamma           = 1,
                            lambda          = 2,
                            min_neighbors   = 6) {
    X <- dslt$getAssay(assay)
    if (is.null(k_neighbors_prop)) k_neighbors_prop <- 80 / nrow(X)
    X_scaled <- scale2(X, center = center, scale = scale)

    cat("LAPLACE, k_neighbors_prop = ", k_neighbors_prop, "snn_threshold_prop = ", snn_threshold_prop, "\n")
    smoothed <- LaplacianAdaptiveDenoizing(
        X_scaled,
        lambda             = lambda,
        k_neighbors_prop   = k_neighbors_prop,
        snn_threshold_prop = snn_threshold_prop,
        gamma              = gamma,
        min_neighbors      = min_neighbors
    )
    dslt$addAssay(new_assay, smoothed)

    dslt <- record_history(dslt, "laplace_denoize",
        params = list(
            assay = assay, new_assay = new_assay,
            center = center, scale = scale,
            k_neighbors_prop = k_neighbors_prop,
            snn_threshold_prop = snn_threshold_prop,
            gamma = gamma, lambda = lambda, min_neighbors = min_neighbors
        ))
    invisible(dslt)
}


# Runs kNN / SNN graph construction + louvain clustering + UMAP on a given assay.
# Adds graphs ("similarity_full","similarity_snn","adjacency_snn","distance_matrix",
# "dist_snn") and embeddings ("louvain_clusters","umap") to dslt under `assay`.
knn_analysis <- function(dslt,
                         assay             = "cGR_smoothed",
                         k_neighbors_prop  = NULL,       # default: 80/nrow
                         snn_threshold_prop = 0.15,
                         gamma             = 1,
                         min_neighbors     = 6,
                         umap_n_neighbors  = 100) {
    X <- dslt$getAssay(assay)
    if (is.null(k_neighbors_prop)) k_neighbors_prop <- 80 / nrow(X)

    res <- runKnnAnalysis(X,
        k_neighbors_prop   = k_neighbors_prop,
        snn_threshold_prop = snn_threshold_prop,
        gamma              = gamma,
        min_neighbors      = min_neighbors)

    dslt$addGraph(assay, names = NULL,
        input = res[c("similarity_full", "similarity_snn",
                      "adjacency_snn", "distance_matrix")])
    dslt$addEmbedding(assay, names = "louvain_clusters",
        input = data.frame(louvain_clusters = res[["louvain_clusters"]]))
    dslt$addGraph(assay, names = "dist_snn",
        input = distance_from_similarity_log(dslt$getGraph(assay, "similarity_snn")))
    dslt$addEmbedding(assay, names = "umap",
        input = uwot::umap(as.dist(dslt$getGraph(assay, "dist_snn")),
                           n_neighbors = umap_n_neighbors))

    dslt <- record_history(dslt, "knn_analysis",
        params = list(
            assay = assay,
            k_neighbors_prop = k_neighbors_prop,
            snn_threshold_prop = snn_threshold_prop,
            gamma = gamma, min_neighbors = min_neighbors,
            umap_n_neighbors = umap_n_neighbors
        ))
    invisible(dslt)
}


# Plots UMAP colored by clusters for a given assay. Side-effect only (prints ggplot).
plot_cluster <- function(dslt,
                         assay           = "cGR_smoothed",
                         umap_name       = "umap",
                         cluster_name    = "louvain_clusters",
                         title           = "Louvain Clusters on UMAP",
                         point_size      = 0.5,
                         legend_size     = 3,
                         colors          = colorscheme_prism_mod.) {
    umap     <- dslt$getEmbedding(assay, umap_name)
    clusters <- as.factor(dslt$getEmbedding(assay, cluster_name))
    p <- ggscatter_colored(umap, clusters,
        ggObj = ggplot(),
        dimnamesXYZ = c("UMAP1", "UMAP2", "Cluster"),
        size = point_size,
        gg_theme = theme_umap_legend,
        colors = colors
    ) +
        guides(color = guide_legend(override.aes = list(size = legend_size))) +
        ggtitle(title) + coord_fixed()
    print(p)
    invisible(p)
}


# Archetypal analysis on a (typically smoothed) assay. Stores alpha as embedding,
# and archetype coordinates as column metadata.
perform_archetypal <- function(dslt,
                               assay        = "cGR_smoothed",
                               center       = FALSE,
                               scale        = TRUE,
                               kappa        = 9,
                               nworkers     = 20,
                               nprojected   = 4,
                               alpha_name   = "archetype_alpha",
                               archetypes_name = "archetypes") {
    X <- as.data.frame(dslt$getAssay(assay))
    X <- scale2(X, center = center, scale = scale)

    arches <- find_params_and_perform_arch(X,
        kappa = kappa, nworkers = nworkers, nprojected = nprojected)

    dslt$addEmbedding(assay, names = alpha_name, input = arches$A)
    dslt$addColumnMetadata(assay, archetypes_name, t(arches$BY))

    dslt <- record_history(dslt, "perform_archetypal",
        params = list(assay = assay, kappa = kappa,
                      nworkers = nworkers, nprojected = nprojected))
    invisible(dslt)
}


# Enforce a common barcode set across all assays (intersection of rownames).
# Generic: uses the dslt's own cleanup() method.
cleanup_dslt <- function(dslt) {
    rows_before <- if (length(dslt$assays)) vapply(dslt$assays, nrow, integer(1)) else integer()
    dslt$cleanup()
    rows_after  <- if (length(dslt$assays)) vapply(dslt$assays, nrow, integer(1)) else integer()

    dslt <- record_history(dslt, "cleanup_dslt",
        params = list(
            n_assays = length(rows_before),
            rows_before_min = if (length(rows_before)) min(rows_before) else NA_integer_,
            rows_after_min  = if (length(rows_after))  min(rows_after)  else NA_integer_
        ))
    invisible(dslt)
}


# Bind arbitrary dslt fields (assays and/or embeddings) into a single matrix.
# Each entry in `fields` is a list with:
#   type = "assay"     → name = assay name
#   type = "embedding" → assay = parent assay, name = embedding name
# Barcodes are intersected across all fields. Optionally scale the result.
create_common_phenotype_table <- function(dslt,
                                          fields,
                                          new_assay    = "common_phenotype",
                                          scale        = TRUE,
                                          center       = FALSE,
                                          center_only  = NULL) {
    parts <- lapply(fields, function(f) {
        if (identical(f$type, "embedding")) {
            as.matrix(dslt$getEmbedding(f$assay, f$name))
        } else {
            as.matrix(dslt$getAssay(f$name))
        }
    })

    common_brcs <- Reduce(intersect, lapply(parts, rownames))
    if (!length(common_brcs))
        stop("create_common_phenotype_table: no common barcodes across fields")

    parts <- lapply(parts, function(m) m[common_brcs, , drop = FALSE])
    combined <- do.call(cbind, parts)

    if (!is.null(center_only)) {
        # center only the specified columns, then scale the whole matrix if requested
        cols_to_center <- intersect(center_only, colnames(combined))
        if (length(cols_to_center)) {
            combined[, cols_to_center] <- scale2(
                combined[, cols_to_center, drop = FALSE],
                center = TRUE, scale = FALSE
            )
        }
        if (scale) {
            combined <- scale2(combined, center = FALSE, scale = TRUE)
        }
    } else if (scale || center) {
        combined <- scale2(combined, center = center, scale = scale)
    }

    dslt$addAssay(new_assay, combined)

    field_labels <- vapply(fields, function(f) {
        if (identical(f$type, "embedding")) paste0(f$assay, "/", f$name) else f$name
    }, character(1))

    dslt <- record_history(dslt, "create_common_phenotype_table",
        params = list(
            fields     = paste(field_labels, collapse = ", "),
            new_assay  = new_assay,
            n_barcodes = nrow(combined),
            n_columns  = ncol(combined),
            scale = scale, center = center,
            center_only = if (!is.null(center_only)) paste(center_only, collapse = ", ") else "none"
        ))
    invisible(dslt)
}


# ================================================================== #
#              (C) VISUALIZATION FUNCTIONS                            #
# ================================================================== #
# Derived from visualize_cGRs_*.R scripts.
# All return ggplot objects invisibly; print for side-effect display.

# Multi-treatment UMAP grid — one panel per treatment column.
plot_multi_umap <- function(dslt,
                            assay        = "cGR_smoothed",
                            umap_name    = "umap",
                            point_size   = 0.2,
                            colors       = rev(RdBl_mod3),
                            symm_quant   = 0.95,
                            show_legend  = FALSE,
                            fixed_coord  = TRUE) {
    umap_coords <- dslt$getEmbedding(assay, umap_name)
    values      <- dslt$getAssay(assay)
    gg_base     <- ggplot()
    if (fixed_coord) gg_base <- gg_base + coord_fixed()
    p <- ggscatter_color_multi(umap_coords, values,
        ggObj = gg_base, size = point_size,
        legend.position = if (show_legend) "right" else "none",
        colors = colors, gg_theme = theme_umap, symmQuant = symm_quant)
    invisible(p)
}


# Bubble heatmap of cluster-summarised scores.
plot_cluster_heatmap <- function(dslt,
                                 assay        = "cGR_smoothed",
                                 cluster_name = "louvain_clusters",
                                 cell_line    = "",
                                 colors       = rev(RdBl_mod3),
                                 symm_quant   = 0.95,
                                 size_range   = c(1, 7),
                                 text_angle_x = 90,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE) {
    values   <- as.data.frame(dslt$getAssay(assay))
    clusters <- dslt$getEmbedding(assay, cluster_name)
    summary  <- summarize_columns(values, clusters)

    p <- ggshape_heatmap(summary, abs(summary),
        size_range   = size_range,
        theme_choice = ggplot2::theme_minimal() +
            theme(plot.title = element_text(size = 16, hjust = 0, vjust = 1)),
        value_label  = "Sensitivity",
        size_label   = "Effect size",
        row_label    = "Cluster",
        column_label = "Treatment",
        title        = paste0(cell_line, " Cluster Analysis of cGR scores"),
        cluster_rows = cluster_rows,
        colorscheme  = colors,
        symmQuant    = symm_quant,
        grid.pars = list(grid.size = 0,
            axis.text.x = element_text(size = 8, angle = 0, hjust = -1),
            grid.color = "#f4f4f4", grid.linetype = "solid"),
        cluster_cols = cluster_cols, text.angle.x = text_angle_x
    ) + coord_fixed()
    invisible(p)
}


# Bubble heatmap of archetype profiles.
plot_archetype_heatmap <- function(dslt,
                                   assay           = "cGR_smoothed",
                                   archetypes_name = "archetypes",
                                   cell_line       = "",
                                   colors          = rev(RdBl_mod3),
                                   symm_quant      = 0.95,
                                   size_range      = c(0.5, 4),
                                   text_angle_x    = 90,
                                   tanh_scale      = 3,
                                   cluster_rows    = TRUE,
                                   cluster_cols    = TRUE) {
    arch_raw <- dslt$getColumnMetadata(assay, archetypes_name)
    scaled   <- tanh(arch_raw / tanh_scale)

    p <- ggshape_heatmap(t(scaled), t(abs(scaled)),
        theme_choice = ggplot2::theme_minimal() +
            theme(plot.title = element_text(size = 18, hjust = 0, vjust = 1)),
        shape_values = 22,
        size_range   = size_range,
        value_label  = "Sensitivity",
        size_label   = "Effect size",
        row_label    = "Archetype",
        column_label = "Treatment",
        title        = paste0(cell_line, " Archetypal Analysis of cGR scores"),
        cluster_rows = cluster_rows,
        colorscheme  = colors,
        symmQuant    = symm_quant,
        grid.pars = list(grid.size = 0, grid.color = "#f4f4f4",
                         grid.linetype = "solid"),
        cluster_cols = cluster_cols, text.angle.x = text_angle_x
    ) + coord_fixed()
    invisible(p)
}


# Density contour plot on UMAP embedding.
plot_density_umap <- function(dslt,
                              assay          = "cGR_smoothed",
                              umap_name      = "umap",
                              title          = "Cell density on UMAP space",
                              adjust         = 0.17,
                              point_size     = 0.4,
                              contour_breaks = seq(0.1, 1, length.out = 6),
                              fixed_coord    = TRUE) {
    umap_coords <- dslt$getEmbedding(assay, umap_name)
    gg_base     <- ggplot()
    if (fixed_coord) gg_base <- gg_base + coord_fixed()
    p <- ggdensity_contours(umap_coords,
        weights        = NULL,
        countour_color = "#333333",
        ggObj          = gg_base,
        adjust         = adjust,
        alpha_line     = 0.36,
        alpha_range    = c(0.1, 0.75),
        size           = 0.25,
        point_size     = point_size,
        point_color    = "#e5e5e5",
        colors         = "#323233",
        gg_theme       = theme_umap,
        contour_breaks = contour_breaks
    ) + ggtitle(title)
    invisible(p)
}


# Archetype alpha weights on UMAP (grayscale per archetype).
plot_archetype_umap <- function(dslt,
                                assay       = "cGR_smoothed",
                                umap_name   = "umap",
                                alpha_name  = "archetype_alpha",
                                size_mult   = 0.17,
                                colors      = Greys,
                                fixed_coord = TRUE) {
    umap_coords <- dslt$getEmbedding(assay, umap_name)
    alpha_vals  <- tanh(dslt$getEmbedding(assay, alpha_name))
    gg_base     <- ggplot()
    if (fixed_coord) gg_base <- gg_base + coord_fixed()
    p <- ggscatter_color_multi(umap_coords,
        ggObj = gg_base, size_mult = size_mult,
        values = alpha_vals, colors = colors,
        legend.position = "none", gg_theme = theme_umap)
    invisible(p)
}


# Assemble the full cluster panel (UMAP + density + heatmap | multi-UMAP).
# Returns the cowplot grob; optionally saves to PDF.
assemble_cluster_panel <- function(dslt,
                                   assay        = "cGR_smoothed",
                                   cell_line    = "",
                                   save_path    = NULL,
                                   width        = 18.1,
                                   height       = 12.25,
                                   ...) {
    p_clust   <- plot_cluster(dslt, assay = assay, title = "Lineage clusters", ...)
    p_density <- plot_density_umap(dslt, assay = assay)
    p_heat    <- plot_cluster_heatmap(dslt, assay = assay, cell_line = cell_line)
    p_multi   <- plot_multi_umap(dslt, assay = assay)

    left  <- cowplot::plot_grid(p_clust, p_density, p_heat,
        nrow = 3, align = "v", rel_heights = c(2, 2, 1.85), axis = "tb")
    panel <- cowplot::plot_grid(left, p_multi,
        nrow = 1, align = "vh", rel_widths = c(2, 3.5), axis = "tb")

    if (!is.null(save_path))
        ggsave(panel, filename = save_path, width = width, height = height)
    invisible(panel)
}


# Assemble the full archetype panel (density + arch-heatmap | arch-UMAP).
# Returns the cowplot grob; optionally saves to PDF.
assemble_archetype_panel <- function(dslt,
                                     assay      = "cGR_smoothed",
                                     cell_line  = "",
                                     save_path  = NULL,
                                     width      = 16.9,
                                     height     = 9.5) {
    p_density <- plot_density_umap(dslt, assay = assay)
    p_heat    <- plot_archetype_heatmap(dslt, assay = assay, cell_line = cell_line)
    p_arch    <- plot_archetype_umap(dslt, assay = assay)

    left  <- cowplot::plot_grid(p_density, p_heat,
        nrow = 2, align = "vh", rel_heights = c(2, 1), axis = "tb")
    panel <- cowplot::plot_grid(left, p_arch,
        nrow = 1, align = "vh", rel_widths = c(2.2, 3), axis = "tb")

    if (!is.null(save_path))
        ggsave(panel, filename = save_path, width = width, height = height)
    invisible(panel)
}


# Ridge distribution of growth counts per time-point.
# Converts the counts matrix to long format and delegates to plot_ridge_distribution().
plot_growth_ridges <- function(dslt,
                               counts_assay = "growth_counts_smooth_mm",
                               cell_line    = "",
                               log_transform = TRUE,
                               log_offset    = 1,
                               save_path     = NULL,
                               width         = 8,
                               height        = 5.5,
                               ...) {
    mat <- as.matrix(dslt$getAssay(counts_assay))
    if (log_transform) mat <- log1p(mat + log_offset)
    df <- to_long_df(mat)
    p  <- plot_ridge_distribution(df,
        title    = paste0(cell_line, " Growth counts distribution"),
        subtitle = counts_assay,
        xlab     = if (log_transform) "log(counts + 1)" else "counts",
        ...)
    if (!is.null(save_path))
        ggsave(p, filename = save_path, width = width, height = height)
    invisible(p)
}


# Growth curves colored by a continuous score (e.g. earliness / lateness).
# Wraps the logic from plot_growth_curves_by_score.R into a reusable function.
plot_growth_curves_by_score <- function(dslt,
                                        curve_assay  = "growth_counts_smooth_mm",
                                        color_assay  = "growth_earliness_score",
                                        cell_line    = "",
                                        color_range  = c(-Inf, Inf),
                                        max_curves   = Inf,
                                        seed         = 1L,
                                        line_alpha   = 0.3,
                                        line_width   = 0.25,
                                        save_path    = NULL,
                                        width        = 8.5,
                                        height       = 5.5) {
    curves <- as.matrix(dslt$getAssay(curve_assay, force = TRUE))
    scores_mat <- as.matrix(dslt$getAssay(color_assay, force = TRUE))
    scores <- setNames(as.numeric(scores_mat[, 1]), rownames(scores_mat))

    common <- intersect(rownames(curves), names(scores))
    curves <- curves[common, , drop = FALSE]
    scores <- scores[common]
    ok <- is.finite(scores) &
          scores >= color_range[1] & scores <= color_range[2]
    curves <- curves[ok, , drop = FALSE]
    scores <- scores[ok]

    if (is.finite(max_curves) && nrow(curves) > max_curves) {
        set.seed(seed)
        idx <- sort(sample.int(nrow(curves), max_curves))
        curves <- curves[idx, , drop = FALSE]
        scores <- scores[idx]
    }

    # extract times
    times <- tryCatch({
        md <- dslt$getColumnMetadata(curve_assay, "time_smooth")
        as.numeric(md[colnames(curves), 1])
    }, error = function(e) seq_len(ncol(curves)))

    plot_df <- data.frame(
        barcode = rep(rownames(curves), each = ncol(curves)),
        time    = rep(times, nrow(curves)),
        value   = as.vector(t(curves)),
        score   = rep(scores, each = ncol(curves)),
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(plot_df,
            ggplot2::aes(x = time, y = value, group = barcode, color = score)) +
        ggplot2::geom_line(alpha = line_alpha, linewidth = line_width) +
        ggplot2::scale_color_gradientn(
            colors = c("#2c7bb6", "#ffffbf", "#d7191c"), name = color_assay) +
        ggplot2::labs(
            title    = paste0(cell_line, " Growth curves"),
            subtitle = paste0(curve_assay, " colored by ", color_assay),
            x = "Time", y = curve_assay) +
        ggplot2::theme_classic(base_size = 11)

    if (!is.null(save_path))
        ggsave(p, filename = save_path, width = width, height = height)
    invisible(p)
}


# ================================================================== #
#                   SAVE / EXPORT                                     #
# ================================================================== #

# Save a dslt to disk (optionally with a layerinfo sidecar).
save_dslt <- function(dslt, path, compress = FALSE, save_layerinfo = TRUE) {
    saveRDS(dslt, path, compress = compress)
    if (save_layerinfo) {
        li <- sub("[.]rds$", "_layerinfo.txt", path, ignore.case = TRUE)
        capture.output(dslt$printLayerNames(), file = li)
    }
    cat("Saved:", path, "\n")
    invisible(path)
}
