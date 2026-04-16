# Pipeline: cGR + Growth phenotype analysis

Minimal pipeline that takes clonal growth-rate (cGR) scores and long-term growth counts for a single cell line, denoises them, builds a joint phenotype space, and produces cluster / archetype / growth visualisations.

cGR input scores are produced by [ClonalResponseAnalysis](https://github.com/YevhenAkimov/ClonalResponseAnalysis).

## Pipeline order (v2)

1. Select active barcodes (fraction threshold)
2. Integrate growth counts
3. Laplacian denoising of cGR
4. Cleanup (intersect barcodes across all assays)
4b. Growth smooth + scores + distances (after cleanup)
5. Common phenotype table (cGR_smoothed + growth scores)
6. kNN / clustering / UMAP on common table
7. Archetypal analysis on common table
8. Visualisation

Growth smoothing runs **after** cleanup so that cGR and growth assays share the same barcode set before scores are computed.

## Data preparation

Bundle the inputs for one cell line into a single RDS so the pipeline has no dependency on Seurat objects or scattered file paths:

```r
cell_line_data <- list(
    result_LMM    = readRDS("Phenomics_scores/result_LMM_H160.rds"),
    growth_counts = readRDS("cGR_paper/data/LT_main/all_lines_growth_counts.rds")[["H160"]]
)
saveRDS(cell_line_data, "pipeline/demo_H160.rds")
```

`result_LMM` is the output of the LMM-based cGR estimation from [ClonalResponseAnalysis](https://github.com/YevhenAkimov/ClonalResponseAnalysis). The DatasetLT built from it contains assays `cGR`, `fraction`, `fraction_control` among others. `growth_counts` is a barcodes-by-timepoints count matrix.

## Running the pipeline

```r
source("pipeline/sources_all.R")   # loads all dependencies in one call
```

### 1. Build DatasetLT and select active barcodes

```r
cell_line <- "H160"
bundle    <- readRDS("pipeline/demo_H160.rds")

dslt <- DatasetLT$new(bundle$result_LMM$result)

# keep barcodes whose max average fraction exceeds a threshold
FRACTION_ASSAYS    <- c("fraction", "fraction_control")
FRACTION_THRESHOLD <- 2.5e-5

fraction_mats <- lapply(FRACTION_ASSAYS, function(a) dslt$getAssay(a))
avg_fraction  <- Reduce("+", fraction_mats) / length(fraction_mats)
active_brcs   <- rownames(avg_fraction)[
    matrixStats::rowMaxs(avg_fraction) > FRACTION_THRESHOLD
]

dslt$setActiveSamples(active_brcs)
```

### 2. Integrate growth counts

```r
growth_counts <- bundle$growth_counts
times         <- parse_times_from_sample_names(colnames(growth_counts))

# stores as assay "growth_counts_longterm"
dslt <- integrate_longterm_growth_curves(dslt,
    counts = growth_counts, times = times,
    assay_name = "growth_counts_longterm"
)
```

### 3. Laplacian denoising of cGR

Reads assay `cGR`, writes new assay `cGR_smoothed`.

```r
dslt <- laplace_denoize(dslt,
    assay              = "cGR",
    new_assay          = "cGR_smoothed",
    scale              = TRUE,
    center             = FALSE,
    gamma              = 2,
    lambda             = 2,
    snn_threshold_prop = 0.15,
    min_neighbors      = 6
)
```

### 4. Cleanup

```r
dslt <- cleanup_dslt(dslt)
```

### 4b. Growth smooth, scores, and distances

Runs after cleanup so that growth and cGR share the same barcode set.

```r
# reads "growth_counts_longterm", creates assays:
#   growth_counts_active, growth_counts_smooth, growth_counts_smooth_mm,
#   growth_gr_smooth, growth_gr_centmass, growth_psds
dslt <- smooth_growth_curves(dslt,
    counts_assay   = "growth_counts_longterm",
    fit_interval   = 3,
    drop_zero_rows = FALSE
)

# reads "growth_counts_smooth_mm" + "growth_gr_centmass",
# creates assays: growth_earliness_score, growth_lateness_score
dslt <- compute_growth_scores(dslt,
    curves_assay   = "growth_counts_smooth_mm",
    centmass_assay = "growth_gr_centmass"
)

# reads "growth_counts_smooth_mm",
# adds graph "growth_sdtw", graph "growth_similarity_snn",
# embedding "growth_mds2"
dslt <- compute_growth_distances(dslt,
    curves_assay = "growth_counts_smooth_mm"
)
```

### 5. Common phenotype table

Combines cGR and growth features into a single matrix. `center_only` centers the growth score columns without centering cGR.

```r
dslt <- create_common_phenotype_table(dslt,
    fields = list(
        list(type = "assay",     name = "cGR_smoothed"),
        list(type = "assay",     name = "growth_earliness_score"),
        list(type = "assay",     name = "growth_lateness_score"),
        list(type = "assay",     name = "growth_psds"),
        list(type = "embedding", assay = "growth_counts_smooth_mm", name = "growth_mds2")
    ),
    new_assay   = "common_phenotype",
    scale       = TRUE,
    center      = FALSE,
    center_only = c("growth_earliness_score", "growth_lateness_score", "growth_psds")
)
```

### 6. Joint neighbourhood analysis and archetypal decomposition

kNN analysis adds graphs (`similarity_full`, `similarity_snn`, `adjacency_snn`, `distance_matrix`, `dist_snn`) and embeddings (`louvain_clusters`, `umap`) under the `common_phenotype` assay. Archetypal analysis adds embedding `archetype_alpha` and column metadata `archetypes`.

```r
dslt <- knn_analysis(dslt,
    assay              = "common_phenotype",
    gamma              = 2,
    snn_threshold_prop = 0.15,
    min_neighbors      = 6,
    umap_n_neighbors   = 100
)

dslt <- perform_archetypal(dslt,
    assay    = "common_phenotype",
    kappa    = 9,
    nworkers = 20
)
```

### 7. Visualisation

```r
assemble_cluster_panel(dslt,   assay = "common_phenotype", cell_line = cell_line)
assemble_archetype_panel(dslt, assay = "common_phenotype", cell_line = cell_line)

plot_growth_ridges(dslt, cell_line = cell_line)
plot_growth_curves_by_score(dslt, cell_line = cell_line, color_assay = "growth_earliness_score")
plot_growth_curves_by_score(dslt, cell_line = cell_line, color_assay = "growth_lateness_score")
```

## File overview

| File | Purpose |
|---|---|
| `sources_all.R` | Single `source()` entry-point that loads all dependencies |
| `pipeline_functions.R` | Step functions (generic and workflow-specific) |
| `pipeline_showcase.R` | Complete runnable script for one cell line |
| `pipeline_orchestrator.R` | Full multi-sample orchestrator (v1) |
| `pipeline_orchestrator_v2.R` | Multi-sample orchestrator with joint phenotype table (v2) |

## Dependencies

R packages: `R6`, `Matrix`, `matrixStats`, `uwot`, `ggplot2`, `cowplot`

GitHub sources (loaded automatically by `sources_all.R`):

- [graphics-R](https://github.com/YevhenAkimov/graphics-R) — plotting helpers and colour palettes
- [DatasetLT](https://github.com/YevhenAkimov/DatasetLT) — lightweight multi-assay container
- [phenomics_scripts](https://github.com/YevhenAkimov/phenomics_scripts) — Laplacian denoising, kNN/SNN, archetypal analysis
- [general_purpose_R](https://github.com/YevhenAkimov/general_purpose_R) — utility functions
- [ClonalResponseAnalysis](https://github.com/YevhenAkimov/ClonalResponseAnalysis) — produces the cGR input scores
