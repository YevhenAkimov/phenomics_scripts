
library(readr)     
library(dplyr)    
library(purrr)    
library(tidyr)  
library(igraph)

#' Scan ICA dimensionalities with ICtest::FOBIasymp and plot p-values
#'
#' For each k (number of non-Gaussian components under H0), this function
#' runs ICtest::FOBIasymp(), collects the test statistic & p-value, and
#' reports results **in terms of the number of Gaussian components g = p - k**,
#' including the null and alternative hypotheses.
#' (FOBIasymp tests if there are p - k Gaussian components; alt: fewer than that). 
#' See ICtest manual FOBIasymp description. [oai_citation:CRAN](https://cran.r-project.org/web/packages/ICtest/ICtest.pdf)
#'
#' @param X           Numeric matrix, rows = observations, cols = variables.
#' @param k_seq       Integers of k to test (non-Gaussian components). Default 0:(p-1).
#' @param type        "S1","S2","S3" (FOBIasymp test type). Default "S2".
#' @param model       "ICA" (default) or "NGCA".
#' @param alpha       Significance level. Default 0.05.
#' @param log10_p     Log-scale p-values on the y-axis. Default TRUE.
#' @param show_plot   Print the plot. Default TRUE.
#' @param verbose     Print progress. Default TRUE.
#'
#' @return list with:
#'   \item{table}{data.frame with k, g = p-k, statistic, p_value, H0, H1, ...}
#'   \item{k_hat}{largest k with p >= alpha}
#'   \item{g_hat}{corresponding number of Gaussian comps = p - k_hat}
#'   \item{alpha}{alpha}
#'   \item{plot}{ggplot object (or NULL)}
#'
#' Scan ICA dimensionalities with ICtest::FOBIasymp and plot p-values
#'
#' For each k (number of non-Gaussian components under H0), this function
#' runs ICtest::FOBIasymp(), collects the test statistic & p-value, and
#' reports results **in terms of the number of Gaussian components g = p - k**,
#' including the null and alternative hypotheses per the ICtest manual
#' (FOBIasymp tests if there are p - k Gaussian components; alt: fewer than that).
#'
#' @param X           Numeric matrix, rows = observations, cols = variables.
#' @param k_seq       Integers of k to test (non-Gaussian components). Default 0:(p-1).
#' @param type        "S1","S2","S3" (FOBIasymp test type). Default "S2".
#' @param model       "ICA" (default) or "NGCA".
#' @param alpha       Significance level. Default 0.05.
#' @param log10_p     Log-scale p-values on the y-axis. Default TRUE.
#' @param show_plot   Print the plot. Default TRUE.
#' @param verbose     Print progress. Default TRUE.
#'
#' @return list with:
#'   \item{table}{data.frame with k, g = p-k, statistic, p_value, H0, H1, ...}
#'   \item{k_hat}{largest k with p >= alpha; if none, returns max(k_seq)}
#'   \item{g_hat}{corresponding number of Gaussian comps = p - k_hat}
#'   \item{alpha}{alpha}
#'   \item{plot}{ggplot object (or NULL)}
#'
scan_ic_k <- function(X,
                      k_seq     = NULL,
                      type      = c("S1", "S2", "S3"),
                      model     = c("ICA", "NGCA"),
                      alpha     = 0.05,
                      log10_p   = TRUE,
                      show_plot = TRUE,
                      verbose   = TRUE) {
  
  if (!requireNamespace("ICtest", quietly = TRUE)) {
    stop("Please install 'ICtest' first: install.packages('ICtest')")
  }
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (any(!is.finite(X))) stop("X contains non-finite values.")
  
  p <- ncol(X)
  if (is.null(k_seq)) k_seq <- 0:(p - 1)
  
  type  <- match.arg(type)
  model <- match.arg(model)
  
  out_list <- vector("list", length(k_seq))
  
  if (verbose) {
    message("Running FOBIasymp() for k in {", paste(k_seq, collapse = ", "), "} ...")
  }
  
  for (i in seq_along(k_seq)) {
    k <- k_seq[i]
    g <- p - k  # number of Gaussian components under H0
    H0 <- sprintf("H0: there are %d Gaussian components (k = %d non-Gaussian)", g, k)
    H1 <- sprintf("H1: there are fewer than %d Gaussian components", g)
    
    fit <- tryCatch(
      ICtest::FOBIasymp(X, k = k, type = type, model = model),
      error = function(e) e
    )
    
    if (inherits(fit, "error")) {
      if (verbose) message("k = ", k, " failed: ", fit$message)
      out_list[[i]] <- data.frame(
        k         = k,
        g         = g,
        statistic = NA_real_,
        p_value   = NA_real_,
        method    = NA_character_,
        data_name = NA_character_,
        H0        = H0,
        H1        = H1,
        stringsAsFactors = FALSE
      )
    } else {
      out_list[[i]] <- data.frame(
        k         = k,
        g         = g,
        statistic = as.numeric(fit$statistic),
        p_value   = as.numeric(fit$p.value),
        method    = if (!is.null(fit$method)) as.character(fit$method) else NA_character_,
        data_name = if (!is.null(fit$data.name)) as.character(fit$data.name) else NA_character_,
        H0        = H0,
        H1        = H1,
        stringsAsFactors = FALSE
      )
    }
  }
  
  tab <- do.call(rbind, out_list)
  
  # ---- k_hat rule (fixed):
  # "largest k with p >= alpha"; if none (e.g., all p == 0), choose max(k_seq)
  ok <- which(!is.na(tab$p_value) & tab$p_value >= alpha)
  if (length(ok) == 0) {
    k_hat <- max(tab$k, na.rm = TRUE)
  } else {
    k_hat <- max(tab$k[ok], na.rm = TRUE)
  }
  g_hat <- p - k_hat
  
  # For plotting on log-scale, avoid log10(0)
  eps <- .Machine$double.xmin
  tab$plot_p <- pmax(tab$p_value, eps)
  
  plt <- NULL
  if (show_plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not installed; using base graphics.")
      y <- if (log10_p) log10(tab$plot_p) else tab$p_value
      ylab <- if (log10_p) expression(log[10](pvalue)) else "p-value"
      plot(tab$g, y, type = "b", pch = 19,
           xlab = "Number of Gaussian components (g = p - k)",
           ylab = ylab,
           main = sprintf("ICtest::FOBIasymp scan (type = %s, model = %s)", type, model))
      abline(h = if (log10_p) log10(alpha) else alpha, col = "red", lty = 2)
      abline(v = g_hat, col = "blue", lty = 3)
    } else {
      ggplot2::theme_set(ggplot2::theme_minimal())
      gg <- ggplot2::ggplot(tab, ggplot2::aes(x = g, y = if (log10_p) plot_p else p_value)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = alpha, linetype = "dashed", colour = "red") +
        ggplot2::geom_vline(xintercept = g_hat, linetype = "dotted", colour = "blue") +
        ggplot2::labs(
          x = "Number of Gaussian components (g = p - k)",
          y = if (log10_p) expression(log[10](p~value)) else "p-value",
          title = sprintf("FOBIasymp scan (type = %s, model = %s)", type, model),
          subtitle = sprintf("alpha = %.3f, chosen g_hat = %d (k_hat = %d)", alpha, g_hat, k_hat)
        )
      
      if (log10_p) {
        gg <- gg + ggplot2::scale_y_log10()
      }
      
      print(gg)
      plt <- gg
    }
  }
  
  invisible(list(
    table  = tab,
    k_hat  = k_hat,
    g_hat  = g_hat,
    alpha  = alpha,
    plot   = plt
  ))
}

distance_from_similarity= function(similarity) { (1 / (MinMax(similarity,0.001,0.999)) - 1) }
distance_from_similarity_log= function(similarity) { (1 / log(MinMax(similarity,0.001,0.999))- 1) } 

                                              
distance_form_similarity= function(similarity) { (1 / (MinMax(similarity,0.001,0.999)) - 1) } ## kept for backward compatibility
distance_form_similarity_log= function(similarity) { (1 / log(MinMax(similarity,0.001,0.999))- 1) } ## for backward compatibility

                                              

AdaptiveKernelDenoizing=function(data,use_ica=T,n.comp=NULL, k_neighbors_prop=0.035, snn_threshold_prop=0.2, gamma=3, make_symmetric_snn=TRUE, min_neighbors=10,center=F,scale=T){
  scaled_inp=scale2(data,center=center,scale=scale)
  
  if (is.null(n.comp)){
    ks=scan_ic_k(scaled_inp)
    ks$plot
    n.comp=ks$k_hat
  }
  if (use_ica){
    inp_data=fastICA::fastICA(scaled_inp, n.comp = n.comp)$S
  } else {
    inp_data=scaled_inp
  }
  
  sim=adaptiveGaussianSNN( inp_data,
                           k_neighbors_prop=k_neighbors_prop,
                           snn_threshold_prop=snn_threshold_prop,
                           gamma=gamma,
                           make_symmetric_snn = TRUE,
                           min_neighbors = min_neighbors) 
  
  unit_similarity=sim$similarity_snn/rowSums(sim$similarity_snn)
  smoothed=unit_similarity %*% scaled_inp 
  
  return(smoothed)
}


runKnnAnalysis=function(data,k_neighbors_prop=0.035,snn_threshold_prop=0.2,gamma=3,make_symmetric_snn=T,graph_snn=T,mode="undirected",weighted=T,calc_dists=F,find_louvain_clusts=T,min_neighbors=8,local_scale_quantile=0.9){
  k_neighbors=k_neighbors_prop*nrow(data)
  snn_threshold=snn_threshold_prop*k_neighbors
  mode=match.arg(mode, c("undirected", "directed")) #
  similObj=adaptiveGaussianSNN(data,k_neighbors,snn_threshold,gamma=gamma,make_symmetric_snn=make_symmetric_snn,min_neighbors=min_neighbors,local_scale_quantile=local_scale_quantile) 
  
  if (graph_snn) {
    similObj=c(similObj,getGraphKNN(similObj$similarity_snn,mode=mode,weighted=weighted,calc_dists=calc_dists,findLouClusts=find_louvain_clusts))
  }  else  {
    similObj=c(similObj,getGraphKNN(similObj$similarity,mode=mode,weighted=weighted,calc_dists=calc_dists,findLouClusts=find_louvain_clusts)) }
  return(similObj)
  
}


getGraphKNN <- function(similarity,mode="undirected",weighted=T,calc_dists=T,findLouClusts=T,resolution=1) {
  
  symmetry_check <- all.equal(similarity, t(similarity))
  cat( "getGraphKNN input is symmetric: ",symmetry_check,"\n")
  
  any_na <- any(is.na(similarity))
  any_nan <- any(is.nan(similarity))
  any_inf <- any(is.infinite(similarity))
  
  
  graph  =  igraph::graph_from_adjacency_matrix(similarity, mode=mode, weighted=weighted )
  
  if (!igraph::is_connected(graph)) {
    message("WARNING, graph is not connected")
  }

  if (calc_dists) {
    
    graph_distances  =  igraph::distances( igraph::graph_from_adjacency_matrix(distance_form_similarity(similarity), mode=mode, weighted=weighted ))
  } else {
    graph_distances=NULL
  }
  
  
  if (findLouClusts) {
    louvain_clusters = igraph::cluster_louvain(graph,resolution=resolution)$membership
    names(louvain_clusters) = rownames(similarity)
  } else {
    louvain_clusters=NULL
  }
  
  
  return(list(graph=graph,graph_distances=graph_distances,louvain_clusters=louvain_clusters))
}


#' Adaptive Gaussian Similarity on a Shared Nearest Neighbour (SNN) Graph
#'
#' @description
#' Builds a Shared Nearest Neighbour (SNN) graph from a distance matrix (or raw
#' feature matrix) and computes an **adaptive** Gaussian similarity kernel with
#' **point-specific (local) bandwidths**. For each node \eqn{i}, the local
#' bandwidth \eqn{\sigma_i} is the \code{local_scale_quantile} quantile (default
#' 0.9) of its SNN neighbour distances. The full similarity can be symmetrised
#' and optionally **sparsified to SNN edges**.
#'
#' You may specify SNN parameters **either** absolutely via
#' \code{k_neighbors} and \code{snn_threshold}, **or** proportionally via
#' \code{k_neighbors_prop} and \code{snn_threshold_prop}:
#' \itemize{
#'   \item If proportions are used, we set
#'     \code{k_neighbors = round(n * k_neighbors_prop)} (clamped to \eqn{\ge 2})
#'     and \code{snn_threshold = round(k_neighbors * snn_threshold_prop)}
#'     (clamped to \eqn{\ge 1}). Here \eqn{n} is the number of rows (samples).
#'   \item **Exactly one** of these pairs may be given; giving both pairs is an error.
#' }
#'
#' @section Methods (scientific insight):
#' Let \eqn{D} be the pairwise distance matrix (raw features are z-scored then
#' \eqn{D} is computed via \code{parallelDist::parDist}). We form an SNN graph
#' using \code{k_neighbors} (size of each k-NN list) and \code{snn_threshold}
#' (minimum number of shared neighbours to keep an edge). For each point
#' \eqn{i}, we define a local scale \eqn{\sigma_i} as the chosen quantile of its
#' SNN neighbour distances. We then compute the adaptive kernel
#' \deqn{
#'   K_{ij} = \exp\!\left( - \frac{d_{ij}^2}{\gamma (\sigma_i^2 + \sigma_j^2)} \right),
#' }
#' where \eqn{\gamma} is a user-supplied global scale. Optionally, a *local
#' connectivity subtraction* is first applied by subtracting, from each row of
#' \eqn{D}, its smallest non-zero distance to enforce a minimum local radius.
#'
#' The function returns both the **full** adaptive similarity
#' (\code{similarity_full}) and the **SNN-sparsified** version
#' (\code{similarity_snn}), plus the logical adjacency and intermediates.
#'
#' @param data A numeric matrix/data.frame (rows = samples) or a \code{dist}
#'   object. If a matrix/data.frame is supplied, it is z-scored before computing
#'   distances.
#' @param k_neighbors Integer (absolute). Size of each k-NN list for building
#'   the SNN graph. **Use either this pair or the proportional pair.**
#' @param snn_threshold Integer (absolute). Minimum number of shared nearest
#'   neighbours (\eqn{kt}) to keep an edge. **Use either this pair or the
#'   proportional pair.**
#' @param k_neighbors_prop Numeric in (0,1]. Proportion of the dataset size to
#'   set \code{k_neighbors = round(n * k_neighbors_prop)}, clamped to \eqn{\ge 2}.
#'   **Use either this pair or the absolute pair.**
#' @param snn_threshold_prop Numeric in (0,1]. Proportion of \code{k_neighbors}
#'   to set \code{snn_threshold = round(k_neighbors * snn_threshold_prop)},
#'   clamped to \eqn{\ge 1}. **Use either this pair or the absolute pair.**
#' @param gamma Numeric > 0. Global scaling factor for the adaptive kernel.
#' @param make_symmetric_snn Logical. If \code{TRUE} (default), enforce symmetry
#'   of the full similarity and of the SNN-sparsified similarity.
#' @param min_neighbors Integer or \code{NULL}. Minimum neighbour count to
#'   enforce per point; if \code{NULL}, defaults to \code{max(round(snn_threshold/5), 2)}.
#'   Points with fewer SNN neighbours are padded using an auxiliary kNN query.

#' @param local_scale_quantile Numeric in (0,1). Quantile of SNN distances used
#'   to estimate \eqn{\sigma_i} (default 0.9).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{raw_data} — the original input.
#'   \item \code{similarity_full} — full adaptive similarity matrix.
#'   \item \code{similarity_snn} — SNN-sparsified adaptive similarity.
#'   \item \code{adjacency_snn} — logical adjacency matrix for SNN edges.
#'   \item \code{distance_matrix} — (optionally corrected) distance matrix.
#'   \item \code{snn} — the \code{dbscan::sNN} result.
#'   \item \code{params} — the effective parameters used.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#'
#' # 1) Absolute specification (k, kt):
#' out_abs <- adaptiveGaussianSNN(
#'   data = X,
#'   k_neighbors = 30,
#'   snn_threshold = 10,
#'   gamma = 2
#' )
#'
#' # 2) Proportional specification (k_prop, kt_prop):
#' out_prop <- adaptiveGaussianSNN(
#'   data = X,
#'   k_neighbors_prop = 0.3,    # k = round(100 * 0.3) = 30
#'   snn_threshold_prop = 0.33, # kt = round(30 * 0.33) = 10
#'   gamma = 2
#' )
#'
#' image(out_abs$similarity_snn)
#' image(out_prop$similarity_snn)
#' }
#'
#' @author Your Name
#' @export
#'
#' @importFrom matrixStats rowMins rowQuantiles
#' @importFrom dbscan sNN kNN
#' @importFrom parallelDist parDist
#'
adaptiveGaussianSNN <- function(data,
                                k_neighbors = NULL,
                                snn_threshold = NULL,
                                k_neighbors_prop = NULL,
                                snn_threshold_prop = NULL,
                                gamma = 2,
                                make_symmetric_snn = TRUE,
                              
                                min_neighbors = 8,
                                
                                local_scale_quantile = 0.9) {
  
  # ---- Packages ----
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Package 'matrixStats' is required but not installed.")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required but not installed.")
  }
  if (!requireNamespace("parallelDist", quietly = TRUE)) {
    stop("Package 'parallelDist' is required but not installed.")
  }
  
  # ---- Basic checks ----
  if (!is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("'gamma' must be a positive scalar.")
  }
  if (!is.numeric(local_scale_quantile) ||
      local_scale_quantile <= 0 || local_scale_quantile >= 1) {
    stop("'local_scale_quantile' must be in (0,1).")
  }
  
  # ---- Resolve (k, kt): absolute vs proportional ----
  abs_pair  <- !is.null(k_neighbors) || !is.null(snn_threshold)
  prop_pair <- !is.null(k_neighbors_prop) || !is.null(snn_threshold_prop)
  
  if (abs_pair && prop_pair) {
    stop("Provide either (k_neighbors & snn_threshold) OR (k_neighbors_prop & snn_threshold_prop), not both.")
  }
  
  # get n (#samples)
  n <- if (inherits(data, "dist")) {
    attr(data, "Size")
  } else {
    nrow(data)
  }
  
  if (prop_pair) {
    if (is.null(k_neighbors_prop) || is.null(snn_threshold_prop)) {
      stop("When using proportional specification, both 'k_neighbors_prop' and 'snn_threshold_prop' must be provided.")
    }
    if (!is.numeric(k_neighbors_prop) || k_neighbors_prop <= 0 || k_neighbors_prop > 1) {
      stop("'k_neighbors_prop' must be in (0,1].")
    }
    if (!is.numeric(snn_threshold_prop) || snn_threshold_prop <= 0 || snn_threshold_prop > 1) {
      stop("'snn_threshold_prop' must be in (0,1].")
    }
    
    k_neighbors  <- max(2L, as.integer(round(n * k_neighbors_prop)))
    snn_threshold <- max(1L, as.integer(round(k_neighbors * snn_threshold_prop)))
    
  } else if (abs_pair) {
    if (is.null(k_neighbors) || is.null(snn_threshold)) {
      stop("When using absolute specification, both 'k_neighbors' and 'snn_threshold' must be provided.")
    }
    if (!is.numeric(k_neighbors) || k_neighbors <= 1) {
      stop("'k_neighbors' must be an integer > 1.")
    }
    if (!is.numeric(snn_threshold) || snn_threshold < 1) {
      stop("'snn_threshold' must be an integer >= 1.")
    }
  } else {
    stop("You must provide either (k_neighbors & snn_threshold) or (k_neighbors_prop & snn_threshold_prop).")
  }
  
  if (snn_threshold > k_neighbors) {
    stop("'snn_threshold' cannot exceed 'k_neighbors'.")
  }
  
  if (is.null(min_neighbors)) {
    min_neighbors <- max(round(snn_threshold / 5), 4)
  }
  if (!is.numeric(min_neighbors) || min_neighbors < 1) {
    stop("'min_neighbors' must be a positive integer.")
  }
  
  raw_data <- data
  
  # ---- Distances ----
  if (!inherits(data, "dist")) {
    data <- scale(as.matrix(data))  # z-score
  }
  
  if (!inherits(data, "dist")) {
    datadist <- parallelDist::parDist(data)
  } else {
    datadist <- data
  }
  
  distance_matrix <- as.matrix(datadist)
  

  
  # ---- Sanity checks ----
  if (any(is.na(as.matrix(datadist)))) {
    stop("Some points have NA in distance matrix. Please check the input data and parameters.")
  }
  
  # ---- SNN construction ----
  snn <- dbscan::sNN(datadist, k = k_neighbors, kt = snn_threshold, jp = FALSE)
  snn_non_na <- rowSums(!is.na(snn[["id"]]))
  message("Median SNNs: ", median(snn_non_na))
  # Ensure at least min_neighbors neighbours
  if (min(snn_non_na) < min_neighbors) {
    warning("Some points have fewer than 'min_neighbors' neighbours; padding using kNN.")
    knn <- dbscan::kNN(datadist, k = min_neighbors)
    
    for (i in seq_along(snn_non_na)) {
      if (snn_non_na[i] < min_neighbors) {
        current_snn_ids  <- snn[["id"]][i, ]
        current_knn_ids  <- knn[["id"]][i, ]
        
        map_idx <- match(current_knn_ids, current_snn_ids)
        
        from_idx <- which(is.na(map_idx))
        to_idx   <- which(is.na(current_snn_ids))[seq_along(from_idx)]
        
        current_snn_ids[to_idx] <- current_knn_ids[from_idx]
        
        snn[["id"]][i, ]  <- current_snn_ids
        knn[["id"]][i, ]  <- current_knn_ids
        snn[["dist"]][i, to_idx] <- knn[["dist"]][i, from_idx]
      }
    }
  }
  
  # ---- Local scales from SNN distances ----

  sigma <- as.numeric(
    matrixStats::rowQuantiles(snn[["dist"]],
                              probs = local_scale_quantile,
                              na.rm = TRUE)
  )
  if (any(sigma == 0))        stop("Some points have zero local scale. Try increasing 'min_neighbors'.")
  if (any(!is.finite(sigma))) stop("Some points have non-finite local scales. Increase the k vs. kt gap.")
  
  sigma2       <- sigma^2
  sum_sigma_sq <- outer(sigma2, sigma2, "+")
  
  # ---- Full adaptive Gaussian similarity ----
  similarity_full <- exp( - (distance_matrix^2) / (gamma * sum_sigma_sq) )
  if (any(is.na(similarity_full))) stop("Some similarity values are NA.")
  
  if (make_symmetric_snn) {
    similarity_full <- 0.5 * (similarity_full + t(similarity_full))
  }
  
  # ---- SNN-sparsified similarity ----
  similarity_snn <- keep_inds_in_square(snn[["id"]], similarity_full)
  adjacency_snn  <- similarity_snn != 0
  
  if (make_symmetric_snn) {
    adjacency_snn  <- adjacency_snn | t(adjacency_snn)
    similarity_snn <- adjacency_snn * similarity_full
  } else {
    similarity_snn <- adjacency_snn * similarity_full
  }
  
  diag(similarity_snn) <- 1
  diag(adjacency_snn)  <- TRUE
  
  message("similarity_snn is symmetric: ", isSymmetric(similarity_snn))
  message("similarity_full is symmetric: ", isSymmetric(similarity_full))
  
  return(list(
    raw_data         = raw_data,
    similarity_full  = similarity_full,
    similarity_snn   = similarity_snn,
    adjacency_snn    = 1*adjacency_snn,
    distance_matrix  = distance_matrix,
    snn              = snn,
    params = list(
      k_neighbors                 = k_neighbors,
      snn_threshold               = snn_threshold,
      k_neighbors_prop            = k_neighbors_prop,
      snn_threshold_prop          = snn_threshold_prop,
      gamma                       = gamma,
      make_symmetric_snn          = make_symmetric_snn,
      min_neighbors               = min_neighbors,
  
      local_scale_quantile        = local_scale_quantile
    )
  ))
}

#' Keep Only SNN Indices in a Square Matrix
#'
#' @description
#' Utility that zeros out all entries in `square_matrix` except those pointed to
#' by the SNN index matrix (typically `snn[['id']]` from `dbscan::sNN`).
#'
#' @param snn_id Integer matrix of indices; rows = points, columns = neighbour IDs.
#' @param square_matrix Square numeric matrix matching the same ordering.
#' @return A square matrix with non-SNN entries set to 0.
#' @keywords internal
#' @noRd
keep_inds_in_square <- function(snn_id, square_matrix) {
  n <- nrow(snn_id)
  out <- matrix(0, nrow(square_matrix), ncol(square_matrix))
  rownames(out) <- rownames(square_matrix)
  colnames(out) <- colnames(square_matrix)
  
  for (i in seq_len(n)) {
    idx <- snn_id[i, ]
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      out[i, idx] <- square_matrix[i, idx]
    }
  }
  out
}


scale2=function(data,center=T,scale=T){
  data=scale(data,center=center,scale=scale)
  attr(data, "scaled:center") <- NULL
  attr(data, "scaled:scale") <- NULL
  return(data)
}


find_params_and_perform_arch= function(input,kappa,nworkers,nprojected=2) { 
  require(archetypal)
  opt_par <- find_pcha_optimal_parameters(
    df         = as.data.frame(input),
    kappas     = kappa,
    #mup2=mup2,mup1=mup1,mdown1=mdown1,mdown2=mdown2,
    nworkers= nworkers
    
  )
  
  aa <- archetypal(
    df           =  as.data.frame(input),
    kappas       = kappa,
    nprojected=nprojected,
    initialrows  = opt_par$sol_initial,
    muAup        = opt_par$mu_up_opt, muAdown = opt_par$mu_down_opt,
    muBup        = opt_par$mu_up_opt, muBdown = opt_par$mu_down_opt,
    verbose      = FALSE,
    nworkers=nworkers
  )
  colnames(aa$A)=paste0("Archtype ", 1:ncol(aa$A) )
  rownames(aa$BY)= paste0("Archtype ", 1:nrow(aa$BY) )
  return(aa)
}



which.colMaxs <- function(mat) {
  apply(mat, 2, which.max)
}

which.colMins <- function(mat) {
  apply(mat, 2, which.min)
}

select_data_rows <- function(data, thr = NULL, select_ratio = 0, ranked = FALSE) {
  # Ensure data is a matrix for consistent processing
  data <- as.matrix(data)
  
  # Input Validation
  if (!is.logical(ranked) || length(ranked) != 1) {
    stop("'ranked' must be a single logical value (TRUE or FALSE).")
  }
  
  if (!ranked && is.null(thr)) {
    stop("'thr' must be provided when 'ranked' is FALSE.")
  }
  
  if (ranked && (!is.numeric(thr) || length(thr) != 1)) {
    stop("'thr' must be a single numeric value when 'ranked' is TRUE.")
  }
  
  if (!is.numeric(select_ratio) || length(select_ratio) != 1 || select_ratio < 0 || select_ratio > 1) {
    stop("'select_ratio' must be a single numeric value between 0 and 1.")
  }
  
  # Initialize logical vector for row selection
  logi <- rep(FALSE, nrow(data))
  
  if (ranked) {
    # Determine the number of top ranks to consider
    top_n <- if (thr < 1) {
      ceiling(nrow(data) * thr)
    } else {
      min(thr, nrow(data))
    }
    
    # Rank the data in descending order for each column
    ranked_data <- apply(-data, 2, rank, ties.method = "min")
    
    # Create a logical matrix where TRUE indicates the element is within the top_n ranks
    logi_matrix <- ranked_data <= top_n
    
    # Calculate the proportion of TRUEs in each row
    prop_true <- rowSums(logi_matrix) / ncol(data)
    
  } else {
    # Create a logical matrix where TRUE indicates the element exceeds the threshold
    logi_matrix <- data > thr
    
    # Calculate the proportion of TRUEs in each row
    prop_true <- rowSums(logi_matrix) / ncol(data)
  }
  
  # Determine row selection based on select_ratio
  if (select_ratio == 0) {
    # Select rows where any element meets the condition
    logi <- prop_true > 0
  } else if (select_ratio == 1) {
    # Select rows where all elements meet the condition
    logi <- prop_true == 1
  } else {
    # Select rows where the proportion of elements meeting the condition is >= select_ratio
    logi <- prop_true >= select_ratio
  }
  
  return(logi)
}


norm_func2=function(data,type="mean",fractions=F,multipl="no") {
  colsums=colSums(data)
  if (is.numeric(multipl)==T) {
    fractions=t(t(data)/colsums)
    data=t(t(fractions)*multipl)
    return(data)
  }
  
  colsums=colSums(data,na.rm=T)
  if (is.numeric(multipl)!=T) {
    if (type=="max")    {data=t(t(data)/colsums)*max(colsums)}
    if (type=="min")    {data=t(t(data)/colsums)*min(colsums)}
    if (type=="mean")    {data=t(t(data)/colsums)*mean(colsums)}
    if (type=="median") {
      require(DESeq2)
      size.fac<-DESeq2::estimateSizeFactorsForMatrix(data)
      data<-t(t(data)/size.fac)
    }
    if (fractions==T) {
      return(get_colFracts(data))
    } else {
      return(data)
    }
    
  }
}

min_max_normalization_cols_matrix <- function(mat, new_min, new_max) {
  col_mins <- colMins(mat,na.rm =T)
  col_maxs <- colMaxs(mat,na.rm =T)
  
  normalized_mat <- ((mat - matrix(col_mins, nrow(mat), ncol(mat), byrow=TRUE)) / 
                       (matrix(col_maxs - col_mins, nrow(mat), ncol(mat), byrow=TRUE))) * 
    (new_max - new_min) + new_min
  
  return(normalized_mat)
}


load_counts_from_files <- function(path,
                                   pattern,
                                   rem_cnts_thr = 0,
                                   header       = NA   # NA = auto-detect
) {
  # ---- libraries -----------------------------------------------------------
  library(readr)     # fast reading of delimited files
  library(dplyr)     # data wrangling; explicit calls as requested
  library(purrr)     # map() / reduce()
  library(tibble)    # column_to_rownames()
  
  # ---- discover files ------------------------------------------------------
  files <- list.files(path       = path,
                      pattern    = pattern,
                      full.names = TRUE)
  if (length(files) == 0)
    stop("No files found: check 'path' and 'pattern'")
  
  # ---- helper to read ONE file --------------------------------------------
  read_one <- function(file) {
    sample_name <- sub("\\..*$", "", basename(file))     # text before 1st dot
    
    # -- decide whether we have a header row ---------------------------------
    if (is.na(header)) {                                # auto-detect
      first_line  <- read_lines(file, n_max = 1)
      has_header  <- grepl("[A-Za-z]", first_line)
    } else {
      has_header  <- isTRUE(header)
    }
    
    # -- read the file -------------------------------------------------------
    dat <- readr::read_csv(
      file,
      col_names = !has_header,               # FALSE ⇒ X1, X2 …
      skip      = if (has_header) 1 else 0,
      show_col_types = FALSE
    )
    
    # Require at least two columns (barcode + count)
    if (ncol(dat) < 2)
      stop(sprintf("File %s has fewer than two columns", basename(file)))
    
    # -- rename by POSITION (dplyr-version-agnostic) -------------------------
    colnames(dat)[1] <- "barcode"
    colnames(dat)[2] <- sample_name
    
    # -- filter low-count barcodes ------------------------------------------
    dat <- dplyr::filter(dat, .data[[sample_name]] >= rem_cnts_thr)
    return(dat)
  }
  
  # ---- read all files & merge ---------------------------------------------
  df_list    <- purrr::map(files, read_one)
  all_counts <- purrr::reduce(df_list, dplyr::full_join, by = "barcode")
  all_counts[is.na(all_counts)] <- 0
  
  # ---- matrix output -------------------------------------------------------
  mat <- all_counts |>
    tibble::column_to_rownames("barcode") |>
    as.matrix()
  
  return(mat)
}

#' Aggregate feature values by cell/lineage mapping:
#' returns unweighted means and (optionally) weighted means
#'
#' The first column of `mapping` is interpreted as the 10x **cell barcode**,
#' the second column as the **lineage barcode**. No renaming is performed.
#' Unweighted aggregation returns the mean across linked lineages per cell/feature.
#' If `weights` are supplied, a weighted mean (sum(w*x)/sum(w)) is also returned.
#'
#' @param values matrix or data.frame.
#'   Rows = lineage barcodes; columns = features. Feature names must be unique and non-empty.
#' @param mapping data.frame.
#'   Two-column mapping: col1 = cell barcode, col2 = lineage barcode. Only these two columns are used.
#'   Duplicate (cell, lineage) pairs are removed (dplyr::distinct()).
#' @param weights matrix or data.frame, optional.
#'   Same labeling as `values` (rownames = lineage barcodes; colnames = features).
#'   Non-finite weights (Inf/-Inf/NaN) are treated as missing.
#'
#' @return list with
#'   - mean_unweighted: data.frame [cells x features]
#'   - weighted_mean:  data.frame [cells x features] (present only if `weights` provided)
#'
#' @export
aggregate_values_by_mapping <- function(values, mapping, weights = NULL) {
  if (is.null(values)) stop("`values` must not be NULL.")
  if (!is.matrix(values) && !is.data.frame(values)) stop("`values` must be a matrix or data.frame.")
  if (is.null(rownames(values))) stop("`values` must have rownames = lineage barcodes.")
  if (is.null(colnames(values))) stop("`values` must have colnames = features.")
  if (any(is.na(colnames(values))) || any(colnames(values) == "")) {
    stop("`values` must have non-missing, non-empty feature names.")
  }
  if (any(duplicated(colnames(values)))) {
    dup <- unique(colnames(values)[duplicated(colnames(values))])
    stop("`values` has duplicated feature names: ", paste(dup, collapse = ", "), ".")
  }
  if (!is.data.frame(mapping) || ncol(mapping) < 2) {
    stop("`mapping` must be a data.frame with at least two columns: first = cell, second = lineage.")
  }
  
  cell_col    <- colnames(mapping)[1] %||% names(mapping)[1]
  lineage_col <- colnames(mapping)[2] %||% names(mapping)[2]
  
  mapping2 <- mapping %>%
    dplyr::select(1:2) %>%
    dplyr::distinct()
  
  feature_order <- colnames(values)
  cell_ids      <- unique(mapping2[[1]])
  
  # Long-form values
  values_long <- values %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = lineage_col) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(lineage_col),
      names_to = "feature",
      values_to = "x"
    )
  
  # Join & unweighted mean
  vx <- mapping2 %>%
    dplyr::inner_join(
      values_long,
      by = lineage_col,
      relationship = "many-to-many"
    )
  
  if (nrow(vx) == 0L) {
    mean_unweighted <- as.data.frame(
      matrix(NA_real_, nrow = length(cell_ids), ncol = length(feature_order),
             dimnames = list(cell_ids, feature_order))
    )
  } else {
    mean_unweighted_long <- vx %>%
      dplyr::group_by(!!rlang::sym(cell_col), feature) %>%
      dplyr::summarise(
        mean_unweighted = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
        .groups = "drop"
      )
    
    mean_unweighted <- mean_unweighted_long %>%
      tidyr::pivot_wider(names_from = "feature", values_from = "mean_unweighted") %>%
      tibble::column_to_rownames(var = cell_col) %>%
      as.data.frame()
    
    # Add missing columns and reorder
    missing_cols <- setdiff(feature_order, colnames(mean_unweighted))
    if (length(missing_cols) > 0) {
      for (mc in missing_cols) mean_unweighted[[mc]] <- NA_real_
    }
    # Add missing rows (cells with no overlapping lineages)
    missing_cells <- setdiff(cell_ids, rownames(mean_unweighted))
    if (length(missing_cells) > 0) {
      add <- as.data.frame(matrix(NA_real_, nrow = length(missing_cells), ncol = ncol(mean_unweighted),
                                  dimnames = list(missing_cells, colnames(mean_unweighted))))
      mean_unweighted <- rbind(mean_unweighted, add)
    }
    mean_unweighted <- mean_unweighted[ cell_ids, feature_order, drop = FALSE ]
  }
  
  result <- list(mean_unweighted = mean_unweighted)
  
  # Optional weighted mean
  if (!is.null(weights)) {
    if (!is.matrix(weights) && !is.data.frame(weights)) {
      stop("`weights` must be a matrix or data.frame when provided.")
    }
    if (is.null(rownames(weights)) || is.null(colnames(weights))) {
      stop("`weights` must have rownames (lineage) and colnames (features).")
    }
    if (any(is.na(colnames(weights))) || any(colnames(weights) == "")) {
      stop("`weights` must have non-missing, non-empty feature names.")
    }
    if (any(duplicated(colnames(weights)))) {
      warning("`weights` has duplicated feature names; only unique feature columns will be matched.")
    }
    
    weights_long <- weights %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = lineage_col) %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(lineage_col),
        names_to = "feature",
        values_to = "w"
      )
    
    vwx <- mapping2 %>%
      dplyr::inner_join(
        values_long,
        by = lineage_col,
        relationship = "many-to-many"
      ) %>%
      dplyr::inner_join(
        weights_long,
        by = c(lineage_col, "feature"),
        relationship = "many-to-many"
      ) %>%
      # drop non-finite or missing weights, and missing x
      dplyr::mutate(w = dplyr::if_else(is.finite(w), w, NA_real_)) %>%
      dplyr::filter(!is.na(x), !is.na(w))
    
    if (any(vwx$w < 0, na.rm = TRUE)) {
      warning("Negative weights detected; interpret weighted means with caution.")
    }
    
    if (nrow(vwx) == 0L) {
      weighted_mean <- as.data.frame(
        matrix(NA_real_, nrow = length(cell_ids), ncol = length(feature_order),
               dimnames = list(cell_ids, feature_order))
      )
    } else {
      weighted_long <- vwx %>%
        dplyr::group_by(!!rlang::sym(cell_col), feature) %>%
        dplyr::summarise(
          num = sum(w * x),
          den = sum(w),
          .groups = "drop"
        ) %>%
        dplyr::mutate(weighted_mean = dplyr::if_else(den > 0, num / den, NA_real_)) %>%
        dplyr::select(dplyr::all_of(c(cell_col, "feature", "weighted_mean")))
      
      weighted_mean <- weighted_long %>%
        tidyr::pivot_wider(names_from = "feature", values_from = "weighted_mean") %>%
        tibble::column_to_rownames(var = cell_col) %>%
        as.data.frame()
      
      # Add missing columns and reorder
      missing_cols_w <- setdiff(feature_order, colnames(weighted_mean))
      if (length(missing_cols_w) > 0) {
        for (mc in missing_cols_w) weighted_mean[[mc]] <- NA_real_
      }
      # Add missing rows
      missing_cells_w <- setdiff(cell_ids, rownames(weighted_mean))
      if (length(missing_cells_w) > 0) {
        add <- as.data.frame(matrix(NA_real_, nrow = length(missing_cells_w), ncol = ncol(weighted_mean),
                                    dimnames = list(missing_cells_w, colnames(weighted_mean))))
        weighted_mean <- rbind(weighted_mean, add)
      }
      weighted_mean <- weighted_mean[ cell_ids, feature_order, drop = FALSE ]
    }
    
    result$weighted_mean <- weighted_mean
  }
  
  return(result)
}

# helper: `%||%`
`%||%` <- function(x, y) if (is.null(x)) y else x




#' Summarize Matrix Rows by Column Groups
#'
#' This function summarizes each row of a matrix by applying a specified summary
#' function to groups of columns defined by a grouping vector. The result is a
#' new matrix where each column corresponds to a unique group, and each entry is
#' the summary statistic of the original row's values within that group.
#'
#' @param mat A numeric matrix to be summarized. Each row represents an observation,
#'            and each column represents a variable.
#' @param groups A vector specifying the group assignment for each column in `mat`.
#'               The length of `groups` must equal the number of columns in `mat`.
#' @param func A summary function to apply to each group within a row (e.g., `sum`, `mean`).
#' @param retain_colnames Logical. If `TRUE`, the row names of `mat` are retained
#'                        in the resulting matrix. Default is `FALSE`.
#'
#' @return A numeric matrix with the same number of rows as `mat` and columns
#'         corresponding to each unique group. Each entry is the result of applying
#'         `func` to the values in the specified group for that row.
#'
#' @examples
#' # Example matrix
#' mat = matrix(1:12, nrow = 3, ncol = 4)
#' colnames(mat) = c("A1", "A2", "B1", "B2")
#'
#' # Define groups
#' groups = c("Group1", "Group1", "Group2", "Group2")
#'
#' # Summarize rows by groups using sum
#' summarized_mat = summarize_matrix_rows_by_column_groups(
#'   mat = mat,
#'   groups = groups,
#'   func = sum,
#'   retain_colnames = TRUE
#' )
#' print(summarized_mat)
#'
#' @export
summarize_matrix_rows_by_column_groups = function(mat, groups, func, retain_colnames = FALSE) {
  
  # Validate that 'mat' is a matrix
  if (!is.matrix(mat)) {
    stop("The 'mat' argument must be a matrix.")
  }
  
  # Validate that 'groups' is a vector with length equal to the number of columns in 'mat'
  if (length(groups) != ncol(mat)) {
    stop("The length of 'groups' must match the number of columns in 'mat'.")
  }
  
  # Identify unique groups
  unique_groups = unique(groups)
  
  # Initialize a matrix to store the summarized results
  # Number of rows remains the same, number of columns equals the number of unique groups
  result = matrix(
    data = NA,
    nrow = nrow(mat),
    ncol = length(unique_groups),
    dimnames = list(if (retain_colnames && !is.null(rownames(mat))) rownames(mat) else NULL,
                    unique_groups)
  )
  
  # Iterate over each row of the matrix and apply the summary function to each group
  for (i in seq_len(nrow(mat))) {
    row_values = mat[i, ]
    summarized_values = sapply(unique_groups, function(g) {
      func(row_values[groups == g])
    })
    result[i, ] = summarized_values
  }
  
  # Retain row names if specified
  if (retain_colnames && !is.null(rownames(mat))) {
    rownames(result) = rownames(mat)
  }
  
  return(result)
}



col_min_except_zero <- function(df_or_mat) {
  m <- as.matrix(df_or_mat)               # works for data.frame or matrix
  res <- apply(m, 2, find_min_except_zero)
  return(res)                             # named numeric vector
}

