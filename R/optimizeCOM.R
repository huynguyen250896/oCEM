#' @title pcaCharts
#'
#' @description Visualizes the screeplot showing the optimal number of principal components.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage function(x, y, optimal = NULL, main = NULL)
#'
#' @param x A list returned by \code{svd()}, typically from the real data.
#' @param y A numeric vector representing the average percent explained variance
#'   of random matrices.
#' @param optimal Optional integer. The selected optimal number of components.
#'   If provided, a vertical line will be drawn at this location.
#' @param main Optional plot title.
#'
#' @return NULL
#' @keywords internal
# NOT EXPORTED

pcaCharts <- function(x, y, optimal = NULL, main = NULL) {

  # Percent explained variance of real data
  x.var <- x$d^2 / sum(x$d^2) * 100

  # Cumulative explained variance of real data
  x.pvar <- cumsum(x.var)

  # Average percent explained variance of random matrices
  x.var1 <- y

  # Save old graphical parameters and restore on exit
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))

  # Left axis range
  y.min <- min(c(x.var, x.var1), na.rm = TRUE)
  y.max <- max(c(x.var, x.var1), na.rm = TRUE)

  # Plot explained variance of observed data and permuted average
  par(mar = c(5, 5, 2, 5))
  par(xaxs = "i", yaxs = "i")

  plot(
    x.var,
    xlab = "k",
    ylab = "% explained by each component",
    axes = FALSE,
    xlim = c(1, length(x.var)),
    ylim = c(y.min, y.max),
    type = "l",
    col = "blue",
    main = main
  )

  lines(x.var1, col = "purple")

  x.ticks <- seq_len(length(x.var))
  axis(1, at = x.ticks, cex.axis = 0.85)

  y.ticks <- pretty(c(y.min, y.max))
  axis(
    2,
    at = y.ticks,
    labels = paste0(y.ticks, "%"),
    las = 1,
    mgp = c(3, 0.75, 0),
    cex.axis = 0.75
  )

  # Draw selected optimal number of components, if provided
  if (!is.null(optimal) && length(optimal) == 1 && !is.na(optimal) && optimal >= 1) {
    abline(v = optimal, col = "red", lwd = 2, lty = 2)
  }

  # Add cumulative explained variance curve on the right axis
  par(new = TRUE)
  par(xaxs = "i", yaxs = "i")

  plot(
    x.pvar,
    xlab = NA,
    ylab = NA,
    axes = FALSE,
    ylim = c(0, 100),
    type = "l",
    col = "green"
  )

  y.ticks1 <- seq(0, 100, 10)
  axis(side = 4, at = y.ticks1, labels = paste0(y.ticks1, "%"), las = 1, cex.axis = 0.75)
  mtext("Cumulative explained variance", side = 4, line = 3)

  # Add legend
  legend(
    "topright",
    legend = c("Observed", "Permuted mean", "Cumulative observed", "Selected k"),
    col = c("blue", "purple", "green", "red"),
    lty = c(1, 1, 1, 2),
    lwd = c(1, 1, 1, 2),
    bty = "n",
    cex = 0.8
  )
}


#' @title optimizeCOM
#'
#' @description Determines the optimal number of principal components using either
#'   the intersection rule or a permutation p-value rule.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage optimizeCOM(
#'   data = NULL,
#'   P = 1000,
#'   standardize = TRUE,
#'   cores = 1,
#'   verbose = TRUE,
#'   type = c("intersection", "perm"),
#'   alpha = 0.05,
#'   adjust.method = "BH",
#'   seed = NULL,
#'   return.details = FALSE
#' )
#'
#' @param data A data frame or matrix. Rows are samples and columns are genes.
#' @param P Positive integer. Number of permutations. Default is 1000.
#' @param standardize Logical. If \code{TRUE}, data are standardized before SVD
#'   and before downstream decomposition.
#' @param cores Positive integer. Number of cores used by \code{parallel::mclapply()}.
#'   Default is 1.
#' @param verbose Logical. Whether to print analysis progress and show scree plot.
#' @param type Character string. Method to determine the optimal number of components.
#'   One of \code{"intersection"} or \code{"perm"}.
#' @param alpha Numeric. Significance threshold used when \code{type = "perm"}.
#' @param adjust.method Character string. Method passed to \code{p.adjust()} when
#'   \code{type = "perm"}. Default is \code{"BH"}.
#' @param seed Optional integer. Random seed for reproducibility.
#' @param return.details Logical. If \code{TRUE}, returns a detailed list. Otherwise
#'   returns only the selected number of components.
#'
#' @return By default, an integer giving the selected number of components.
#'   If \code{return.details = TRUE}, returns a list with additional outputs.
#'
#' @examples
#' \dontrun{
#' optimizeCOM(data = exp, type = "intersection")
#' optimizeCOM(data = exp, type = "perm", alpha = 0.05, seed = 250896)
#' }
#'
#' @export

optimizeCOM <- function(data = NULL,
                        P = 1000,
                        standardize = TRUE,
                        cores = 1,
                        verbose = TRUE,
                        type = c("intersection", "perm"),
                        alpha = 0.05,
                        adjust.method = "BH",
                        seed = 250896,
                        return.details = FALSE) {

  # ---------------------------------------------------------------------------
  # 1. Input checks
  # ---------------------------------------------------------------------------

  if (missing(data) || is.null(data)) {
    stop("Error: gene expression data are missing.\n")
  }

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Error: data must be a matrix or data.frame.\n")
  }

  data <- as.matrix(data)

  if (!is.numeric(data)) {
    stop("Error: data must be numeric.\n")
  }

  if (nrow(data) < 2 || ncol(data) < 2) {
    stop("Error: data must have at least 2 rows and 2 columns.\n")
  }

  if (anyNA(data)) {
    stop("Error: data contains missing values. Please remove or impute them first.\n")
  }

  if (!is.logical(standardize) || length(standardize) != 1 || is.na(standardize)) {
    stop("Error: standardize must be a single TRUE or FALSE.\n")
  }

  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("Error: verbose must be a single TRUE or FALSE.\n")
  }

  if (!is.logical(return.details) || length(return.details) != 1 || is.na(return.details)) {
    stop("Error: return.details must be a single TRUE or FALSE.\n")
  }

  if (!is.numeric(P) || length(P) != 1 || is.na(P) || P < 1) {
    stop("Error: P must be a positive integer.\n")
  }

  if (!is.numeric(cores) || length(cores) != 1 || is.na(cores) || cores < 1) {
    stop("Error: cores must be a positive integer.\n")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Error: alpha must be a single number between 0 and 1.\n")
  }

  P <- as.integer(P)
  cores <- as.integer(cores)
  type <- match.arg(type)

  # ---------------------------------------------------------------------------
  # Optional seed for reproducibility
  # ---------------------------------------------------------------------------

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
      stop("Error: seed must be NULL or a single numeric value.\n")
    }
    seed <- as.integer(seed)
    set.seed(seed)
  }

  # ---------------------------------------------------------------------------
  # Logging helper
  # ---------------------------------------------------------------------------

  now <- Sys.time()

  mlog <- if (!verbose) {
    function(...) {}
  } else {
    function(...) {
      message(...)
      flush.console()
    }
  }

  mlog("Starting optimizeCOM...")

  # ---------------------------------------------------------------------------
  # SVD on the observed data
  # ---------------------------------------------------------------------------

  # If requested, standardize data before SVD
  if (isTRUE(standardize)) {
    svd_obj <- svd(scale(data))
  } else {
    svd_obj <- svd(data)
  }

  # Percent explained variance of observed data
  x.var <- svd_obj$d^2 / sum(svd_obj$d^2) * 100

  # ---------------------------------------------------------------------------
  # 5. Build permutation reference distribution
  # ---------------------------------------------------------------------------

  # For each permutation:
  #   - independently permute entries within each sample (row)
  #   - compute singular values
  #   - convert them into percent explained variance
  #
  # Result:
  #   ranMat is a matrix of dimension:
  #   number_of_components x P
  #
  # Each column corresponds to one permutation.
  mlog("Running permutations...")

  ranMat.list <- parallel::mclapply(seq_len(P), mc.cores = cores, FUN = function(p) {

    # Generate a random matrix by independently permuting values in each row
    permutated_matrix <- t(apply(data, 1, sample))

    # Standardize permuted matrix if requested
    if (isTRUE(standardize)) {
      d_perm <- svd(scale(permutated_matrix))$d
    } else {
      d_perm <- svd(permutated_matrix)$d
    }

    # Return percent explained variance of this permuted matrix
    d_perm^2 / sum(d_perm^2) * 100
  })

  ranMat <- do.call(cbind, ranMat.list)

  # Average reference curve across permutations
  average <- rowMeans(ranMat)

  # ---------------------------------------------------------------------------
  # 6. Choose optimal number of components
  # ---------------------------------------------------------------------------

  # Two supported strategies:
  #
  # A. "intersection"
  #    Find where observed explained variance curve intersects the average
  #    permutation curve.
  #
  # B. "perm"
  #    For each component j, compute permutation p-value:
  #      p_j = (1 + #{perm_var_j >= obs_var_j}) / (P + 1)
  #    Then adjust p-values and select the largest j that remains significant.

  pvals <- NULL
  pvals.adj <- NULL
  optimal.raw <- NULL

  if (identical(type, "intersection")) {

    mlog("Selecting components using intersection rule...")

    # Find indices where the logical relation changes between consecutive points
    # Example: observed > average changes from TRUE to FALSE
    intersect.points <- which(diff(x.var > average) != 0)

    if (length(intersect.points) == 0) {
      # If no formal intersection exists, fallback to counting how many observed
      # components exceed the permutation mean
      optimal1 <- sum(x.var > average)

      # Ensure at least 1 component is selected
      if (optimal1 < 1) {
        optimal1 <- 1L
      }

      optimal.raw <- optimal1

    } else {
      # Compute local slopes around each candidate intersection
      x1.slopes <- x.var[intersect.points + 1] - x.var[intersect.points]
      x2.slopes <- average[intersect.points + 1] - average[intersect.points]

      # Avoid division by zero in pathological cases
      valid <- which((x1.slopes - x2.slopes) != 0)

      if (length(valid) == 0) {
        optimal1 <- sum(x.var > average)
        if (optimal1 < 1) {
          optimal1 <- 1L
        }
        optimal.raw <- optimal1
      } else {
        intersect.points <- intersect.points[valid]
        x1.slopes <- x1.slopes[valid]
        x2.slopes <- x2.slopes[valid]

        # Linear interpolation of intersection location
        optimal.raw <- intersect.points +
          ((average[intersect.points] - x.var[intersect.points]) / (x1.slopes - x2.slopes))

        # Choose the first intersection
        optimal1 <- min(ceiling(optimal.raw))

        if (!is.finite(optimal1) || is.na(optimal1) || optimal1 < 1) {
          optimal1 <- 1L
        }
      }
    }

  } else if (identical(type, "perm")) {

    mlog("Selecting components using permutation p-values...")

    # Permutation p-value for each component
    pvals <- vapply(seq_along(x.var), function(j) {
      (1 + sum(ranMat[j, ] >= x.var[j])) / (P + 1)
    }, numeric(1))

    # Multiple testing correction
    pvals.adj <- p.adjust(pvals, method = adjust.method)

    # Significant components
    sig.idx <- which(pvals.adj < alpha)

    if (length(sig.idx) == 0) {
      # If none are significant, select 1 to avoid downstream crash
      # You may change this to 0 if your package logic wants that behavior
      optimal1 <- 1L
    } else {
      optimal1 <- max(sig.idx)
    }

    optimal.raw <- optimal1
  }

  # Clamp optimal1 to valid range
  max.comp <- min(nrow(data), ncol(data))
  optimal1 <- as.integer(max(1, min(optimal1, max.comp)))

  # ---------------------------------------------------------------------------
  # 7. Plot scree plot if requested
  # ---------------------------------------------------------------------------

  if (isTRUE(verbose)) {
    plot.title <- if (identical(type, "intersection")) {
      "Optimal number of components by intersection"
    } else {
      "Optimal number of components by permutation test"
    }

    pcaCharts(svd_obj, average, optimal = optimal1, main = plot.title)
  }
    # ---------------------------------------------------------------------------
  # 8. Report selected number of components
  # ---------------------------------------------------------------------------

  cat(">> oCEM suggests choosing the optimal number of components is:", optimal1, "\n")

  # ---------------------------------------------------------------------------
  # 9. Final timing message
  # ---------------------------------------------------------------------------

  timediff <- Sys.time() - now
  mlog("Done in ", timediff, " ", units(timediff), ".")

  # ---------------------------------------------------------------------------
  # 10. Return result
  # ---------------------------------------------------------------------------

  if (isTRUE(return.details)) {
    return(list(
      optimal_ncomp = optimal1,
      method = type,
      seed = seed,
      standardize = standardize,
      observed_variance = x.var,
      permuted_variance = ranMat,
      permuted_mean_variance = average,
      raw_intersection = optimal.raw,
      pvals = pvals,
      pvals_adjusted = pvals.adj
    ))
  } else {
    return(optimal1)
  }
}



#   # ---------------------------------------------------------------------------
#   # 8. Prepare data for ICA and IPCA
#   # ---------------------------------------------------------------------------

#   # Standardize data before decomposition if requested
#   data.proc <- data
#   if (isTRUE(standardize)) {
#     data.proc <- scale(data.proc)
#   }

#   # Re-apply seed before ICA for reproducibility
#   if (!is.null(seed)) {
#     set.seed(seed)
#   }

#   # ---------------------------------------------------------------------------
#   # 9. Run ICA
#   # ---------------------------------------------------------------------------

#   mlog("Running ICA...")

#   ICA <- fastICA::fastICA(
#     t(data.proc),
#     n.comp = optimal1,
#     fun = "logcosh",
#     alpha = 1,
#     method = "C",
#     row.norm = FALSE,
#     maxit = 1000,
#     tol = 0.0001
#   )

#   rownames(ICA$S) <- colnames(data.proc)

#   # Signatures of genes from ICA
#   kurt_ICA <- ICA$S

#   # ---------------------------------------------------------------------------
#   # 10. Run IPCA
#   # ---------------------------------------------------------------------------

#   mlog("Running IPCA...")

#   IPCA <- mixOmics::ipca(
#     data.proc,
#     ncomp = optimal1,
#     fun = "logcosh",
#     mode = "deflation"
#   )

#   # Signatures of genes from IPCA
#   kurt_IPCA <- IPCA[["loadings"]][["X"]]

#   # ---------------------------------------------------------------------------
#   # 11. Compute kurtosis for ICA and IPCA components
#   # ---------------------------------------------------------------------------

#   n.comp.final <- min(ncol(kurt_ICA), ncol(kurt_IPCA))

#   cc2 <- data.frame(
#     component = paste("Huy", seq_len(n.comp.final)),
#     ICA = numeric(n.comp.final),
#     IPCA = numeric(n.comp.final),
#     stringsAsFactors = FALSE
#   )

#   for (i in seq_len(n.comp.final)) {
#     cc2$ICA[i] <- e1071::kurtosis(kurt_ICA[, i])
#     cc2$IPCA[i] <- e1071::kurtosis(kurt_IPCA[, i])
#   }

#   # ---------------------------------------------------------------------------
#   # 12. Message based on kurtosis rule
#   # ---------------------------------------------------------------------------

#   cat(">> oCEM suggests choosing the optimal number of components is:", optimal1, "\n")

#   if (all(cc2$ICA < 3) && any(cc2$IPCA > 3)) {
#     cat(">> oCEM also suggests using IPCA-FDR for your case.\n")
#   }

#   if (all(cc2$ICA < 3) && all(cc2$IPCA < 3)) {
#     cat(">> oCEM cannot suggest which method should be selected. Please use a more stringent approach to make the best decision for your case.\n")
#   }

#   if (any(cc2$ICA > 3) && all(cc2$IPCA < 3)) {
#     cat(">> oCEM also suggests using ICA for your case.\n")
#   }

#   if (any(cc2$ICA > 3) && any(cc2$IPCA > 3)) {
#     cat(">> Both ICA and IPCA-FDR are appropriate for your case. Please use a more stringent approach to make the best decision.\n")
#   }

#   # ---------------------------------------------------------------------------
#   # 13. Final timing message
#   # ---------------------------------------------------------------------------

#   timediff <- Sys.time() - now
#   mlog("Done in ", timediff, " ", units(timediff), ".")

#   # ---------------------------------------------------------------------------
#   # 14. Return result
#   # ---------------------------------------------------------------------------

#   if (isTRUE(return.details)) {
#     return(list(
#       optimal_ncomp = optimal1,
#       method = type,
#       seed = seed,
#       standardize = standardize,
#       observed_variance = x.var,
#       permuted_variance = ranMat,
#       permuted_mean_variance = average,
#       raw_intersection = optimal.raw,
#       pvals = pvals,
#       pvals_adjusted = pvals.adj,
#       kurtosis = cc2,
#       ICA = ICA,
#       IPCA = IPCA
#     ))
#   } else {
#     return(optimal1)
#   }
# }