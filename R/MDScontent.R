#' Multidimensional scaling (MDS) Map for Content Validity Ratings (Item–Trait Correspondence)
#'
#' @description
#' `MDScontent()` builds a **conceptual map** (MDS) from judges' ratings of the
#' correspondence between each item and each trait (attribute). The function:
#' (1) aggregates judges' ratings into an **Items × Traits** profile matrix,
#' (2) computes **item-to-item dissimilarities** (default: Euclidean),
#' (3) obtains a **2D MDS configuration** (default: non-metric, with STRESS),
#' (4) computes **trait centroids** in the MDS space (theoretical via `key`, or empirical),
#' (5) optionally returns an **item–trait distance matrix** (`Dmatrix = TRUE`),
#' and (6) plots either the main map, a biplot, or both.
#'
#' @param data Numeric `matrix` or `data.frame` of size J × (I*T): rows are judges;
#'   each column corresponds to a specific (item, trait) pair rating.
#' @param item Vector of length `ncol(data)` indicating the **item id** for each column.
#' @param trait Vector of length `ncol(data)` indicating the **trait id** for each column.
#' @param key Optional vector of length I (number of unique items). Theoretical trait
#'   assignment for each item (used to define centroids). If `NULL`, items are assigned
#'   empirically to the trait with the largest profile score.
#' @param score Aggregation of judges into the Items × Traits profile matrix:
#'   `"mean"` (default), `"median"`, or `"p_ge"` (proportion of ratings >= `cut`).
#' @param cut Integer threshold for `score = "p_ge"`. Default is 4.
#' @param distance Item-to-item dissimilarity computed from profiles:
#'   `"euclid"` (default), `"cor"` (1 - correlation), or `"cosine"`.
#' @param mds MDS type: `"nonmetric"` (default, returns STRESS) or `"classic"`.
#' @param k Dimensionality for MDS solution. Default 2 (plots assume k = 2).
#' @param centroid How to compute centroid coordinates per trait:
#'   `"mean"` (default) or `"median"` (coordinate-wise).
#' @param display Plot type: `"items"` (main map), `"biplot"`, or `"both"`.
#' @param Dmatrix Logical; if `TRUE`, returns the item–trait distance matrix Δ (I × T).
#' @param label.items Logical; label items in plot(s). Default `TRUE`.
#' @param label.traits Logical; label trait centroids in plot(s). Default `TRUE`.
#' @param ... Additional arguments forwarded to the MDS routine:
#'   - for nonmetric: `trymax`, `maxit`, `tol`
#'
#' @return A list with elements:
#' \describe{
#'   \item{profile}{Items × Traits matrix used to build distances.}
#'   \item{dist_items}{`dist` object of item-to-item dissimilarities.}
#'   \item{coords_items}{data.frame of MDS coordinates for items.}
#'   \item{centroids}{data.frame of centroid coordinates for traits and centroid type.}
#'   \item{fit}{list with `mds`, `k`, `stress` (if nonmetric), and `gof` (R^2 of distances).}
#'   \item{Dmatrix}{(Optional) matrix Δ of item–trait distances in MDS space.}
#' }
#'
#' @details
#' ## Conceptual rationale
#' The function provides a **geometric visualization** of content validity structure.
#' Ratings from judges are first aggregated into an Items × Traits profile matrix.
#' Item-to-item dissimilarities are computed (default: Euclidean distance),
#' and a multidimensional scaling (MDS) solution is obtained (default: non-metric).
#'
#' Trait centroids are computed in the MDS space using either:
#' - Theoretical assignment (`key`), or
#' - Empirical assignment (trait with maximum profile value per item).
#'
#' The resulting map represents:
#' - Items as points,
#' - Traits as centroids,
#' - Optional segments from each item to its assigned trait centroid,
#' - Optional trait vectors (biplot) indicating gradients of trait correspondence.
#'
#' ## Interpretation
#' This function is intended primarily for **visual diagnostic purposes**.
#' It does not replace quantitative content validity coefficients (e.g., Aiken's V,
#' CVI, SVAL, or Hit-based approaches), but complements them by examining
#' structural coherence.
#'
#' The map may be interpreted as follows:
#'
#' - **Item clustering**: Items targeting the same trait should appear spatially close.
#'
#' - **Distance to centroid**: Items located near their theoretical trait centroid
#'   indicate strong conceptual specificity.
#'
#' - **Ambiguous items**: Items located between centroids or far from their expected
#'   centroid may reflect conceptual overlap or insufficient specificity.
#'
#' - **Trait redundancy**: If trait centroids are very close or vectors point in
#'   similar directions, this may indicate conceptual overlap between traits.
#'
#' - **STRESS and GOF**: These indices refer to the geometric quality of the
#'   low-dimensional representation and should not be interpreted as evidence
#'   of content validity per se.
#'
#' Therefore, `MDScontent()` should be interpreted as a structural visualization tool
#' that complements coefficient-based evidence.
#'
#' @examples
#' ## Example 1
#'#  This simulated data has the following structure: 12 judges evaluated the fit of 6
#'#  items in  4 attributes. That is: 12 judges × 6 items × 4 traits matrix. In the administration
#'#  of the validity survey, the items were presented in 12 rows, and each item
#'#  was evaluated for its correspondence to four attributes. The database should
#'#  be structured as follows:
#'#
#'#  item item1.judge1 item1.judge2 item1.judge3 item1.judge4 item2.judge1 item2.judge2 ...
#'#
#'#  For this example:
#'
#'#  Item1 (Fit to Trait 1)
#'#  Item2 (Fit to Trait 1)
#'#  Item3 (Fit to Trait 2)
#'#  Item4 (Fit to Trait 2)
#'#  Item5 (Approximate fit to Trait 3)
#'#  Item6 (Approximate fit to Trait 4)
#'#
#'#   Example2 <- matrix(c(
#'#     5,2,2,1,  5,1,2,1,  2,5,1,2,  1,5,2,1,  2,1,5,2,  1,2,2,5,
#'#     5,1,2,2,  4,2,1,1,  1,5,2,1,  2,5,1,1,  1,2,5,2,  2,1,1,5,
#'#     4,2,1,1,  5,1,2,1,  2,4,1,2,  1,5,2,2,  2,1,5,1,  1,2,2,5,
#'#     5,1,2,1,  5,2,1,1,  1,5,2,1,  2,4,1,2,  1,2,5,2,  2,1,1,5,
#'#     4,2,1,2,  5,1,2,1,  2,5,1,1,  1,5,2,1,  2,1,4,2,  1,2,2,4,
#'#     5,1,2,1,  4,2,1,2,  1,5,2,1,  2,5,1,1,  1,2,5,1,  2,1,2,5,
#'#     5,2,1,1,  5,1,2,1,  2,4,1,2,  1,5,2,1,  2,1,5,2,  1,2,1,5,
#'#     4,1,2,1,  5,2,1,1,  1,5,2,2,  2,4,1,1,  1,2,5,1,  2,1,2,5,
#'#     5,2,1,1,  4,1,2,2,  2,5,1,1,  1,5,2,1,  2,1,5,2,  1,2,1,5,
#'#     5,1,2,2,  5,2,1,1,  1,5,2,1,  2,5,1,2,  1,2,5,1,  2,1,2,5,
#'#     4,2,1,1,  5,1,2,1,  2,4,1,2,  1,5,2,1,  2,1,4,2,  1,2,1,5,
#'#     5,1,2,1,  4,2,1,1,  1,5,2,1,  2,4,1,2,  1,2,5,1,  2,1,2,5),
#'#     nrow = 12,
#'#     byrow = TRUE)
#'
#'## Short cut for grouping items
#'rep(1:6, each = 4)
#'
#'## Short cut for grouping traits
#'#'rep(1:4, times = 6)
#'
#'## theorethical correspondence
#'c(1,1,2,2,3,4)
#'
#'MDScontent(data = dat2,
#'item = rep(1:6, each = 4),
#'trait = rep(1:4, times = 6),
#'key = key <- c(1,1,2,2,3,4),
#'score = "mean",
#'distance = "euclid",
#'mds = "nonmetric",
#'display = "items")
#'
#' @export
MDScontent <- function(
    data,
    item,
    trait,
    key = NULL,
    score = c("mean", "median", "p_ge"),
    cut = 4,
    distance = c("euclid", "cor", "cosine"),
    mds = c("nonmetric", "classic"),
    k = 2,
    centroid = c("mean", "median"),
    display = c("items", "biplot", "both"),
    Dmatrix = FALSE,
    label.items = TRUE,
    label.traits = TRUE,
    ...
) {
  # -----------------------
  # Minimal validations
  # -----------------------
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data)) stop("`data` must be a numeric matrix/data.frame.")
  if (ncol(data) < 2) stop("`data` must have at least 2 columns (>= 2 item×trait pairs).")

  if (missing(item) || missing(trait)) stop("Both `item` and `trait` must be provided.")
  if (length(item) != ncol(data)) stop("`item` must have length equal to ncol(data).")
  if (length(trait) != ncol(data)) stop("`trait` must have length equal to ncol(data).")
  if (anyNA(item) || anyNA(trait)) stop("`item` and `trait` must not contain NA.")

  score <- match.arg(score)
  distance <- match.arg(distance)
  mds <- match.arg(mds)
  centroid <- match.arg(centroid)
  display <- match.arg(display)

  if (!is.numeric(cut) || length(cut) != 1) stop("`cut` must be a single numeric value.")
  if (!is.numeric(k) || length(k) != 1 || k < 2) stop("`k` must be an integer >= 2 (plotting is implemented for 2D).")
  k <- as.integer(k)

  # -----------------------
  # Build Items × Traits profile matrix X
  # -----------------------
  item_levels <- unique(item)
  trait_levels <- unique(trait)
  I <- length(item_levels)
  Tt <- length(trait_levels)

  item_index <- match(item, item_levels)
  trait_index <- match(trait, trait_levels)

  X <- matrix(NA_real_, nrow = I, ncol = Tt,
              dimnames = list(paste0("Item", item_levels),
                              paste0("Trait", trait_levels)))

  agg_fun <- switch(
    score,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    p_ge = function(x) mean(x >= cut, na.rm = TRUE)
  )

  for (i in seq_len(I)) {
    for (t in seq_len(Tt)) {
      cols <- which(item_index == i & trait_index == t)
      if (length(cols) == 0) next
      v <- as.numeric(data[, cols, drop = FALSE])
      X[i, t] <- agg_fun(v)
    }
  }

  if (anyNA(X)) {
    stop("Some (item, trait) cells are missing (NA) in the profile matrix. ",
         "Ensure every item is rated on every trait (or provide complete columns).")
  }

  # -----------------------
  # Compute item-to-item dissimilarities (Fixed and cleaned)
  # -----------------------
  dist_items <- switch(
    distance,
    euclid = stats::dist(X, method = "euclidean"),
    cor = {
      # Correlation between item profiles (rows of X) -> cor(t(X))
      Ci <- stats::cor(t(X), method = "pearson", use = "pairwise.complete.obs")
      as.dist(pmax(0, 1 - Ci))
    },
    cosine = {
      # Cosine distance: 1 - cosine similarity
      norms <- sqrt(rowSums(X^2))
      # Protection against division by zero
      norms[norms == 0] <- 1
      Xn <- X / norms
      S <- Xn %*% t(Xn)
      as.dist(pmax(0, 1 - S))
    }
  )

  # -----------------------
  # MDS
  # -----------------------
  coords <- NULL
  stress <- NA_real_
  gof <- NA_real_

  if (mds == "classic") {
    # Use eig = TRUE to get standard GOF
    cmd_fit <- stats::cmdscale(dist_items, k = k, eig = TRUE)
    coords <- cmd_fit$points
    # GOF for classic MDS (proportion of variance explained)
    if (any(cmd_fit$eig > 0)) {
      gof <- sum(cmd_fit$eig[1:k]) / sum(cmd_fit$eig[cmd_fit$eig > 0])
    } else {
      gof <- NA
    }
  } else {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for nonmetric MDS (isoMDS). Please install it.")
    }
    dots <- list(...)
    trymax <- if (!is.null(dots$trymax)) dots$trymax else 50
    maxit  <- if (!is.null(dots$maxit))  dots$maxit  else 200
    tol    <- if (!is.null(dots$tol))    dots$tol    else 1e-3

    init <- stats::cmdscale(dist_items, k = k)
    fit_mds <- MASS::isoMDS(dist_items, y = init, k = k,
                            maxit = maxit, tol = tol, trace = FALSE)
    coords <- fit_mds$points
    stress <- fit_mds$stress

    if (!is.na(stress) && stress < 0.01) {
      warning("STRESS is very close to zero (", signif(stress, 3),
              "). Data may be over-structured.", call. = FALSE)
    }

    # GOF as R^2 between original and reproduced distances
    d2 <- stats::dist(coords[, 1:2, drop = FALSE])
    gof <- suppressWarnings(stats::cor(as.vector(dist_items), as.vector(d2), use = "complete.obs")^2)
  }

  coords <- as.matrix(coords)
  colnames(coords) <- paste0("Dim", seq_len(ncol(coords)))
  rownames(coords) <- rownames(X)

  # -----------------------
  # Trait assignment for centroid computation
  # -----------------------
  if (is.null(key)) {
    row_max <- apply(X, 1, max, na.rm = TRUE)
    ties <- rowSums(X == row_max, na.rm = TRUE) > 1
    if (any(ties)) {
      warning(sum(ties), " item(s) had tied maximum trait scores. Assigning to the first tied trait.", call. = FALSE)
    }
    key_idx <- max.col(X, ties.method = "first")
    key_type <- "empirical"
  } else {
    if (length(key) != I) stop("`key` must have length equal to number of unique items.")
    key_idx <- match(key, trait_levels)
    if (anyNA(key_idx)) {
      stop("`key` contains values not present in `trait` levels.")
    }
    key_type <- "theoretical"
  }

  # -----------------------
  # Centroids per trait in MDS space
  # -----------------------
  Ccoords <- matrix(NA_real_, nrow = Tt, ncol = 2,
                    dimnames = list(paste0("Trait", trait_levels), c("Dim1", "Dim2")))
  for (t in seq_len(Tt)) {
    ids <- which(key_idx == t)
    if (length(ids) == 0) {
      warning("Trait ", trait_levels[t], " has no assigned items. Centroid cannot be computed.", call. = FALSE)
      next
    }
    if (centroid == "mean") {
      Ccoords[t, ] <- colMeans(coords[ids, 1:2, drop = FALSE])
    } else {
      Ccoords[t, ] <- apply(coords[ids, 1:2, drop = FALSE], 2, stats::median)
    }
  }

  centroids_df <- data.frame(
    trait = trait_levels,
    Dim1 = Ccoords[, 1],
    Dim2 = Ccoords[, 2],
    type = key_type,
    stringsAsFactors = FALSE
  )

  # -----------------------
  # Dmatrix: item–trait distances in the MDS plane
  # -----------------------
  Dmat <- NULL
  if (isTRUE(Dmatrix)) {
    Dmat <- matrix(NA_real_, nrow = I, ncol = Tt,
                   dimnames = list(rownames(coords), paste0("Trait", trait_levels)))
    for (i in seq_len(I)) {
      for (t in seq_len(Tt)) {
        if (anyNA(Ccoords[t, ])) next
        Dmat[i, t] <- sqrt(sum((coords[i, 1:2] - Ccoords[t, ])^2))
      }
    }
  }

  fit <- list(mds = mds, k = k, stress = stress, gof = gof,
              score = score, distance = distance, centroid = centroid,
              centroid_type = key_type)

  # -----------------------
  # Plotting (Aesthetically Upgraded)
  # -----------------------
  # Palette setup
  trait_colors <- rainbow(Tt, s = 0.6, v = 0.85)
  names(trait_colors) <- trait_levels
  item_colors <- trait_colors[as.character(trait_levels[key_idx])]

  # Pre-calculate Biplot Vectors
  V <- matrix(NA_real_, nrow = Tt, ncol = 2)
  for (t in seq_len(Tt)) {
    V[t, 1] <- suppressWarnings(stats::cor(X[, t], coords[, 1], use = "pairwise.complete.obs"))
    V[t, 2] <- suppressWarnings(stats::cor(X[, t], coords[, 2], use = "pairwise.complete.obs"))
  }
  V[is.na(V)] <- 0

  # Scaling for arrows (proportional to axis span, preserves angles)
  rx <- diff(range(coords[, 1])); ry <- diff(range(coords[, 2]))
  target_len <- 0.7 * min(rx, ry)
  lens <- sqrt(rowSums(V^2))
  s <- if (max(lens) > 0) target_len / max(lens) else 1
  V2 <- V * s

  .plot_items <- function(overlay_biplot = FALSE) {
    x <- coords[, 1]; y <- coords[, 2]

    # Adjust limits if overlaying biplot
    xlim <- range(x); ylim <- range(y)
    if (overlay_biplot) {
      xlim <- range(c(xlim, V2[, 1]))
      ylim <- range(c(ylim, V2[, 2]))
    }
    xlim <- xlim + c(-0.15, 0.15) * diff(xlim)
    ylim <- ylim + c(-0.15, 0.15) * diff(ylim)

    plot(x, y, type = "n", xlab = "Dimension 1", ylab = "Dimension 2",
         xlim = xlim, ylim = ylim, main = "MDS Content Map",
         bty = "l", las = 1, col.axis = "gray40", col.lab = "gray20")
    grid(col = "gray90", lty = 1)

    # Draw convex hulls and segments
    for (t in seq_len(Tt)) {
      ids <- which(key_idx == t)
      if (length(ids) == 0 || anyNA(Ccoords[t, ])) next

      segments(x0 = x[ids], y0 = y[ids], x1 = Ccoords[t, 1], y1 = Ccoords[t, 2],
               col = adjustcolor(trait_colors[t], alpha.f = 0.4), lwd = 1.5)

      if (length(ids) >= 3) {
        hull <- chull(coords[ids, 1:2])
        polygon(coords[ids[hull], 1], coords[ids[hull], 2],
                col = adjustcolor(trait_colors[t], alpha.f = 0.1),
                border = adjustcolor(trait_colors[t], alpha.f = 0.3), lwd = 1)
      }
    }

    # Plot items
    points(x, y, pch = 21, bg = adjustcolor(item_colors, alpha.f = 0.8), col = "gray20", cex = 1.5)
    if (isTRUE(label.items)) text(x, y, labels = rownames(coords), pos = 3, cex = 0.8, col = "gray30")

    # Plot centroids
    valid_centroid <- !is.na(Ccoords[, 1])
    points(Ccoords[valid_centroid, 1], Ccoords[valid_centroid, 2], pch = 22,
           bg = trait_colors[valid_centroid], col = "gray20", cex = 1.8, lwd = 1.5)
    if (isTRUE(label.traits)) {
      text(Ccoords[valid_centroid, 1], Ccoords[valid_centroid, 2],
           labels = rownames(Ccoords)[valid_centroid], pos = 1, font = 2, cex = 0.9)
    }

    if (overlay_biplot) {
      arrows(0, 0, V2[, 1], V2[, 2], length = 0.1, col = trait_colors, lwd = 2)
      if (isTRUE(label.traits)) {
        text(V2[, 1], V2[, 2], labels = colnames(X), pos = 4, font = 2, cex = 0.9, col = trait_colors)
      }
    }

    mtext(paste0("Method: ", mds, " | Distance: ", distance, " | STRESS: ", signif(stress, 3)),
          side = 3, line = 0.5, cex = 0.8, col = "gray40")
    legend("topright", legend = rownames(Ccoords), pch = 22, pt.bg = trait_colors,
           col = "gray20", bty = "n", title = "Traits", cex = 0.9)
  }

  .plot_biplot <- function() {
    x <- coords[, 1]; y <- coords[, 2]

    xlim <- range(c(x, V2[, 1])) + c(-0.15, 0.15) * diff(range(c(x, V2[, 1])))
    ylim <- range(c(y, V2[, 2])) + c(-0.15, 0.15) * diff(range(c(y, V2[, 2])))

    plot(x, y, type = "n", xlab = "Dimension 1", ylab = "Dimension 2",
         xlim = xlim, ylim = ylim, main = "MDS Biplot",
         bty = "l", las = 1, col.axis = "gray40", col.lab = "gray20")
    grid(col = "gray90", lty = 1)

    points(x, y, pch = 21, bg = adjustcolor("gray50", alpha.f = 0.7), col = "gray20", cex = 1.2)
    if (isTRUE(label.items)) text(x, y, labels = rownames(coords), pos = 3, cex = 0.8, col = "gray30")

    arrows(0, 0, V2[, 1], V2[, 2], length = 0.1, col = trait_colors, lwd = 2)
    if (isTRUE(label.traits)) {
      text(V2[, 1], V2[, 2], labels = colnames(X), pos = 4, font = 2, cex = 0.9, col = trait_colors)
    }

    mtext(paste0("Method: ", mds, " | Distance: ", distance),
          side = 3, line = 0.5, cex = 0.8, col = "gray40")
    legend("topright", legend = colnames(X), pch = 22, pt.bg = trait_colors,
           col = "gray20", bty = "n", title = "Traits", cex = 0.9)
  }

  if (display == "items") {
    .plot_items(overlay_biplot = FALSE)
  } else if (display == "biplot") {
    .plot_biplot()
  } else if (display == "both") {
    .plot_items(overlay_biplot = TRUE)
  }

  # -----------------------
  # Return
  # -----------------------
  coords_df <- data.frame(
    item = item_levels,
    Dim1 = coords[, 1],
    Dim2 = coords[, 2],
    assigned_trait = trait_levels[key_idx],
    assignment_type = key_type,
    stringsAsFactors = FALSE
  )

  out <- list(
    profile = X,
    dist_items = dist_items,
    coords_items = coords_df,
    centroids = centroids_df,
    fit = fit,
    Dmatrix = Dmat
  )
  class(out) <- c("MDScontent", class(out))
  out
}

# Custom print method for clean console output
#' @export
print.MDScontent <- function(x, ...) {
  cat("--- MDS Content Analysis ---\n")
  cat("MDS Method   :", x$fit$mds, "\n")
  cat("Dimensions   :", x$fit$k, "\n")
  if (x$fit$mds == "nonmetric") {
    cat("STRESS       :", format(x$fit$stress, digits = 4), "\n")
  }
  cat("GOF (R^2)    :", format(x$fit$gof, digits = 4), "\n")
  cat("Score        :", x$fit$score, "\n")
  cat("Distance     :", x$fit$distance, "\n")
  cat("Centroid type:", x$fit$centroid_type, "(", x$fit$centroid, ")\n")
  cat("----------------------------\n")
  invisible(x)
}
