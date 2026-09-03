#' Hinkin–Tracey Content Validity Indices with Confidence Intervals
#'
#' @description
#' Computes Hinkin and Tracey (1999) content validity indices for multiple items
#' using the standardized operationalizations and empirical benchmarks proposed by
#' Colquitt et al. (2019). For each item, the function computes:
#' \itemize{
#'   \item Descriptive means of ratings across all evaluated constructs.
#'   \item \code{htc}: Hinkin–Tracey correspondence index based on the mean rating of the
#'         target construct.
#'   \item \code{htd}: Hinkin–Tracey distinctiveness index based on the rescaled mean
#'         difference between ratings of the target construct and orbiting constructs.
#'   \item \code{htd} computed pairwise between the target construct and each individual
#'         orbiting construct.
#' }
#'
#' Asymmetric confidence intervals are constructed using Wilson's score interval method
#' for proportions for \code{htc} (Penfield & Miller, 2004), and Willink's (2005) asymmetric
#' interval incorporating the third moment and Cornish-Fisher expansion for \code{htd},
#' eliminating out-of-bounds limits and providing realistic coverage even under small expert sample sizes.
#'
#' @param data A data frame or matrix in wide format, where each column
#'   corresponds to an item–construct rating. Column names must follow the pattern
#'   \code{"item.construct"} (e.g., \code{"item1.c1"}, \code{"item1.c2"}).
#' @param items Character vector with the base names of the items (e.g.,
#'   \code{c("item1", "item2")}).
#' @param constructs Character vector with construct labels (e.g.,
#'   \code{c("c1", "c2", "c3")}).
#' @param key Named character vector mapping each item to its target construct.
#'   Names must match \code{items} and values must belong to \code{constructs}
#'   (e.g., \code{c(item1 = "c2", item2 = "c1")}).
#' @param anchors Integer. Number of response scale options (e.g., 5 or 7). Ratings
#'   are assumed to range from 1 to \code{anchors}.
#' @param ci Logical. If \code{TRUE}, asymmetric confidence intervals are
#'   computed for \code{htc} (Wilson score) and \code{htd} (Willink Cornish-Fisher). Default is \code{FALSE}.
#' @param conf.level Confidence level for intervals (e.g., \code{0.90}, \code{0.95}).
#'   Default is \code{0.95}.
#' @param na.rm Logical. If \code{TRUE} (default), missing values are removed pair-wise or
#'   list-wise depending on the sub-index calculation.
#' @param nd Integer. Number of decimal places for rounding output values. Default is \code{3}.
#'
#' @details
#' \strong{Empirical Benchmarks (Colquitt et al., 2019)}
#'
#' Based on empirical decile distributions across extensive content validation studies,
#' Colquitt et al. (2019) suggest the following evaluation criteria:
#' \itemize{
#'   \item \strong{Definitional Correspondence (\code{htc}):}
#'     \itemize{
#'       \item Strong: \eqn{\ge 0.86}
#'       \item Moderate: \eqn{0.78} to \eqn{0.85}
#'       \item Weak: \eqn{\le 0.77}
#'     }
#'   \item \strong{Definitional Distinctiveness (\code{htd}):}
#'     \itemize{
#'       \item Strong: \eqn{\ge 0.27}
#'       \item Moderate: \eqn{0.21} to \eqn{0.26}
#'       \item Weak: \eqn{\le 0.20}
#'     }
#' }
#'
#' \strong{Rationale and Construction of Confidence Intervals}
#'
#' Standard Studentized \emph{t}-intervals often fail in content validity tasks because rating
#' distributions near scale bounds produce zero sample variance or skewness. To overcome this:
#' \itemize{
#'   \item \strong{\code{htc} Transformation:} Computed via Wilson's score interval (Penfield & Miller, 2004)
#'         after mapping the raw target mean to proportion space \eqn{p \in [0, 1]}.
#'   \item \strong{\code{htd} Transformation:} Computed via Willink's (2005) asymmetric interval method,
#'         which corrects for sample skewness using the third central moment and Cornish-Fisher expansion,
#'         with post-hoc truncation to the natural domain \eqn{[-1, 1]}.
#' }
#'
#' \strong{Methodological Note on Comparisons:}
#' Confidence intervals should be used to evaluate the precision of individual coefficients or to
#' compare items within the \emph{same} metric type (e.g., comparing \code{htc} between Item 1 and Item 2).
#' Comparing an \code{htc} interval directly against an \code{htd} interval is methodologically invalid,
#' as they capture fundamentally different theoretical constructs and variance structures.
#'
#' @return A list with three data frames:
#' \describe{
#'   \item{\code{Item.descriptive}}{Item names, target constructs, judge counts (\code{nj}), and mean ratings.}
#'   \item{\code{Item.criteria}}{Global \code{htc} and \code{htd} indices with corresponding confidence intervals.}
#'   \item{\code{Pairwise.criteria}}{Pairwise \code{htd} indices comparing target vs. each orbiting construct.}
#' }
#'
#' @references
#' Colquitt, J. A., Sabey, T. B., Rodell, J. B., & Hill, E. T. (2019). Content validation guidelines:
#' Evaluation criteria for definitional correspondence and definitional distinctiveness.
#' \emph{Journal of Applied Psychology, 104}(10), 1243–1265.
#'
#' Hinkin, T. R., & Tracey, J. B. (1999). An analysis of variance approach to content validation.
#' \emph{Organizational Research Methods, 2}(2), 175–186.
#'
#' Penfield, R. D., & Miller, J. M. (2004). Improving content validation studies using an asymmetric
#' confidence interval for the mean of expert ratings. \emph{Applied Measurement in Education, 17}(4), 359–370.
#'
#' Willink, R. (2005). A confidence interval for a mean with a correction for skewness.
#' \emph{Metrika}, 61(2), 159–173.
#'
#' @export
ColquittHT <- function(
    data,
    items,
    constructs,
    key,
    anchors,
    ci         = FALSE,
    conf.level = 0.95,
    na.rm      = TRUE,
    nd         = 3
) {
  # --- Basic argument validation ---
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (!is.character(items) || length(items) == 0) {
    stop("'items' must be a non-empty character vector.")
  }

  if (!is.character(constructs) || length(constructs) < 2) {
    stop("'constructs' must be a character vector with at least two constructs.")
  }

  if (missing(key) || is.null(key)) {
    stop("'key' must be provided as a named character vector mapping items to target constructs.")
  }

  if (is.null(names(key)) || any(names(key) == "")) {
    stop("'key' must have names corresponding to item names.")
  }

  if (!all(items %in% names(key))) {
    stop("All 'items' must appear as names in 'key'.")
  }

  if (!all(key[items] %in% constructs)) {
    stop("All target constructs in 'key' must be included in 'constructs'.")
  }

  if (!is.numeric(anchors) || length(anchors) != 1 || anchors <= 1) {
    stop("'anchors' must be a numeric value > 1 (e.g., 5 or 7).")
  }

  if (!is.logical(ci) || length(ci) != 1) {
    stop("'ci' must be a single logical value (TRUE/FALSE).")
  }

  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a numeric value between 0 and 1.")
  }

  # --- Internal Helper: Wilson Score Confidence Interval for Proportions (HTC) ---
  wilson_ci <- function(p, n, conf.level) {
    if (is.na(p) || is.na(n) || n <= 0) {
      return(list(lwr = NA_real_, upr = NA_real_))
    }
    p <- max(0, min(1, p))

    alpha <- 1 - conf.level
    z <- stats::qnorm(1 - alpha / 2)
    z2 <- z^2

    denominator <- 1 + z2 / n
    center <- (p + z2 / (2 * n)) / denominator
    spread <- (z * sqrt((p * (1 - p) / n) + (z2 / (4 * n^2)))) / denominator

    lwr <- max(0, center - spread)
    upr <- min(1, center + spread)

    list(lwr = lwr, upr = upr)
  }

  # --- Internal Helper: Willink Asymmetric Confidence Interval for HTD (2005) ---
  willink_ci_htd <- function(x, conf.level) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 3) {
      return(list(lwr = NA_real_, upr = NA_real_))
    }

    x_bar <- mean(x)
    s <- stats::sd(x)
    alpha <- (1 - conf.level) / 2

    if (s == 0) {
      return(list(lwr = x_bar, upr = x_bar))
    }

    # Unbiased third central moment
    mu3_hat <- (n / ((n - 1) * (n - 2))) * sum((x - x_bar)^3)

    # Scaled skewness coefficient (Willink, 2005)
    a <- (mu3_hat / (s^3)) / (6 * sqrt(n))

    t_low <- stats::qt(alpha, df = n - 1)
    t_high <- stats::qt(1 - alpha, df = n - 1)

    G <- function(r) {
      if (abs(a) < 1e-6) return(r)
      inner <- 1 + 6 * a * (r - a)
      if (inner < 0) inner <- 0
      (1 / (2 * a)) * (inner^(1/3) - 1)
    }

    lwr <- x_bar - G(t_high) * (s / sqrt(n))
    upr <- x_bar - G(t_low) * (s / sqrt(n))

    # Truncate to natural bounds of htd [-1, 1]
    lwr <- max(-1, lwr)
    upr <- min(1, upr)

    list(lwr = lwr, upr = upr)
  }

  item_desc_list <- list()
  item_crit_list <- list()
  pairwise_list  <- list()

  idx_desc <- 1L
  idx_crit <- 1L
  idx_pair <- 1L

  # --- Loop over items ---
  for (item in items) {

    target <- key[[item]]

    # Construct standard column names for the current item
    item_cols <- paste0(item, ".", constructs)

    # Verify column existence in input dataset
    missing_cols <- setdiff(item_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop(
        "For item '", item, "', the following columns are missing in 'data': ",
        paste(missing_cols, collapse = ", ")
      )
    }

    # Extract construct rating vectors
    ratings_list <- lapply(item_cols, function(nm) data[[nm]])
    names(ratings_list) <- constructs

    target_vec <- ratings_list[[target]]

    # Determine effective judge sample size on target construct
    nj <- if (na.rm) sum(!is.na(target_vec)) else length(target_vec)

    # --- 1. Descriptive statistics row ---
    means_c <- sapply(
      ratings_list,
      function(v) if (na.rm) mean(v, na.rm = TRUE) else mean(v)
    )

    desc_row <- data.frame(
      item   = item,
      target = target,
      nj     = nj,
      t(means_c),
      check.names = FALSE
    )
    mean_col_names <- paste0("M.", constructs)
    colnames(desc_row)[(ncol(desc_row) - length(constructs) + 1):ncol(desc_row)] <- mean_col_names

    item_desc_list[[idx_desc]] <- desc_row
    idx_desc <- idx_desc + 1L

    # --- 2. Global Criteria: HTC and HTD ---
    mean_target <- means_c[target]
    htc         <- mean_target / anchors

    htc_lci <- htc_uci <- NA_real_
    if (ci && nj > 0) {
      p_target <- (mean_target - 1) / (anchors - 1)
      ci_target <- wilson_ci(p = p_target, n = nj, conf.level = conf.level)

      mean_lwr <- ci_target$lwr * (anchors - 1) + 1
      mean_upr <- ci_target$upr * (anchors - 1) + 1
      htc_lci  <- mean_lwr / anchors
      htc_uci  <- mean_upr / anchors
    }

    # Global HTD computation (judge-level mean difference across orbiting constructs)
    orbiting_constructs <- setdiff(constructs, target)
    D_vec <- rep(NA_real_, length(target_vec))

    for (i in seq_along(target_vec)) {
      t_i <- target_vec[i]
      if (is.na(t_i)) next
      orbit_vals <- vapply(
        orbiting_constructs,
        function(cc) ratings_list[[cc]][i],
        FUN.VALUE = NA_real_
      )
      if (na.rm) {
        orbit_vals <- orbit_vals[!is.na(orbit_vals)]
      }
      if (length(orbit_vals) == 0) next
      D_vec[i] <- t_i - mean(orbit_vals)
    }

    D_mean <- if (na.rm) mean(D_vec, na.rm = TRUE) else mean(D_vec)
    htd    <- D_mean / (anchors - 1)

    htd_lci <- htd_uci <- NA_real_
    if (ci) {
      judge_htd <- D_vec / (anchors - 1)
      ci_htd <- willink_ci_htd(judge_htd, conf.level = conf.level)
      htd_lci <- ci_htd$lwr
      htd_uci <- ci_htd$upr
    }

    crit_row <- data.frame(
      item    = item,
      htc     = htc,
      htc.lci = htc_lci,
      htc.uci = htc_uci,
      htd     = htd,
      htd.lci = htd_lci,
      htd.uci = htd_uci,
      check.names = FALSE
    )

    item_crit_list[[idx_crit]] <- crit_row
    idx_crit <- idx_crit + 1L

    # --- 3. Pairwise Criteria: HTD per orbiting construct ---
    for (cc in orbiting_constructs) {
      vec_c <- ratings_list[[cc]]

      if (na.rm) {
        idx_valid <- !is.na(target_vec) & !is.na(vec_c)
        d_c <- target_vec[idx_valid] - vec_c[idx_valid]
      } else {
        d_c <- target_vec - vec_c
      }

      D_c_mean <- if (na.rm) mean(d_c, na.rm = TRUE) else mean(d_c)
      htd_c    <- D_c_mean / (anchors - 1)

      htd_c_lci <- htd_c_uci <- NA_real_
      if (ci) {
        judge_htd_c <- d_c / (anchors - 1)
        ci_dc <- willink_ci_htd(judge_htd_c, conf.level = conf.level)
        htd_c_lci <- ci_dc$lwr
        htd_c_uci <- ci_dc$upr
      }

      pair_row <- data.frame(
        item     = item,
        target   = target,
        orbiting = cc,
        htd      = htd_c,
        htd.lci  = htd_c_lci,
        htd.uci  = htd_c_uci,
        check.names = FALSE
      )

      pairwise_list[[idx_pair]] <- pair_row
      idx_pair <- idx_pair + 1L
    }
  }

  # Combine rows across items
  Item.descriptive  <- do.call(rbind, item_desc_list)
  Item.criteria     <- do.call(rbind, item_crit_list)
  Pairwise.criteria <- do.call(rbind, pairwise_list)

  # Round numeric output columns
  round_numeric_df <- function(df, digits) {
    is_num <- vapply(df, is.numeric, logical(1))
    df[is_num] <- lapply(df[is_num], round, digits = digits)
    df
  }

  Item.descriptive  <- round_numeric_df(Item.descriptive, nd)
  Item.criteria     <- round_numeric_df(Item.criteria, nd)
  Pairwise.criteria <- round_numeric_df(Pairwise.criteria, nd)

  list(
    Item.descriptive  = Item.descriptive,
    Item.criteria     = Item.criteria,
    Pairwise.criteria = Pairwise.criteria
  )
}
