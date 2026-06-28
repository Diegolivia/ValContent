#' Hinkin–Tracey Content Validity Indices for Multiple Items
#'
#' @description
#' Computes Hinkin and Tracey (1999) content validity indices for multiple items
#' when judges rate how well each item corresponds to several construct
#' definitions in the same judgement task. For each item, the function returns:
#' \itemize{
#'   \item Descriptive means of ratings for each construct.
#'   \item \code{htc}: correspondence index based on the mean rating of the
#'         target construct.
#'   \item \code{htd}: distinctiveness index based on the mean rating
#'         difference between the target construct and all orbiting constructs.
#'   \item \code{htd} computed pairwise between the target construct and each
#'         orbiting construct.
#' }
#'
#' Ratings are assumed to be given on a discrete Likert-type scale from
#' 1 to \code{anchors}. Confidence intervals for the indices are based on
#' a t-based interval for the mean, truncated to the possible rating range,
#' and then rescaled to the \code{htc}/\code{htd} metrics.
#'
#' @param data A data frame or matrix in wide format, where each column
#'   corresponds to one item–construct combination. Columns are expected
#'   to follow the pattern \code{"item.construct"}, e.g.,
#'   \code{"item1.c1"}, \code{"item1.c2"}, \code{"item1.c3"}.
#' @param items Character vector with the base names of the items to be
#'   analyzed (i.e., the prefixes before the dot), e.g.
#'   \code{c("item1", "item2", "item3")}.
#' @param constructs Character vector with the construct labels (suffixes
#'   after the dot), e.g. \code{c("c1", "c2", "c3")}. The function will
#'   search for columns \code{paste0(item, ".", construct)} for every
#'   combination of \code{items} and \code{constructs}.
#' @param key Named character vector mapping each item to its target
#'   construct. The names must match \code{items}, and the values must be
#'   one of the \code{constructs}. For example:
#'   \code{c(item1 = "c2", item2 = "c1", item3 = "c3")}.
#' @param anchors Integer. Number of response options in the rating scale
#'   (e.g., 5 or 7). Ratings are assumed to range from 1 to \code{anchors}.
#' @param ci Logical. If \code{TRUE}, asymmetric confidence intervals are
#'   computed for \code{htc}, \code{htd} and the pairwise \code{htd} indices.
#'   Default is \code{FALSE}.
#' @param conf.level Confidence level for the confidence intervals (e.g.,
#'   \code{0.90}, \code{0.95}, \code{0.99}). Used only if \code{ci = TRUE}.
#'   Default is \code{0.95}.
#' @param na.rm Logical. If \code{TRUE} (default), missing values are removed
#'   when computing means and differences. If \code{FALSE}, missing values
#'   will propagate \code{NA} in the corresponding statistics.
#' @param nd Integer. Number of decimal places to use when rounding
#'   the numeric results. Default is \code{3}.
#'
#' @details
#' For each item, the function extracts the ratings for each construct from
#' columns named \code{"item.construct"}. The target construct is specified
#' via the \code{key} argument.
#'
#' \strong{Item-level descriptive statistics}
#'
#' For each item, \code{$Item.descriptive} reports:
#' \itemize{
#'   \item \code{item}: item name (base).
#'   \item \code{target}: target construct for that item.
#'   \item \code{nj}: number of judges with a non-missing rating on the
#'         target construct.
#'   \item \code{M.<construct>}: mean rating for each construct.
#' }
#'
#' \strong{Item-level criteria}
#'
#' For each item, \code{$Item.criteria} reports:
#' \itemize{
#'   \item \code{htc}: the correspondence index, defined as the mean rating
#'         on the target construct divided by \code{anchors}.
#'   \item \code{htd}: the distinctiveness index, defined as:
#'         \eqn{\text{htd} = \bar{D} / (a - 1)}, where \eqn{\bar{D}} is the
#'         mean (over judges) of the average difference between the rating
#'         on the target construct and the ratings on all orbiting
#'         constructs, and \eqn{a} is \code{anchors}.
#'   \item If \code{ci = TRUE}, lower and upper confidence limits
#'         (\code{htc.lci}, \code{htc.uci}, \code{htd.lci}, \code{htd.uci})
#'         based on a t-interval for the mean truncated to the possible
#'         rating range and then rescaled.
#' }
#'
#' \strong{Pairwise distinctiveness}
#'
#' For each item and each orbiting construct (i.e., each construct other
#' than the target), \code{$Pairwise.criteria} reports:
#' \itemize{
#'   \item \code{item}: item name.
#'   \item \code{target}: target construct for that item.
#'   \item \code{orbiting}: the orbiting construct.
#'   \item \code{htd}: a pairwise distinctiveness index defined as:
#'         \eqn{\text{htd}_{c} = \bar{D}_c / (a - 1)}, where \eqn{\bar{D}_c}
#'         is the mean (over judges) of the difference between the rating on
#'         the target construct and the rating on the orbiting construct
#'         \eqn{c}.
#'   \item If \code{ci = TRUE}, \code{htd.lci} and \code{htd.uci} as
#'         truncated t-based confidence intervals rescaled to the
#'         \code{htd} metric.
#' }
#'
#' This function is intended for the Hinkin and Tracey (1999) method of
#' content validation, where judges rate how well each item corresponds
#' to each construct definition. The same wide format can also be useful
#' for implementing criteria inspired by Colquitt et al. (2019).
#'
#' Ratings are assumed to be coded from 1 to \code{anchors}. If another
#' coding scheme is used, users should recode the ratings before calling
#' the function.
#'
#' @return
#' A list with three components:
#' \describe{
#'   \item{\code{Item.descriptive}}{A data frame with one row per item,
#'     containing item-level descriptive statistics:
#'     \code{item}, \code{target}, \code{nj}, and the construct means
#'     (\code{M.<construct>} columns).}
#'   \item{\code{Item.criteria}}{A data frame with one row per item,
#'     containing the global Hinkin–Tracey indices:
#'     \code{item}, \code{htc}, \code{htc.lci}, \code{htc.uci},
#'     \code{htd}, \code{htd.lci}, \code{htd.uci}. If \code{ci = FALSE},
#'     the confidence interval columns are filled with \code{NA}.}
#'   \item{\code{Pairwise.criteria}}{A data frame in long format, with one
#'     row per item–orbiting construct combination, containing:
#'     \code{item}, \code{target}, \code{orbiting}, \code{htd},
#'     \code{htd.lci}, \code{htd.uci}. If \code{ci = FALSE}, the confidence
#'     interval columns are filled with \code{NA}.}
#' }
#'
#' @references
#' Hinkin, T. R., & Tracey, J. B. (1999). An analysis of variance approach
#' to content validation. \emph{Organizational Research Methods, 2}(2), 175–186.
#'
#' Colquitt, J. A., Sabey, T. B., Rodell, J. B., & Hill, E. T. (2019).
#' Content validation guidelines: Evaluation criteria for definitional
#' correspondence and definitional distinctiveness. \emph{Journal of Applied
#' Psychology, 104}(10), 1243–1265.
#'
#' Penfield, R. D., & Miller, J. M. (2004). Improving content validation
#' studies using an asymmetric confidence interval for the mean of expert
#' ratings. \emph{Applied Measurement in Education, 17}(4), 359–370.
#'
#' @examples
#' \dontrun{
#' ## Example structure (toy data)
#' set.seed(123)
#' HTmult.data <- data.frame(
#'   item1.c1 = sample(1:5, 50, replace = TRUE),
#'   item1.c2 = sample(1:5, 50, replace = TRUE),
#'   item1.c3 = sample(1:5, 50, replace = TRUE),
#'   item2.c1 = sample(1:5, 50, replace = TRUE),
#'   item2.c2 = sample(1:5, 50, replace = TRUE),
#'   item2.c3 = sample(1:5, 50, replace = TRUE)
#' )
#'
#'HTmult(
#'   data       = dat,
#'   items      = c("item1", "item2"),
#'   constructs = c("c1", "c2", "c3"),
#'   key        = c(item1 = "c2", item2 = "c3"),
#'   anchors    = 5,
#'   ci         = TRUE,
#'   conf.level = 0.95)
#' }
#'
#' @export
HTmult <- function(
    data,
    items,
    constructs,
    key,
    anchors,
    ci         = FALSE,
    conf.level = 0.95,
    na.rm      = TRUE,
    nd     = 3
) {
  # --- Basic checks and setup ---
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

  # Helper: t-based CI for mean with truncation to [lower, upper]
  mean_ci_truncated <- function(x, lower, upper, conf.level) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n == 0) {
      return(list(mean = NA_real_, lwr = NA_real_, upr = NA_real_))
    }
    m <- mean(x)
    if (n == 1 || is.na(sd(x)) || sd(x) == 0) {
      # No variability or single observation: CI collapses to the truncated mean
      lwr <- max(min(m, upper), lower)
      upr <- lwr
    } else {
      alpha <- 1 - conf.level
      se <- sd(x) / sqrt(n)
      tcrit <- stats::qt(1 - alpha / 2, df = n - 1)
      lwr <- m - tcrit * se
      upr <- m + tcrit * se
      # truncate to possible range
      lwr <- max(lwr, lower)
      upr <- min(upr, upper)
    }
    list(mean = m, lwr = lwr, upr = upr)
  }

  item_desc_list   <- list()
  item_crit_list   <- list()
  pairwise_list    <- list()

  idx_desc <- 1L
  idx_crit <- 1L
  idx_pair <- 1L

  # --- Loop over items ---
  for (item in items) {

    target <- key[[item]]

    # Build column names for this item
    item_cols <- paste0(item, ".", constructs)

    # Check that all required columns exist
    missing_cols <- setdiff(item_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop(
        "For item '", item, "', the following columns are missing in 'data': ",
        paste(missing_cols, collapse = ", ")
      )
    }

    # Extract ratings per construct
    ratings_list <- lapply(item_cols, function(nm) data[[nm]])
    names(ratings_list) <- constructs

    # Target ratings
    target_vec <- ratings_list[[target]]

    # Number of judges (non-missing in target)
    if (na.rm) {
      nj <- sum(!is.na(target_vec))
    } else {
      nj <- length(target_vec)
    }

    # --- Item.descriptive row ---
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
    # Rename construct mean columns to M.<construct>
    mean_col_names <- paste0("M.", constructs)
    colnames(desc_row)[(ncol(desc_row) - length(constructs) + 1):ncol(desc_row)] <- mean_col_names

    item_desc_list[[idx_desc]] <- desc_row
    idx_desc <- idx_desc + 1L

    # --- Item.criteria row: HTC and HTD global ---
    # HTC
    mean_target <- means_c[target]
    htc         <- mean_target / anchors

    htc_lci <- htc_uci <- NA_real_
    if (ci) {
      # Ratings assumed 1..anchors
      ci_target <- mean_ci_truncated(
        x      = target_vec,
        lower  = 1,
        upper  = anchors,
        conf.level = conf.level
      )
      if (!is.na(ci_target$lwr)) {
        htc_lci <- ci_target$lwr / anchors
        htc_uci <- ci_target$upr / anchors
      }
    }

    # HTD global: average difference target - orbiting (over constructs),
    # then averaged over judges and normalized by (anchors - 1)
    orbiting_constructs <- setdiff(constructs, target)
    # For each judge, compute mean difference across orbiting constructs
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

    # Range of D: [-(anchors-1), +(anchors-1)]
    D_mean <- if (na.rm) mean(D_vec, na.rm = TRUE) else mean(D_vec)
    htd    <- D_mean / (anchors - 1)

    htd_lci <- htd_uci <- NA_real_
    if (ci) {
      ci_D <- mean_ci_truncated(
        x      = D_vec,
        lower  = -(anchors - 1),
        upper  =  (anchors - 1),
        conf.level = conf.level
      )
      if (!is.na(ci_D$lwr)) {
        htd_lci <- ci_D$lwr / (anchors - 1)
        htd_uci <- ci_D$upr / (anchors - 1)
      }
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

    # --- Pairwise.criteria: HTD per orbiting construct ---
    for (cc in orbiting_constructs) {
      vec_c <- ratings_list[[cc]]

      # Differences per judge for this pair
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
        ci_dc <- mean_ci_truncated(
          x      = d_c,
          lower  = -(anchors - 1),
          upper  =  (anchors - 1),
          conf.level = conf.level
        )
        if (!is.na(ci_dc$lwr)) {
          htd_c_lci <- ci_dc$lwr / (anchors - 1)
          htd_c_uci <- ci_dc$upr / (anchors - 1)
        }
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

  # Bind results
  Item.descriptive <- do.call(rbind, item_desc_list)
  Item.criteria    <- do.call(rbind, item_crit_list)
  Pairwise.criteria <- do.call(rbind, pairwise_list)

  # Round numeric columns
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
