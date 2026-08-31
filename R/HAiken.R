#' Coefficient of Homogeneity of Response (Aiken's H)
#'
#' @description
#' Calculates Aiken's H coefficient of homogeneity for each item and, optionally,
#' an overall H (total) across all items. The function uses bootstrap resampling
#' of judges to obtain confidence intervals.
#'
#' @param data A data frame with judges in columns and items in rows.
#' @param ncat Number of response categories.
#' @param conf.level Confidence level for the intervals (e.g., .90, .95).
#' @param na.rm Logical. If TRUE, rows with missing values are removed.
#' @param overall Logical. If TRUE (default), the overall H (total) is added as the last row.
#' @param B Integer. Number of bootstrap resamples (default = 1000).
#' @param ci.type Character: "logit" (default, recommended), "perc" (percentile), or "norm" (normal approximation).
#'        The "logit" method transforms the bootstrap H values to the logit scale,
#'        computes percentiles there, and back-transforms, ensuring the CI lies within [0,1].
#'
#' @details
#' The procedure for obtaining confidence intervals with ci.type = "logit" is as follows:
#' \enumerate{
#'   \item Compute the point estimate of H for each item and for the total (if overall = TRUE).
#'   \item Generate B bootstrap samples by resampling the judges (columns) with replacement.
#'         For each bootstrap sample, recompute H for each item and the total.
#'   \item Apply the logit transformation to each bootstrap H value:
#'         \deqn{L = \log(H / (1 - H))}.
#'         To avoid infinities when H = 0 or H = 1, a small constant \eqn{\varepsilon = 10^{-6}}
#'         is added (or subtracted) so that \eqn{H^* = \max(\varepsilon, \min(1-\varepsilon, H))}.
#'   \item Obtain the percentiles of the logit-transformed bootstrap distribution:
#'         \eqn{L_{\text{inf}} = \text{percentile}_{\alpha/2}(L)}, \eqn{L_{\text{sup}} = \text{percentile}_{1-\alpha/2}(L)}.
#'   \item Back-transform to the original scale:
#'         \eqn{H_{\text{inf}} = \exp(L_{\text{inf}}) / (1 + \exp(L_{\text{inf}}))},
#'         \eqn{H_{\text{sup}} = \exp(L_{\text{sup}}) / (1 + \exp(L_{\text{sup}}))}.
#' }
#' This ensures that the confidence interval respects the [0,1] bounds and is asymmetric when appropriate.
#'
#' The methods "perc" and "norm" use the bootstrap percentiles or normal approximation directly on H,
#' which may produce limits outside [0,1] in extreme cases; they are provided for comparison but are not recommended.
#'
#' @return A data frame with columns:
#'   \item{Item}{Item number, or "Total" for the overall coefficient.}
#'   \item{H}{Aiken's H coefficient.}
#'   \item{lwr.ci}{Lower bound of the confidence interval.}
#'   \item{upr.ci}{Upper bound of the confidence interval.}
#'   \item{n.jueces}{Number of judges (same for all rows).}
#'
#' @references
#' Aiken, L. R. (1980). Content validity and reliability of single items or questionnaires.
#'   \emph{Educational and Psychological Measurement, 40}, 955-959.
#'
#' Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings.
#'   \emph{Educational and Psychological Measurement, 45}, 131-142.
#'
#' @examples
#' \dontrun{
#' # Sample data: 8 items, 18 judges, ratings 1-4
#' datos <- data.frame(t(file1[, c("claA6", "claA18", "claA2", "claA16",
#'                                 "claA4", "claA12", "claA9", "claA21")]))
#' Haiken(datos, ncat = 4, conf.level = .90, overall = TRUE, B = 1000, ci.type = "logit")
#' }
#'
#' @export
Haiken <- function(data, ncat, conf.level = .95, na.rm = FALSE,
                   overall = TRUE, B = 1000,
                   ci.type = c("logit", "perc", "norm")) {

  ci.type <- match.arg(ci.type)

  if (!is.data.frame(data)) data <- as.data.frame(data)

  if (!na.rm && anyNA(data)) {
    stop("Missing values found. Set na.rm=TRUE or handle them first.")
  }
  if (na.rm) data <- na.omit(data)

  n_items <- nrow(data)
  m <- ncol(data)
  if (n_items < 1 || m < 2) stop("Need at least 1 item and 2 judges.")

  delta <- ifelse(m %% 2 == 0, 1, 0)
  denom_item <- (ncat - 1) * (m^2 - delta)
  denom_total <- (ncat - 1) * n_items * (m^2 - delta)

  # --- Internal function to compute H (individual and total) ---
  calc_H <- function(mat) {
    ni <- nrow(mat)
    mi <- ncol(mat)
    delta_i <- ifelse(mi %% 2 == 0, 1, 0)
    denom_i <- (ncat - 1) * (mi^2 - delta_i)
    denom_total_i <- (ncat - 1) * ni * (mi^2 - delta_i)

    H_ind <- numeric(ni)
    sum_diffs_total <- 0
    all_valid <- TRUE

    for (i in 1:ni) {
      # Convert row to numeric robustly
      row_i <- as.numeric(mat[i, ])
      # If any NA, try to convert from factor/character
      if (anyNA(row_i)) {
        row_i <- as.numeric(as.character(mat[i, ]))
      }
      # If still any NA, impute with median of non-NA values (or skip)
      if (anyNA(row_i)) {
        row_i[is.na(row_i)] <- median(row_i, na.rm = TRUE)
      }

      # Compute pairwise absolute differences
      diffs <- tryCatch(
        combn(row_i, 2, function(p) abs(p[1] - p[2])),
        error = function(e) rep(NA, choose(mi, 2))
      )

      if (anyNA(diffs)) {
        all_valid <- FALSE
        H_ind[i] <- NA
        sum_diffs_total <- NA
        break
      } else {
        sum_diffs <- sum(diffs, na.rm = TRUE)
        sum_diffs_total <- sum_diffs_total + sum_diffs
        H_ind[i] <- 1 - (4 * sum_diffs) / denom_i
      }
    }

    H_total <- if (anyNA(H_ind) || !all_valid) {
      NA
    } else {
      1 - (4 * sum_diffs_total) / denom_total_i
    }

    list(H_individual = H_ind, H_total = H_total, valid = all_valid)
  }

  # --- Point estimates ---
  est <- calc_H(data)
  H_ind <- est$H_individual
  H_total <- est$H_total

  if (anyNA(H_ind)) stop("Point estimates contain NA. Check data and ncat.")
  if (overall && is.na(H_total)) stop("Total H is NA. Check data.")

  # --- Bootstrap over judges (columns) with resampling until valid ---
  boot_H_ind <- matrix(NA, nrow = B, ncol = n_items)
  boot_H_total <- numeric(B)

  for (b in 1:B) {
    valid_boot <- FALSE
    attempts <- 0
    while (!valid_boot && attempts < 10) {
      idx_col <- sample(1:m, size = m, replace = TRUE)
      boot_data <- data[, idx_col, drop = FALSE]
      boot_est <- calc_H(boot_data)
      if (boot_est$valid && !anyNA(boot_est$H_individual) && !is.na(boot_est$H_total)) {
        valid_boot <- TRUE
        boot_H_ind[b, ] <- boot_est$H_individual
        boot_H_total[b] <- boot_est$H_total
      }
      attempts <- attempts + 1
    }
    if (!valid_boot) {
      # Fallback: use point estimates
      boot_H_ind[b, ] <- H_ind
      boot_H_total[b] <- H_total
    }
  }

  # --- Helper to compute confidence intervals (robust) ---
  get_ci <- function(boot_vals, point_est, conf.level, ci.type) {
    boot_vals <- boot_vals[!is.na(boot_vals)]
    if (length(boot_vals) == 0) {
      return(c(lwr = point_est, upr = point_est))
    }
    if (length(unique(boot_vals)) == 1) {
      return(c(lwr = point_est, upr = point_est))
    }

    alpha <- 1 - conf.level

    if (ci.type == "perc") {
      lwr <- quantile(boot_vals, probs = alpha/2, na.rm = TRUE, names = FALSE)
      upr <- quantile(boot_vals, probs = 1 - alpha/2, na.rm = TRUE, names = FALSE)
      if (is.na(lwr)) lwr <- min(boot_vals)
      if (is.na(upr)) upr <- max(boot_vals)
      return(c(lwr = lwr, upr = upr))

    } else if (ci.type == "norm") {
      se <- sd(boot_vals, na.rm = TRUE)
      z <- qnorm(1 - alpha/2)
      lwr <- point_est - z * se
      upr <- point_est + z * se
      return(c(lwr = lwr, upr = upr))

    } else { # "logit" (default)
      eps <- 1e-6
      vals <- boot_vals
      vals[vals <= 0] <- eps
      vals[vals >= 1] <- 1 - eps

      logit_vals <- log(vals / (1 - vals))
      lwr_logit <- quantile(logit_vals, probs = alpha/2, na.rm = TRUE, names = FALSE)
      upr_logit <- quantile(logit_vals, probs = 1 - alpha/2, na.rm = TRUE, names = FALSE)
      if (is.na(lwr_logit)) lwr_logit <- min(logit_vals)
      if (is.na(upr_logit)) upr_logit <- max(logit_vals)

      lwr <- exp(lwr_logit) / (1 + exp(lwr_logit))
      upr <- exp(upr_logit) / (1 + exp(upr_logit))
      return(c(lwr = lwr, upr = upr))
    }
  }

  # --- Compute CIs for each item ---
  lwr_ind <- numeric(n_items)
  upr_ind <- numeric(n_items)
  for (i in 1:n_items) {
    ci <- get_ci(boot_H_ind[, i], H_ind[i], conf.level, ci.type)
    lwr_ind[i] <- ci["lwr"]
    upr_ind[i] <- ci["upr"]
  }

  # --- CI for the total H (if overall = TRUE) ---
  if (overall) {
    ci_total <- get_ci(boot_H_total, H_total, conf.level, ci.type)
    lwr_total <- ci_total["lwr"]
    upr_total <- ci_total["upr"]
  } else {
    lwr_total <- NA
    upr_total <- NA
  }

  # --- Build results data.frame ---
  resultados <- data.frame(
    Item = 1:n_items,
    H = round(H_ind, 3),
    lwr.ci = round(lwr_ind, 3),
    upr.ci = round(upr_ind, 3),
    n.jueces = m,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (overall) {
    fila_total <- data.frame(
      Item = "Total",
      H = round(H_total, 3),
      lwr.ci = if (is.na(lwr_total)) NA else round(lwr_total, 3),
      upr.ci = if (is.na(upr_total)) NA else round(upr_total, 3),
      n.jueces = m,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    resultados <- rbind(resultados, fila_total)
  }

  attr(resultados, "ci.type") <- ci.type
  attr(resultados, "B") <- B
  attr(resultados, "conf.level") <- conf.level

  return(resultados)
}
