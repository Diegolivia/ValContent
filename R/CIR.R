#' MOVER-R Confidence Interval for the Ratio of Two Bounded Proportion Estimates
#'
#' Computes the confidence interval for the ratio of two independent bounded
#' proportion estimates (e.g., Aiken's V coefficients) using the MOVER-R
#' (Method of Variance Estimates Recovery for Ratios) closed-form procedure.
#'
#' @param group1 A \code{data.frame} containing the point estimates and their
#'   confidence limits for the first group. Must include columns specified in
#'   \code{coef.col}, \code{lwr.col}, \code{upr.col}, and an \code{"Item"}
#'   column for matching.
#' @param group2 A \code{data.frame} containing the point estimates and their
#'   confidence limits for the second group. Same column requirements as
#'   \code{group1}.
#' @param coef.col Character string indicating the column name for the point
#'   estimates.
#' @param lwr.col Character string indicating the column name for the lower
#'   confidence limits.
#' @param upr.col Character string indicating the column name for the upper
#'   confidence limits.
#' @param na.rm Logical. If \code{TRUE}, rows containing missing values in the
#'   relevant columns are removed with a warning. If \code{FALSE}, missing
#'   values trigger an error.
#'
#' @return A \code{data.frame} with the following columns:
#'   \item{Item}{Common item identifiers present in both groups.}
#'   \item{R}{The estimated ratio (coefficient of group1 divided by
#'     coefficient of group2).}
#'   \item{lwr.ci}{The lower bound of the MOVER-R confidence interval.}
#'   \item{upr.ci}{The upper bound of the MOVER-R confidence interval.}
#'
#' @details
#' The MOVER-R interval is computed using the closed-form solution proposed by
#' Donner and Zou (2012). Given two independent estimates \eqn{\hat{\theta}_1}
#' and \eqn{\hat{\theta}_2} with corresponding confidence limits
#' \eqn{(l_1, u_1)} and \eqn{(l_2, u_2)}, the interval for the ratio
#' \eqn{R = \theta_1 / \theta_2} is obtained as:
#'
#' \deqn{L = \frac{m - \sqrt{m^2 - A \cdot D}}{D}, \quad
#'       U = \frac{m + \sqrt{m^2 - B \cdot C}}{C}}
#'
#' where \eqn{m = \hat{\theta}_1 \hat{\theta}_2},
#' \eqn{A = l_1(2\hat{\theta}_1 - l_1)},
#' \eqn{B = u_1(2\hat{\theta}_1 - u_1)},
#' \eqn{C = l_2(2\hat{\theta}_2 - l_2)}, and
#' \eqn{D = u_2(2\hat{\theta}_2 - u_2)}.
#'
#' For bounded proportion estimates (e.g., Aiken's V, which lies in [0, 1]),
#' it is common for the lower confidence limit to equal exactly 0 or the upper
#' limit to equal exactly 1. This yields \eqn{C = 0} or \eqn{D = 0},
#' respectively, causing a singular division. To maintain numerical stability,
#' the function applies a small perturbation (\eqn{\epsilon = 10^{-10}}) to
#' any zero-valued auxiliary terms \eqn{C} or \eqn{D}. A warning is issued
#' whenever such a correction is applied. Users should consider this
#' adjustment carefully and, if necessary, resort to alternative methods
#' (e.g., bootstrap) for items with boundary confidence limits.
#'
#' @references
#' Donner, A., & Zou, G. Y. (2012). Closed-form confidence intervals for
#' functions of the normal mean and standard deviation. \emph{Statistical
#' Methods in Medical Research}, 21(4), 347-359.
#'
#' Zou, G., Donner, A., & Qiu, S. (2025). MOVER-R for Confidence Intervals
#' of Ratios. In N. Balakrishnan, T. Colton, B. Everitt, W. Piegorsch,
#' F. Ruggeri and J.L. Teugels (Eds.), \emph{Wiley StatsRef: Statistics
#' Reference Online}. doi:10.1002/9781118445112.stat08085
#'
#' @examples
#' \dontrun{
#' # Suppose 'group1' and 'group2' contain columns 'est', 'low', and 'high'
#' res <- CIR(group1, group2, coef.col = "est", lwr.col = "low", upr.col = "high")
#' print(res)
#' }
#'
#' @export
CIR <- function(group1, group2, coef.col, lwr.col, upr.col, na.rm = FALSE) {

  # ----------------------------------------------------------------------
  # 1. Column validation
  # ----------------------------------------------------------------------

  # Ensure that the required columns exist in both data frames
  if (!all(c(coef.col, lwr.col, upr.col) %in% colnames(group1))) {
    stop("'group1' must contain the columns specified in 'coef.col', 'lwr.col', and 'upr.col'.")
  }
  if (!all(c(coef.col, lwr.col, upr.col) %in% colnames(group2))) {
    stop("'group2' must contain the columns specified in 'coef.col', 'lwr.col', and 'upr.col'.")
  }

  # Identify the relevant columns (including the matching key "Item")
  cols_g1 <- intersect(c(coef.col, lwr.col, upr.col, "Item"), colnames(group1))
  cols_g2 <- intersect(c(coef.col, lwr.col, upr.col, "Item"), colnames(group2))

  # ----------------------------------------------------------------------
  # 2. Handling of missing values
  # ----------------------------------------------------------------------

  has_na_g1 <- any(is.na(group1[, cols_g1, drop = FALSE]))
  has_na_g2 <- any(is.na(group2[, cols_g2, drop = FALSE]))

  if (has_na_g1 || has_na_g2) {
    if (!na.rm) {
      stop("Missing values detected. Use 'na.omit()' first or set 'na.rm = TRUE'.")
    } else {
      # Subset to the estimation columns (excluding "Item") for complete.cases
      est_cols_g1 <- intersect(c(coef.col, lwr.col, upr.col), colnames(group1))
      est_cols_g2 <- intersect(c(coef.col, lwr.col, upr.col), colnames(group2))

      keep_g1 <- complete.cases(group1[, est_cols_g1, drop = FALSE])
      keep_g2 <- complete.cases(group2[, est_cols_g2, drop = FALSE])

      nrem1 <- sum(!keep_g1)
      nrem2 <- sum(!keep_g2)

      group1 <- group1[keep_g1, , drop = FALSE]
      group2 <- group2[keep_g2, , drop = FALSE]

      warning(sprintf(
        "na.rm = TRUE: removed %d rows with missing values from 'group1' and %d from 'group2'.",
        nrem1, nrem2
      ))
    }
  }

  # ----------------------------------------------------------------------
  # 3. Matching items across groups
  # ----------------------------------------------------------------------

  common_items <- intersect(group1$Item, group2$Item)

  if (length(common_items) < nrow(group1) || length(common_items) < nrow(group2)) {
    warning(
      "Some items in 'group1' and 'group2' do not match. ",
      "Only the common items will be processed."
    )
  }

  # Retain only the common items and sort for deterministic merging
  group1 <- group1[group1$Item %in% common_items, ]
  group2 <- group2[group2$Item %in% common_items, ]
  group1 <- group1[order(group1$Item), ]
  group2 <- group2[order(group2$Item), ]

  # ----------------------------------------------------------------------
  # 4. Combine the two data frames
  # ----------------------------------------------------------------------

  combined <- merge(group1, group2, by = "Item", suffixes = c("_g1", "_g2"))

  # Pre-allocate the results data frame
  results <- data.frame(
    Item   = combined$Item,
    R      = NA_real_,
    lwr.ci = NA_real_,
    upr.ci = NA_real_,
    stringsAsFactors = FALSE
  )

  # Construct column names for easier access within the loop
  coef_g1 <- paste0(coef.col, "_g1")
  coef_g2 <- paste0(coef.col, "_g2")
  lwr_g1  <- paste0(lwr.col, "_g1")
  upr_g1  <- paste0(upr.col, "_g1")
  lwr_g2  <- paste0(lwr.col, "_g2")
  upr_g2  <- paste0(upr.col, "_g2")

  # ----------------------------------------------------------------------
  # 5. MOVER-R interval computation (Donner & Zou, 2012)
  # ----------------------------------------------------------------------

  # Numerical tolerance for detecting boundary cases
  EPS <- 1e-10

  for (i in seq_len(nrow(combined))) {

    # Retrieve the estimates and their confidence limits
    th1 <- combined[[coef_g1]][i]
    th2 <- combined[[coef_g2]][i]
    l1  <- combined[[lwr_g1]][i]
    u1  <- combined[[upr_g1]][i]
    l2  <- combined[[lwr_g2]][i]
    u2  <- combined[[upr_g2]][i]

    # Point estimate of the ratio
    R <- th1 / th2

    # Auxiliary quantities for MOVER-R
    # A = l1 * (2*th1 - l1)
    # B = u1 * (2*th1 - u1)
    # C = l2 * (2*th2 - l2)
    # D = u2 * (2*th2 - u2)
    A <- l1 * (2 * th1 - l1)
    B <- u1 * (2 * th1 - u1)
    C <- l2 * (2 * th2 - l2)
    D <- u2 * (2 * th2 - u2)

    # --------------------------------------------------------------------
    # Small correction for singularities (C or D equal to zero)
    # This frequently occurs when the confidence limits of a bounded
    # proportion (e.g., Aiken's V) are exactly 0 or 1.
    # --------------------------------------------------------------------
    if (C <= EPS) {
      warning(sprintf(
        "Item '%s': auxiliary term C = %.10f (l2 = %.4f, th2 = %.4f). ",
        "Setting C to %.1e to avoid division by zero.",
        combined$Item[i], C, l2, th2, EPS
      ))
      C <- EPS
    }

    if (D <= EPS) {
      warning(sprintf(
        "Item '%s': auxiliary term D = %.10f (u2 = %.4f, th2 = %.4f). ",
        "Setting D to %.1e to avoid division by zero.",
        combined$Item[i], D, u2, th2, EPS
      ))
      D <- EPS
    }

    # Product m = theta1 * theta2
    m <- th1 * th2

    # -------------------- Lower bound (L) --------------------
    # L = (m - sqrt(m^2 - A*D)) / D, provided that A*D <= m^2
    rad_low <- m^2 - A * D

    # Correct for minuscule negative radicands due to rounding
    if (rad_low < 0 && abs(rad_low) < 1e-12) {
      rad_low <- 0
    }

    if (rad_low < 0) {
      warning(sprintf(
        "Item '%s': negative radicand (%.10f) for the lower bound. Setting lwr.ci to NA.",
        combined$Item[i], rad_low
      ))
      L <- NA_real_
    } else {
      L <- (m - sqrt(rad_low)) / D
      # If the lower bound is trivially negative (due to rounding), set to 0
      if (!is.na(L) && L < 0 && abs(L) < 1e-10) L <- 0
    }

    # -------------------- Upper bound (U) --------------------
    # U = (m + sqrt(m^2 - B*C)) / C, provided that B*C <= m^2
    rad_up <- m^2 - B * C

    if (rad_up < 0 && abs(rad_up) < 1e-12) {
      rad_up <- 0
    }

    if (rad_up < 0) {
      warning(sprintf(
        "Item '%s': negative radicand (%.10f) for the upper bound. Setting upr.ci to NA.",
        combined$Item[i], rad_up
      ))
      U <- NA_real_
    } else {
      U <- (m + sqrt(rad_up)) / C
    }

    # --------------------------------------------------------------------
    # Store the results with rounding to three decimal places
    # --------------------------------------------------------------------
    results$R[i]      <- round(R, 3)
    results$lwr.ci[i] <- round(L, 3)
    results$upr.ci[i] <- round(U, 3)
  }

  return(results)
}
