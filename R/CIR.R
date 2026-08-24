#' MOVER-R Confidence Interval for the Ratio of Content Validity Coefficients
#'
#' Computes the confidence interval for the ratio of two independent
#' content validity coefficients  using the MOVER-R (Method of Variance Estimates Recovery for Ratios)
#' closed-form procedure.
#'
#' @param group1 data.frame containing estimates and CIs for the first group.
#'   Must include columns specified in \code{coef.col}, \code{lwr.col},
#'   \code{upr.col}, and an \code{"Item"} column for matching.
#' @param group2 data.frame containing estimates and CIs for the second group.
#'   Same column requirements as \code{group1}.
#' @param coef.col Character string. Name of the column with the point estimates.
#' @param lwr.col Character string. Name of the column with the lower bounds.
#' @param upr.col Character string. Name of the column with the upper bounds.
#' @param na.rm Logical. If \code{TRUE}, rows with missing values in the relevant
#'   columns are removed with a warning. If \code{FALSE}, missing values
#'   trigger an error.
#'
#' @return A data.frame with columns:
#'   \item{Item}{Common item names between both groups.}
#'   \item{R}{Estimated ratio (coefficient group1 / coefficient group2).}
#'   \item{lwr.ci}{Lower bound of the MOVER-R confidence interval.}
#'   \item{upr.ci}{Upper bound of the MOVER-R confidence interval.}
#'
#' @references
#' Donner, A., & Zou, G. Y. (2012). Closed-form confidence intervals for
#' functions of the normal mean and standard deviation. \emph{Statistical
#' Methods in Medical Research}, 21(4), 347-359.
#'
#' Zou, G., Donner, A., & Qiu, S. (2025). MOVER-R for Confidence Intervals
#' of Ratios. In N. Balakrishnan, T. Colton, B. Everitt, W. Piegorsch,
#' F. Ruggeri and J.L. Teugels (Eds.), \emph{Wiley StatsRef: Statistics
#' Reference Online}. https://doi.org/10.1002/9781118445112.stat08085
#'
#' @examples
#' \dontrun{
#' # Assuming you have data frames 'g1' and 'g2' with columns 'est', 'low', 'high'
#' result <- CIR(g1, g2, coef.col = "est", lwr.col = "low", upr.col = "high")
#' print(result)
#' }
#' @export
CIR <- function(group1, group2, coef.col, lwr.col, upr.col, na.rm = FALSE) {

  # ------ 1. Validación de columnas ------
  if (!all(c(coef.col, lwr.col, upr.col) %in% colnames(group1))) {
    stop("group1 must contain the columns specified in coef.col, lwr.col, and upr.col")
  }
  if (!all(c(coef.col, lwr.col, upr.col) %in% colnames(group2))) {
    stop("group2 must contain the columns specified in coef.col, lwr.col, and upr.col")
  }

  # Columns of interest (including "Item" for the merge)
  cols_g1 <- intersect(c(coef.col, lwr.col, upr.col, "Item"), colnames(group1))
  cols_g2 <- intersect(c(coef.col, lwr.col, upr.col, "Item"), colnames(group2))

  # ------ 2. Handling Missing Values (NA) ------
  has_na_g1 <- any(is.na(group1[, cols_g1, drop = FALSE]))
  has_na_g2 <- any(is.na(group2[, cols_g2, drop = FALSE]))

  if (has_na_g1 || has_na_g2) {
    if (!na.rm) {
      stop("Missing values detected. Use na.omit() first or set na.rm=TRUE")
    } else {

      # Only the "Estimate" and "CI" columns (not "Item") for complete cases
      est_cols_g1 <- intersect(c(coef.col, lwr.col, upr.col), colnames(group1))
      est_cols_g2 <- intersect(c(coef.col, lwr.col, upr.col), colnames(group2))

      keep_g1 <- complete.cases(group1[, est_cols_g1, drop = FALSE])
      keep_g2 <- complete.cases(group2[, est_cols_g2, drop = FALSE])

      nrem1 <- sum(!keep_g1)
      nrem2 <- sum(!keep_g2)

      group1 <- group1[keep_g1, , drop = FALSE]
      group2 <- group2[keep_g2, , drop = FALSE]

      warning(sprintf("na.rm=TRUE: removed %d rows with missing values from group1 and %d from group2",
                      nrem1, nrem2))
    }
  }

  # ------ 3. Intersection of Items ------
  common_items <- intersect(group1$Item, group2$Item)
  if (length(common_items) < nrow(group1) || length(common_items) < nrow(group2)) {
    warning("Some items do not match between group1 and group2. Only the common items will be processed..")
  }

  group1 <- group1[group1$Item %in% common_items, ]
  group2 <- group2[group2$Item %in% common_items, ]

  # Sort by Item to ensure a match (optional, but safe)
  group1 <- group1[order(group1$Item), ]
  group2 <- group2[order(group2$Item), ]

  # ------ 4. Mixing (merge) ------
  combined <- merge(group1, group2, by = "Item", suffixes = c("_g1", "_g2"))

  # Pre-assign results
  results <- data.frame(
    Item = combined$Item,
    R    = NA_real_,
    lwr.ci = NA_real_,
    upr.ci = NA_real_,
    stringsAsFactors = FALSE
  )

  # Column names in the combined data.frame
  coef_g1 <- paste0(coef.col, "_g1")
  coef_g2 <- paste0(coef.col, "_g2")
  lwr_g1  <- paste0(lwr.col, "_g1")
  upr_g1  <- paste0(upr.col, "_g1")
  lwr_g2  <- paste0(lwr.col, "_g2")
  upr_g2  <- paste0(upr.col, "_g2")

  # ------ 5. Calculation of the MOVER-R Interval (Donner & Zou, 2012) ------
  for (i in 1:nrow(combined)) {

    # Extraer valores
    th1 <- combined[[coef_g1]][i]
    th2 <- combined[[coef_g2]][i]
    l1  <- combined[[lwr_g1]][i]
    u1  <- combined[[upr_g1]][i]
    l2  <- combined[[lwr_g2]][i]
    u2  <- combined[[upr_g2]][i]

    # Point ratio
    R <- th1 / th2

    # --- MOVER-R Components ---
    # A = l1*(2*th1 - l1)  ;  B = u1*(2*th1 - u1)
    # C = l2*(2*th2 - l2)  ;  D = u2*(2*th2 - u2)
    A <- l1 * (2 * th1 - l1)
    B <- u1 * (2 * th1 - u1)
    C <- l2 * (2 * th2 - l2)
    D <- u2 * (2 * th2 - u2)

    # Product m = theta1 * theta2
    m <- th1 * th2

    # Initialize limits as NA
    L <- NA_real_
    U <- NA_real_

    # --- Lower limit ---
    # D must be greater than 0, and the radicand must be greater than or equal to 0.
    if (D > 0) {
      rad_low <- m^2 - A * D
      if (rad_low < 0) {
        # If the radicand is slightly negative due to rounding, we set it to 0
        if (abs(rad_low) < 1e-12) rad_low <- 0
      }
      if (rad_low >= 0) {
        L <- (m - sqrt(rad_low)) / D
        # By the definition of Aiken's V (ratio), it should not be negative,
        # but if it is, we set it to 0 (conservative approach).
        if (L < 0 && L > -1e-10) L <- 0
      } else {
        warning(sprintf("Negative value found for the lower bound in Item '%s'. NA is assigned.",
                        combined$Item[i]))
      }
    } else {
      warning(sprintf("D <= 0 for the item '%s'. The lower bound cannot be calculated (D = %f).",
                      combined$Item[i], D))
    }

    # --- Upper limit ---
    # C > 0 and the square root is >= 0 are required
    if (C > 0) {
      rad_up <- m^2 - B * C
      if (rad_up < 0) {
        if (abs(rad_up) < 1e-12) rad_up <- 0
      }
      if (rad_up >= 0) {
        U <- (m + sqrt(rad_up)) / C
      } else {
        warning(sprintf("Negative value found for the upper limit in Item '%s'. NA is assigned.",
                        combined$Item[i]))
      }
    } else {
      warning(sprintf("C <= 0 for the item '%s'. The upper bound cannot be calculated (C = %f).",
                      combined$Item[i], C))
    }

    # Assign results (rounded to 3 decimal places)
    results$R[i]      <- round(R, 3)
    results$lwr.ci[i] <- round(L, 3)
    results$upr.ci[i] <- round(U, 3)
  }

  return(results)
}
