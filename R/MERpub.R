#' @title Mean of expert ratings from published summaries (MERpub)
#' @description Calculates the asymmetric score confidence interval (Penfield, 2003)
#'              for the mean of expert ratings, using only the mean rating and the
#'              number of raters per item (commonly reported in published studies).
#' @param M Numeric vector of mean ratings for each item.
#' @param n Integer vector of number of raters for each item (same length as M).
#' @param ncat Number of possible values or categories used in rating (e.g., 5 for a
#'             1-to-5 scale, 7 for a 0-to-6 scale).
#' @param start Minimum possible value of the rating scale (0 or 1).
#' @param conf.level Confidence level for confidence intervals (e.g., .90, .95, .99).
#' @param na.rm Logical. If FALSE (default) the function stops when missing values
#'              are detected in 'M' or 'n'. If TRUE, missing values are removed
#'              pairwise (i.e., rows with NA in either vector are dropped).
#'
#' @return
#' dataframe with columns: Item (sequential identifier), MER (mean expert rating),
#' lwr.ci (lower limit of the confidence interval), upr.ci (upper limit).
#'
#' @details
#' This function implements the Wilson score method generalized by Penfield (2003)
#' for constructing asymmetric confidence intervals around the mean of a rating
#' scale item. It is specifically designed for use with published summary data
#' (i.e., when only means and sample sizes are available), making it ideal for
#' meta-analytic purposes or secondary analyses of content validity studies.
#'
#' The score confidence interval does not require the standard deviation of the
#' ratings; it derives the standard error from the mean and the number of raters,
#' assuming the data follow a bounded binomial‑like distribution. This approach
#' is more appropriate for content validity ratings, which typically involve small
#' numbers of experts (often < 10) and few response categories.
#'
#' For an item with mean M and n raters, the interval is computed as:
#' \itemize{
#'   \item Compute p = M/k if start = 0, or p = (M-1)/k if start = 1, where k = ncat - 1.
#'   \item Obtain π_L and π_U via Wilson score limits for the proportion p.
#'   \item Compute LCL = M - z * sqrt(k * π_L * (1 - π_L) / n)
#'   \item Compute UCL = M + z * sqrt(k * π_U * (1 - π_U) / n)
#' }
#' The resulting limits are truncated to the possible scale range [start, start + ncat - 1].
#'
#' @references
#' Penfield, R. D. (2003). A score method of constructing asymmetric confidence
#' intervals for the mean of a rating scale item. Psychological Methods, 8(2), 149-163.
#' \url{https://doi.org/10.1037/1082-989x.8.2.149}
#'
#' Penfield, R. D., & Miller, J. M. (2004). Improving content validation studies
#' using an asymmetric confidence interval for the mean of expert ratings.
#' Applied Measurement in Education, 17(4), 359–370.
#' \url{https://doi.org/10.1207/s15324818ame1704_2}
#'
#' Miller, J. M., & Penfield, R. D. (2005). Using the score method to construct
#' asymmetric confidence intervals: An SAS program for content validation in
#' scale development. Behavior Research Methods, 37, 450-452.
#' \url{https://doi.org/10.3758/BF03192713}
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical
#' inference. Journal of the American Statistical Association, 22, 209-212.
#' \url{https://doi.org/10.2307/2276774}
#'
#' @seealso
#' \code{\link{MER}} for the version that works with raw ratings (dataframe of judges).
#'
#' @examples
#' \dontest{
#' # Example with vectors directly
#' medias <- c(4.2, 3.8, 4.5)
#' ns     <- c(8, 10, 7)
#'
#' # Compute 90% asymmetric confidence intervals (scale from 1 to 5)
#' MERpub(M = medias, n = ns, ncat = 5, start = 1, conf.level = .90)
#'
#' # For a scale from 0 to 4 (start = 0)
#' MERpub(M = medias, n = ns, ncat = 5, start = 0, conf.level = .95)
#'
#' # Using columns from a data frame (with dplyr or base R)
#' df <- data.frame(item = c("A","B","C"), mean = c(4.1,3.9,4.6), raters = c(9,8,10))
#' MERpub(M = df$mean, n = df$raters, ncat = 5, start = 1)
#' }
#'
#' @export
MERpub <- function(M, n, ncat, start, conf.level = 0.90, na.rm = FALSE) {

  # Check that M and n are numeric vectors of the same length
  if (!is.numeric(M)) stop("'M' must be a numeric vector")
  if (!is.numeric(n)) stop("'n' must be a numeric vector")
  if (length(M) != length(n)) stop("'M' and 'n' must have the same length")

  # Handle missing values
  if (!na.rm) {
    if (any(is.na(M)) || any(is.na(n))) {
      stop("Missing values detected in 'M' or 'n'. Use na.rm=TRUE to remove them.")
    }
  } else {
    # Pairwise removal
    valid <- !is.na(M) & !is.na(n)
    M <- M[valid]
    n <- n[valid]
    if (length(M) == 0) stop("No valid observations after removing NAs")
  }

  # Other checks
  if (!is.numeric(ncat) || ncat <= 1) stop("'ncat' must be a number greater than 1")
  if (!is.numeric(start) || !start %in% c(0, 1)) stop("'start' must be 0 or 1")
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be between 0 and 1")
  }

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  k <- ncat - 1
  min_val <- start
  max_val <- start + ncat - 1

  # Compute for each pair (M, n)
  results <- mapply(function(mean_val, n_rater) {
    if (n_rater < 2) {
      warning("Item has n < 2; confidence interval may be unreliable.")
    }

    # Compute p
    if (start == 0) {
      p <- mean_val / k
    } else {
      p <- (mean_val - 1) / k
    }
    p <- max(1e-10, min(1 - 1e-10, p))

    term1 <- 2 * p * n_rater * k + z^2
    term2 <- z * sqrt(4 * n_rater * k * p * (1 - p) + z^2)
    denom <- 2 * (n_rater * k + z^2)

    pi_L <- (term1 - term2) / denom
    pi_U <- (term1 + term2) / denom

    LCL <- mean_val - z * sqrt(k * pi_L * (1 - pi_L) / n_rater)
    UCL <- mean_val + z * sqrt(k * pi_U * (1 - pi_U) / n_rater)

    LCL <- max(LCL, min_val)
    UCL <- min(UCL, max_val)

    c(MER = round(mean_val, 3),
      lwr.ci = round(LCL, 3),
      upr.ci = round(UCL, 3))
  }, M, n, SIMPLIFY = FALSE)

  # Combine results into a data frame
  results_df <- as.data.frame(do.call(rbind, results))
  results_df$Item <- seq_len(nrow(results_df))
  results_df <- results_df[, c("Item", "MER", "lwr.ci", "upr.ci")]

  return(results_df)
}
