#' @title CVIpub: Confidence Intervals for CVI or CVI.R published
#'
#' @description
#' Calculates Wilson asymmetric confidence intervals for published Content Validity Index (CVI),
#' with optional correction for chance agreement based on Polit et al. (2007).
#'
#' When \code{correct = TRUE}, the function computes the adjusted coefficient \code{CVI.R}
#' and applies the Wilson method to its absolute value. For computational stability,
#' CVI values of exactly 1.0 or 0.0 are replaced by 0.9999 and 0.0001 respectively.
#' If an extreme negative value of CVI.R is detected (e.g., from CVI ≈ 0), the value is
#' truncated to 0 and the confidence interval is not estimated.
#'
#' @param cvi Numeric vector of CVI values (between 0 and 1).
#' @param n Integer or numeric vector. Number of judges per item.
#' @param conf.level Confidence level for Wilson interval (e.g., 0.90, 0.95, 0.99). Default is 0.95.
#' @param item.names Optional vector of item names. If NULL, defaults to Item1, Item2, etc.
#' @param correct Logical. If TRUE, calculates CVI.R (Polit correction) and its Wilson CI.
#'
#' @return A data.frame with:
#' \itemize{
#'   \item If \code{correct = FALSE}: CVI, lower and upper Wilson confidence bounds.
#'   \item If \code{correct = TRUE}: CVI.R, lower and upper Wilson confidence bounds (or NA if truncated).
#' }
#'
#' @examples
#' \donttest{
#' CVIpub(cvi = c(.90, .8, .9, .9),
#' n = c(10, 10, 10, 10),
#' conf.level = .90)
#' }
#' 
#' @export
CVIpub <- function(cvi, n, conf.level = 0.95, item.names = NULL, correct = FALSE) {
  # Internal Wilson CI function
  get_wilson_CI <- function(p_hat, n, conf.level) {
    crit <- qnorm(1 - (1 - conf.level) / 2)
    omega <- n / (n + crit^2)
    A <- p_hat + crit^2 / (2 * n)
    B <- crit * sqrt(p_hat * (1 - p_hat) / n + crit^2 / (4 * n^2))
    CI <- c(lower = omega * (A - B), upper = omega * (A + B))
    return(CI)
  }
  
  # Validation
  if (any(cvi < 0 | cvi > 1)) stop("All CVI values must be between 0 and 1.")
  if (length(n) == 1) n <- rep(n, length(cvi))
  if (length(n) != length(cvi)) stop("Length of 'n' must match length of 'cvi'.")
  if (is.null(item.names)) item.names <- paste0("Item", seq_along(cvi))
  
  # Replace exact 0 and 1 to avoid instability
  cvi.safe <- ifelse(cvi == 1, 0.9999, ifelse(cvi == 0, 0.0001, cvi))
  
  # Output vectors
  est <- lwr <- upr <- numeric(length(cvi))
  truncation_reported <- FALSE
  
  for (i in seq_along(cvi)) {
    if (correct) {
      A <- round(cvi.safe[i] * n[i])
      
      # Cálculo estable usando log-factorial
      log_pc <- lchoose(n[i], A) + n[i] * log(0.5)
      Pc <- exp(log_pc)
      
      cvir <- (cvi.safe[i] - Pc) / (1 - Pc)
      
      # Check for extreme negative or undefined values
      if (!is.finite(cvir) || cvir < 0) {
        cvir <- 0
        lwr[i] <- upr[i] <- NA
        if (!truncation_reported) {
          message("Extreme negative CVI.R value detected (CVI ≈ 0). Truncated to 0; confidence interval not estimated.")
          truncation_reported <- TRUE
        }
      } else {
        abs_cvir <- abs(cvir)
        CI <- get_wilson_CI(abs_cvir, n[i], conf.level)
        lwr[i] <- if (cvir < 0) -CI["upper"] else CI["lower"]
        upr[i] <- if (cvir < 0) -CI["lower"] else CI["upper"]
      }
      
      est[i] <- cvir
    } else {
      est[i] <- cvi[i]
      CI <- get_wilson_CI(cvi[i], n[i], conf.level)
      lwr[i] <- CI["lower"]
      upr[i] <- CI["upper"]
    }
  }
  
  # Assemble output
  result <- data.frame(
    Item = item.names,
    coef = round(est, 3),
    lwr.ci = round(lwr, 3),
    upr.ci = round(upr, 3)
  )
  
  names(result)[2] <- if (correct) "CVI.R" else "CVI"
  return(result)
}
