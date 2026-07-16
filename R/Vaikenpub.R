#' @title Vaikenpub: Confidence Intervals for Published Aiken's V
#'
#' @description
#' Calculates Wilson asymmetric confidence intervals for reported Aiken's V coefficients,
#' assuming a fixed number of raters and known scale range. Based on the method of Penfield & Giacobbi (2004),
#' which treats Aiken's V as a proportion.
#'
#' @param V Numeric vector of published Aiken's V values (between 0 and 1).
#' @param n Integer. Number of judges or raters who rated each item.
#' @param lo Numeric. Minimum value of the rating scale.
#' @param hi Numeric. Maximum value of the rating scale.
#' @param conf.level Confidence level for the Wilson interval (default = 0.95).
#' @param item.names Optional. Vector of item names. If NULL, defaults to Item1, Item2, etc.
#'
#' @return A data.frame with item names, V values, and lower and upper confidence intervals.
#'
#' @references
#' Penfield, R. D., & Giacobbi, P. R., Jr. (2004). Applying a score confidence interval to Aiken’s item content-relevance index. 
#' \emph{Measurement in Physical Education and Exercise Science, 8}(4), 213–225.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. 
#' \emph{Journal of the American Statistical Association, 22}, 209–212.
#'
#' @examples
#' \donttest{
#' Vaikenpub(
#'   V = c(0.90, 0.75, 0.60),
#'   n = 10,
#'   lo = 1, hi = 5,
#'   conf.level = 0.95,
#'   item.names = c("Item A", "Item B", "Item C"))
#' }
#'
#' @export
Vaikenpub <- function(V, n, lo, hi, conf.level = 0.90, item.names = NULL) {
  # Wilson interval internal
  get_wilson_CI <- function(p_hat, N, conf.level) {
    crit <- qnorm(1 - (1 - conf.level) / 2)
    omega <- N / (N + crit^2)
    A <- p_hat + crit^2 / (2 * N)
    B <- crit * sqrt(p_hat * (1 - p_hat) / N + crit^2 / (4 * N^2))
    c(lower = omega * (A - B), upper = omega * (A + B))
  }
  
  if (any(V < 0 | V > 1)) stop("All V values must be between 0 and 1.")
  if (is.null(item.names)) item.names <- paste0("Item", seq_along(V))
  
  k <- hi - lo
  N <- n * k
  
  lwr <- upr <- numeric(length(V))
  
  for (i in seq_along(V)) {
    CI <- get_wilson_CI(V[i], N, conf.level)
    lwr[i] <- CI["lower"]
    upr[i] <- CI["upper"]
  }
  
  data.frame(
    Item = item.names,
    V = round(V, 3),
    lwr.ci = round(lwr, 3),
    upr.ci = round(upr, 3)
  )
}
