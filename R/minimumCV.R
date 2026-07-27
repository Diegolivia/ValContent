#' Minimum Empirical Content Validity Coefficient Required for a Target Lower Bound
#'
#' Determines the minimum observed coefficient (e.g., Aiken's V, CVI, CVC) that a
#' panel of expert judges must achieve so that the lower bound of the Wilson
#' score confidence interval reaches a specified target value. This function
#' solves the inverse problem of the Wilson interval exactly, avoiding the
#' inaccuracies of asymptotic (Wald-type) approximations, which are particularly
#' problematic for small panels and proportions near the boundaries.
#'
#' @param target.value A numeric value in (0, 1) representing the desired lower
#'   bound of the confidence interval for the content validity coefficient. In
#'   practice, this is often set to 0.70 or 0.60 depending on the stringency of
#'   the validation context.
#' @param conf.level A numeric value in (0, 1) specifying the confidence level
#'   for the interval. Common values are 0.90 (conservative) or 0.80 (more
#'   liberal) in content validity research.
#' @param Nsize An integer vector specifying the panel sizes (i.e., number of
#'   expert judges) to evaluate. Defaults to `5:25`.
#' @param alternative A character string, either `"two.sided"` or `"one.sided"`.
#'   Specifies whether to use a two-sided or a one-sided (lower) confidence
#'   interval. For guaranteeing that the coefficient is at least `target.value`,
#'   a one-sided interval is conceptually appropriate and statistically less
#'   conservative. Defaults to `"two.sided"` for consistency with conventional
#'   reporting, though `"one.sided"` is recommended for sample size planning.
#' @param plot A logical value. If `TRUE`, generates a scatter plot displaying
#'   the required minimum coefficient as a function of panel size. Feasible
#'   solutions are shown as blue points; infeasible combinations (where the
#'   target cannot be achieved even with perfect agreement) are marked as red
#'   crosses. When `plot = TRUE`, the results data frame is also printed to the
#'   console, so you see both the graph and the numerical table. Defaults to
#'   `FALSE`.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{N}: The panel size evaluated.
#'   \item \code{Min.coef}: The minimum empirical coefficient required. If the
#'     target is unattainable for that panel size, the value is \code{NA}.
#'   \item \code{Factible}: A logical indicator (\code{TRUE/FALSE}) denoting
#'     whether the target lower bound is mathematically achievable given the
#'     panel size and confidence level.
#' }
#' The data frame is always returned (visibly) so that it can be assigned to a
#' variable for further analysis.
#'
#' @details
#' The Wilson score interval for an observed proportion \eqn{\hat{p}} has a lower
#' bound \eqn{L} given by:
#' \deqn{
#'   L = \frac{\hat{p} + \frac{z^2}{2n} - z\sqrt{\frac{\hat{p}(1-\hat{p})}{n} + \frac{z^2}{4n^2}}}{1 + \frac{z^2}{n}}
#' }
#' where \eqn{z} is the appropriate quantile of the standard normal distribution
#' (determined by \code{conf.level} and \code{alternative}), and \eqn{n} is the
#' number of expert judges.
#'
#' Setting \eqn{L = T} (where \eqn{T} is the \code{target.value}) and solving
#' for \eqn{\hat{p}} yields the closed-form solution:
#' \deqn{
#'   \hat{p}_{min} = T + z \cdot \sqrt{\frac{T(1-T)}{n}}
#' }
#'
#' This solution is valid only if \eqn{\hat{p}_{min} \le 1}. This condition is
#' equivalent to:
#' \deqn{
#'   T \le \frac{1}{1 + z^2/n}
#' }
#' When this condition is violated, even a perfect empirical coefficient
#' (\eqn{\hat{p} = 1}) would produce a confidence interval whose lower bound
#' falls below \eqn{T}. The function flags these cases as infeasible.
#'
#' @references
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#' statistical inference. \emph{Journal of the American Statistical Association},
#' 22(158), 209-212.
#'
#' Aiken, L. R. (1985). Three coefficients for analyzing the reliability and
#' validity of ratings. \emph{Educational and Psychological Measurement}, 45(1),
#' 131-142.
#'
#' @examples
#' \donttest{
#' # Basic usage: typical thresholds for content validity
#' minimumCV(target.value = 0.70, conf.level = 0.90,
#'           Nsize = c(5, 10, 15, 20))
#'           }
#'
#' \donttest{
#' # With one-sided interval (recommended for planning)
#' minimumCV(target.value = 0.60, conf.level = 0.80,
#'           Nsize = 3:15,
#'           alternative = "one.sided")
#'           }
#'
#' # Plot and table
#'  \donttest{
#'  minimumCV(target.value = 0.75, conf.level = 0.95,
#'           Nsize = 3:30, alternative = "one.sided", plot = TRUE)
#'           }
#'
#' @export
minimumCV <- function(target.value,
                      conf.level = 0.95,
                      Nsize = 5:25,
                      alternative = c("two.sided", "one.sided"),
                      plot = FALSE) {

  # --- Argument validation ---
  alternative <- match.arg(alternative)

  if (target.value <= 0 || target.value >= 1) {
    stop("The 'target.value' must be strictly between 0 and 1.")
  }
  if (conf.level <= 0 || conf.level >= 1) {
    stop("The 'conf.level' must be strictly between 0 and 1.")
  }
  if (!is.numeric(Nsize) || any(Nsize < 1)) {
    stop("'Nsize' must be a vector of positive integers.")
  }
  if (!is.logical(plot)) {
    stop("'plot' must be logical (TRUE or FALSE).")
  }

  # --- Compute the standard normal quantile ---
  if (alternative == "one.sided") {
    z <- qnorm(conf.level)
  } else { # two.sided
    z <- qnorm(1 - (1 - conf.level) / 2)
  }

  # --- Pre-allocate results ---
  resultados <- data.frame(
    N        = Nsize,
    Min.coef = as.numeric(NA),
    Factible = as.logical(NA)
  )

  # --- Main loop over panel sizes ---
  for (i in seq_along(Nsize)) {
    n <- Nsize[i]
    z2_over_n <- (z^2) / n

    if (target.value <= 1 / (1 + z2_over_n)) {
      p_min <- target.value + z * sqrt(target.value * (1 - target.value) / n)
      resultados$Min.coef[i] <- round(p_min, 4)
      resultados$Factible[i] <- TRUE
    } else {
      resultados$Min.coef[i] <- NA
      resultados$Factible[i] <- FALSE
    }
  }

  # --- Plotting (if requested) ---
  if (plot) {
    feasible_data <- resultados[resultados$Factible == TRUE, ]
    infeasible_data <- resultados[resultados$Factible == FALSE, ]

    y_max <- max(c(1, feasible_data$Min.coef), na.rm = TRUE)
    y_min <- min(c(target.value, feasible_data$Min.coef), na.rm = TRUE)
    y_range <- c(max(0, y_min - 0.05), min(1, y_max + 0.05))

    if (nrow(feasible_data) == 0) {
      plot(1, type = "n",
           xlim = range(Nsize),
           ylim = c(0, 1),
           xlab = "Number of Expert Judges (n)",
           ylab = "Minimum Required Empirical Coefficient",
           main = "Sample Size Planning for Content Validity (Wilson Inversion)")
      legend("topright", legend = c("Target Bound", "Infeasible"),
             lty = c(2, NA), pch = c(NA, 4), col = c("gray30", "red"), bty = "n")
      abline(h = target.value, lty = 2, col = "gray30")
      points(infeasible_data$N,
             rep(0.5, nrow(infeasible_data)),
             pch = 4, col = "red", cex = 1.2)
      text(infeasible_data$N,
           rep(0.5, nrow(infeasible_data)),
           labels = "Unattainable", pos = 3, cex = 0.7, col = "red")
    } else {
      plot(feasible_data$N,
           feasible_data$Min.coef,
           pch = 19,
           col = "steelblue",
           xlim = range(Nsize),
           ylim = y_range,
           xlab = "Number of Expert Judges (n)",
           ylab = "Minimum Required Empirical Coefficient",
           main = "Sample Size Planning for Content Validity (Wilson Inversion)")

      if (nrow(infeasible_data) > 0) {
        y_infeasible <- min(1, y_max + 0.02)
        points(infeasible_data$N,
               rep(y_infeasible, nrow(infeasible_data)),
               pch = 4, col = "red", cex = 1.2)
        text(infeasible_data$N,
             rep(y_infeasible, nrow(infeasible_data)),
             labels = "Unattainable", pos = 3, cex = 0.7, col = "red")
      }

      abline(h = target.value, lty = 2, col = "gray30")
      legend("topright",
             legend = c("Feasible (Wilson Inversion)", "Target Bound", "Infeasible"),
             pch = c(19, NA, 4),
             lty = c(NA, 2, NA),
             col = c("steelblue", "gray30", "red"),
             bty = "n")
    }
  }

  # --- Return results data frame (visible, so it prints) ---
  return(resultados)
}
