#' @title Aiken's V coefficient
#' @description Calculate the V coefficient, known as Aiken's V.
#'
#' @param data dataframe, with the columns assigned to each rater, and the rows assigned to each evaluated item.
#' @param min minimum possible rating value
#' @param max maximum possible rating value
#' @param conf.level confidence level for confidence intervals (ex., .90, .95, .99)
#' @param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#' @param overall Logical. If TRUE, adds a row with the overall V coefficient.
#' @param overall.method Character. Method for the overall V:
#'        \itemize{
#'          \item \code{"global"}: treats the entire matrix as a single "super-item" (all ratings pooled).
#'          \item \code{"Aiken"}: computes V for each judge (based on all items) and then averages them.
#'        }
#' @param overall.ci Character. Method for the overall confidence interval:
#'        \itemize{
#'          \item \code{"Wilson"}: Wilson score interval (Penfield & Giacobbi, 2004).
#'          \item \code{"MER"}: score interval for the mean of ratings, then transformed to V (Penfield, 2003).
#'        }
#'        Note: When \code{overall.method = "Aiken"}, only \code{"Wilson"} is used (others ignored).
#'
#' @return dataframe with V coefficients for all items, and confidence intervals.
#'         If overall = TRUE, an extra row with the overall index is included.
#'
#' @details
#' The overall V provides a global estimate of content validity for the whole instrument.
#' Two aggregation methods are available:
#' \itemize{
#'   \item \code{"global"}: all ratings are pooled (treating the matrix as one super-item).
#'         The V is the mean of all transformed scores.
#'   \item \code{"Aiken"}: V is computed for each judge (using their ratings across all items)
#'         and then averaged. This follows the large-sample procedure described by Aiken (1985).
#' }
#' For confidence intervals, the Wilson score method (Penfield & Giacobbi, 2004) is available for
#' both approaches. For the \code{"global"} method, the MER method (Penfield, 2003) is also offered,
#' which constructs the CI on the mean of ratings and then transforms to V.
#'
#' **Overall V and its confidence interval**
#'
#' When `overall = TRUE`, the function treats the entire matrix of ratings (all items × all judges)
#' as a single "super-item". The overall V coefficient is computed as the mean of all transformed
#' scores `(rating - min) / (max - min)`, which is equivalent to `(mean(all_ratings) - min) / (max - min)`.
#' This provides a global estimate of content validity for the whole instrument.
#' Two methods are available to construct the asymmetric confidence interval for this overall V:
#'
#'  - **`overall_method = "Wilson"`** (default): Applies the Wilson score interval directly to the overall
#'  proportion V, following the logic of Penfield & Giacobbi (2004). The effective sample size
#'  is `n_judges × (max - min)`, respecting the independence of judges. This method is the
#'  most direct extension of the standard item-level V confidence interval to the global index, and
#'  is consistent with the original proposal by Aiken and later developments.
#'
#'  - **`overall_method = "mer"`**: Implements the score confidence interval for the **mean** of
#'  the ratings (MER; Penfield, 2003; Penfield & Miller, 2004) in the original scale, and then
#'  transforms the lower and upper bounds to the V metric. This approach first computes an
#'  asymmetric CI for the mean `M` of all ratings, and then applies the linear
#'  transformation `(CI - min) / (max - min)`. By using the mean as the primary parameter, this
#'  method explicitly models the variability among judges and does not rely on the binomial expansion
#'  used in the Wilson method.
#'
#'  Both methods yield a V total that is identical in point estimate, but the confidence intervals may
#'  differ slightly. In either case, it is recommended to report the method used and,
#'  if possible, to provide both intervals in supplementary materials for transparency.
#'
#' @references
#' Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings.
#'     Educational and Psychological Measurement, 45, 131-142.
#' Penfield, R. D. (2003). A score method of constructing asymmetric confidence intervals
#'     for the mean of a rating scale item. Psychological Methods, 8(2), 149-163.
#' Penfield, R. D. & Giacobbi, P. R., Jr. (2004) Applying a score confidence interval
#'     to Aiken’s item content-relevance index. Measurement in Physical Education and Exercise Science, 8(4), 213-225.
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference.
#'     Journal of the American Statistical Association, 22, 209-212.
#'
#' @export
#'
#' @examples
#' data2Tst <- data.frame(
#'   J1 = c(4, 1, 1, 1, 4),
#'   J2 = c(4, 1, 2, 2, 3),
#'   J3 = c(4, 1, 3, 3, 5),
#'   J4 = c(4, 1, 4, 5, 5),
#'   J5 = c(4, 1, 5, 5, 5),
#'   J6 = c(4, 1, 3, 5, 5))
#'
#' # Original: item-level V with Wilson CI
#' Vaiken(data2Tst, min = 1, max = 5, conf.level = .90)
#'
#' # Overall V with global pooling, Wilson CI
#' Vaiken(data2Tst, min = 1, max = 5, overall = TRUE, overall.method = "global", overall.ci = "Wilson")
#'
#' # Overall V with Aiken's method (average of judge V's), Wilson CI
#' Vaiken(data2Tst, min = 1, max = 5, overall = TRUE, overall.method = "Aiken")
Vaiken <- function(data, min, max, conf.level = 0.95, na.rm = FALSE,
                   overall = FALSE,
                   overall.method = c("global", "Aiken"),
                   overall.ci = c("Wilson", "MER")) {

  # ---- Missing values ----
  if (!na.rm) {
    if (any(is.na(data))) {
      stop("Missing values detected. Use na.omit() first or set na.rm=TRUE.")
    }
  } else {
    data <- na.omit(data)
  }

  # ---- Input validation ----
  if (!is.data.frame(data)) stop("'data' must be a data.frame.")
  if (!is.numeric(as.matrix(data))) stop("'data' must contain only numeric values.")
  if (!is.numeric(min) || !is.numeric(max)) stop("'min' and 'max' must be numeric.")
  if (min >= max) stop("'min' must be less than 'max'.")
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be between 0 and 1.")
  if (any(data < min | data > max))
    stop("All ratings must be between 'min' and 'max'.")

  overall.method <- match.arg(overall.method)
  overall.ci <- match.arg(overall.ci)

  # If overall.method == "Aiken", force overall.ci = "Wilson" (with a message)
  if (overall.method == "Aiken" && overall.ci != "Wilson") {
    warning("For overall.method = 'Aiken', only 'Wilson' CI is available. Switching to 'Wilson'.")
    overall.ci <- "Wilson"
  }

  n_jueces <- ncol(data)
  n_items <- nrow(data)
  z <- qnorm(1 - (1 - conf.level) / 2)
  k <- max - min   # range of the scale

  # ---- Internal function for item-level V and Wilson CI (unchanged) ----
  calcular_valores <- function(puntajes) {
    M <- mean(puntajes)
    V <- (M - min) / k
    nk <- length(puntajes) * k   # n_jueces * k
    term1 <- 2 * nk * V + z^2
    term2 <- z * sqrt(4 * nk * V * (1 - V) + z^2)
    denom <- 2 * (nk + z^2)
    lwr.ci <- max(0, min(1, (term1 - term2) / denom))
    upr.ci <- max(0, min(1, (term1 + term2) / denom))
    return(c(V = V, lwr.ci = lwr.ci, upr.ci = upr.ci))
  }

  resultados <- t(apply(data, 1, calcular_valores))
  resultados_df <- data.frame(
    Item = 1:n_items,
    V = round(resultados[, "V"], 3),
    lwr.ci = round(resultados[, "lwr.ci"], 3),
    upr.ci = round(resultados[, "upr.ci"], 3)
  )

  # ---- Overall calculation (if requested) ----
  if (overall) {

    if (overall.method == "global") {
      # ---- Global super-item approach (pool all ratings) ----
      all_scores <- as.vector(as.matrix(data))
      V_total <- mean((all_scores - min) / k)   # mean of transformed scores

      if (overall.ci == "Wilson") {
        # Wilson CI directly on V_total, using n_jueces * k as effective N
        nk_total <- n_jueces * k
        term1_t <- 2 * nk_total * V_total + z^2
        term2_t <- z * sqrt(4 * nk_total * V_total * (1 - V_total) + z^2)
        denom_t <- 2 * (nk_total + z^2)
        lwr_total <- max(0, min(1, (term1_t - term2_t) / denom_t))
        upr_total <- max(0, min(1, (term1_t + term2_t) / denom_t))
        etiqueta <- "Total (global, Wilson)"

      } else { # overall.ci == "MER"
        # MER: Wilson CI on mean of ratings, then transform to V
        M_total <- mean(all_scores)
        p <- V_total   # same as (M - min)/k
        p <- max(1e-10, min(1 - 1e-10, p))

        n <- n_jueces
        term1_mer <- 2 * p * n * k + z^2
        term2_mer <- z * sqrt(4 * n * k * p * (1 - p) + z^2)
        denom_mer <- 2 * (n * k + z^2)
        pi_L <- (term1_mer - term2_mer) / denom_mer
        pi_U <- (term1_mer + term2_mer) / denom_mer

        LCL_mean <- M_total - z * sqrt(k * pi_L * (1 - pi_L) / n)
        UCL_mean <- M_total + z * sqrt(k * pi_U * (1 - pi_U) / n)
        LCL_mean <- max(min, min(max, LCL_mean))
        UCL_mean <- max(min, min(max, UCL_mean))

        lwr_total <- max(0, min(1, (LCL_mean - min) / k))
        upr_total <- max(0, min(1, (UCL_mean - min) / k))
        etiqueta <- "Total (global, MER)"
      }

    } else { # overall.method == "Aiken"
      # ---- Aiken's method: average of V's computed for each judge ----
      V_por_juez <- apply(data, 2, function(col) {
        (mean(col) - min) / k
      })
      V_total <- mean(V_por_juez)

      # Wilson CI on the mean V (treating it as a proportion)
      # Effective N = n_jueces * k
      nk_total <- n_jueces * k
      term1_t <- 2 * nk_total * V_total + z^2
      term2_t <- z * sqrt(4 * nk_total * V_total * (1 - V_total) + z^2)
      denom_t <- 2 * (nk_total + z^2)
      lwr_total <- max(0, min(1, (term1_t - term2_t) / denom_t))
      upr_total <- max(0, min(1, (term1_t + term2_t) / denom_t))
      etiqueta <- "Total (Aiken, Wilson)"
    }

    # Add overall row to results
    fila_total <- data.frame(
      Item = etiqueta,
      V = round(V_total, 3),
      lwr.ci = round(lwr_total, 3),
      upr.ci = round(upr_total, 3)
    )
    resultados_df <- rbind(resultados_df, fila_total)
  }

  return(resultados_df)
}
