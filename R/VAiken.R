#' @title Aiken's V coefficient
#' @description Calculate the V coefficient, known as Aiken's V.
#'
#' @param data dataframe, with the columns assigned to each rater, and the rows assigned to each evaluated item.
#' @param min minimum possible rating value
#' @param max maximum possible rating value
#' @param conf.level confidence level for confidence intervals (ex., .90, .95, .99)
#' @param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#'
#'@return
#'dataframe with V coefficients for all items analyzed, and confidence intervals
#'
#'@details
#'Calculate the V coefficient (Aiken, 1980, 1985), with the formula of Penfield & Giacobbi (2004). It also calculates asymmetric confidence intervals (Wilson, 1927; Penfield & Giacobbi, 2004). The results should be complemented by an estimator of variability or inter-judge agreement. The function uses the modified formula presented by Penfield & Giacobbi (2004).
#'This function substantially improves on Vaiken (Merino, & Livia, 2009) because it calculates for multiple items and any confidence level.
#'
#'Note: The function has not yet been prepared to resolve missing values, so the user must remove or impute any missing values.
#'
#'
#'@references
#'Aiken, L. R. (1980). Content validity and reliability of single items or questionnaires. \emph{Educational and. Psychological Measurement, 40}, 955-959. \doi{10.1177/001316448004000419}
#'
#'Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. \emph{Educational and Psychological Measurement, 45}, 131-142. \doi{10.1177/0013164485451012}
#'
#'Merino, C., & Livia, J. (2009). Intervalos de confianza asimetricos para el indice de validez de contenido: un programa Visual Basic para la V de Aiken. \emph{Anales de Psicologia, 25}(1), 169-171. \url{https://revistas.um.es/analesps/article/view/71631}
#'
#'Penfield, R. D. & Giacobbi, P. R., Jr. (2004) Applying a score confidence interval to Aiken’s item content-relevance index. \emph{Measurement in Physical Education and Exercise Science, 8}(4), 213-225. \doi{10.1207/s15327841mpee0804_3}
#'
#'Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. \emph{Journal of the American Statistical Association, 22}, 209-212. \doi{10.2307/2276774}
#'
#'@seealso
#'\code{\link[PropCIs:scoreci]{PropCIs::scoreci}} for score method confidence interval
#'
#' @author
#' Diego Livia-Ortiz (\email{diegolivia@hotmail.com})
#' Cesar Merino-Soto (\email{sikayax@yahoo.com.ar})
#'
#' @export
#'
#'@examples
#'
#'### Example 2
#'data2Tst <- data.frame(
#'J1 = c(4, 1, 1, 1, 4),
#'J2 = c(4, 1, 2, 2, 3),
#'J3 = c(4, 1, 3, 3, 5),
#'J4 = c(4, 1, 4, 5, 5),
#'J5 = c(4, 1, 5, 5, 5),
#'J6 = c(4, 1, 3, 5, 5))
#'Vaiken(data = data2Tst, min = 1, max = 5, conf.level = .90)

Vaiken <- function(data, min, max, conf.level = 0.95, na.rm = FALSE) {

  # Detection of Missing Values
  if (!na.rm) {
    if (any(is.na(data))) {
      stop("Missing values detected. Use na.omit() first, or set na.rm=TRUE.")
    }
  } else {
    data <- na.omit(data)
  }

  # Validation: data
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data.frame.")
  }
  if (!is.numeric(as.matrix(data))) {
    stop("The 'data' field must contain only numeric values.")
  }

  # Validation: min, max
  if (!is.numeric(min) || !is.numeric(max)) {
    stop("'min' and 'max' must be numeric values.")
  }
  if (min >= max) {
    stop("'min' must be less than 'max'.")
  }

  # Validation: confidence level
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be between 0 and 1 (excluding the extremes).")
  }

  # Verification: NA
  if (any(is.na(data))) {
    stop("'data' contains NA values. Please clean the data before using the function.")
  }

  # Verification: valid range of data
  if (any(data < min | data > max)) {
    stop("All scores in 'data' must be between 'min' and 'max'.")
  }

  # N judges
  n <- ncol(data)

  # Critic value Z for confidence level
  z <- qnorm(1 - (1 - conf.level) / 2)

  # Calculate V Aiken and confidence interval
  calcular_valores <- function(puntajes) {
    M <- mean(puntajes) # Average of the scores
    V <- (M - min) / (max - min) # Calculation of Aiken's V

    # lower and upper limits, method Wilson's score
    nk <- n * (max - min)
    term1 <- 2 * nk * V + z^2
    term2 <- z * sqrt(4 * nk * V * (1 - V) + z^2)
    denom <- 2 * (nk + z^2)

    lwr.ci <- (term1 - term2) / denom
    upr.ci <- (term1 + term2) / denom

    return(c(V = V, lwr.ci = lwr.ci, upr.ci = upr.ci))
  }

  # Aiken V for every item (row in the data frame)
  resultados <- t(apply(data, 1, calcular_valores))

  # Data frame for the output
  resultados_df <- data.frame(
    Item = 1:nrow(data),
    V = round(resultados[, "V"], 3),
    lwr.ci = round(resultados[, "lwr.ci"], 3),
    upr.ci = round(resultados[, "upr.ci"], 3)
  )

  return(resultados_df)
}
