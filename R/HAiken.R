#'@title Coefficient of homogeneity of response
#'@description Calculate the coefficient of homogeneity of response for each item (Aiken, 1980, 1985).
#'
#'@param data dataframe, with the columns assigned to each judge, and the rows assigned to each evaluated item.
#'@param ncat number of response categories or options used in the rating
#'@param conf.level confidence level for the confidence intervals (eg., .90, .95, .99)
#'@param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#' 
#'@return
#'dataframe with H coefficients for all items analyzed, and their confidence intervals.
#'
#'@details
#'Compute the H coefficient (Aiken, 1980, 1985) to estimate the homogeneity of response of the judges/scorers to the items.
#'To maintain consistency with the methods usually associated with content validity, 'HAiken' is proposed as an option.
#''HAiken' also compute asymmetric confidence intervals use the method of Wilson (1927). and adapted by Penfield and Giacobbi (2004) for Aiken's V coefficient.
#'The H coefficient, or equivalent coefficients, should complement the results of the content validity coefficients.
#'Other methods for estimating judges' agreement or homogeneity of response may also be useful.
#'Note: The function has not yet been prepared to resolve missing values, so the user must remove or impute any missing values.
#'
#'@references
#'Aiken, L. R. (1980). Content validity and reliability of single items or questionnaires. Educational and. Psychological Measurement, 40, 955-959. https://doi.org/10.1177/001316448004000419
#'
#'Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. Educational and Psychological Measurement, 45, 131-142. https://doi.org/10.1177/0013164485451012
#'
#'Penfield, R. D., & Miller, J. M. (2004). Improving Content Validation Studies Using an Asymmetric Confidence Interval for the Mean of Expert Ratings. Applied Measurement in Education, 17(4), 359–370. https://doi.org/10.1207/s15324818ame1704_2
#'
#'Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association, 22, 209-212. https://doi.org/10.2307/2276774
#'
#'@seealso
#'\code{\link[PropCIs:scoreci]{PropCIs::scoreci}} for score method confidence interval
#'
#'@author
#'Cesar Merino-Soto (\email{sikayax@yahoo.cam.ar})
#'

#'@examples
#'### Example 1 --------------
#'
#'#Load data
#'Ej2 <- data.frame(
#'  j1 = c(4, 1, 1, 1, 4),
#'  j2 = c(4, 1, 2, 2, 3),
#'  j3 = c(4, 1, 3, 3, 5),
#'  j4 = c(4, 1, 4, 5, 5),
#'  j5 = c(4, 1, 5, 5, 5),
#'  j6 = c(4, 1, 3, 5, 5)
#')
#'
#'# Run HAiken
#'Haiken(Ej2, ncat = 5, conf.level = .90)
#'
#'### Example 2 ----------------
#'# In a dataframe where the rows are the items and the columns are the raters,
#'# H can be calculated for the raters (columns) by simply transposing the data
#'# and entering it as a data frame.
#'
#'Haiken(as.data.frame(t(Ej2)), ncat = 5, conf.level = .90)
#'
#'@export
Haiken <- function(data, ncat, conf.level, na.rm = FALSE) {

  # Detección de valores perdidos
  if (!na.rm) {
    if (any(is.na(data))) {
      stop("Valores perdidos detectados. Usa na.omit() primero o establece na.rm=TRUE.")
    }
  } else {
    data <- na.omit(data)
  }

  # Numero de filas
  num_filas <- nrow(data)

  # Calcular j segun si ncol(data) es par o impar
  j <- ifelse(ncol(data) %% 2 == 0, 0, 1)

  # data frame para almacenar resultados
  resultados <- data.frame(
    Item = 1:num_filas,
    H = numeric(num_filas),
    lwr.ci = numeric(num_filas),
    upr.ci = numeric(num_filas),
    n.subj = numeric(num_filas)
  )

  # Calcular las diferencias por pares (diferencias absolutas) y los intervalos de confianza para cada fila
  for (i in 1:num_filas) {
    fila_actual <- data[i, ]

    # Calcular las diferencias por pares como diferencias absolutas y desempaquetarlas
    diff_abs <- combn(colnames(fila_actual), 2, function(pair) {
      col1 <- fila_actual[pair[1]]
      col2 <- fila_actual[pair[2]]
      diff_abs <- abs(col1 - col2)
      return(diff_abs)
    }, simplify = FALSE)  # Utilizar simplify = FALSE para obtener una lista

    # Desempaquetar las diferencias absolutas de la lista
    diff_abs_unlist <- unlist(diff_abs)

    # Calcular el resultado para la fila actual usando ncat, diff_abs_unlist y j
    resultado_actual <- 1 - (4 * sum(diff_abs_unlist) / ((ncat - 1) * (ncol(data)^2) - j))

    get_wilson_CI <- function(x, n, conf.level) {
      n <- n
      p_hat <- x
      SE_hat_sq <- p_hat * (1 - p_hat) / n
      crit <- qnorm(1 - conf.level / 2)
      omega <- n / (n + crit^2)
      A <- p_hat + crit^2 / (2 * n)
      B <- crit * sqrt(SE_hat_sq + crit^2 / (4 * n^2))
      CI <- c('lower' = omega * (A - B),
              'upper' = omega * (A + B))
      return(CI)
    }
    # Calcular el intervalo de confianza con la funcion get_wilson_CI y el nivel de confianza alpha
    wilson_interval <- get_wilson_CI(resultado_actual, ncol(data), conf.level)

    # Almacenar los resultados en el data frame
    resultados[i, "H"] <- resultado_actual
    resultados[i, "lwr.ci"] <- wilson_interval['lower']
    resultados[i, "upr.ci"] <- wilson_interval['upper']
    resultados[i, "n.subj"] <- num_filas

  }

  return(round(resultados, 3))
}

