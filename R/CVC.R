#'@title CVC
#'@description Calculates the content validity coefficient (CVC; Hernandez-Nieto, 2002). CVC makes an adjustment for random response.
#'@param data dataframe, with n columns (judges or scorers), and k rows (evaluated items).
#'@param max maximum possible rating value used.
#'@param conf.level confidence level for confidence intervals (eg., .90, .95, .99).
#'@param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#' 
#'@return dataframe with CVC coefficients and confidence intervals.
#'
#'@details
#'This function calculates the content validity coefficient CVC (Hernandez-Nieto, 2002). Asymmetric confidence intervals are also calculated (Wilson, 1927; Penfield & Giacobbi, 2004). CVC' is the second coefficient that adjusts for possible random response of the raters, while another proposal for the CVI coefficient was created by Polit, & Beck (2007).
#'
#'Note: The function has not yet been prepared to resolve missing values, so the user must remove or impute any missing values.
#'
#'@references
#'Hernandez-Nieto, R. A. (2002). Contributions to Statistical Analysis. Merida, Venezuela: Universidad de Los Andes.
#'Penfield, R. D. & Giacobbi, P. R., Jr. (2004) Applying a score confidence interval to Aiken's item content-relevance index. Measurement in Physical Education and Exercise Science, 8(4), 213-225. https://doi.org/10.1207/s15327841mpee0804_3
#'Polit DF, Beck CT, Owen SV.  Is the CVI an acceptable indicator of content validity?  Appraisal and recommendations. Research in Nursing & Health. 2007;30(4):459-67. https://doi.org/10.1002/nur.20199
#'Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association, 22, 209-212. https://doi.org/10.2307/2276774
#'
#'@seealso
#'\code{\link[PropCIs:scoreci]{PropCIs::scoreci}} for score method confidence interval
#'
#' @author
#' Cesar Merino-Soto (\email{sikayax@yahoo.com.ar})
#'
#' @export
#'
#'@examples
#'### Example 1
#'### Fictitious data Ej1: 15 judges, and 15 items evaluated by the judges
#'
#'Ej1 <- data.frame(
#'  J1 = c(5, 5, 6, 6, 6, 6, 6, 6, 5, 6, 5, 6, 3, 4, 4),
#'  J2 = c(5, 2, 6, 6, 6, 6, 6, 5, 4, 6, 5, 6, 4, 4, 3),
#'  J3 = c(5, 5, 6, 6, 6, 5, 6, 5, 4, 5, 5, 6, 5, 3, 5),
#'  J4 = c(5, 5, 6, 6, 6, 6, 6, 6, 5, 6, 3, 6, 5, 3, 5),
#'  J5 = c(5, 5, 6, 6, 6, 6, 6, 6, 6, 5, 5, 6, 3, 4, 5),
#'  J6 = c(5, 2, 6, 6, 6, 6, 6, 5, 6, 6, 3, 6, 4, 4, 4),
#'  J7 = c(2, 4, 5, 6, 6, 6, 6, 5, 6, 6, 5, 6, 5, 3, 6),
#'  J8 = c(5, 5, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 5, 3, 5),
#'  J9 = c(4, 5, 5, 6, 5, 4, 6, 5, 5, 6, 4, 6, 5, 4, 4),
#'  J10 = c(4, 2, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 3),
#'  J11 = c(5, 4, 6, 6, 6, 6, 5, 6, 6, 6, 6, 6, 5, 4, 4),
#'  J12 = c(5, 4, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 5, 4, 4),
#'  J13 = c(2, 4, 6, 6, 6, 6, 5, 6, 4, 2, 4, 6, 3, 3, 4),
#'  J14 = c(5, 5, 6, 6, 6, 6, 6, 5, 6, 6, 6, 6, 5, 4, 5),
#'  J15 = c(5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 3, 4))
#'
#'## Run CVC.
#'CVC(Ej1,max = 6, conf.level = .90)
#'
#'### Example 2
#'### Random data
#'## Aritificial data
#'random_data <- data.frame(replicate(n = 6,sample(1:5,12,rep=TRUE)))
#'
#'## Run CVC
#'CVC(random_data,max = 5, conf.level = .90)
#'


# Funcion principal CVC
CVC <- function(data, max, conf.level, na.rm = FALSE) {
  # Verificar si el data.frame contiene solo valores numericos
  if (!all(sapply(data, is.numeric))) {
    stop("El data.frame debe contener solo valores numericos.")
  }

  # Detección de valores perdidos
  if (!na.rm) {
    if (any(is.na(data))) {
      stop("Valores perdidos detectados. Usa na.omit() primero o establece na.rm=TRUE.")
    }
  } else {
    data <- na.omit(data)
  }

  # Numero de jueces (columnas)
  num_jueces <- ncol(data)

  # Calcular Pe segun la formula dada
  Pe <- (1 / num_jueces) ^ num_jueces

  # Calcular la media de cada item (fila)
  medias_items <- rowMeans(data)

  # Calcular el CVC para cada item
  CVC_items <- round((medias_items - Pe) / max, 3)

  # Calcular los intervalos de confianza de Wilson para cada CVC

  get_wilson_CI <- function(x, n, conf.level) {
    p_hat <- x
    SE_hat_sq <- p_hat * (1 - p_hat) / n
    crit <- qnorm(1 - (1 - conf.level) / 2)
    omega <- n / (n + crit^2)
    A <- p_hat + crit^2 / (2 * n)
    B <- crit * sqrt(SE_hat_sq + crit^2 / (4 * n^2))
    CI <- c('lower' = omega * (A - B),
            'upper' = omega * (A + B))
    return(CI)
  }
  intervalos_CI <- t(sapply(CVC_items, get_wilson_CI, n = num_jueces, conf.level = conf.level))

  # Crear un data.frame con los resultados
  resultado_df <- data.frame(
    Item = 1:nrow(data),
    CVC = round(CVC_items, 3),
    lwr.ci = round(intervalos_CI[, "lower"],3),
    upr.ci = round(intervalos_CI[, "upper"],3)
  )

  return(resultado_df)
}

