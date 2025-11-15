#'@title Content Validity Index Revised (with random agreement adjustment)
#'@description Calculate the CVI coefficient with random agreement adjustment (Polit & Beck, 2007).
#'@param data dataframe, with the columns assigned to each judge, and the rows assigned to each evaluated item.
#'@param cut Specific cut-off point in the response rating, from which the item is considered valid (or relevant, clear, etc.).
#'@param conf.level confidence level for asymmetric confidence intervals (ex., .90, .95, .99)
#'
#'@return
#'dataframe with `CVIR` coefficients for all items analyzed, and confidence intervals.
#'
#'@details
#'The `CVIR` function calculates the CVI, with adjustment for chance agreement (Polit, Beck, & Owen, 2007).
#'As with the 'CVI' function, in 'CVIR' items can be rated with a wide range of point scales (Polit, Beck, & Owen, 2007; Polit, & Beck, 2006). This rating scale can be 3-point, or usually 4- or 5-point (Lynn, 1986; Waltz & Bausell, 1981).
#'The `CVIR` cutoff point, for identifying valid items, is generally based on the last two top ratings; for example, on a relevance scale of 1 to 4, the cutoff point may be 3 (Beck & Gable, 2001; Grant & Davis, 1997). 'CVIR' uses 'cut' to dichotomize the judges' responses and transform them into the CVI.
#'The asymmetric confidence intervals use Wilson's (1927) approach, as applied to Aiken's V coefficient (Penfield, & Giacobbi, 2004).
#'The results should be complemented by an estimator of inter-judge variability or agreement.
#'
#'Note: The function has not yet been prepared to resolve missing values, so the user must remove or impute any missing values.
#'
#'@references
#'Davis, L.L. (1992). Instrument review: Getting the most from your panel of experts. Applied Nursing Research, 5, 194-197. https://doi.org/10.1016/S0897-1897(05)80008-4
#'
#'Grant, J.S., & Davis, L.T. (1997). Selection and use of content experts in instrument development. Research in Nursing & Health, 20, 269-274. https://doi.org/10.1002/(sici)1098-240x(199706)20:3<269::aid-nur9>3.0.co;2-g
#'
#'Lynn, M.R. (1986). Determination and quantification of content validity. Nursing Research, 35, 382-385. https://doi.org/10.1097/00006199-198611000-00017
#'
#'Martuza, V.R. (1977). Applying norm-referenced and criterion-referenced measurement in education. Boston: Allyn & Bacon
#'
#'Penfield, R. D. & Giacobbi, P. R., Jr. (2004) Applying a score confidence interval to Aiken's item content-relevance index. Measurement in Physical Education and Exercise Science, 8(4), 213-225. https://doi.org/10.1207/s15327841mpee0804_3
#'
#'Polit, D. F., & Beck, C. T. (2006). The content validity index: are you sure you know what's being reported? Critique and recommendations. Research in nursing & health, 29(5), 489-497. https://doi.org/10.1002/nur.20147
#'
#'Polit, D.F., Beck, C.T. and Owen, S.V. (2007), Is the CVI an acceptable indicator of content validity? Appraisal and recommendations. Res. Nurs. Health, 30: 459-467. https://doi.org/10.1002/nur.20199
#'
#'Waltz, C.F., & Bausell, R.B. (1981). Nursing research: Design, statistics, and computer analysis. Philadelphia: F. A. Davis.
#'
#'Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association, 22, 209-212. doi: 10.2307/2276774
#'
#'@seealso
#'\code{\link[PropCIs:scoreci]{PropCIs::scoreci}} for score method confidence interval
#'\code{\link[ValCont:CVI]{ValCont:CVI}} for CVI
#'
#'@examples
#'## Load data
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
#'## Run CVIR
#'
#'CVIR(data = Ej1, cut = 4, conf.level = .90)
#'
#'@author
#'Cesar Merino-Soto (\email{sikayax@yahoo.cam.ar})
#'
#'@export
CVIR <- function(data, cut, conf.level) {
  # FunciOn interna para calcular el intervalo de confianza de Wilson
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
  # FunciOn interna para calcular P_C usando la fOrmula combinatoria
  calculate_pc <- function(N, A) {
    if (A > N || A < 0) {
      return(0) # Evita errores si A es mayor que N o menor que 0
    }
    # Aplicamos la fOrmula combinatoria
    Pc <- (factorial(N) / (factorial(A) * factorial(N - A))) * (0.5^N)
    return(Pc)
  }
  # Verificar si los datos son numEricos
  if (!all(sapply(data, is.numeric))) {
    stop("Todas las columnas deben contener datos numEricos.")
  }
  # NUmero total de jueces
  num_jueces <- ncol(data)
  # Inicializar listas para almacenar resultados
  cvipb_vals <- numeric(nrow(data))
  lwr_cis <- numeric(nrow(data))
  upr_cis <- numeric(nrow(data))
  # Calcular CVI ajustado para cada Item
  for (i in 1:nrow(data)) {
    Na <- sum(data[i, ] >= cut)  # Jueces en acuerdo
    Nna <- num_jueces - Na       # Jueces en desacuerdo
    N <- Na + Nna                # Total de jueces
    CVI_val <- Na / N
    # Calcular P_C usando la nueva fOrmula combinatoria dentro de la funciOn
    P_C <- calculate_pc(N, Na)
    # Calcular kappa de Lynn
    cvipb_vals[i] <- (CVI_val - P_C) / (1 - P_C)
    # Calcular intervalo de confianza de Wilson para kappa usando el valor absoluto
    abs_kappa <- abs(cvipb_vals[i])
    CI_kappa <- get_wilson_CI(abs_kappa, N, conf.level)
    # Si kappa es negativo, devolver el signo negativo a los lImites del CI
    if (!is.na(cvipb_vals[i]) && cvipb_vals[i] < 0) {
      lwr_cis[i] <- -CI_kappa['upper']  # Intercambia el lImite superior e inferior
      upr_cis[i] <- -CI_kappa['lower']
    } else {
      lwr_cis[i] <- CI_kappa['lower']
      upr_cis[i] <- CI_kappa['upper']
    }
  }
  # Crear un data frame con los resultados y redondear a 3 decimales
  resultado <- data.frame(
    Item = 1:nrow(data),
    CVIR = round(cvipb_vals, 3),
    lwr.ci = round(lwr_cis, 3),
    upr.ci = round(upr_cis, 3)
  )
  return(resultado)
}
