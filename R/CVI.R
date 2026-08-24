#'@title Content Validity Index (CVI)
#'@description Calculate the Content Validity Index (CVI), without adjustment for random agreement.
#'@param data dataframe, with the columns assigned to each judge, and the rows assigned to each assessed item.
#'@param cut Specific cut-off point, from which the item is considered valid (or relevant, clear, etc.)
#'@param conf.level confidence level for asymmetric confidence intervals (ex., .90, .95, .99)
#'@param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#'@return
#'dataframe with CVI for all items analyzed, and confidence intervals.
#'
#'@details
#'Calculate the CVI coefficient (Martuza, 1977; Lynn, 1986), and asymmetric confidence intervals (Wilson, 1927; Penfield & Giacobbi, 2004). The results should be complemented by an estimator of inter-judge variability or agreement.
#'The content validity index (CVI) was first proposed by Martuza (1977), but frequently associated with other researchers (Polit, & Beck, 2006), e.g., Lynn (1986), Davis (1992), among others.
#'The \strong{CVI} function calculates the CVI without adjustment for chance agreement (Polit, & Beck, 2007).
#'Items can be rated according to traditional recommendations (Polit, Beck, & Owen, 2007; Polit, & Beck, 2006), based on an ordinal scale between 3 to 5 points (Lynn, 1986), or on a 4-point scale (Lynn, 1986; Waltz & Bausell, 1981).
#'The usual labels to evaluate the relevance for each item are: not relevant, somewhat relevant, etc. (Davis, 1992).
#'The usual CVI cut-off point for identifying valid from invalid items is generally in the top two ratings (3 or higher on a relevance scale of 1 to 4; Beck & Gable, 2001; Grant & Davis, 1997). 'cut' dichotomizes the judges' responses to calculate CVI.
#'The \strong{CVI} function can be used for rated items with any rating range, and any chosen cut point ('cut').
#'Asymmetric confidence intervals use Wilson's (1927) approach, as used for Aiken's V coefficient (Penfield, & Giacobbi, 2004).
#'
#'@references
#'Davis, L.L. (1992). Instrument review: Getting the most from your panel of experts. \emph{Applied Nursing Research, 5}, 194-197. https://doi.org/10.1016/S0897-1897(05)80008-4
#'
#'Grant, J.S., & Davis, L.T. (1997). Selection and use of content experts in instrument development. \emph{Research in Nursing & Health, 20}, 269-274. https://doi.org/10.1002/(sici)1098-240x(199706)20:3<269::aid-nur9>3.0.co;2-g
#'
#'Lynn, M.R. (1986). Determination and quantification of content validity. \emph{Nursing Research, 35}, 382-385.
#'
#'Martuza, V.R. (1977). \emph{Applying norm-referenced and criterion-referenced measurement in education}. Boston: Allyn & Bacon
#'
#'Penfield, R. D. & Giacobbi, P. R., Jr. (2004) Applying a score confidence interval to Aiken's item content-relevance index. \emph{Measurement in Physical Education and Exercise Science, 8}(4), 213-225. \doi{10.1207/s15327841mpee0804_3}
#'
#'Polit, D. F., & Beck, C. T. (2006). The content validity index: are you sure you know what's being reported? Critique and recommendations. \emph{Research in Nursing & Health, 29}(5), 489-497. \doi{10.1002/nur.20147}
#'
#'Polit, D.F., Beck, C.T. and Owen, S.V. (2007), Is the CVI an acceptable indicator of content validity? Appraisal and recommendations. \emph{Research in Nursing & Health, 30}, 459-467. \doi{10.1002/nur.20199}
#'
#'Waltz, C.F., & Bausell, R.B. (1981). \emph{Nursing research: Design, statistics, and computer analysis}. Philadelphia: F. A. Davis.
#'
#'Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. \emph{Journal of the American Statistical Association, 22}, 209-212. \doi{10.2307/2276774}
#'
#'@seealso
#'\code{\link[PropCIs:scoreci]{PropCIs::scoreci}} for score method confidence interval
#'
#'@examples
#'
#'# Load data
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
#'# Run CVI
#'CVI(data = Ej1, cut = 4, conf.level = .90)
#'
#'@export
CVI <- function(data, cut, conf.level, na.rm = FALSE) {
  # Built-in function for calculating the Wilson confidence interval
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

  # Check whether the data is numeric
  if (!all(sapply(data, is.numeric))) {
    stop("All columns must contain numeric data.")
  }

  # Detection of Missing Values
  if (!na.rm) {
    if (any(is.na(data))) {
      stop("Missing values detected. Use na.omit() first, or set na.rm=TRUE.")
    }
  } else {
    data <- na.omit(data)
  }

  # Total number of judges
  num_jueces <- ncol(data)

  # Initialize lists to store results
  CVI_vals <- numeric(nrow(data))
  lwr_cis <- numeric(nrow(data))
  upp_cis <- numeric(nrow(data))

  # Calculate CVI, agreement percentage, and Wilson's CI for each item
  for (i in 1:nrow(data)) {
    Na <- sum(data[i, ] >= cut)
    Nna <- num_jueces - Na
    N <- Na + Nna

    CVI_vals[i] <- Na / N

    # Calculate the Wilson confidence interval
    CI <- get_wilson_CI(CVI_vals[i], N, conf.level)
    lwr_cis[i] <- CI['lower']
    upp_cis[i] <- CI['upper']
  }

  # Data frame with the results, rounded to 3 decimal places
  resultado <- data.frame(
    Item = 1:nrow(data),
    CVI = round(CVI_vals, 3),
    lwr.ci = round(lwr_cis, 3),
    upr.ci = round(upp_cis, 3)
  )

  return(resultado)
}

