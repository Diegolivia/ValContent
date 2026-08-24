#'@title  Mean of expert ratings (MER)
#'@description Calculate the Mean of expert ratings (MER), with asymmetric confidence intervals.
#'@param data dataframe, with the columns assigned to each judge, and the rows assigned to each evaluated item.
#'@param ncat Number of possible values or categories used in rating
#'@param start Minimum possible value (0 or 1)
#'@param conf.level Confidence level for confidence intervals (ex., .90, .95, .99)
#'@param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#'
#'@return
#'dataframe with MERs for all items analyzed, and confidence intervals.
#'
#'@details
#'Calculate the average rating of the judges for each item, based on the proposal of Penfield & Miller (2004), and asymmetric confidence intervals (Wilson, 1927; Penfield, 2003; Penfield & Miller, 2004). MER' is a modification of the two syntax previous (Merino-Soto and Livia-Segovia, 2022; Penfield, & Miller, 2004).
#'Due to the usual small number of participants in content validity studies, asymmetric confidence intervals may be an optimal approach to inferentially assess item validity.
#'The results should be supplemented with an estimator of inter-judge variability or agreement.
#'
#'@references
#'Aiken, L. R. (1980). Content validity and reliability of single items or questionnaires. \emph{Educational and. Psychological Measurement, 40}, 955-959. \doi{10.1177/001316448004000419}
#'
#'Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. \emph{Educational and Psychological Measurement, 45}, 131-142. \doi{10.1177/0013164485451012}
#'
#'Merino, C., & Livia, J. (2009). Intervalos de confianza asimetricos para el indice de validez de contenido: un programa Visual Basic para la V de Aiken. \emph{Anales de Psicologia, 25}(1), 169-171. \url{https://revistas.um.es/analesps/article/view/71631}
#'
#'Miller, J. M., & Penfield, R. D. (2005). Using the score method to construct asymmetric confidence intervals: An SAS program for content validation in scale development. \emph{Behavior Research Methods, 37}, 450-452. \doi{10.3758/BF03192713}
#'
#'Penfield, R. D. (2003). A score method of constructing asymmetric confidence intervals for the mean of a rating scale item. \emph{Psychological methods, 8}(2), 149-163. \doi{10.1037/1082-989x.8.2.149}
#'
#'Penfield, R. D. & Giacobbi, P. R., Jr. (2004) Applying a score confidence interval to Aiken’s item content-relevance index. \emph{Measurement in Physical Education and Exercise Science, 8}(4), 213-225. \doi{10.1207/s15327841mpee0804_3}
#'
#'Penfield, R. D., & Miller, J. M. (2004). Improving Content Validation Studies Using an Asymmetric Confidence Interval for the Mean of Expert Ratings. \emph{Applied Measurement in Education, 17}(4), 359–370. \doi{10.1207/s15324818ame1704_2}
#'
#'Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. \emph{Journal of the American Statistical Association, 22}, 209-212. \doi{10.2307/2276774}
#'
#'@seealso
#'\code{\link[ValCont:Haiken]{ValCont::HAiken}} for a homogeneity coefficient
#'
#'@examples
#'
#'## Example
#'
#'#Load data
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
#'## Run MER
#'
#'MER(data = Ej1, ncat = 6, start = 1, conf.level = .90)
#'
#'@export
MER <- function(data, ncat, start, conf.level = 0.90, na.rm = FALSE) {

  # Missing values detection
  if (!na.rm) {
    if (any(is.na(data))) {
      stop("Missing values detected. Use na.omit() first or set na.rm=TRUE.")
    }
  } else {
    data <- na.omit(data)
  }

  # Initial checks
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!is.numeric(ncat) || ncat <= 1) stop("'ncat' must be a number greater than 1")
  if (!is.numeric(start) || !start %in% c(0, 1)) stop("'start' must be 0 or 1")
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be between 0 and 1")
  }

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)          # normal quantile (not t)
  k <- ncat - 1                      # range of the scale (0 to k)
  min_val <- start
  max_val <- start + ncat - 1

  # Function to compute CI for one item (vector of ratings)
  compute_item <- function(item_ratings) {
    M <- mean(item_ratings)
    n <- length(item_ratings)

    # Compute p depending on start
    if (start == 0) {
      p <- M / k
    } else { # start == 1
      p <- (M - 1) / k
    }

    # Avoid p exactly 0 or 1 for numerical stability (though formulas handle it)
    p <- max(1e-10, min(1 - 1e-10, p))

    # Wilson score limits for proportion
    term1 <- 2 * p * n * k + z^2
    term2 <- z * sqrt(4 * n * k * p * (1 - p) + z^2)
    denom <- 2 * (n * k + z^2)

    pi_L <- (term1 - term2) / denom
    pi_U <- (term1 + term2) / denom

    # Confidence limits for the mean
    LCL <- M - z * sqrt(k * pi_L * (1 - pi_L) / n)
    UCL <- M + z * sqrt(k * pi_U * (1 - pi_U) / n)

    # Truncate to the possible scale range (robustness)
    LCL <- max(LCL, min_val)
    UCL <- min(UCL, max_val)

    c(MER = round(M, 3),
      lwr.ci = round(LCL, 3),
      upr.ci = round(UCL, 3))
  }

  results <- apply(data, 1, compute_item)
  results_df <- as.data.frame(t(results))
  results_df$Item <- seq_len(nrow(results_df))
  results_df <- results_df[, c("Item", "MER", "lwr.ci", "upr.ci")]

  return(results_df)
}
