#' @title Substantive Validity for single Items
#' @description This function calculates substantive coefficients ('psa' and 'svc'; Anderson & Gerbing, 1991), including asymmetric confidence intervals using the Wilson method.
#' @param item Numerical vector with the responses to an item. Each response represents a chosen construct, usually forming a multinomial variable.
#' @param conf.level Confidence level for asymmetric confidence intervals (e.g., .90, .95, .99).
#'
#' @return
#' A list of two data frames with the psa and svc results for each category (i.e., construct evaluated), together with asymmetric confidence intervals.
#'
#' @details
#' The function analyzes substantive validity using the approach of Anderson and
#' Gerbing (1991). Two coefficients are estimated: proportion of substantive agreement
#' ('psa') and the substantive validity coefficient ('svc'). Both coefficients are
#' calculated for each construct in order to compare the target construct with the rest.
#' Each 'psa' and 'svc' is supplemented with asymmetric confidence intervals
#' (Penfield & Giacobbi, 2004; Wilson, 1927).
#'
#' The input vector for 'item' can come from a 'data.frame', where rows are
#' the participants (experiential judges or expert judges), and the values represent
#' the chosen construct. Values may be: 1, 2, 3, 4 or other numerical codes (e.g., 2, 4, 10, 6).
#'
#' Note: This function does not yet handle missing values. Users should remove or impute missing data before using it.
#'
#' @references
#' Anderson, J. C., & Gerbing, D. W. (1991). Predicting the performance of measures in
#' a confirmatory factor analysis with a pretest assessment of their substantive
#' validities. *Journal of Applied Psychology*, 76(5), 732–740. https://doi.org/10.1037/0021-9010.76.5.732
#'
#' Cabedo-Peris, J., Merino-Soto, C., Chans, G.M., & Marti-Vilar, M. (2024).
#'  Exploring the Loss Aversion Scale’s psychometric properties in Spain.
#'  *Scientific Reports*, 14, 15756. https://doi.org/10.1038/s41598-024-66695-6
#'
#' Penfield, R. D., & Giacobbi, P. R., Jr. (2004). Applying a score confidence interval
#' to Aiken’s item content-relevance index. *Measurement in Physical Education and
#' Exercise Science*, 8(4), 213–225. https://doi.org/10.1207/s15327841mpee0804_3
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical
#' inference. *Journal of the American Statistical Association*, 22, 209–212.
#' https://doi.org/10.2307/2276774
#'
#' @seealso
#' \code{\link[PropCIs:scoreci]{scoreci}} for score method confidence interval
#' \code{\link[ValCont:SVALmult]{SVALmult}} for multi-item analysis. \cr
#' \code{\link[ValCont:SVALplot]{SVALplot}} for visualizing results of SVAL analyses.
#'
#' @examples
#' # Example 1
#' gais <- data.frame(
#'   gais2 = c(4,4,4,4,4,2,3,2,4,4,3,4,4,3,4,4,2,1,2,3,4,3,4,4,2,3,4,4,4),
#'   gais3 = c(4,4,4,4,4,1,4,2,1,4,3,4,4,3,4,4,4,1,2,3,3,3,2,4,4,3,4,4,4),
#'   gais4 = c(4,4,4,4,4,1,1,1,1,2,4,4,1,1,4,4,2,1,2,3,1,3,4,4,2,4,2,4,4),
#'   gais5 = c(4,4,4,2,4,1,2,1,1,4,2,1,4,4,4,4,2,1,2,3,3,3,4,1,1,4,1,3,4)
#' )
#' SVALsingle(item = gais$gais2, conf.level = 0.95)
#'
#' @author
#' Cesar Merino-Soto (\email{sikayax@yahoo.com.ar})
#'
#' @export
SVALsingle <- function(item, conf.level = 0.95)
{
  # Wilson score interval function
  scoreci <- function(x, n, conf.level) {
    zalpha <- abs(qnorm((1 - conf.level) / 2))
    phat <- x / n
    bound <- (zalpha * sqrt((phat * (1 - phat) + (zalpha^2) / (4 * n)) / n)) /
      (1 + (zalpha^2) / n)
    midpnt <- (phat + (zalpha^2) / (2 * n)) / (1 + (zalpha^2) / n)

    up.ci <- round(midpnt + bound, digits = 3)
    lwr.ci <- round(midpnt - bound, digits = 3)

    return(data.frame(lwr.ci = lwr.ci, up.ci = up.ci))
  }

  # Frequency table
  freq_table <- as.data.frame(table(item))
  colnames(freq_table) <- c("Cat", "Freq")
  freq_table$Cat <- as.numeric(as.character(freq_table$Cat))
  Ntotal <- sum(freq_table$Freq)

  # PSA
  psa_results <- freq_table
  psa_results$psa <- round(freq_table$Freq / Ntotal, 3)
  psa_ci <- scoreci(freq_table$Freq, Ntotal, conf.level)
  psa_results <- cbind(psa_results, psa_ci)

  # SVC
  svc_results <- freq_table
  svc_results$svc <- round((freq_table$Freq - (Ntotal - freq_table$Freq)) / Ntotal, 3)
  svc_abs <- abs(svc_results$svc)
  svc_ci <- scoreci(svc_abs * Ntotal, Ntotal, conf.level)
  svc_ci <- svc_ci * sign(svc_results$svc)
  svc_results <- cbind(svc_results, svc_ci)

  return(list("PSA" = psa_results, "SVC" = svc_results))
}
