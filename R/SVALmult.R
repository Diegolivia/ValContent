#' @title Substantive Validity for an set of Items
#' @description
#' For a data frame of items, this function calculates two substantive validity
#' coefficients—\code{psa} (proportion of substantive agreement) and \code{svc}
#' (substantive validity coefficient; Anderson & Gerbing, 1991)—along with asymmetric
#' confidence intervals using the Wilson method.
#'
#' @param data A data frame with item responses. Each response represents a chosen construct, typically forming a multinomial variable.
#' @param columns A vector of column indices or names to be analyzed.
#' @param conf.level Confidence level for asymmetric confidence intervals (e.g., .90, .95, .99).
#'
#' @return
#' A named list where each element corresponds to one of the analyzed items. Each element contains two data frames:
#' \itemize{
#'   \item \code{psa_results}: A data frame with the following columns for each category (construct):
#'     \describe{
#'       \item{\code{Cat}}{The category value (construct code).}
#'       \item{\code{Freq}}{The frequency of responses in that category.}
#'       \item{\code{psa}}{The proportion of substantive agreement.}
#'       \item{\code{lwr.ci}}{Lower bound of the Wilson confidence interval for \code{psa}.}
#'       \item{\code{up.ci}}{Upper bound of the Wilson confidence interval for \code{psa}.}
#'     }
#'
#'   \item \code{svc_results}: A data frame with the following columns for each category:
#'     \describe{
#'       \item{\code{Cat}}{The category value (construct code).}
#'       \item{\code{Freq}}{The frequency of responses in that category.}
#'       \item{\code{svc}}{The substantive validity coefficient.}
#'       \item{\code{lwr.ci}}{Lower bound of the Wilson confidence interval for \code{svc}.}
#'       \item{\code{up.ci}}{Upper bound of the Wilson confidence interval for \code{svc}.}
#'     }
#' }
#'
#' @details
#' The function analyzes substantive validity using the approach of Anderson and Gerbing (1991). It estimates:
#' \itemize{
#'   \item \strong{psa}: Proportion of responses in favor of the target category.
#'   \item \strong{svc}: Substantive validity coefficient, calculated by comparing the frequency of the target category with the average of the non-target categories.
#' }
#'
#' Both coefficients are calculated for each construct to compare the target construct with the others.
#' Asymmetric confidence intervals are computed using the Wilson score method (Penfield & Giacobbi, 2004; Wilson, 1927).
#' Examples of use can be founded in Cabedo et al. (2024), and Merino-Soto et al. (2021).
#'
#' Each column in \code{data} should contain numerical responses representing constructs
#'  (e.g., 1, 2, 3, 4). These values do not need to be sequential but must be distinct.
#'
#' \strong{Note}: The function does not handle missing values. Users must impute or remove NAs before using this function.
#'
#' \strong{Interpretation of Negative SVC Values}: Negative values of the \code{svc} coefficient are valid and interpretable. They indicate that the selected category was chosen less frequently than the average of the non-target categories, suggesting low substantive agreement for that construct.
#'
#' @references
#' Anderson, J. C., & Gerbing, D. W. (1991). Predicting the performance of measures in a confirmatory factor analysis with a pretest assessment of their substantive validities. \emph{Journal of Applied Psychology, 76}(5), 732–740. https://doi.org/10.1037/0021-9010.76.5.732
#'
#' Cabedo-Peris, J., Merino-Soto, C., Chans, G.M., & Marti-Vilar, M. (2024). Exploring the Loss Aversion Scale’s psychometric properties in Spain. \emph{Scientific Reports, 14}, 15756. https://doi.org/10.1038/s41598-024-66695-6
#'
#' Merino-Soto, C., Calderón-De la Cruz, G., Gil-Monte, P., & Juárez-García, A. (2021). Substantive validity within the framework of content validity: Application in the Workload Scale. \emph{Revista Argentina de Ciencias del Comportamiento, 13}(1), 81–92. \url{https://revistas.unc.edu.ar/index.php/racc/article/view/20547/33426}
#'
#' Penfield, R. D., & Giacobbi, P. R. Jr. (2004). Applying a score confidence interval to Aiken’s item content-relevance index. \emph{Measurement in Physical Education and Exercise Science, 8}(4), 213–225. https://doi.org/10.1207/s15327841mpee0804_3
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. \emph{Journal of the American Statistical Association, 22}, 209–212. https://doi.org/10.2307/2276774
#'
#' @seealso
#' \code{\link[PropCIs]{scoreci}} for score-based confidence intervals. \cr
#' \code{\link[ValContent:SVALsingle]{SVALsingle}} for computing substantive validity on individual items. \cr
#' \code{\link[ValContent:SVALplot]{SVALplot}} for visualizing the results of SVALmult.
#'
#' @examples
#' data.gais <- data.frame(
#'  gais2 = c(4,4,4,4,4,2,3,2,4,4,3,4,4,3,4,4,2,1,2,3,4,3,4,4,2,3,4,4,4),
#'  gais3 = c(4,4,4,4,4,1,4,2,1,4,3,4,4,3,4,4,4,1,2,3,3,3,2,4,4,3,4,4,4),
#'  gais4 = c(4,4,4,4,4,1,1,1,1,2,4,4,1,1,4,4,2,1,2,3,1,3,4,4,2,4,2,4,4),
#'  gais5 = c(4,4,4,2,4,1,2,1,1,4,2,1,4,4,4,4,2,1,2,3,3,3,4,1,1,4,1,3,4))
#'
#' SVALmult(data = data.gais, columns = c("gais2", "gais3"), conf.level = .90)
#' SVALmult(data = data.gais, columns = c(1,2,3), conf.level = .90)
#'
#' @author
#' César Merino-Soto (\email{sikayax@yahoo.com.ar})
#'
#' @export
SVALmult <- function(data, columns, conf.level = 0.95) {
  # Subfunction: Wilson CI
  scoreci <- function(x, n, conf.level) {
    zalpha <- abs(qnorm((1 - conf.level) / 2))
    phat <- x / n
    bound <- (zalpha * sqrt((phat * (1 - phat) + (zalpha^2) / (4 * n)) / n)) /
      (1 + (zalpha^2) / n)
    midpnt <- (phat + (zalpha^2) / (2 * n)) / (1 + (zalpha^2) / n)
    data.frame(lwr.ci = round(midpnt - bound, 3), up.ci = round(midpnt + bound, 3))
  }

  # Subfunction: core calculation per item
  vsub_sing_internal <- function(item, conf.level) {
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

    list(psa_results = psa_results, svc_results = svc_results)
  }

  if (is.numeric(columns)) {
    columns <- colnames(data)[columns]
  }

  results <- lapply(columns, function(col) {
    item <- data[[col]]
    vsub_sing_internal(item, conf.level)
  })
  names(results) <- columns
  return(results)
}
