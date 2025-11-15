#' @title Critical value for CVR, based on Baghestani et al., (2019)
#' @description Calculates critical values for Lawshe's content validity ratio (CVR) using the Bayesian method
#' proposed by Baghestani et al. (2019), based on Beta priors and posterior probabilities.
#' @param num_jueces Number of judges who provided their ratings (must be >= 2).
#' @param prior Type of prior distribution: "Jeffreys" (default) or "Uniform".
#' @param alpha Significance level for the critical value calculation (default: 0.05).
#' @param tails Distribution tails: "one" for one-tailed test (default) or "two" for two-tailed test.
#'
#' @return
#' A data frame with numeric results and six columns:
#' "Method": The method used (Baghestani).
#' - "Prior": The type of prior used ("Jeffreys" or "Uniform").
#' - "Alpha": The chosen significance level.
#' - "Tails": The type of test ("one" or "two").
#' - "CritVal": The calculated critical value for CVR.
#' - "MinJudges": The minimum number of judges needed for the critical CVR.
#'
#' @details
#' Lawshe's content validity ratio (CVR; ) requires critical values for interpretation, based on the
#' alpha level chosen (0.05 or 0.01) and the number of judges. `CVRcut.Bag` estimates this critical
#' value based on the Bayesian statistic (Baghestani et al., 2017). This approach can maximize the
#' accuracy of the critical values for CVR, due to the usual small sample size to obtain evidence
#' of content validity. This method requires key information about the possible prior distribution.
#' Baghestani et al., (2017) uses noninformative priors for CVR, based on: a) Jeffrey (1935),
#' Beta(1/2, 1/2); or b) Berger (2013)'s uniform, Beta(1,1).
#'
#' @references
#' Baghestani, A. R., Ahmadi, F., Tanha, A., & Meshkat, M. (2019). Bayesian Critical Values for Lawshe's
#' Content Validity Ratio. Measurement and Evaluation in Counseling and Development, 52(1), 69-73.
#' https://doi.org/10.1080/07481756.2017.1308227
#'
#' Lawshe, C. H. (1975). A quantitative approach to content validity. Personnel psychology, 28, 563-575.
#' https://doi.org/10.1111/j.1744-6570.1975.tb01393.x
#'
#' @seealso
#'\code{\link[ValCont:CVR]{ValCont::CVR}}
#'\code{\link[ValCont:CVRcut]{ValCont::CVRcut}}
#' @examples
#'
#' ### Example 1 ----------
#' #' # N judges = 45, uniform prior, alpha = .05, one-tailed test
#' CVRcut.Bag(num_jueces = 45, prior = "Uniform", alpha = .05, tails = "one")
#'
#' ### Example 2 ----------
#' # N judges = 45, Jeffreys prior, alpha = .05, one-tailed test
#' CVRcut.Bag(num_jueces = 45, prior = "Jeffreys", alpha = .05, tails = "one")
#'
#' ### Example 3 ----------
#' # N judges = 45, Jeffreys prior, alpha = .01, two-tailed test
#' CVRcut.Bag(num_jueces = 45, prior = "Jeffreys", alpha = .01, tails = "two")
#'
#' @author
#'Cesar Merino-Soto (\email{sikayax@yahoo.cam.ar})
#'
#' @export
CVRcut.Bag <- function(num_jueces, prior = "Jeffreys", alpha = 0.05, tails = "one") {
  # Validar numero de jueces
  if (num_jueces < 2) {
    stop("El numero de jueces debe ser al menos 2 para un calculo valido.")
  }

  # Ajustar alpha si es una prueba de dos colas
  if (tails == "two") {
    alpha <- alpha / 2
  } else if (tails != "one") {
    stop("El argumento tails debe ser one o two")
  }

  # Parametros de los priors
  if (prior == "Jeffreys") {
    a <- 0.5
    b <- 0.5
  } else if (prior == "Uniform") {
    a <- 1
    b <- 1
  } else {
    stop("El prior debe ser 'Jeffreys' o 'Uniform'.")
  }

  # Calcular n_e critico
  calcular_ne <- function(N, alpha, a, b) {
    for (ne in 0:N) {
      # Probabilidad acumulada con pbeta
      p_null <- pbeta(0.5, ne + a, N - ne + b)
      if (p_null < alpha) {
        return(ne)
      }
    }
    return(NA)  # Si no encuentra un valor
  }

  min_jueces <- calcular_ne(num_jueces, alpha, a, b)

  # Calcular CVR critico
  if (!is.na(min_jueces)) {
    cutoff <- (2 * min_jueces - num_jueces) / num_jueces
  } else {
    cutoff <- NA
  }

  # Resultado
  result <- data.frame(
    Method = "Baghestani",
    Prior = prior,
    Alpha = alpha,
    Tails = tails,
    CritVal = cutoff,
    MinJudges = min_jueces
  )

  return(result)
}
