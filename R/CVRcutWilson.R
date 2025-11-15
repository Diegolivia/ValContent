#' @title Critical value for CVR, based on Wilson et al. (2012)
#' @description Calculates critical values for Lawshe's content validity ratio (CVR) based
#' on the method proposed by Wilson et al. (2012).
#' @param num_jueces Number of judges who provided their ratings (must be >= 2).
#' @param alpha Significance level for the critical value calculation (default: 0.05).
#' @param tails Distribution tails: "one" for one-tailed test (default) or "two" for
#' two-tailed test.
#'
#' @return
#' A data frame with numeric results and the following columns:
#' - "Method": The method used ("Wilson et al. (2012)").
#' - "Alpha": The chosen significance level.
#' - "Tails": The type of test ("one" or "two").
#' - "CritVal": The calculated critical value for CVR (rounded to 3 decimals).
#' - "MinJudges": The minimum number of judges needed for the critical CVR.
#'
#' @details
#' The content validity ratio (CVR; Lawshe, 1975) requires critical values for interpretation,
#' based on the alpha level chosen (e.g., 0.05 or 0.01) and the number of judges (N).
#' `CVRcut.Wilson` calculates the critical values using a normal approximation to the binomial
#' distribution, as recalculated by Wilson et al. (2012).
#'
#' @references
#' Lawshe, C. H. (1975). A quantitative approach to content validity. Personnel psychology, 28, 563-575.
#' https://doi.org/10.1111/j.1744-6570.1975.tb01393.x
#'
#' Wilson, F. R., Pan, W., & Schumsky, D. A. (2012). Recalculation of the critical values for Lawshe's
#' content validity ratio. *Measurement and Evaluation in Counseling and Development, 45*, 197-210.
#' https://doi.org/10.1177/0748175612440286
#'
#' @seealso
#'\code{\link[ValCont:CVR]{ValCont::CVR}}
#'\code{\link[ValCont:CVRcut]{ValCont::CVRcut}}
#'
#' @examples
#' CVRcut.Wilson(num_jueces = 45, alpha = .05, tails = "one")
#' CVRcut.Wilson(num_jueces = 45, alpha = .01, tails = "one")
#' @export
CVRcut.Wilson <- function(num_jueces, alpha = 0.05, tails = "one") {
  # Validar numero de jueces
  if (num_jueces < 2) {
    stop("El numero de jueces debe ser al menos 2 para un calculo valido.")
  }

  # Validar argumento tails
  tails <- tolower(tails)
  if (!tails %in% c("one", "two")) {
    stop("El argumento 'tails' debe ser 'one' o 'two'.")
  }

  # Determinar z_alpha segun el tipo de prueba
  if (tails == "one") {
    z_alpha <- stats::qnorm(1 - alpha)
  } else {
    z_alpha <- stats::qnorm(1 - alpha / 2)
  }

  # Calcular CVR crítico y numero minimo de jueces
  cutoff <- round(z_alpha / sqrt(num_jueces), 3)
  min_jueces <- ceiling(z_alpha * sqrt(num_jueces / 2) + num_jueces / 2)

  # Resultado
  result <- data.frame(
    Method = "Wilson et al. (2012)",
    Alpha = alpha,
    Tails = tails,
    CritVal = cutoff,
    MinJudges = min_jueces
  )

  return(result)
}
