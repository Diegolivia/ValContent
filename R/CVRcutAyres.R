#' @title Critical value for CVR, based on Ayres & Scally (2014)
#' @description Calculates critical values for Lawshe’s content validity ratio (CVR) using the method
#' proposed by Ayres & Scally (2014), based on exact binomial probabilities.
#' @param num_jueces Number of judges who provided their ratings (must be >= 2).
#' @param alpha Significance level for the critical value calculation (default: 0.05).
#' @param tails Distribution tails: "one" for one-tailed test (default) or "two" for two-tailed test.
#'
#' @return
#' A data frame with numeric results and four columns:
#' - "Method": The method used (Ayres & Scally).
#' - "Alpha": The chosen significance level.
#' - "Tails": The type of test ("one" or "two").
#' - "CritVal": The calculated critical value for CVR.
#' - "MinJudges": The minimum number of judges needed for the critical CVR.
#'
#' @details
#' Content validity ratio (CVR; Lawshe, 1975) requires critical values for interpretation,
#' based on the alpha level chosen (0.05 or 0.01) and the number of judges. `CVRcut.Ayres`
#' estimates this critical value based on the Exact Binomial Probabilities method, an approach
#' apparently appropriate for small samples, as is usual in content validity studies,
#' and for proportions approaching 0 or 1. The function allows for both one-tailed and two-tailed
#' tests.
#'
#' @references
#' Ayres, C., & Scally, A. J. (2014). Critical values for Lawshe's content validity ratio. Measurement and Evaluation in Counseling and Development, 47, 79–86. https://doi.org/10.1177/0748175613513808
#'
#' Lawshe, C. H. (1975). A quantitative approach to content validity. Personnel psychology, 28, 563–575. https://doi.org/10.1111/j.1744-6570.1975.tb01393.x
#'
#'@seealso
#'\code{\link[ValContent:CVR]{ValContent::CVR}}
#'\code{\link[ValContent:CVRcut]{ValContent::CVRcut}}
#' @examples
#'
#' ### Example 1 ----------
#' # N judges = 10, alpha = .05
#' # Expected output:
#' CVRcut.Ayres(num_jueces = 10, alpha = .05)
#'
#' ### Example 2 ----------
#' # N judges = 45, alpha = .05
#' # (Expected output: Critical value and minimum number of judges for CVR, in alpha .05.)
#' CVRcut.Ayres(num_jueces = 45, alpha = .05)
#'
#' ### Example 3 ----------
#' # N judges = 45, alpha = .01
#' # (Expected output: Lower alpha (.01), higher critical values.)
#' CVRcut.Ayres(num_jueces = 45, alpha = .01)
#'
#' @export
CVRcut.Ayres <- function(num_jueces, alpha = 0.05, tails = "one") {
  # Validaciones
  if (!is.numeric(num_jueces) || num_jueces < 2 || num_jueces != as.integer(num_jueces)) {
    stop("El numero de jueces debe ser un entero mayor o igual a 2.")
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("El valor de 'alpha' debe estar entre 0 y 1 (exclusivo).")
  }

  # Guardar valor original de alpha
  alpha_original <- alpha

  # Ajustar alpha si es una prueba de dos colas
  if (tails == "two") {
    alpha <- alpha / 2
  } else if (tails != "one") {
    stop("El argumento 'tails' debe ser 'one' o 'two'.")
  }

  # Funcion interna para calcular minimo n_e
  calcular_ne <- function(N, alpha) {
    for (ne in 0:N) {
      p_cum <- 1 - pbinom(ne - 1, size = N, prob = 0.5)
      if (p_cum <= alpha) {
        return(ne)
      }
    }
    return(NA)
  }

  # Calcular valores
  min_jueces <- calcular_ne(num_jueces, alpha)
  cutoff <- if (!is.na(min_jueces)) (2 * min_jueces - num_jueces) / num_jueces else NA

  # Resultado
  result <- data.frame(
    Method = "Ayres & Scally (2014)",
    Alpha = alpha_original,
    Tails = tails,
    CritVal = cutoff,
    MinJudges = min_jueces
  )

  return(result)
}
