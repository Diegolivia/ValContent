#' @title CVRcut: CVR Cutoff Table Generator
#' @description Builds a table of critical values for Lawshe's Content Validity Ratio (CVR), using one of three methods: Wilson (confidence intervals), Ayres, or Baghestani (Bayesian).
#'
#' @keywords validity CVR content-analysis
#'
#' @param N_min Integer. Minimum number of judges (must be < `N_max`).
#' @param N_max Integer. Maximum number of judges.
#' @param method Character. Method to compute the CVR cutoff values. Options are: `"Wilson"`, `"Ayres"`, or `"Bag"`.
#' @param alpha Numeric. Significance level (used only in `"Wilson"` and `"Bag"` methods). Default is 0.05.
#' @param tails Character. Tail type for confidence interval calculation (only for `"Wilson"` method). Options are `"one"` or `"two"`. Default is `"one"`.
#' @param prior Character. Prior distribution to use in the `"Bag"` method (e.g., `"Jeffreys"`). Default is `"Jeffreys"`.
#' @param interpolacion Character. Reserved for future use. Default is `"none"`.
#'
#' @return
#' A data frame with columns indicating the number of judges (`N`) and the corresponding CVR critical values.
#' The structure of the output depends on the method used.
#'
#' @details
#' This function serves to generate interpretation tables for the CVR coefficient
#' according to Lawshe (1975), using different statistical approaches:
#' - \strong{Wilson method}: Calculates confidence intervals for proportions.
#' - \strong{Ayres method}: Uses the original proposal by Ayres (1972).
#' - \strong{Bag method}: Bayesian method proposed by Baghestani (1993), with prior selection.
#'
#' @references
#' Lawshe, C. H. (1975). A quantitative approach to content validity. *Personnel Psychology*, 28(4), 563–575. https://doi.org/10.1111/j.1744-6570.1975.tb01393.x
#'
#' Ayre, C., & Scally, A. J. (2014). Critical Values for Lawshe’s Content Validity Ratio: Revisiting the Original Methods of Calculation. *Measurement and Evaluation in Counseling and Development*, 47(1), 79–86. https://doi.org/10.1177/0748175613513808
#'
#' Baghestani, A. R., Ahmadi, F., Tanha, A., & Meshkat, M. (2017). Bayesian Critical Values for Lawshe’s Content Validity Ratio. *Measurement and Evaluation in Counseling and Development*, 52(1), 69–73. https://doi.org/10.1080/07481756.2017.1308227
#'
#' Wilson, F. R., Pan, W., & Schumsky, D. A. (2012). Recalculation of the Critical Values for Lawshe’s Content Validity Ratio. *Measurement and Evaluation in Counseling and Development*, 45(3), 197–210. https://doi.org/10.1177/0748175612440286
#'
#' @seealso [CVRcut.Wilson()], [CVRcut.Ayres()], [CVRcut.Bag()]
#'
#' @examples
#' CVRcut(N_min = 5, N_max = 15, method = "Wilson", alpha = 0.05, tails = "two")
#'
#' @export
CVRcut <- function(N_min, N_max, method, alpha = 0.05, tails = "one",
                   prior = "Jeffreys", interpolacion = "none") {
  # Validar metodo
  if (!method %in% c("Wilson", "Ayres", "Bag")) {
    stop("Metodo no valido. Use 'Wilson', 'Ayres' o 'Bag'.")
  }

  # Validar rango de N
  if (!is.numeric(N_min) || !is.numeric(N_max) || N_min < 1 || N_max < 1 || N_min >= N_max) {
    stop("'N_min' debe ser menor que 'N_max' y ambos mayores que 0.")
  }

  # Validar alpha
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' debe ser un numero entre 0 y 1.")
  }

  # Validar tails
  if (!tails %in% c("one", "two")) {
    stop("'tails' debe ser 'one' o 'two'.")
  }

  # Crear un data.frame vacio para los resultados
  resultados <- data.frame()

  # Iterar sobre el rango de N
  for (N in N_min:N_max) {
    if (method == "Wilson") {
      res <- CVRcut.Wilson(N, alpha = alpha, tails = tails)
    } else if (method == "Ayres") {
      res <- CVRcut.Ayres(N)
    } else if (method == "Bag") {
      res <- CVRcut.Bag(N, prior = prior, alpha = alpha)
    }
    res$N <- N
    resultados <- rbind(resultados, res)
  }

  return(resultados)
}
