#'@title Confidence interval for difference: single items
#'@description Calculate confidence intervals for the difference of specific content validity coefficients on a single item, based on Method of Variance Estimates Recovery (MOVER).
#'@param coef1 numerical value obtained for coefficient 1
#'@param coef2 numerical value obtained for coefficient 2
#'@param lci1 Lower limit for coefficient 1
#'@param uci1 Upper limit for coefficient 1
#'@param lci2 Lower limit for coefficient 2
#'@param uci2 Upper limit for coefficient 2
#'
#'@return
#'dataframe with Confidence intervals for the difference of the analyzed item, and confidence intervals.
#'
#'@details
#' `CID` uses Method of Variance Estimates Recovery (MOVER; Zou, & Donner, 2008). Because data produced by judges'
#'judgments tend to be asymmetrically distributed (if the item is rated, on a scale of 1 to 5, as predominantly
#'valid then its values will be > 3), MOVER is appropriate for non-normal distributions.
#'The application of MOVER for content validity coefficient was initially published by Merino-Soto (2018) for the
#'difference between V coefficients (Aiken, 1980, 1985). Later, Merino-Soto (2023) extended this approach for
#'Aiken's V by adding a point estimator of the difference, based on the standardized difference between proportions.
#'There is \code{\link[ValCont:CID]{ValCont::CID}} for this method for multiple items, however `CIDsingle` function compares the
#'results of two specific content validity coefficients with their respective confidence intervals. For example, the
#'user can compare the coefficient calculated on his data (however, this coefficient must also have its confidence
#'intervals calculated), with another coefficient obtained from a previous study.
#'The compared content validity coefficients obtained should be of the same type, and the estimated confidence intervals
#'for these coefficients should also come from the same level; for example, at .95 or .90.
#'Singer (2010) observed that at extremely low values (e.g., proportions near .0), the coverage of this method is not
#'as good. In the context of comparing content validity coefficients, treated as proportions, it is rare to find such
#'low coefficients (and their confidence intervals). However, unless the items have been constructed extremely poorly
#'in content, it is unlikely to obtain content validity coefficients close to zero.
#'
#'@references
#'Aiken, L. R. (1980). Content validity and reliability of single items or questionnaires. \emph{Educational and. Psychological Measurement, 40}, 955-959. \doi{10.1177/001316448004000419}
#'
#'Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. \emph{Educational and Psychological Measurement, 45}, 131-142. \doi{10.1177/0013164485451012}
#'
#'Merino-Soto, C. (2018) Confidence interval for difference between coefficients of content validity (Aiken's V): a SPSS syntax. \emph{Anales de  Psicologia, 34}(3), 587-590. \doi{10.6018/analesps.34.3.283481}
#'
#'Merino-Soto, C. (2023). Coeficientes V de Aiken: diferencias en los juicios de validez de contenido. \emph{MHSalud, 20}(1), 23-32. \doi{10.15359/mhs.20-1.3}
#'
#'Singer, J. (2010). Construction of confidence limits about effect measures: A general approach, by G. Y. Zou and A. Donner, Statistics in Medicine 2008; 27:1693-1702. \emph{Statistics in Medicine, 29}(16), 1757–1759. \doi{10.1002/sim.3887}
#'
#'Zou, G.Y. and Donner, A. (2008) Construction of confidence limits about effect measures: a general approach. \emph{Statistics in Medicine, 27}, 1693–1702. \doi{10.1002/sim.3095}
#'
#'@seealso
#'\code{\link[ratesci:moverci]{ratesci::moverci}} for MOVER method ofr ratios
#'\code{\link[PropCIs:scoreci]{PropCIs::scoreci}} for score method confidence interval
#'\code{\link[ValCont:CID]{ValCont::CID}} for this method for multiple items
#'
#'@examples
#'
#'# Example
#'CIDsingle(
#'  coef1 = 0.85, coef2 = 0.78,
#'  lci1 = 0.80, uci1 = 0.90,
#'  lci2 = 0.73, uci2 = 0.83
#')
#'
#'@author
#'Cesar Merino-Soto (\email{sikayax@yahoo.com.ar})
#'
#'
#'@export
CIDsingle <- function(coef1, coef2, lci1, uci1, lci2, uci2) {
  # Verify that all arguments are numeric
  if (!all(sapply(list(coef1, coef2, lci1, uci1, lci2, uci2), is.numeric))) {
    stop("All arguments must be numeric values.")
  }

  # Validar que los limites sean coherentes
  if (!(lci1 <= coef1 && coef1 <= uci1)) stop("The coefficient 1 is not within its confidence interval.")
  if (!(lci2 <= coef2 && coef2 <= uci2)) stop("The coefficient 2 is not within its confidence interval.")

  # Calculate the difference between the coefficients
  difference <- coef1 - coef2

  # Calculate the confidence interval for the difference using MOVER
  lower_diff <- difference - sqrt((coef1 - lci1)^2 + (uci2 - coef2)^2)
  upper_diff <- difference + sqrt((uci1 - coef1)^2 + (coef2 - lci2)^2)

  # Return the results
  return(data.frame(
    Difference = round(difference, 3),
    lwr.ci = round(lower_diff, 3),
    upr.ci = round(upper_diff, 3)
  ))
}
