#'@title Confidence interval for diference of content validity coefficients
#'@description Calculates confidence interval for the difference of content validity coefficients, based on Method of Variance Estimates Recovery (MOVER).
#'
#'@param group1 Output dataframe 1 of one of the 'ValCont' functions to obtain content validity coefficients.
#'@param group2 Output Output dataframe 2 of one of the 'ValCont' functions to obtain content validity coefficients.
#'@param coef.col Name of the column in the dataframe storing the calculated coefficient
#'@param lwr.col Name of the column in the dataframe that stores the lower bound of the confidence interval
#'@param upr.col Name of the column in the dataframe that stores the upper limit of the confidence interval
#'@param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#'@return
#'dataframe with four columns: label of the items, difference between the coefficients, the upper and upper limit of the confidence interval of the difference.
#'
#'@details
#'`CID` uses Method of Variance Estimates Recovery (MOVER; Zou, & Donner, 2008).
#'Because data produced by judges' judgments tend to be asymmetrically distributed (if the item is rated, on a scale of 1 to 5, as predominantly valid then its values will be > 3), MOVER is appropriate for non-normal distributions.
#'MOVER depends on the quality or precision of the confidence intervals calculated for the coefficients in each group.
#'The application of MOVER for content validity coefficient was initially published by Merino-Soto (2018) for the difference between V coefficients (Aiken, 1980, 1985). Later, Merino-Soto (2023) extended this approach for Aiken's V by adding a point estimator of the difference, based on the standardized difference between proportions
#'The compared content validity coefficients obtained should be of the same type, and the estimated confidence intervals for these coefficients should also come from the same level; for example, at .95 or .90.
#'Singer (2010) observed that at extremely low values (e.g., proportions near .0), the coverage of this method is not as good. In the context of comparing content validity coefficients, treated as proportions, it is rare to find such low coefficients (and their confidence intervals). Unless the items are extremely poor in content.
#'
#'Note: The function has not yet been prepared to resolve missing values, so the user must remove or impute any missing values.
#'
#'@references
#'Aiken, L. R. (1980). Content validity and reliability of single items or questionnaires. Educational and. Psychological Measurement, 40, 955-959. https://doi.org/10.1177/001316448004000419
#'
#'Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. Educational and Psychological Measurement, 45, 131-142. https://doi.org/10.1177/0013164485451012
#'
#'Merino-Soto, C. (2018) Confidence interval for difference between coefficients of content validity (Aiken's V): a SPSS syntax. Anales de  Psicologia, 34(3), 587-590. https://doi.org/10.6018/analesps.34.3.283481
#'
#'Merino-Soto, C. (2023). Coeficientes V de Aiken: diferencias en los juicios de validez de contenido. MHSalud, 20(1), 23-32. https://doi.org/10.15359/mhs.20-1.3
#'
#'Singer, J. (2010). Construction of confidence limits about effect measures: A general approach, by G. Y. Zou and A. Donner, Statistics in Medicine 2008; 27:1693-1702. Statistics in Medicine, 29(16), 1757–1759. https://doi.org/10.1002/sim.3887
#'
#'Zou, G.Y. and Donner, A. (2008) Construction of confidence limits about effect measures: a general approach. Stat. Med. 27, 1693–1702. https://doi.org//10.1002/sim.3095
#'
#'@author
#'Cesar Merino-Soto (\email{sikayax@yahoo.cam.ar})
#'
#'@seealso
#'\code{\link[ratesci:moverci]{ratesci::moverci}}
#'\code{\link[ValCont:CIDsingle]{ValCont::CIDsingle}}
#'
#'@examples
#'
#'### Example 1 -----------
#'
#'## Group 1 output
#'
#'Vgroup1 <- data.frame(
#'  V = c(0.85, 0.78, 0.90),
#'  lwr.ci = c(0.80, 0.75, 0.88),
#'  upr.ci = c(0.90, 0.82, 0.92)
#')
#'rownames(Vgroup1) <- c("Item1", "Item2", "Item3")
#'
#'## Group 2 output
#'
#'Vgroup2 <- data.frame(
#'  V = c(0.80, 0.75),
#'  lwr.ci = c(0.76, 0.70),
#'  upr.ci = c(0.84, 0.78)
#')
#'rownames(Vgroup2) <- c("Item1", "Item2")
#'
#'## Run
#'CID(Vgroup1, Vgroup2,"V","lwr.ci","upr.ci")
#'
#'### Example 2 -----------
#'
#'## Random data (Low ratings): 11 items (rows), 4 raters (columns)
#'
#'random_data2 <- data.frame(
#'  juez1 = sample(1:2, 11, replace = TRUE),
#'  juez2 = sample(2:3, 11, replace = TRUE),
#'  juez3 = sample(1:2, 11, replace = TRUE),
#'  juez4 = sample(1:3, 11, replace = TRUE))
#'
#'## Random data (High ratings): 10 items (rows), 6 raters (columns)
#'
#'random_data3 <- data.frame(
#'  obs1 = sample(5:7, 10, replace = TRUE),
#'  obs2 = sample(4:7, 10, replace = TRUE),
#'  obs3 = sample(4:7, 10, replace = TRUE),
#'  obs4 = sample(5:7, 10, replace = TRUE),
#'  obs5 = sample(5:7, 10, replace = TRUE),
#'  obs6 = sample(6:7, 10, replace = TRUE))
#'
#'## Saving results from Aiken's V analysis, for each group
#'
#'output.V1 <- Vaiken(data = random_data2, min = 1, max = 7, conf.level = .90)
#'output.V2 <- Vaiken(data = random_data3, min = 1, max = 7, conf.level = .90)
#'
#'## Run CID
#'
#'CID(group1 = output.V2,
#'    group2 = output.V1,
#'    coef.col = "V",
#'    lwr.col = "lwr.ci",
#'    upr.col = "upr.ci")
#'
#'@export
CID <- function(group1, group2, coef.col = "coef", lwr.col = "lwr.ci", upr.col = "upr.ci", na.rm = FALSE) {
  # Validar que los argumentos son data.frames
  if (!is.data.frame(group1) || !is.data.frame(group2)) {
    stop("Ambos argumentos 'group1' y 'group2' deben ser data.frames.")
  }

  # Detección de valores perdidos
  if (!na.rm) {
    if (any(is.na(group1)) || any(is.na(group2))) {
      stop("Valores perdidos detectados. Usa na.omit() primero o establece na.rm=TRUE.")
    }
  } else {
    group1 <- na.omit(group1)
    group2 <- na.omit(group2)
  }

  # Validar que las columnas necesarias existen en los data.frames
  required_cols <- c(coef.col, lwr.col, upr.col)
  if (!all(required_cols %in% colnames(group1))) {
    stop("El data.frame 'group1' no contiene las columnas necesarias.")
  }
  if (!all(required_cols %in% colnames(group2))) {
    stop("El data.frame 'group2' no contiene las columnas necesarias.")
  }

  # Combinar los datos para realizar las comparaciones
  combined_data <- merge(
    group1,
    group2,
    by = "row.names",
    suffixes = c("_1", "_2")
  )
  rownames(combined_data) <- combined_data$Row.names
  combined_data$Row.names <- NULL

  # Verificar si hubo exclusiones debido a items ausentes en un grupo
  if (nrow(combined_data) < nrow(group1) || nrow(combined_data) < nrow(group2)) {
    warning("Algunos items no tienen correspondencia entre los grupos y han sido excluidos de la comparacion.")
  }

  # Calcular la diferencia de los coeficientes y los intervalos de confianza
  combined_data <- within(combined_data, {
    Delta <- round(get(paste0(coef.col, "_1")) - get(paste0(coef.col, "_2")), 3)
    lwr.ci <- round(
      Delta - sqrt((get(paste0(coef.col, "_1")) - get(paste0(lwr.col, "_1")))^2 +
                     (get(paste0(upr.col, "_2")) - get(paste0(coef.col, "_2")))^2), 3)
    upr.ci <- round(
      Delta + sqrt((get(paste0(upr.col, "_1")) - get(paste0(coef.col, "_1")))^2 +
                     (get(paste0(coef.col, "_2")) - get(paste0(lwr.col, "_2")))^2), 3)
  })

  # Seleccionar las columnas finales para el reporte y reiniciar nombres de filas
  result <- combined_data[, c("Delta", "lwr.ci", "upr.ci")]
  row.names(result) <- NULL

  return(result)
}
