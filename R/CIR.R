#'@title Confidence Intervals of Ratios of content validity coefficients
#'@description Calculates confidence interval for the ratio of two content validity coefficients, based on method of variance recovery for ratios (MOVER-R).
#'
#'@param group1 Output dataframe 1 of one of the 'ValCont' functions to obtain content validity coefficients.
#'@param group2 Output Output dataframe 2 of one of the 'ValCont' functions to obtain content validity coefficients.
#'@param coef.col Name of the column in the dataframe storing the calculated coefficient
#'@param lwr.col Name of the column in the dataframe that stores the lower bound of the confidence interval
#'@param upr.col Name of the column in the dataframe that stores the upper limit of the confidence interval
#'@param na.rm Logical. If FALSE (default) the function stops when missing values are detected.
#'              If TRUE rows with missing values in the relevant columns are removed before processing.
#'
#'@return
#'dataframe with four columns: label of the items, ratio between the coefficients, the upper and upper limit of the confidence interval of the difference.
#'
#'@details
#''CIR' uses method of variance recovery for ratios (MOVER-R; Zou, Donner, & Qiu, 2025), as a general approach for two non-normal quantities.
#'Because data produced by judges' judgments tend to be asymmetrically distributed (if the item is rated, on a scale of 1 to 5, as predominantly valid then its values will be > 3), MOVER-R is appropriate for non-normal distributions.
#'MOVER-R depends on the quality or precision of the confidence intervals calculated for the coefficients in each group.
#'The compared content validity coefficients obtained should be of the same type, and the estimated confidence intervals for these coefficients should also come from the same level; for example, at .95 or .90.
#'If the two dataframes have different numbers of evaluated items (i.e., different numbers of rows), 'CIR' function matches the commonly labeled items, assuming they are the same items.
#'Note: The function has not yet been prepared to resolve missing values, so the user must remove or impute any missing values.
#'
#'@references
#'Zou, G., Donner, A. & Qiu, S. (2025). MOVER-R for Confidence Intervals of Ratios. In  N. Balakrishnan, T. Colton, B. Everitt, W. Piegorsch, F. Ruggeri and J.L. Teugels (Eds.), Wiley StatsRef: Statistics Reference Online. https://doi.org/10.1002/9781118445112.stat08085
#'
#' @author
#'Cesar Merino-Soto (\email{sikayax@yahoo.cam.ar})
#'
#'@seealso
#'\code{\link[ratesci:moverci]{ratesci::moverci}}
#'\code{\link[ValCont:CID]{ValCont::CID}}
#'
#'@examples
#'
#'### Example 1 -----------
#'
#'## Random data (Low ratings): 11 items (rows), 4 raters (columns)
#'Data1 <- data.frame(
#'  juez1 = sample(1:2, 11, replace = TRUE),
#'  juez2 = sample(2:3, 11, replace = TRUE),
#'  juez3 = sample(1:2, 11, replace = TRUE),
#'  juez4 = sample(1:3, 11, replace = TRUE))
#'
#'## Random data (High ratings): 10 items (rows), 6 raters (columns)
#'Data2 <- data.frame(
#'  obs1 = sample(5:7, 10, replace = TRUE),
#'  obs2 = sample(4:7, 10, replace = TRUE),
#'  obs3 = sample(4:7, 10, replace = TRUE),
#'  obs4 = sample(5:7, 10, replace = TRUE),
#'  obs5 = sample(5:7, 10, replace = TRUE),
#'  obs6 = sample(6:7, 10, replace = TRUE))
#'
#'## Saving results from Aiken's V analysis, for each group
#'group1 <- Vaiken(data = Data1, min = 1, max = 7, conf.level = .90)
#'group2 <- Vaiken(data = Data2, min = 1, max = 7, conf.level = .90)
#'
#'CIR(group1 = group1,
#'    group2 = group2,
#'    coef.col = "V",
#'    lwr.col = "lwr.ci",
#'    upr.col = "upr.ci")
#'
#'
#' @export
CIR <- function(group1, group2, coef.col, lwr.col, upr.col, na.rm = FALSE) {
  # Validaciones iniciales
  if (!all(c(coef.col, lwr.col, upr.col) %in% colnames(group1))) {
    stop("group1 debe contener las columnas especificadas en coef.col, lwr.col y upr.col")
  }
  if (!all(c(coef.col, lwr.col, upr.col) %in% colnames(group2))) {
    stop("group2 debe contener las columnas especificadas en coef.col, lwr.col y upr.col")
  }

    # Detectar valores perdidos en las columnas relevantes
  cols_g1 <- intersect(c(coef.col, lwr.col, upr.col, "Item"), colnames(group1))
  cols_g2 <- intersect(c(coef.col, lwr.col, upr.col, "Item"), colnames(group2))

  has_na_g1 <- any(is.na(group1[, cols_g1, drop = FALSE]))
  has_na_g2 <- any(is.na(group2[, cols_g2, drop = FALSE]))

  if (has_na_g1 || has_na_g2) {
    if (!na.rm) {
      stop("Missing values detected. Use na.omit() first or set na.rm=TRUE")
    } else {
      # Eliminar filas con NA en las columnas de coef, lwr o upr (no es necesario Item para complete.cases)
      keep_g1 <- complete.cases(group1[, intersect(c(coef.col, lwr.col, upr.col), colnames(group1)), drop = FALSE])
      keep_g2 <- complete.cases(group2[, intersect(c(coef.col, lwr.col, upr.col), colnames(group2)), drop = FALSE])
      nrem1 <- sum(!keep_g1)
      nrem2 <- sum(!keep_g2)
      group1 <- group1[keep_g1, , drop = FALSE]
      group2 <- group2[keep_g2, , drop = FALSE]
      warning(sprintf("na.rm=TRUE: removed %d rows with missing values from group1 and %d from group2", nrem1, nrem2))
    }
  }

  # Identificar items comunes
  common_items <- intersect(group1$Item, group2$Item)

  if (length(common_items) < nrow(group1) || length(common_items) < nrow(group2)) {
    warning("Algunos items no coinciden entre group1 y group2. Se procesaran unicamente los items comunes.")
  }

  # Filtrar los data.frames por items comunes
  group1 <- group1[group1$Item %in% common_items, ]
  group2 <- group2[group2$Item %in% common_items, ]

  # Combinar ambos data.frames por Item
  combined <- merge(group1, group2, by = "Item", suffixes = c("_g1", "_g2"))

  # Inicializar las columnas de salida
  results <- data.frame(
    Item = combined$Item,
    R = NA,
    lwr.ci = NA,
    upr.ci = NA
  )

  # Calcular R y los intervalos de confianza para cada item
  for (i in 1:nrow(combined)) {
    coef_g1 <- combined[[paste0(coef.col, "_g1")]][i]
    coef_g2 <- combined[[paste0(coef.col, "_g2")]][i]
    lwr_g1 <- combined[[paste0(lwr.col, "_g1")]][i]
    upr_g1 <- combined[[paste0(upr.col, "_g1")]][i]
    lwr_g2 <- combined[[paste0(lwr.col, "_g2")]][i]
    upr_g2 <- combined[[paste0(upr.col, "_g2")]][i]

    R <- coef_g1 / coef_g2
    var_g1 <- (upr_g1 - lwr_g1)^2 / 4
    var_g2 <- (upr_g2 - lwr_g2)^2 / 4
    lwr.ci <- max(0, R - sqrt(var_g1 + var_g2))
    upr.ci <- R + sqrt(var_g1 + var_g2)

    # Guardar los resultados
    results$R[i] <- round(R, 3)
    results$lwr.ci[i] <- round(lwr.ci, 3)
    results$upr.ci[i] <- round(upr.ci,3)
  }

  return(results)
}

