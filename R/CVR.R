#'@title Content Validity Ratio
#'@description Calculate the CVR coefficient (content validity ratio)
#' 
#'@param data dataframe, with the columns assigned to each judge, and the rows assigned to each evaluated item.
#'
#'@return
#'dataframe with CVR coefficients for all analyzed items
#'
#'@details
#'The CVR function calculates the content validity ratio (CVR; Lashe, 1975). The CVR function analyses data in which the items were rated in three categories: Not important (1), useful but not essential (2), essential (3). CVR dichotomizes each rating, and yields the proportion of judges who considered the item essential.
#'Under this model, the usual cut-off point for identifying essential items in the CVR should be located at the highest possible rating, i.e., on a scale of 1 to 3, the cut-off point is 3.
#'
#'Note: The function has not yet been prepared to resolve missing values, so the user must remove or impute any missing values.
#'
#'@references
#'Lawshe, C. H. (1975). A quantitative approach to content validity. Personnel psychology, 28, 563–575. https://doi.org/10.1111/j.1744-6570.1975.tb01393.x
#'
#'Martuza, V.R. (1977). Applying norm-referenced and criterion-referenced measurement in education. Boston: Allyn & Bacon
#'
#'@author
#' Cesar Merino-Soto (\email{sikayax@yahoo.com.ar})
#'
#'@examples
#'
#'## Artificial data
#'CVRdata <- data.frame(j1 = c(3,3,3,3,3,3,3,2,2,2),
#'j2 = c(3,2,3,2,2,2,3,3,3,3),
#'j3 = c(3,2,3,2,2,2,3,3,3,2),
#'j4 = c(3,2,3,3,2,3,3,3,3,3),
#'j5 = c(3,2,3,3,2,2,2,2,3,2))
#'
#'## Run CVR
#'CVR(data = CVRdata)
#'
#'@export
CVR <- function(data) {
# Transformar las respuestas: 3 -> 1 (essential), otros -> 0
data_binaria <- as.data.frame(t(apply(data, 1, function(row) {
  ifelse(row == 3, 1, 0)
})))

# Crear un data.frame para almacenar los resultados
resultados <- data.frame(
  Item = seq_len(nrow(data_binaria)),
  CVR = numeric(nrow(data_binaria))
)

# Iterar sobre las filas (items)
for (i in seq_len(nrow(data_binaria))) {
  # Contar el numero de respuestas "essential" (ne) y total (N)
  ne <- sum(data_binaria[i, ])  # Respuestas 1 ("essential")
  N <- ncol(data_binaria)       # Total de respuestas (jueces)

  # Calcular CVR para el item
  CVR <- (ne - (N / 2)) / (N / 2)

  # Almacenar los resultados
  resultados$CVR[i] <- round(CVR, 3)
}

# Retornar el data.frame con los resultados
return(resultados)
}
