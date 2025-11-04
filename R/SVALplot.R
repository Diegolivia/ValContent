#' @title Plot Substantive Validity Coefficients
#' @description
#' Visualizes the results from \code{SVALmult()} for a specific item and coefficient
#' type (\code{"svc"} or \code{"psa"}), including confidence intervals. Optionally,
#' a target construct category can be highlighted and labeled.
#'
#' @param results A list returned by \code{SVALmult()}, containing \code{psa_results} and \code{svc_results} per item.
#' @param item Character. Name of the item to be plotted (must match one of the names in \code{results}).
#' @param type Character. Coefficient to plot. Either \code{"svc"} or \code{"psa"}. Defaults to \code{"svc"}.
#' @param target Optional. Numeric value indicating the target construct category. It will be marked in the plot.
#' @param labels Optional. Character vector to replace category codes with meaningful labels (must match the number of categories).
#'
#' @return
#' A \code{ggplot2} object displaying the selected coefficient values with confidence intervals.
#'
#' @details
#' This function allows plotting substantive validity coefficients obtained from \code{SVALmult()}.
#'  If a target construct category is provided via the \code{target} argument,
#'  it will be highlighted. You can also provide a vector of labels to improve interpretability of the category axis.
#'
#' @seealso
#'\code{\link[ValContent:SVALmult]{ValContent::SVALmult}} to compute substantive validity coefficients across multiple items. \cr
#'\code{\link[ValContent:SVALsingle]{ValContent::SVALsingle}} for single-item analysis.
#'
#' @examples
#' \dontrun{
#' results <- SVALmult(data = data.gais, columns = c("gais2", "gais3"))
#' SVALplot(results, item = "gais3", type = "svc", target = 4)
#' SVALplot(results, item = "gais3", type = "psa", labels = c("Anger", "Fear", "Joy", "Sadness"))}
#'
#' @author
#' César Merino-Soto (\email{sikayax@yahoo.com.ar})
#'
#' @export
SVALplot <- function(results, item, type = "svc", target = NULL, labels = NULL) {
  if (!item %in% names(results)) stop("Item not found in results.")
  if (!type %in% c("svc", "psa")) stop("Type must be either 'svc' or 'psa'.")

  df <- results[[item]][[paste0(type, "_results")]]

  # Crear variable auxiliar para labels del eje X
  df$CategoryLabel <- as.character(df$Cat)

  # Reemplazar por labels si el usuario los proporciona
  if (!is.null(labels)) {
    if (length(labels) != nrow(df)) stop("Length of 'labels' must match number of categories.")
    df$CategoryLabel <- labels
  }

  # Añadir asterisco al target si se indica
  if (!is.null(target)) {
    df$CategoryLabel[df$Cat == target] <- paste0(df$CategoryLabel[df$Cat == target], "*")
  }

  ggplot(df, aes(x = CategoryLabel, y = .data[[type]])) +
    geom_point(size = 3, color = "#0072B2") +
    geom_errorbar(
      aes(ymin = lwr.ci, ymax = up.ci), width = 0.2, color = "#0072B2"
    ) +
    labs(
      title = paste("Substantive Validity -", toupper(type), "-", item),
      x = "Construct Category", y = toupper(type)
    ) +
    theme_minimal()
}

