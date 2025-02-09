MER <- function(data, ncat, start, conf.level = 0.90) {
  # Verificaciones iniciales
  if (!is.data.frame(data)) stop("'data' debe ser un data.frame")
  if (!is.numeric(ncat) || ncat <= 1) stop("'ncat' debe ser un numero mayor que 1")
  if (!is.numeric(start) || !start %in% c(0, 1)) stop("'start' debe ser 0 o 1")
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' debe estar entre 0 y 1")
  }

  # Ajuste del nivel de confianza
  alpha <- 1 - conf.level
  
  # C�lculo de la media y del intervalo de confianza para cada �tem
  results <- apply(data, 1, function(item_ratings) {
    mean_rating <- mean(item_ratings) # Media del �tem
    n <- length(item_ratings)        # N�mero de jueces
    sd_rating <- sd(item_ratings)    # Desviaci�n est�ndar
    
    # C�lculo de los l�mites del IC asim�trico
    error_margin <- qt(1 - alpha / 2, df = n - 1) * (sd_rating / sqrt(n))
    lwr <- max(mean_rating - error_margin, start) # L�mite inferior ajustado al m�nimo posible
    upr <- min(mean_rating + error_margin, start + ncat - 1) # L�mite superior ajustado al m�ximo posible
    
    c(MER = round(mean_rating, 3), 
      lwr.ci = round(lwr, 3), 
      upr.ci = round(upr, 3))
  })
  
  # Transformar los resultados en un data.frame
  results_df <- as.data.frame(t(results))
  results_df$Item <- seq_len(nrow(results_df)) # Agregar la columna 'Item'
  results_df <- results_df[, c("Item", "MER", "lwr.ci", "upr.ci")] # Reordenar columnas
  
  return(results_df)
}
