library(sf)
# Poner los mapas de sedimentos en una cuadrícula
par(mfrow = c(4, 4), mar = c(3, 4, 3, 2)) 

# Extraer y ordenar valores de ma
ma_values <- as.numeric(sub(".*_(\\d+)Ma.*", "\\1", basename(complete_paths)))
order_indices <- order(ma_values)

for (i in order_indices) {
  file_name <- basename(complete_paths[i])
  ma_value <- ma_values[i]
  bbox <- st_bbox(shp_list[[i]])
  xlim <- c(-150, 150)
  ylim <- c(-90, 90)
  
  # Filtrar fósiles para la edad actual
  fossils_ma <- fossil_data[fossil_data$ma == ma_value, ]
  
  # Convertir a objeto espacial sf
  if (nrow(fossils_ma) > 0) {
    fossils_sf <- st_as_sf(fossils_ma, coords = c("LONG", "LAT"), crs = 4326)
  }


  # Plot del mapa de sedimentos
  plot(st_geometry(shp_list[[i]]),
       main = paste(ma_value, "Ma"),
       col = "black",
       xlim = xlim, ylim = ylim,
       axes = FALSE)
  
  # Agregar ejes y caja
  axis(1, at = seq(-150, 150, by = 50), labels = seq(-150, 150, by = 50), las = 1)
  axis(2, at = c(-50, 0, 50), labels = c("-50", "0", "50"), las = 2)
  box()
  
  # Agregar puntos de fósiles si existen
  if (nrow(fossils_ma) > 0) {
    points(st_coordinates(fossils_sf), col = "red", pch = 19, cex = 0.2) 
  }
}
