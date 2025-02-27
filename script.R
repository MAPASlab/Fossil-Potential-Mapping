## METADATA ===============================================================
## Description: 
## scripts for identifying regions with high fossil preservation potential using 
## climate data (HadCM3 model (Valdes et al., 2017)), fossil records (The NOW Community (2024) and Paleobiology Database (2024)), and sediment zones (Chorlton, 2007). 
## We apply the Köppen-Geiger climate reclassification function (Galván et al., 2023; https://github.com/MAPASlab/KoppenGeiger_inR) to map biomes across 15 geological periods. 
## 
## R version: 4.2.2 for Windows
## Date: 2024-10-18 16:59:26
## License: GPL3 (Marta Please check)
## Author: Marta Matamala, Sara Varela, Oskar Hagen.
##=======================================================================##

library(terra)
library(rgplates)
library(sf)
library(maps)
library(dplyr)
source("support_functions.R")

########################
#DATA PREPARATION
# 1) WORLD DATA

# Unzip the data file (only the first time)
untar("./data/scotese/Scotese_temp_precip_Ceno.tar.nc", exdir="./data/scotese")

# LIST AND LOAD TEMPERATURE DATA
# List temperature files
list_temp.nc <- list.files("./data/scotese/formatted_data/tfkea", pattern = ".temp.nc")
list_prec.nc <- list.files("./data/scotese/formatted_data/tfkea", pattern = ".precip.nc")

# Create a raster stack from the temperature files
stack_temp <- rotate (rast(file.path("./data/scotese/formatted_data/tfkea", list_temp.nc)))
stack_prec <- rotate (rast(file.path("./data/scotese/formatted_data/tfkea", list_prec.nc)))
names(stack_temp)<- list_temp.nc
names(stack_prec)<- list_prec.nc


# dir <- "./data/scotese/formatted_data"
# setwd(dir)
scotese_times <- data.frame(
  "name" = list.files("./data/scotese/formatted_data"),
  "ma" = c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
)

# Create a list of month abbreviations
month_ID <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

# Initialize lists for temperature and precipitation data
alldata_temp <- list()
alldata_prec <- list()


# Process data for each time period
dir <- "./data/scotese/formatted_data"
for (i in 1:nrow(scotese_times)) {
  # i <- 1
  time_ID <- scotese_times[i, 1]
  my_result_temp <- monthly_stacks(time_ID, "1_5m_temp", dir)
  my_result_prec <- monthly_stacks(time_ID, "precip", dir)
  alldata_temp[[i]] <- my_result_temp
  alldata_prec[[i]] <- my_result_prec
}

# Name the lists
ma <- c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
names(alldata_temp) <- paste("temp_", ma, sep = "")
names(alldata_prec) <- paste("prec_", ma, sep = "")
myears<- paste (ma, "Ma", sep="_")

# create stack of world climatic data for temp and prec with average prec and average temp
world_temp<- rast()
world_prec<- rast()

for (i in 1:15){
  temp <-rast(alldata_temp [[i]])
  world_temp<- c(world_temp, mean (temp))
  names (world_temp)[i]<- myears [i]
  prec <- rast(alldata_prec [[i]])
  world_prec<- c(world_prec, sum(prec))
  names (world_prec) [i] <- myears [i]
} 

#plot (world_temp, world_prec, xlim=c(-40, 40), ylim=c(0,6000), main=ma)

########################
# 2) SEDIMENTS

sediment_files <- list.files("./data/possible_fossil_reconstructed_dissolved",  pattern = "\\.shp$")
# create stack of sediment data for temp and prec with average prec and average temp
allextract_temp <- list()
allextract_prec <- list()


for (i in (1:length(sediment_files))) { 
  shapefile <- read_sf(file.path("./data/possible_fossil_reconstructed_dissolved", paste0("possible_fossil_recons_", ma [i], "Ma_dissolved.shp")))
  shapefile<- vect (shapefile)
  extracted_data_temp <- terra::extract (world_temp [[i]], shapefile, xy = TRUE, ID = FALSE)
  extracted_data_prec <- terra::extract (world_prec [[i]], shapefile, xy = TRUE, ID = FALSE)
  allextract_temp[[paste0(ma [i], "Ma")]] <- extracted_data_temp
  allextract_prec[[paste0(ma [i], "Ma")]] <- extracted_data_prec
}


for (i in (1:length(sediment_files))) {
  rast_temp_name <- paste0("rast", ma[i], "_temp")
  assign(rast_temp_name, rast (na.omit(allextract_temp[[i]][, c(2, 3, 1)])))
  
  rast_prec_name <- paste0("rast", ma[i], "_prec")
  assign(rast_prec_name, rast(na.omit(allextract_prec[[i]][, c(2, 3, 1)])))
}

# str (allextract_prec)

########################
# 3) FOSSILS

setwd("../DATA")
fossil_data <- read.csv("fossil_data/pbdb_data_27-02-25processed.csv")

#str(fossil_data)
#head(fossil_data)
#summary(fossil_data)

# define the period limits (to reclassify the fossils within our working periods)
periods <- setNames(lapply(ma, function(x) {
  if (x==0){c(0,0.01)}
  else if (x==3){c(x - 0.2, x + 0.2)}
  else { c(x - 1, x + 1)}
}), as.character(ma))

periods <- list(
  "0" = c(0, 0.01),          
  "3" = c(2.8, 3.2),         
  "11" = c(10, 12),          
  "15" = c(14, 16),          
  "20" = c(19, 21),          
  "26" = c(25, 27),          
  "31" = c(30, 32),          
  "36" = c(35, 37),          
  "40" = c(39, 41),          
  "45" = c(44, 46),          
  "52" = c(51, 53),          
  "56" = c(55, 57),          
  "61" = c(60, 62),          
  "66" = c(65, 67),          
  "69" = c(68, 70))

# reclassify the fossils to the period (ma) that corresponds to our conditions

fossil_data$ma <- NA
for (i in 1:15) {
  fossil_data$ma[fossil_data$MIN_AGE >= min(periods[[i]]) & fossil_data$MAX_AGE <= max(periods[[i]])] <- names(periods[i])
  fossil_data$ma[fossil_data$MIN_AGE <= max(periods[[i]]) & fossil_data$MAX_AGE >= min(periods[[i]])] <- names(periods[i])
}


# convert fossil to work
fossil_data<- as.data.frame (fossil_data)
str (fossil_data)
fossil_data$ma<- as.numeric (as.character (fossil_data$ma))
fossil_data_old<- fossil_data
fossil_data<- fossil_data_old [complete.cases (fossil_data_old), ]

# apply palaeorotations (rgplates) to all fossil occurances (palaeolat and palaeolon)
fossil_data$paleolong<- NA
fossil_data$paleolat<- NA
names (fossil_data)


for (i in sort (unique (fossil_data$ma))){
  fossil_data [fossil_data$ma==i, 9:10]<- reconstruct(fossil_data [fossil_data$ma==i, 4:3], 
                                                      age = i, model = "PALEOMAP")
}

# add prec and temp data to each fossil as a function of paleolat and paleolon
sort (unique (fossil_data$ma))
# world_temp 
fossil_data$prec<- NA
fossil_data$temp<- NA
tiempos <- sort (unique (fossil_data$ma))


for (i in 1:length(tiempos)){
  # create a mask for time
  tiempos_mask <- fossil_data$ma==tiempos [i]
  fossil_data$temp [tiempos_mask]<- extract (world_temp [[i]], fossil_data[tiempos_mask ,9:10 ])[, 2]
  fossil_data$prec [tiempos_mask]<- extract (world_prec[[i]], fossil_data[tiempos_mask ,9:10 ])[, 2]
}

########################
# PLOTS

# plot the sediments for the 15 periods

dir <- "possible_fossil_reconstructed_dissolved"
complete_paths <- file.path(dir, sediment_files)
shp_list <- lapply(complete_paths, st_read) # load shp into a list

par(mfrow = c(4, 4), mar = c(3, 4, 3, 2)) # margins

for (i in seq_along(shp_list)) {
  file_name <- basename(complete_paths[i])
  ma_value <- as.numeric(sub(".*_(\\d+)Ma.*", "\\1", file_name)) # name according to age (Ma)
  bbox <- st_bbox(shp_list[[i]])
  xlim <- c(-180, 180)
  ylim <- c(-90, 90)    
  
  plot(st_geometry(shp_list[[i]]), 
       main = paste(ma_value, "Ma"), 
       col = "black", 
       xlim = xlim, ylim = ylim, 
       axes = FALSE, asp = 1)
  
  axis(1, at = seq(-150, 150, by = 50), labels = seq(-150, 150, by = 50), las = 1)
  axis(2, at = c(-50, 0, 50), labels = c("-50", "0", "50"), las = 2) 
  box()
}


# plot T-P diagram with world, sediments and fossils data
tiempos

par(mfrow = c(3, 5))
for (i in 1:14){
  plot (world_temp[[i]], world_prec[[i]], xlim=c(-40, 40), ylim=c(0,6000), main=ma[[i]])
  points (allextract_temp[[i]] [,1], allextract_prec[[i]][,1], col=2, pch=16)
  #plot(world_temp[[i]], world_prec[[i]], add=T)
  points(fossil_data$temp [fossil_data$ma == tiempos [i]], fossil_data$prec[fossil_data$ma == tiempos [i]], col=3, pch=16)
}


# boxplot temp all periods
# TODO add legend to this plot
par(mfrow = c(3, 5))
for (i in 1:14){
  boxplot(as.vector (world_temp[[i]]), allextract_temp[[i]][,1],        
           fossil_data$temp [fossil_data$ma == tiempos [i]], 
           main=ma[[i]], na.rm=T, ylim=c(-40, 40), names = c("Wrld", "Sedi", "Foss"))
}


# boxplot prec all periods
# TODO add legend to this plot
par(mfrow = c(3, 5))
for (i in 1:14){
  boxplot (as.vector (world_prec[[i]]), allextract_prec[[i]][,1],
           fossil_data$prec [fossil_data$ma == tiempos [i]], 
           main=ma[[i]], na.rm=T, ylim=c(0, 1500), names = c("Wrld", "Sedi", "Foss"))
}

# evolution of boxplot medians over time
par(mfrow = c(2, 1))

# TEMPERATURE MEDIANS
# vector for storing medians
median_allextract <- numeric(14)
median_fossil <- numeric(14)

# calculate medians for each period
for (i in 1:14) {
  median_allextract[i] <- median(allextract_temp[[i]][,1], na.rm = TRUE)
  median_fossil[i] <- median(fossil_data$temp[fossil_data$ma == tiempos[i]], na.rm = TRUE)
}

# plot
plot(tiempos, median_allextract, type = "l", col = "blue", lty = 1, ylim = range(c(median_allextract, median_fossil), na.rm = TRUE),
     ylab = "Temperature median", xlab = "Tiempo (Ma)", axes = FALSE)
lines(tiempos, median_fossil, col = "red", pch = 16, lty = 1)

axis(1) # add axis separately
axis(2)

abline(h = 0, col = "black") # add line to value 0

legend("bottomright", legend = c("allextract_temp", "fossil_data"), col = c("blue", "red"), lty = 1, bty = "n")

# PRECIPITATION MEDIANS
# vector for storing medians
median_allextract_prec <- numeric(14)
median_fossil_prec <- numeric(14)

# calculate medians for each period
for (i in 1:14) {
  median_allextract_prec[i] <- median(allextract_prec[[i]][,1], na.rm = TRUE)
  median_fossil_prec[i] <- median(fossil_data$prec[fossil_data$ma == tiempos[i]], na.rm = TRUE)
}

# plot
plot(tiempos, median_allextract_prec, type = "l", col = "blue", lty = 1, ylim = range(c(median_allextract_prec, median_fossil_prec), na.rm = TRUE),
     ylab = "Precipitation median", xlab = "Tiempo (Ma)", axes = FALSE)
lines(tiempos, median_fossil_prec, col = "red", lty = 1)

axis(1) # add axis separately
axis(2)

abline(h = 0, col = "black") # add line to value 0

legend("bottomright", legend = c("allextract_prec", "fossil_data"), col = c("blue", "red"), lty = 1, bty = "n")

########################
#BIOMES RECLASSIFICATION

# apply the kg_class function to all periods
biomas<-list ()
for (i in 1:length(ma)){
  t<- as.array(rast(alldata_temp[[i]]))
  p<- as.array(rast(alldata_prec[[i]]))
  biomas [[i]]<-kg_reclass (Temp=t, Prec=p, type = "broadclass")
}

# create a stack of biomes for all periods (projected in the WGS84)
biomes_stack <- NULL
for (i in 1:length(ma)){
  biomes <- rast (biomas [[i]], extent=ext(-181.875, 178.125, -91.25, 91.25), 
                  crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
  biomes_stack <- c (biomes_stack, biomes)
}
# warning because the first raster is empty and ignored

# assign names to each layer of the stack according to age (in millions of years) and plot
names (biomes_stack) <- ma

# plot biomes_stack
biome_colors <- c("#55A868", "#E2A76F", "#4C72B0", "#8172B3", "#CCCCCC")
biome_values <- c(1, 2, 3, 4, 5)
biome_labels <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")

# create an auxiliary SpatRaster containing all possible biome values
template <- biomes_stack[[1]]  # take the first raster as a template
template[] <- NA  # fill with NA

biomes_stack <- rast(biomes_stack)
biomes_stack_fixed <- rast(biomes_stack)
for (i in 1:nlyr(biomes_stack)) {
  layer <- biomes_stack[[i]]
  for (value in biome_values) {
    if (!value %in% values(layer, na.rm = TRUE)) {
      layer[1] <- value
    }
  }
  biomes_stack_fixed[[i]] <- layer
}

par(mfrow = c(4, 4), mar = c(2, 2, 3, 1)/2) 

for (i in 1:nlyr(biomes_stack_fixed)) {
  plot(biomes_stack_fixed[[i]], col = biome_colors, legend = FALSE, 
       main = paste(names(biomes_stack_fixed)[i], "Ma"))
}
# plot an empty box at position 16 for legend
plot.new()
legend("center", legend = biome_labels, fill = biome_colors, bty = "n", cex = 1.2)


# Rasterise the reconstructed sediment shapefiles for each period
sed_ras <- rast()
for (i in 1:15) {
  shapefile <- read_sf(file.path("possible_fossil_reconstructed_dissolved", paste0("possible_fossil_recons_", ma [i], "Ma_dissolved.shp")))
  ras <- rasterize (x = shapefile, y = biomes_stack[[1]])
  sed_ras <- c (sed_ras, ras)
}
# warning because the first raster is empty and ignored


########################
#AREA TOTAL BIOMES VS SEDIMENTS BIOMES

# Define the path to the directory containing the reconstructed sediment data.
setwd("possible_fossil_reconstructed_dissolved")

# Calculate the area of the cells in 100 million km2
area_rast <- cellSize(biomes_stack[[i]]) / 100000000
par(mfrow = c(1, 1))
plot(area_rast) 

# Initialise a matrix to store results from the biome and sediment area (of 15 columns for the 15 periods with 15 rows in which only 10 will be filled, 5 for the biomes without sediment filter 1-5, and 5 for the biomes with sediment filter 11-15)
res <- matrix(0, 15, 15)

# Loop to quantify the area of sediments and biomes in different time periods
for (i in 2:15) {
  shapefile <- read_sf(paste0("possible_fossil_recons_", ma[i], "Ma_dissolved.shp"))
  sed_raster <- rasterize(shapefile, area_rast)
  sed_raster[is.na(sed_raster[])] <- 0
  sed_raster <- sed_raster * 10
  
  # Sum the biome and sediment rasters and convert them into a dataframe
  all_sed_biomes <- biomes_stack[[i]] + sed_raster
  all_sed_biomes_area <- c(all_sed_biomes, area_rast)
  res_area <- as.data.frame(all_sed_biomes_area)
  res_area <- res_area[complete.cases(res_area), 1:2]
  
  # Aggregate the total area of each biome and store it in the results matrix
  res_final <- aggregate((res_area$area) ~ res_area[, 1], FUN = sum)
  res[res_final[, 1], i] <- res_final[, 2]
}

# Initialise lists to store areas of different biomes in different time periods.
tropical <- c()
arid <- c()
temperate <- c()
cold <- c()
polar <- c()

# Loop to calculate the area of each type of biome
for (i in 2:15) {
  tropicals <- res[i, 1] + res[i, 6]
  arids <- res[i, 2] + res[i, 7]
  temperates <- res[i, 3] + res[i, 8]
  colds <- res[i, 4] + res[i, 9]
  
  # Check if index 10 is not NA to calculate polar biomes
  if (!is.na(res[i, 10])) {
    polars <- res[i, 5] + res[i, 10]
  }
  
  # We store the results in the respective lists
  tropical <- c(tropicals, tropical)
  arid <- c(arids, arid)
  temperate <- c(temperates, temperate)
  cold <- c(colds, cold)
  polar <- c(polars, polar)
}

# Define the names of the biomes for visualisation
biome_names <- c("tropical", "arid", "temperate", "cold", "polar")

# Calculate the sediments area (already stored in res[11:15, ])
total_area <- res[1:5, ]
sediments_area <- res[11:15, ]

# Plot
par(mfrow = c(1, 5), mai = c(1, 1, 2, 0)/2)

# Loop to create bar charts for each biome
for (i in 1:5) {
  # Bar chart of total area
  barplot(total_area[i, ],
          names.arg = ma,  # Time period labels on the x-axis
          col = "gray",
          border = NA,
          main = paste(biome_names[i], "biome"),
          xlab = "Time (Ma)",
          ylab = "Area (thousands of millions of km²)",
          ylim = c(0, 600000)
  )
  
  # Bar chart of the sediment area
  barplot(
    sediments_area[i, ],
    names.arg = ma,
    col = "red",
    border = NA,
    add = TRUE,
    ylim = c(0, 600000)
  )
  
  # Add a legend only in the last graph (when i == 5)
  if (i == 5) {
    legend(
      "topright",
      legend = c("Total Area", "Sediment Area"),
      fill = c("gray", "red")
    )
  }
}

########################
#% OF NON SAMPLED BIOME AREA

# Biome, period and colour labels
biome_names <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
ma <- c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
biome_colors <- c("#55A868", "#E2A76F", "#4C72B0", "#8172B3", "#CCCCCC")

# Calculate percentage of unsampled area
percentage_unsampled <- ((total_area - sediments_area) / total_area) * 100

# Replace negative values with NA (optional, to handle invalid data)
percentage_unsampled[percentage_unsampled < 0] <- NA
percentage_unsampled[is.nan(percentage_unsampled)] <- 0
percentage_unsampled[is.na(percentage_unsampled)] <- 0

# Plot
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
plot(
  ma, percentage_unsampled[1, ], type = "l", col = biome_colors[1], lwd = 2,
  xlab = "Time (Ma)", ylab = "% of unsampled biome area",
  ylim = c(60,100), xaxt = "n"
)
axis(1, at = ma, labels = ma) # Customise X-axis with ma labels

# Add lines for the other biomes
for (i in 2:nrow(percentage_unsampled)) {
  lines(ma, percentage_unsampled[i, ], col = biome_colors[i], lwd = 2)
}

# Add legend
legend(
  "bottomright", legend = biome_names, col = biome_colors,
  lwd = 2, bty = "n"
)

