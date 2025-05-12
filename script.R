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
library(tidyr)
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
 
setwd("./data/fossils")
# fossil_data <- read.csv("pbdb_data_27-02-25processed.csv")
# 
# str(fossil_data)
# head(fossil_data)
# #summary(fossil_data)
# 
# # define the period limits (to reclassify the fossils within our working periods)
# periods <- setNames(lapply(ma, function(x) {
#   if (x==0){c(0,0.01)}
#   else if (x==3){c(x - 0.2, x + 0.2)}
#   else { c(x - 1, x + 1)}
# }), as.character(ma))
# 
# periods <- list(
#   "0" = c(0, 0.01),          
#   "3" = c(2.8, 3.2),         
#   "11" = c(10, 12),          
#   "15" = c(14, 16),          
#   "20" = c(19, 21),          
#   "26" = c(25, 27),          
#   "31" = c(30, 32),          
#   "36" = c(35, 37),          
#   "40" = c(39, 41),          
#   "45" = c(44, 46),          
#   "52" = c(51, 53),          
#   "56" = c(55, 57),          
#   "61" = c(60, 62),          
#   "66" = c(65, 67),          
#   "69" = c(68, 70))
# 
# # reclassify the fossils to the period (ma) that corresponds to our conditions
# 
# fossil_data$ma <- NA
# for (i in 1:15) {
#   fossil_data$ma[fossil_data$MIN_AGE >= min(periods[[i]]) & fossil_data$MAX_AGE <= max(periods[[i]])] <- names(periods[i])
#   fossil_data$ma[fossil_data$MIN_AGE <= max(periods[[i]]) & fossil_data$MAX_AGE >= min(periods[[i]])] <- names(periods[i])
# }
# 
# 
# # convert fossil to work
# fossil_data<- as.data.frame (fossil_data)
# str (fossil_data)
# fossil_data$ma<- as.numeric (as.character (fossil_data$ma))
# fossil_data_old<- fossil_data
# fossil_data<- fossil_data_old [complete.cases (fossil_data_old), ]
# 
# # apply palaeorotations (rgplates) to all fossil occurances (palaeolat and palaeolon)
# fossil_data$paleolong<- NA
# fossil_data$paleolat<- NA
# names (fossil_data)
# 
# 
# for (i in sort (unique (fossil_data$ma))){
#   fossil_data [fossil_data$ma==i, 8:9]<- reconstruct(fossil_data [fossil_data$ma==i, 4:3], 
#                                                       age = i, model = "PALEOMAP")
# }
# 
# write.csv(fossil_data, "fossil_data_rotated.csv", row.names=TRUE)

# Upload CSV file
fossil_data <- read.csv("fossil_data_rotated.csv", row.names = 1)
fossil_data <- read.csv("../fossils/fossil_data_rotated.csv", row.names = 1)


# add prec and temp data to each fossil as a function of paleolat and paleolon
sort (unique (fossil_data$ma))
# world_temp 
fossil_data$prec<- NA
fossil_data$temp<- NA
tiempos <- sort (unique (fossil_data$ma))


for (i in 1:length(tiempos)){
  # create a mask for time
  tiempos_mask <- fossil_data$ma==tiempos [i]
  fossil_data$temp [tiempos_mask]<- terra::extract (world_temp [[i]], fossil_data[tiempos_mask, 8:9 ])[, 2]
  fossil_data$prec [tiempos_mask]<- terra::extract (world_prec[[i]], fossil_data[tiempos_mask, 8:9 ])[, 2]
}

########################
# PLOTS

# 1) plot the sediments for the 15 periods
setwd("../possible_fossil_reconstructed_dissolved")
dir <- "./"
complete_paths <- file.path(dir, sediment_files)
shp_list <- lapply(complete_paths, st_read) # load shp into a list

par(mfrow = c(4, 4), mar = c(3, 4, 3, 2)) # margins

# extract and sort ma values
ma_values <- as.numeric(sub(".*_(\\d+)Ma.*", "\\1", basename(complete_paths)))
order_indices <- order(ma_values)

for (i in order_indices) {
  file_name <- basename(complete_paths[i])
  ma_value <- ma_values[i]
  bbox <- st_bbox(shp_list[[i]])
  xlim <- c(-150, 150)
  ylim <- c(-90, 90)
  
  plot(st_geometry(shp_list[[i]]),
       main = paste(ma_value, "Ma"),
       col = "black",
       xlim = xlim, ylim = ylim,
       axes = FALSE)
  
  axis(1, at = seq(-150, 150, by = 50), labels = seq(-150, 150, by = 50), las = 1)
  axis(2, at = c(-50, 0, 50), labels = c("-50", "0", "50"), las = 2)
  box()
}

# 2) plot T-P diagram with world, sediments and fossils data
par(mfrow = c(3, 5))
for (i in 1:14){
  plot (world_temp[[i]], world_prec[[i]], xlim=c(-40, 40), ylim=c(0,6000), main=ma[[i]])
  points (allextract_temp[[i]] [,1], allextract_prec[[i]][,1], col=2, pch=16)
  #plot(world_temp[[i]], world_prec[[i]], add=T)
  points(fossil_data$temp [fossil_data$ma == tiempos [i]], fossil_data$prec[fossil_data$ma == tiempos [i]], col=3, pch=16)
}


# 3) boxplots temp 
# 3.1) all periods (for separated)
par(mfrow = c(3, 5))
for (i in 1:14){
  boxplot(as.vector (world_temp[[i]]), allextract_temp[[i]][,1],        
           fossil_data$temp [fossil_data$ma == tiempos [i]], 
           main=ma[[i]], na.rm=T, ylim=c(-40, 40), names = c("Wrld", "Sedi", "Foss"))
}

# 3.2) all periods in the same timeline (separated)
#WORLD CLIMATE
par(mfrow = c(1, 1))
clima_model_temp<- values (world_temp)
colnames (clima_model_temp)<- tiempos
clima_model_temp <- clima_model_temp[, rev(order(as.numeric(colnames(clima_model_temp))))]
#boxplot (clima_model_temp[,2:15], col="lightgrey", border="black")

#SEDIMENTS
sedi<- allextract_temp
clima_sed_temp<- list ()
for (i in 1:15){
  clima_sed_temp [[i]]<- sedi[[i]][,1] 
}
clima_sed_temp<- as.data.frame(do.call(cbind, clima_sed_temp)) # warning message
names (clima_sed_temp)<- tiempos
clima_sed_temp <- clima_sed_temp[rev(order(as.numeric(names(clima_sed_temp))))]
#boxplot (clima_sed_temp[2:15], col="lightgrey", border="black")

#FOSSILS
fossil_data_temp<-data.frame(ma=fossil_data$ma, temp=fossil_data$temp)
fossil_data_temp<- subset(fossil_data_temp, ma <= 66) #without 69 Ma data
fossil_data_temp$ma <- factor(fossil_data_temp$ma, levels = rev(sort(unique(fossil_data_temp$ma)))) #order ma in descending order
#boxplot(temp ~ ma, data = fossil_data_temp, col="lightgrey", border="black")

#Plot
par(mfrow = c(3, 1))
boxplot (clima_model_temp[,2:15], col="lightgrey", ylim=c(-40, 45), xlab = "Time (Ma)", ylab = "Temperature (ºC)", main="World climate", outline=FALSE)
boxplot (clima_sed_temp[2:15], col="#332288", ylim=c(-40, 45), xlab = "Time (Ma)", ylab = "Temperature (ºC)", main="Sediments", outline=FALSE)
boxplot(temp ~ ma, data = fossil_data_temp, col="#AA4499", ylim=c(-40, 45), xlab = "Time (Ma)", ylab = "Temperature (ºC)", main="Fossils", outline=FALSE)

# 3.3) all periods in the same timeline (together)
#WORLD CLIMATE
clima_model_temp2 <- as.data.frame(clima_model_temp[, 2:15])

cl_mo_t <- data.frame()
for (i in colnames(clima_model_temp2)) {
  v <- data.frame(value = clima_model_temp2[, i], year = i, source = "model")
  cl_mo_t <- rbind(cl_mo_t, v)
}
#str(cl_mo_t)

#SEDIMENTS
clima_sed_temp2 <- clima_sed_temp[, 2:15]

cl_se_t <- data.frame()
for (i in colnames(clima_sed_temp2)) {
  v <- data.frame(value = clima_sed_temp2[, i], year = i, source = "sediment")
  cl_se_t <- rbind(cl_se_t, v)
}
#str(cl_se_t)

# FOSSILS
df_fossil <- fossil_data_temp     
ma <- sort(unique(df_fossil$ma))  

le_fossil <- c()
for (i in ma) {
  da <- df_fossil[which (df_fossil$ma == i), ]
  le_fossil <- c(le_fossil, nrow(da))
}

# Create a data frame with a column 'n' ranging from 1 to max occurrences
dat_fossil <- data.frame(n = 1:max(le_fossil))

for (i in ma) {
  da <- df_fossil[which (df_fossil$ma == i), 2]
  length(da) <- max(le_fossil)
  dat_fossil <- cbind (dat_fossil, da)
}

fossil_data_temp2 <- dat_fossil[, -1]
colnames(fossil_data_temp2) <- ma
#str(fossil_data_prec2)
#dim(fossil_data_prec2)

fo_t <- data.frame()
for (i in colnames(fossil_data_temp2)) {
  v <- data.frame(value = fossil_data_temp2[, i], year = rep(i, nrow(fossil_data_temp2)), source = rep ("fossil", nrow (fossil_data_temp2)))
  fo_t <- rbind(fo_t, v)
}
#str(cl_fo_pr)


#Put the 3 data together
all_temp2 <- rbind (cl_mo_t, cl_se_t, fo_t)
#str(all_prec2)
all_temp2$year <- factor(all_temp2$year, levels = unique (all_temp2$year))
all_temp2$source <- factor(all_temp2$source, levels = c ("model", "sediment", "fossil"))

# vector for plot spacing
v <- c()
for (i in seq (1, 54, 4)) {
  v <- c (v, i:(i+2))
}

# vector for ma (axis titles)
vl <- c()
for (i in rev(ma[1:14])) {
  j <- c("", i, "")
  vl <- c(vl, j )
}

#Plot
cols = c ("model" = "lightgrey", "sediment" = "#332288" , "fossil" = "#AA4499")

par(mfrow = c(1, 1))
boxplot (value ~ source + year, data=all_temp2, at = v,
         names = rev(vl),col=cols, ylim=c(-45, 45), xlab = "Time (Ma)", ylab = "Temperature (ºC)", main="ALL", outline=FALSE)
legend("bottomleft", fill = cols, legend = c("model", "sediment", "fossil"), horiz = T, bty = "n")

# 4) boxplots prec 
# 4.1) all periods (for separated)
par(mfrow = c(3, 5))
for (i in 1:14){
  boxplot (as.vector (world_prec[[i]]), allextract_prec[[i]][,1],
           fossil_data$prec [fossil_data$ma == tiempos [i]], 
           main=ma[[i]], na.rm=T, ylim=c(0, 1500), names = c("Wrld", "Sedi", "Foss"))
}

# 4.2) all periods in the same timeline
#WORLD CLIMATE
clima_model_prec<- values (world_prec)
colnames (clima_model_prec)<- tiempos
clima_model_prec <- clima_model_prec[, rev(order(as.numeric(colnames(clima_model_prec))))]
#boxplot (clima_model_prec, col="lightgrey", border="black")

#SEDIMENTS
preci<- allextract_prec
clima_sed_prec<- list ()
for (i in 1:15){
  clima_sed_prec [[i]]<- preci[[i]][,1] 
}
clima_sed_prec<- as.data.frame(do.call(cbind, clima_sed_prec)) #warning message
names (clima_sed_prec)<- tiempos
clima_sed_prec <- clima_sed_prec[rev(order(as.numeric(names(clima_sed_prec))))]
#boxplot (clima_sed_prec[2:15], col="lightgrey", border="black")

#FOSSILS
fossil_data_prec<-data.frame(ma=fossil_data$ma, prec=fossil_data$prec)
fossil_data_prec<- subset(fossil_data_prec, ma <= 66) #without 69 Ma data
fossil_data_prec$ma <- factor(fossil_data_prec$ma, levels = rev(sort(unique(fossil_data_prec$ma)))) #order ma in descending order
#boxplot(prec ~ ma, data = fossil_data_prec, col="lightgrey", border="black")

#plot
par(mfrow = c(3, 1))
boxplot (clima_model_prec[,2:15],col="lightgrey", ylim=c(0, 3000), xlab = "Time (Ma)", ylab = "Precipitation (mm/year)", main="World climate", outline=FALSE) 
boxplot (clima_sed_prec[2:15],col="#332288", ylim=c(0, 3000), xlab = "Time (Ma)", ylab = "Precipitation (mm/year)", main="Sediments", outline=FALSE) 
boxplot (prec~ma, data=fossil_data_prec, col= "#AA4499", ylim=c(0, 3000), xlab = "Time (Ma)", ylab = "Precipitation (mm/year)", main="Fossils", outline=FALSE)

# 4.3) all periods in the same timeline (together)
#WORLD CLIMATE
clima_model_prec2 <- as.data.frame(clima_model_prec[, 2:15])

cl_mo_pr <- data.frame()
for (i in colnames(clima_model_prec2)) {
  v <- data.frame(value = clima_model_prec2[, i], year = rep(i, nrow(clima_model_prec2)), source = rep ("model", nrow (clima_model_prec2)))
  cl_mo_pr <- rbind(cl_mo_pr, v)
}
#str(cl_mo_pr)

#SEDIMENTS
clima_sed_prec2 <- clima_sed_prec[, 2:15]

cl_se_pr <- data.frame()
for (i in colnames(clima_sed_prec2)) {
  v <- data.frame(value = clima_sed_prec2[, i], year = rep(i, nrow(clima_sed_prec2)), source = rep ("sediment", nrow (clima_sed_prec2)))
  cl_se_pr <- rbind(cl_se_pr, v)
}
#str(cl_se_pr)

# FOSSILS
df_fossil <- fossil_data_prec      
ma <- sort(unique(df_fossil$ma))  

le_fossil <- c()
for (i in ma) {
  da <- df_fossil[which (df_fossil$ma == i), ]
  le_fossil <- c(le_fossil, nrow(da))
}

dat_fossil <- data.frame(n = 1:max(le_fossil))
for (i in ma) {
  da <- df_fossil[which (df_fossil$ma == i), 2]
  length(da) <- max(le_fossil)
  dat_fossil <- cbind (dat_fossil, da)
}
fossil_data_prec2 <- dat_fossil[, -1]
colnames(fossil_data_prec2) <- ma
#str(fossil_data_prec2)
#dim(fossil_data_prec2)

cl_fo_pr <- data.frame()
for (i in colnames(fossil_data_prec2)) {
  v <- data.frame(value = fossil_data_prec2[, i], year = i, source = "fossil")
  cl_fo_pr <- rbind(cl_fo_pr, v)
}
#str(cl_fo_pr)


#Put the 3 data together
all_prec2 <- rbind (cl_mo_pr, cl_se_pr, cl_fo_pr)
#str(all_prec2)
all_prec2$year <- factor(all_prec2$year, levels = unique (all_prec2$year))
all_prec2$source <- factor(all_prec2$source, levels = c ("model", "sediment", "fossil"))

#Plot
cols = c ("model" = "lightgrey", "sediment" = "#332288" , "fossil" = "#AA4499")

v <- c()
for (i in seq (1, 54, 4)) {
  v <- c (v, i:(i+2))
}
v

vl <- c()
for (i in rev(ma[1:14])) {
  j <- c("", i, "")
  vl <- c(vl, j )
}
vl

par(mfrow = c(1, 1))
boxplot (value ~ source + year, data=all_prec2, 
         at = v,
         names = rev(vl), 
         xaxs = FALSE,
         col=cols, ylim=c(0, 3000), xlab = "Time (Ma)", ylab = "Precipitation (mm/year)", main="ALL", outline=FALSE)
legend("top", fill = cols, legend = c("model", "sediment", "fossil"), horiz = T, bty = "n")

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
names (biomes_stack) <- rev(ma)

biome_colors <- c("#117733", "#DDCC77", "#CC6677", "#88CCEE", "#888888")
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

# 1) plot biomes_stack
par(mfrow = c(4, 4)) 

for (i in nlyr(biomes_stack_fixed):1) {
  plot(biomes_stack_fixed[[i]], col = biome_colors, legend = FALSE, 
       main = paste(names(biomes_stack_fixed)[i], "Ma"))
}
# plot an empty box at position 16 for legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = biome_labels, fill = biome_colors, bty = "n", cex = 1.2)


# Rasterise the reconstructed sediment shapefiles for each period
setwd("../possible_fossil_reconstructed_dissolved")
## tiempos is the vector with all the periods (69 ma included)
sed_ras <- rast()
for (i in 1:15) {
  shapefile <- read_sf(paste0("possible_fossil_recons_", tiempos [i], "Ma_dissolved.shp"))
  ras <- rasterize (x = shapefile, y = biomes_stack[[1]])
  sed_ras <- c (sed_ras, ras)
}
# warning because the first raster is empty and ignored

# 2) plot 3 maps (1 for each layer used: sediments, biomes, fossils) for 1 example time (45 Ma) for material and methods
par(mfrow = c(3, 3)) 

#BIOMES
plot(biomes_stack_fixed[[10]], 
     col = biome_colors, legend = FALSE, 
     main = "World climatic zones")
# plot an empty box at position 16 for legend
#plot.new()
#legend("center", legend = biome_labels, fill = biome_colors, bty = "n", cex = 1.2)

#SEDIMENTS
plot(biomes_stack_fixed[[10]], main = "Sediments",
     col = "grey", legend = FALSE) #grey raster as a basis
plot(st_geometry(shp_list[[10]]), 
     col = "#332288", border = "#332288", add=T)

#FOSSILS
plot(biomes_stack_fixed[[10]], main = "Fossil data",
     col = "grey", legend = FALSE) #grey raster as a basis
#plot(st_geometry(shp_list[[10]]), col = "black", add=T) #sediment cover
# Add fossils as red points
fossil_filtered <- subset(fossil_data, ma == 45)
points(fossil_filtered$paleolong, fossil_filtered$paleolat, col =  "#AA4499", pch = 19, cex = 0.2)

#######################
#AREA TOTAL BIOMES VS SEDIMENTS BIOMES

# Define the path to the directory containing the reconstructed sediment data.
setwd("../possible_fossil_reconstructed_dissolved")

# Calculate the area of the cells in 100 million km2
area_rast <- cellSize(biomes_stack[[1]]) / 100000000
par(mfrow = c(1, 1))
plot(area_rast) 

# Initialise a matrix to store results from the biome and sediment area (of 15 columns for the 15 periods with 15 rows in which only 10 will be filled, 5 for the biomes without sediment filter 1-5, and 5 for the biomes with sediment filter 11-15)
# res <- matrix(0, 15, 15)
res <- matrix(0, 15, 14)

# Loop to quantify the area of sediments and biomes in different time periods
for (i in 1:14) {
  shapefile <- read_sf(paste0("possible_fossil_recons_", tiempos[i], "Ma_dissolved.shp"))
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

# Initialise vectors to store areas of different biomes in different time periods.
tropical <- c()
arid <- c()
temperate <- c()
cold <- c()
polar <- c()

# Loop to calculate the area of each type of biome
for (i in 1:14) {
  tropicals <- res[1, i] + res[11, i]
  arids <- res[2, i] + res[12, i]
  temperates <- res[3, i] + res[13, i]
  colds <- res[4, i] + res[14, i]
  
  # Check if index 10 is not NA to calculate polar biomes
  if (!is.na(res[5, i])) {
    polars <- res[5, i] + res[15, i]
  } else {
    polars <-0
  }
  
  # We store the results in the respective vector
  tropical <- c(tropical, tropicals)
  arid <- c(arid, arids)
  temperate <- c(temperate, temperates)
  cold <- c(cold, colds)
  polar <- c(polar, polars)
}

# Define the names of the biomes for visualisation
biome_names <- c("tropical", "arid", "temperate", "cold", "polar")

# Calculate the sediments area (already stored in res[11:15, ])
total_area <- res[1:5, ] + res[11:15, ]
sediments_area <- res[11:15, ]

# Plot
par(mfrow = c(1, 5), mar = c(5, 5, 4, 2) )

# Loop to create bar charts for each biome
for (i in 1:5) {
  # Bar chart of total area
  barplot(rev(total_area[i, 1:14]),
          names.arg = ma[1:14],  # Time period labels on the x-axis
          col = "gray",
          border = NA,
          main = paste(biome_names[i], "biome"),
          xlab = "Time (Ma)",
          ylab = "Area (thousands of millions of km²)",
          ylim = c(0, 600000)
  )
  
  # Bar chart of the sediment area
  barplot(
    rev(sediments_area[i, 1:14]),
    names.arg = ma[1:14],
    col = "#332288",
    border = NA,
    add = TRUE,
    ylim = c(0, 600000)
  )
  
  # Add a legend only in the last graph (when i == 5)
  if (i == 5) {
    legend(
      "topright",
      legend = c("Total Area", "Sediment Area"),
      fill = c("gray", "#332288")
    )
  }
}

percentage_sediments <- (sediments_area / total_area) * 100
row_means_sediments <- rowMeans(percentage_sediments, na.rm = TRUE)

########################
#% OF NON SAMPLED BIOME AREA

# Biome, period and colour labels
biome_names <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
# ma <- c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
ma <- c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66)
biome_colors <- c("#117733", "#DDCC77", "#CC6677", "#88CCEE", "#888888")

# Calculate percentage of unsampled area
percentage_unsampled <- ((total_area - sediments_area) / total_area) * 100

# Replace NAN with NA
percentage_unsampled[is.nan(percentage_unsampled)] <- NA

# OPTION 1: Plot with polar line
par(mfrow = c(1, 1))
matplot(
  rev(ma), t(percentage_unsampled[, ncol(percentage_unsampled):1]), type = "l", 
  lty = 1, lwd = 2, col = biome_colors, 
  xlab = "Time (Ma)", ylab = "% of unsampled biome area",
  ylim = c(35, 100), xaxt = "n", xlim = rev(range(ma)) # Force x-axis inversion
)

# Add inverted x-axis labels
axis(1, at = rev(ma), labels = rev(ma))

# Add legend
legend(
  "bottomleft", legend = biome_names, col = biome_colors,
  lwd = 2, bty = "n"
)

# OPTION 2: Plot without polar line, for better visualization of lines for the other biomes by removing the last row (i.e. polar data) of the df
percentage_unsampled2 <- percentage_unsampled[-5,]
par(mfrow = c(1, 1))
matplot(
  rev(ma), t(percentage_unsampled2[, ncol(percentage_unsampled2):1]), type = "l", 
  lty = 1, lwd = 2, col = biome_colors, 
  xlab = "Time (Ma)", ylab = "% of unsampled biome area",
  ylim = c(70, 100), xaxt = "n", xlim = rev(range(ma)) # Force x-axis inversion
)

# Add inverted x-axis labels
axis(1, at = rev(ma), labels = rev(ma))

# Add legend
legend(
  "bottomleft", legend = biome_names, col = biome_colors,
  lwd = 2, bty = "n"
)

########################
#ACCUMULATED BARPLOTS OF %
biome_names <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
biome_colors <- c("#117733", "#DDCC77", "#CC6677", "#88CCEE", "#888888")

#climate distribution across time (%)
# total_area_climate <- res[1:5, ] + res[11:15, ]
total_area_climate <- total_area
sum_area_climate<- colSums(total_area_climate)
perc_total_climate<- round (sweep(total_area_climate, 2, sum_area_climate, `/`)*100, 2)

colSums(perc_total_climate)
colnames (perc_total_climate)<- tiempos[1:14]
row.names(perc_total_climate)<- c(biome_names)
perc_total_climate_rev <- perc_total_climate[, rev(colnames(perc_total_climate))] # Reverse x-axis

par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
barplot(perc_total_climate_rev, col=biome_colors, border="white", xlab="Time", ylab="%", main="World climate")
legend("topright", legend=biome_names, fill=biome_colors, border="black", cex=0.8)

#SEDIMENTS distribution across time (%)
# total_area_sediments <- res[11:15, ]
total_area_sediments <-sediments_area 
sum_area_sediments<- colSums(total_area_sediments)
perc_total_sediments<- round (sweep(total_area_sediments, 2, sum_area_sediments, `/`)*100, 2)

colnames (perc_total_sediments)<- tiempos[1:14]
row.names(perc_total_sediments)<- c(biome_names)
perc_total_sediments_rev <- perc_total_sediments[, rev(colnames(perc_total_sediments))] # Reverse x-axis

barplot(perc_total_sediments_rev, col=biome_colors, border="white", xlab="Time", ylab="%", main="Sediments")
legend("topright", legend=biome_names, fill=biome_colors, border="black", cex=0.8)

#FOSSIL RECORD distribution across time (%)
for (i in 1:14) {
  # Create a mask for the current time period
  tiempos_mask <- fossil_data$ma == tiempos[i]
  fossil_data$biomes[tiempos_mask] <- terra::extract(biomes_stack_fixed[[i]], fossil_data[tiempos_mask, 8:9])[, 2]
}

# Create df with time period (ma) and biomes
fossils <- data.frame(ma = fossil_data$ma, biomes = fossil_data$biomes)

# Remove rows with missing values
fossils <- fossils[complete.cases(fossils), ]

# Count the number of fossils per biome and time period
df_fossils <- fossils %>%
  group_by(biomes, ma) %>%
  summarise(n())

# Rename the third column to "n"
colnames(df_fossils)[3] <- "n"

# Convert the data to a wide format (time periods as columns)
df_wide <- df_fossils %>%
  ungroup() %>% 
  pivot_wider(names_from = ma, values_from = n, values_fill = 0)  

colnames(df_wide)

# Calculate the total number of fossils per time period
sum_area_fossils <- colSums(df_wide)

# Calculate % of fossils per biome and time period
perc_total_fossils <- round(sweep(df_wide, 2, sum_area_fossils, `/`) * 100, 2)

# Exclude the 1st column (biomes)
perc_total_fossils <- as.matrix(perc_total_fossils[, 2:15])
row.names(perc_total_fossils)<- c(biome_names)
perc_total_fossils_rev <- perc_total_fossils[, rev(colnames(perc_total_fossils))] # Reverse x-axis

barplot(perc_total_fossils_rev, col=biome_colors, border="white", xlab="Time", ylab="%", main="Fossil data")
legend("topright", legend=biome_names, fill=biome_colors, border="black", cex=0.8)

########################
# DELTA PLOT OF %

# % SEDIMENTS VS % WORLD_CLIMATE
# Calculate delta
delta_climased<- round ((perc_total_sediments_rev - perc_total_climate_rev)/perc_total_climate_rev *100, 2)

# Plot
par(mfrow = c(2, 1), mar = c(3, 4, 3, 2))
plot (colnames (delta_climased), rev(delta_climased [1,]), type="l", ylim=c(-100, 500), col="#117733", xlab="Time (Ma)", ylab="DELTA (% sediments - % climate)", xaxt="n", lwd = 3)
lines (colnames (delta_climased), rev(delta_climased [2,]), col="#DDCC77", lwd = 3)
lines (colnames (delta_climased), rev(delta_climased [3,]), col="#CC6677", lwd = 3)
lines (colnames (delta_climased), rev(delta_climased [4,]), col="#88CCEE", lwd = 3)
lines (colnames (delta_climased), rev(delta_climased [5,]), col="#888888", lwd = 3)
abline (h=0, col="black", lty=1)
axis(1, at=rev(colnames(delta_climased)), labels=rev(ma[1:14])) 
legend("top", legend=c(biome_labels), col=c("#117733", "#DDCC77", "#CC6677", "#88CCEE", "#888888"), lwd=2, bty="n", ncol= 5)

# % FOSSILS VS % WORLD_CLIMATE
# Calculate delta
delta_climafoss<- round ((perc_total_fossils_rev - perc_total_climate_rev)/perc_total_climate_rev *100, 2)

# Plot
plot (colnames (delta_climased), rev(delta_climafoss [1,]), type="l", ylim=c(-100, 500), col="#117733", xlab="Time (Ma)", ylab="DELTA (% fossils - % climate)", xaxt="n", lwd = 3)
lines (colnames (delta_climased), rev(delta_climafoss [2,]), col="#DDCC77", lwd = 3)
lines (colnames (delta_climased), rev(delta_climafoss [3,]), col="#CC6677", lwd = 3)
lines (colnames (delta_climased), rev(delta_climafoss [4,]), col="#88CCEE", lwd = 3)
lines (colnames (delta_climased), rev(delta_climafoss [5,]), col="#888888", lwd = 3)
abline (h=0, col="black", lty=1)
axis(1, at=rev(colnames(delta_climased)[2:15]), labels=rev(ma[1:14])) 
#legend("top", legend=c(biome_labels), col=c("#55A868", "#E2A76F", "#4C72B0", "#8172B3", "#CCCCCC"), lwd=2, bty="n", ncol=5)

