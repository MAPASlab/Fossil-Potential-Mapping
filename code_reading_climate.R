## METADATA ===============================================================
## Description: 
## scripts for identifying regions with high fossil preservation potential using 
## climate data (Scotese model), fossil records, and sediment zones. 
## We apply Köppen-Geiger climate classification to map biomes across 15 geological periods. 
## 
## R version: 4.2.2 for Windows
## Date: 2024-10-18 16:59:26
## License: GPL3 (Marta Please check)
## Author: Marta Matamala, Sarava Varela, Oskar Hagen.
##=======================================================================##

library(terra)
library(rgplates)
library(sf)
library (maps)
source("support_functions.R")

########################
#DATA PREPARATION
# 1) WORLD DATA

# Unzip the data file (only the first time)
untar("./data/scotese/Scotese_temp_precip_Ceno.tar.nc", exdir="./data/scotese")


### VARIABLES ------------
# years in millions
ma <- c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
# Create a list of month abbreviations
month_ID <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
dir_scot <- "./data/scotese/formatted_data"
dir_tfkea <- file.path(dir_scot, "tfkea")


# LIST AND LOAD TEMPERATURE DATA
# List temperature files
list_temp.nc <- list.files(dir_tfkea, pattern = ".temp.nc")
list_prec.nc <- list.files(dir_tfkea, pattern = ".precip.nc")

# Set the directory to access the temperature files
# setwd("C:/Users/Usuario/OneDrive - campus.udg.edu/Escriptori/Oskar_sediments/scotese/formatted_data/tfkea")

# Create a raster stack from the temperature files
stack_temp <- rotate (rast(file.path(dir_tfkea, list_temp.nc)))
stack_prec <- rotate (rast(file.path(dir_tfkea, list_prec.nc)))
names(stack_temp)<- list_temp.nc
names(stack_prec)<- list_prec.nc


# dir <- "./data/scotese/formatted_data"
# setwd(dir)
scotese_times <- data.frame(
  "name" = list.files("./data/scotese/formatted_data"),
  "ma" = c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
)



# Initialize lists for temperature and precipitation data
alldata_temp <- list()
alldata_prec <- list()



# Process data for each time period

for (i in 1:nrow(scotese_times)) {
  # i <- 1
  time_ID <- scotese_times[i, 1]
  my_result_temp <- monthly_stacks(time_ID, "1_5m_temp", dir)
  my_result_prec <- monthly_stacks(time_ID, "precip", dir)
  alldata_temp[[i]] <- my_result_temp
  alldata_prec[[i]] <- my_result_prec
}

# Name the lists
names(alldata_temp) <- paste("temp_", ma, sep = "")
names(alldata_prec) <- paste("prec_", ma, sep = "")
myears<- paste (ma, "Ma", sep="_")

# create stack of world climatic data for temp and prec with average prec and average temp
world_temp<- NULL
world_prec<- NULL
for (i in 1:length(ma)){
  # i <- 1
  world_temp<- c(world_temp, mean (alldata_temp [[i]]))
  names(world_temp)[i]<- myears [i]
  world_prec<- c(world_prec, sum(alldata_prec [[i]]))
  names(world_prec) [i] <- myears [i]
} 


#plot (world_temp, world_prec, xlim=c(-40, 40), ylim=c(0,6000), main=ma)

########################
# 2) SEDIMENTS

# setwd("C:/Users/Usuario/OneDrive - campus.udg.edu/Escriptori/Oskar_sediments/reconstructed_sediment")
sediment_files <- list.files("./data/reconstructed_sediment",  pattern = "\\.shp$")
# create stack of sediment data for temp and prec with average prec and average temp
allextract_temp <- list()
allextract_prec <- list()

ommited_sediments <- 1 # assuming ommited sediments are allways the present ones

for (i in (1:length(sediment_files))+ommited_sediments) { # ask Marta, why starting at 2 here?
  shapefile <- read_sf(file.path("./data/reconstructed_sediment", paste0("sediment_recons_", ma [i], "Ma.shp")))
  extracted_data_temp <- extract (world_temp [[i]], shapefile, xy = TRUE, ID = FALSE)
  extracted_data_prec <- extract (world_prec [[i]], shapefile, xy = TRUE, ID = FALSE)
  allextract_temp[[paste0(ma [i], "Ma")]] <- extracted_data_temp
  allextract_prec[[paste0(ma [i], "Ma")]] <- extracted_data_prec
}
# ask Marta, does this resolution makes sense? 1km²?
#plot(rast(allextract_temp[[14]][,c(2,3,1)]))

## beware that here if you put 15 layers of sediment, you have to remove the minus 1 from the "i"!! and start the loop at 1.

for (i in (1:length(sediment_files))+ommited_sediments) {
  # i <- 2
  #print(i)
  rast_temp_name <- paste0("rast", ma[i], "_temp")
  assign(rast_temp_name, rast (na.omit(allextract_temp[[i-ommited_sediments]][, c(2, 3, 1)])))
  
  rast_prec_name <- paste0("rast", ma[i], "_prec")
  assign(rast_prec_name, rast(na.omit(allextract_prec[[i-ommited_sediments]][, c(2, 3, 1)])))
}

# str (allextract_prec)
########################
# 3) FOSSILS

fossil_data <- read.csv("./data/now_database/now_database_09-09-24processed.csv")

str(fossil_data)
head(fossil_data)
summary(fossil_data)

# define the period limits (to reclassify the fossils within our working periods)
periods <- setNames(lapply(ma, function(x) {
  if (x==0){c(0,0.01)}
  else if (x==3){c(x - 0.2, x + 0.2)}
  else { c(x - 1, x + 1)}
}), as.character(ma))

# reclassify the fossils to the period (ma) that corresponds to our conditions

# This is mess restrictive...
fossil_data$ma<- NA
for (i in 1:length(ma)){
  fossil_data$ma[fossil_data$MIN_AGE>=min(periods [[i]]) & fossil_data$MAX_AGE<= max(periods[[i]])]<- names(periods [i])
  fossil_data$ma[fossil_data$MIN_AGE>=min(periods [[i]]) & fossil_data$MIN_AGE<= max(periods[[i]])]<- names(periods [i])
  fossil_data$ma[fossil_data$MAX_AGE>=min(periods [[i]]) & fossil_data$MAX_AGE<= max(periods[[i]])]<- names(periods [i])
}

# # this is more restrictive.....
# fossil_data$ma_alternative<- NA
# for (i in 1:length(ma)) {
#   condition <- fossil_data$MIN_AGE >= min(periods[[i]]) & fossil_data$MAX_AGE <= max(periods[[i]])
#   fossil_data$ma_alternative[condition] <- names(periods[i])
# }
# 

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

# plot T-P diagram with world, sediments and fossils data.
# sediment has no time 0 (i.e. 14 scenarios). world has time 0 (i.e. 15 scenarios). "i" 's are unmatched.
tiempos

par(mfrow = c(3, 5))
for (i in 2:14){
  plot (world_temp[[i]], world_prec[[i]], xlim=c(-40, 40), ylim=c(0,6000), main=ma[[i]])
  points (allextract_temp[[i-1]] [,1], allextract_prec[[i-1]][,1], col=2, pch=16)
  plot(world_temp[[i]], world_prec[[i]], add=T)
  points(fossil_data$temp [fossil_data$ma == tiempos [i]], fossil_data$prec[fossil_data$ma == tiempos [i]], col="darkblue", pch=1, cex=1)
}

# boxplot temp all periods
dev.off()
# TODO add legend to this plot
par(mfrow = c(3, 5))
for (i in 2:14){
  boxplot(as.vector (world_temp[[i]]), allextract_temp[[i-1]][,1],        
          fossil_data$temp [fossil_data$ma == tiempos [i]], 
          main=ma[[i]], na.rm=T, ylim=c(-40, 40))
  
}

# boxplot prec all periods
# TODO add legend to this plot
par(mfrow = c(3, 5))
for (i in 2:14){
  boxplot (as.vector (world_prec[[i]]), allextract_prec[[i-1]][,1],
           fossil_data$prec [fossil_data$ma == tiempos [i]], 
           main=ma[[i]], na.rm=T, ylim=c(0, 1500))
  
}

########################
#BIOMES RECLASSIFICATION

# Convert temperature and precipitation data from the first period into arrays
t0<- as.array(alldata_temp[[1]])
str (t0)
p0<- as.array(alldata_prec[[1]])

# TEST: First period climate classification using the kg_reclass function 
# Returns the biome (broad climate category) for each cell in the spatial matrix
tiempo0<- kg_reclass (Temp=t0, 
                      Prec=p0,type="broadclass")

# Apply the kg_class function to all periods
biomas<-list ()
for (i in 1:length(ma)){
  t<- as.array(alldata_temp[[i]])
  p<- as.array(alldata_prec[[i]])
  biomas [[i]]<-kg_reclass (Temp=t, Prec=p,type="broadclass")
}

# TEST: Project the first period in the WGS84 coordinate system and calculate cellsizes
r<- alldata_temp[[1]]
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
area_raster<- cellSize (r)
dev.off() # TODO, is this necessary?
plot (r [[1]])
map (add=T)
area<- values (area_raster)
area_raster

dim (area)
plot (area_raster)

# Create a stack of biomes for all periods (projected in the WGS84)
biomes_stack <- NULL
for (i in 1:length(ma)){
  biomes <- rast (biomas [[i]], extent=ext(-181.875, 178.125, -91.25, 91.25), 
                  crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
  biomes_stack <- c (biomes_stack, biomes)
}
# warning because the first raster is empty and ignored

# Assign names to each layer of the stack according to age (in millions of years) and plot
names (biomes_stack) <- ma
plot (biomes_stack)
# 15 warnings because first rasters are empty and ignored

# Rasterise the reconstructed sediment shapefiles for each period
sed_ras <- rast()
for (i in 2:15) {
  shapefile <- read_sf(file.path("./data/reconstructed_sediment", paste0("sediment_recons_", ma [i], "Ma.shp")))
  ras <- rasterize (x = shapefile, y = biomes_stack[[1]])
  sed_ras <- c (sed_ras, ras)
}
# warning because the first raster is empty and ignored


# 1)Calculate the areas occupied by each biome within the fossil sediment zones
biomes_areas <- data.frame(ma = ma [-1], tropical = NA, arid = NA, 
                           temperate = NA, cold = NA, polar = NA)
biomes_areas_percent <- data.frame(ma = ma [-1], tropical = NA, arid = NA, 
                                   temperate = NA, cold = NA, polar = NA)

for (i in 2:15) {
  tabla <- data.frame (values (biomes_stack[[i]]),
                       values(area_raster)/1000000, #Areas of each cell in millions of km²
                       values(sed_ras[[i-1]]))
  names(tabla) <- c("biome", "area", "sediment")
  
  vec <- c()
  for (j in 1:5) {
    area <- sum (tabla[which(tabla$biome == j & tabla$sediment == 1), "area"])
    vec <- c (vec, area)
  }
  
  # Save areas in the data frame
  biomes_areas [i-1, 2:6 ] <- vec
  biomes_areas_percent [i-1, 2:6 ] <- round (vec * 100 / sum (vec), 2)
  
}


# 2)Calculate the areas occupied by each biome around the world
biomes_areas_mundo <- data.frame (ma = ma , tropical = NA, arid = NA, 
                                  temperate = NA, cold = NA, polar = NA)
biomes_areas_percent_mundo <- data.frame (ma = ma, tropical = NA, arid = NA, 
                                          temperate = NA, cold = NA, polar = NA)

for (i in 1:15) {
  tabla <- data.frame (values (biomes_stack[[i]]),
                       values (area_raster)/1000000) #Areas of each cell in millions of km²
  names(tabla) <- c("biome", "area")
  
  vec <- c()
  for (j in 1:5) {
    area <- sum (tabla[which(tabla$biome == j), "area"])
    vec <- c (vec, area)
  }
  
  # Save areas in the data frame
  biomes_areas_mundo [i, 2:6 ] <- vec
  biomes_areas_percent_mundo [i, 2:6 ] <- round (vec * 100 / sum (vec), 2)
  
}


#PROPORTIONS
# Proportion of areas occupied by each biome in sediment areas to the world total
proportion_by_biome <- round (biomes_areas [,-1] * 100 / biomes_areas_mundo [-1, -1], 2)
proportion_by_biome$ma <- ma[-1]

# Plot proportions
# think of a type of graph that represents well the proportions of the 5 biomes between the areas occupied by each biome in sediment zones with respect to the world total
plot (1:5, proportion_by_biome [1,1:5])
points (1:5, proportion_by_biome [2,1:5])



