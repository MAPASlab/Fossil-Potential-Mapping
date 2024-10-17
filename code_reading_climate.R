library(terra)
library(rgplates)
library(sf)

########################
#DATA PREPARATION
# 1) WORLD DATA

# Set the working directory
#setwd("C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/scotese")
# Unzip the data file (only the first time)
#untar("Scotese_temp_precip_Ceno.tar.nc")

# LIST AND LOAD TEMPERATURE DATA
# List temperature files
list_temp.nc <- list.files("C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/scotese/formatted_data/tfkea", pattern = ".temp.nc")
list_prec.nc <- list.files("C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/scotese/formatted_data/tfkea", pattern = ".precip.nc")

# Set the directory to access the temperature files
setwd("C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/scotese/formatted_data/tfkea")

# Create a raster stack from the temperature files
stack_temp <- rotate (rast(list_temp.nc))
stack_prec <- rotate (rast(list_prec.nc))
names (stack_temp)<- list_temp.nc
names (stack_prec)<- list_prec.nc


dir <- "C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/scotese/formatted_data"
setwd(dir)
scotese_times <- data.frame(
  "name" = list.files(dir),
  "ma" = c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
)

# Create a list of month abbreviations
month_ID <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

# Initialize lists for temperature and precipitation data
alldata_temp <- list()
alldata_prec <- list()

# Function for processing temperature and precipitation data for all periods

#time_ID=integer, numbers codifying periods, valid numbers = c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)
# var_name = "1_5m_temp", "precip" 

monthly_stacks <- function(time_ID, var_name, dir) {
  setwd(file.path(dir, time_ID))
  land_mask <- rotate(rast(paste(time_ID, ".qrparm.mask.nc", sep = "")))
  envar_stack <- rast()
  
  for (i in 1:length(month_ID)) {
    file_name <- paste(time_ID, "a.pdcl", month_ID[i], "_", var_name, ".nc", sep = "")
    if (file.exists(file_name)) {
      envar_rast <- rotate(rast(file_name))
      
      if (var_name == "1_5m_temp") {
        envar_rast <- round(envar_rast - 273.15, 1)  # Convert from Kelvin to Celsius
      } else if (var_name == "precip") {
        envar_rast <- round(envar_rast * 60 * 60 * 24 * 30, 0)  # Convert to mm/month
      }
      
      envar_masked <- terra::mask(envar_rast, land_mask[[1]], maskvalues = 0)
      names (envar_masked)<- paste (time_ID, var_name, month_ID[i], sep="_")
      envar_stack <- c(envar_stack, envar_masked)
    }
  }
  
  return(envar_stack)
}


# Process data for each time period
dir <- "C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/scotese/formatted_data"
for (i in 1:nrow(scotese_times)) {
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
  world_temp<- c(world_temp, mean (alldata_temp [[i]]))
  names (world_temp)[i]<- myears [i]
  world_prec<- c(world_prec, sum(alldata_prec [[i]]))
  names (world_prec) [i] <- myears [i]
} 


#plot (world_temp, world_prec, xlim=c(-40, 40), ylim=c(0,6000), main=ma)

########################
# 2) SEDIMENTS

setwd("C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/reconstructed_sediment")
list.files ()
# create stack of sediment data for temp and prec with average prec and average temp
allextract_temp <- list()
allextract_prec <- list()

for (i in 2:15) {
  shapefile <- read_sf(paste0("sediment_recons_", ma [i], "Ma.shp"))
  extracted_data_temp <- extract (world_temp [[i]], shapefile, xy = TRUE, ID = FALSE)
  extracted_data_prec <- extract (world_prec [[i]], shapefile, xy = TRUE, ID = FALSE)
  allextract_temp[[paste0(ma [i], "Ma")]] <- extracted_data_temp
  allextract_prec[[paste0(ma [i], "Ma")]] <- extracted_data_prec
}

## beware that here if you put 15 layers of sediment, you have to remove the minus 1 from the "i"!! and start the loop at 1.

for (i in 2:15) {
  rast_temp_name <- paste0("rast", ma[i], "_temp")
  assign(rast_temp_name, rast (na.omit(allextract_temp[[i-1]][, c(2, 3, 1)])))
  
  rast_prec_name <- paste0("rast", ma[i], "_prec")
  assign(rast_prec_name, rast(na.omit(allextract_prec[[i-1]][, c(2, 3, 1)])))
}

str (allextract_prec)
########################
# 3) FOSSILS

fossil_data <- read.csv("C:/Users/Usuario/OneDrive - Universidade de Vigo/PhD/SEDIMENTOS/DATA/now_database/now_database_09-09-24processed.csv")

str(fossil_data)
head(fossil_data)
summary(fossil_data)

# define the period limits (to reclassify the fossils within our working periods)
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

periods

# reclassify the fossils to the period (ma) that corresponds to our conditions
fossil_data$ma<- NA
for (i in 1:15){
  fossil_data$ma [fossil_data$MIN_AGE>=min (periods [[i]]) & fossil_data$MAX_AGE<= max (periods[[i]])]<- names (periods [i])
  fossil_data$ma [fossil_data$MIN_AGE>=min (periods [[i]]) & fossil_data$MIN_AGE<= max (periods[[i]])]<- names (periods [i])
  fossil_data$ma  [fossil_data$MAX_AGE>=min (periods [[i]]) & fossil_data$MAX_AGE<= max (periods[[i]])]<- names (periods [i])
}

# convert fossil to work
fossil_data<- as.data.frame (fossil_data)
str (fossil_data)
fossil_data$ma<- as.numeric (as.character (fossil_data$ma))
fossil_data_old<- fossil_data
fossil_data<- fossil_data_old [complete.cases (fossil_data_old), ]

# apply palaeorotations (rgplates) to all fossils (palaeolat and palaeolon)
fossil_data$paleolong<- NA
fossil_data$paleolat<- NA
names (fossil_data)

for (i in sort (unique (fossil_data$ma))){
  fossil_data [fossil_data$ma==i, 9:10]<- reconstruct(fossil_data [fossil_data$ma==i, 4:3], 
                                                      age = i, model = "PALEOMAP")
}

# add prec and temp data to each fossil as a function of paleolat and paleolon
sort (unique (fossil_data$ma))
world_temp 
fossil_data$prec<- NA
fossil_data$temp<- NA
tiempos<- sort (unique (fossil_data$ma))
for (i in 1:14){
  fossil_data$temp [fossil_data$ma==tiempos [i]]<- extract (world_temp [[i]], fossil_data [fossil_data$ma==tiempos[i] ,9:10 ])[, 2]
  fossil_data$prec [fossil_data$ma==tiempos [i]]<- extract (world_prec[[i]], fossil_data [fossil_data$ma==tiempos [i] ,9:10 ])[, 2]
}

########################
# PLOTS

# plot T-P diagram with world, sediments and fossils data.
# sediment has no time 0 (has 14 scenarios). world has time 0 (has 15 scenarios). "i" 's are unmatched.
tiempos

par(mfrow = c(3, 5))
for (i in 2:14){
  plot (world_temp[[i]], world_prec[[i]], xlim=c(-40, 40), ylim=c(0,6000), main=ma[[i]])
  points (allextract_temp[[i-1]] [,1], allextract_prec[[i-1]][,1], col=2, pch=16)
  points(fossil_data$temp [fossil_data$ma == tiempos [i]], fossil_data$prec[fossil_data$ma == tiempos [i]], col=3, pch=16)
}

# boxplot temp all periods
dev.off()
par(mfrow = c(3, 5))
for (i in 2:14){
  boxplot (as.vector (world_temp[[i]]), allextract_temp[[i-1]][,1],        
           fossil_data$temp [fossil_data$ma == tiempos [i]], 
           main=ma[[i]], na.rm=T, ylim=c(-40, 40))
  
}

# boxplot prec all periods
par(mfrow = c(3, 5))
for (i in 2:14){
  boxplot (as.vector (world_prec[[i]]), allextract_prec[[i-1]][,1],
           fossil_data$prec [fossil_data$ma == tiempos [i]], 
           main=ma[[i]], na.rm=T, ylim=c(0, 1500))
  
}

########################
#BIOMES RECLASSIFICATION

#  "kg_reclass" function  (S. Galván)

# R functions based on Beck et al. (2018) "Koppengeiger" function on Matlab.
# "Temp" and "Prec" arguments represent temperature and precipitation regimenes, 
# and should be provided as three dimentional arrays with the third dimension 
# representing time (12 months).
# Temp data should be provided in degrees Celsius, and Prep data in mm/months.

# "type" argument has to be "class" (for results indicating the numeric indentifier
# of the specific class 1-30) or "broadclass" (for results indicating the numeric 
# identifier of the broad class (1-5).
library (maps)


kg_reclass <- function(Temp, Prec, type) {
  
  if (identical(dim(Temp), dim(Prec))==F){
    stop("Data matrices have not the same dimensions")
  }
  
  if (any(Prec<0, na.rm = T)==T) {
    stop("Precipitation data can not include negative values") #Added na.rm
  }
  
  if (type != "class" && type != "broadclass") {
    stop("The type argument provided does not exist")
  }
  
  Temp[is.na(Temp)] <- -99 #Added
  Prec[is.na(Prec)] <- -99 #Added
  
  T_ONDJFM <- apply(Temp[, , c(1,2,3,10,11,12)], c(1,2), mean) 
  T_AMJJAS <- apply(Temp[, , c(4, 5, 6, 7, 8, 9)], c(1,2), mean)
  tmp <- T_AMJJAS>T_ONDJFM
  SUM_SEL <- array(as.logical(0), dim(Temp))
  SUM_SEL[,, c(4, 5, 6, 7, 8, 9)] <- rep(tmp, 6)
  SUM_SEL[,, c(10, 11, 12, 1, 2, 3)] <- rep(!tmp,6)
  rm(tmp)
  
  Pw <- apply(Prec*!SUM_SEL, c(1,2), sum) 
  Ps <- apply(Prec*SUM_SEL, c(1,2), sum) 
  
  Pdry <- apply (Prec, c(1,2), min) 
  
  tmp <- SUM_SEL
  tmp[tmp==0] = NA
  Psdry <- apply(Prec*tmp, c(1,2), min, na.rm = TRUE) 
  Pswet <- apply(Prec*tmp, c(1,2), max, na.rm = TRUE) 
  
  tmp <- !SUM_SEL
  tmp[tmp==0] = NA
  Pwdry <- apply(Prec*tmp, c(1,2), min, na.rm = TRUE) 
  Pwwet <- apply(Prec*tmp, c(1,2), max, na.rm = TRUE) 
  
  MAT <- apply(Temp, c(1,2), mean)
  MAP <- apply(Prec, c(1,2), sum) 
  Tmon10 <- apply(Temp > 10, c(1,2), sum) 
  Thot <- apply(Temp, c(1,2), max)
  Tcold <- apply(Temp, c(1,2), min)
  
  Pthresh <- 2*MAT+14 #where temp = -99, this is -184
  Pthresh[Pw>Ps*2.333] <- 2*MAT[Pw>Ps*2.333]       
  Pthresh[Ps>Pw*2.333] <- 2*MAT[Ps>Pw*2.333]+28 
  
  B <- MAP < 10*Pthresh 
  BW <- B & MAP < 5*Pthresh 
  BWh <- BW & MAT >= 18
  BWk <- BW & MAT < 18
  BS <- B & MAP >= 5*Pthresh
  BSh <- BS & MAT >= 18
  BSk <- BS & MAT < 18
  
  A <- Tcold >= 18 & !B 
  Af <- A & Pdry >= 60
  Am <- A & !Af & Pdry >= 100-MAP/25
  Aw <- A & !Af & Pdry < 100-MAP/25
  
  C <- Thot > 10 & Tcold > 0 & Tcold < 18 & !B 
  Cs <- C & Psdry<40 & Psdry<Pwwet/3
  Cw <- C & Pwdry<Pswet/10
  overlap <- Cs & Cw
  Cs[overlap & Ps>Pw] <- 0
  Cw[overlap & Ps<=Pw] <- 0
  Csa <- Cs & Thot >= 22
  Csb <- Cs & !Csa & Tmon10 >= 4
  Csc <- Cs & !Csa & !Csb & Tmon10>=1 & Tmon10<4
  Cwa <- Cw & Thot >= 22
  Cwb <- Cw & !Cwa & Tmon10 >= 4
  Cwc <- Cw & !Cwa & !Cwb & Tmon10>=1 & Tmon10<4
  Cf <- C & !Cs & !Cw
  Cfa <- Cf & Thot >= 22
  Cfb <- Cf & !Cfa & Tmon10 >= 4
  Cfc <- Cf & !Cfa & !Cfb & Tmon10>=1 & Tmon10<4
  
  D <- Thot>10 & Tcold<=0 & !B   
  Ds <- D & Psdry<40 & Psdry<Pwwet/3
  Dw <- D & Pwdry<Pswet/10
  overlap <- Ds & Dw
  Ds[overlap & Ps>Pw] <- 0
  Dw[overlap & Ps<=Pw] <- 0
  Dsa <- Ds & Thot>=22
  Dsb <- Ds & !Dsa & Tmon10 >= 4
  Dsd <- Ds & !Dsa & !Dsb & Tcold < (-38) 
  Dsc <- Ds & !Dsa & !Dsb & !Dsd
  
  Dwa <- Dw & Thot>=22
  Dwb <- Dw & !Dwa & Tmon10 >= 4
  Dwd <- Dw & !Dwa & !Dwb & Tcold < (-38)
  Dwc <- Dw & !Dwa & !Dwb & !Dwd
  Df <- D & !Ds & !Dw
  Dfa <- Df & Thot>=22
  Dfb <- Df & !Dfa & Tmon10 >= 4
  Dfd <- Df & !Dfa & !Dfb & Tcold < (-38)
  Dfc <- Df & !Dfa & !Dfb & !Dfd
  
  E <- Thot <= 10 & Thot > (-90) & !B    #Added & Thot > (-90)
  ET <- E & Thot>0
  EF <- E & Thot<=0 & Thot> (-90) # Added & Thot> (-90)
  
  
  if (type == "class") {
    Class <- list(Af, Am, Aw, BWh, BWk, BSh, BSk, Csa, Csb, Csc, Cwa, Cwb,
                  Cwc, Cfa, Cfb, Cfc, Dsa, Dsb, Dsc, Dsd, Dwa, Dwb, Dwc, Dwd, Dfa,
                  Dfb, Dfc, Dfd, ET, EF)
    Class_cont <- array(0, dim(Temp[,,1]))
    for (i in 1:30){
      Class_cont[Class[[i]]==1] <- i
    }
    Class_cont[Class_cont == 0] <- NA #Added
    return(Class_cont)
  } 
  if (type == "broadclass") {
    Broadclass <- list(A, B, C, D, E)
    BroadClass_cont <- array(0, dim(Temp[,,1]))
    for (i in 1:5){
      BroadClass_cont[Broadclass[[i]]==1] <- i
    }
    BroadClass_cont[BroadClass_cont == 0] <- NA #Added 
    return(BroadClass_cont)
  }   
}

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
for (i in 1:15){
  t<- as.array(alldata_temp[[i]])
  p<- as.array(alldata_prec[[i]])
  biomas [[i]]<-kg_reclass (Temp=t, Prec=p,type="broadclass")
}

# TEST: Project the first period in the WGS84 coordinate system and calculate cellsizes
r<- alldata_temp[[1]]
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
area_raster<- cellSize (r)
dev.off()
plot (r [[1]])
map (add=T)
area<- values (area_raster)
area_raster

dim (area)
plot (area_raster)

# Create a stack of biomes for all periods (projected in the WGS84)
biomes_stack <- rast ()
for (i in 1:15){
  biomes <- rast (biomas [[i]], extent=ext(-181.875, 178.125, -91.25, 91.25), 
                  crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
  biomes_stack <- c (biomes_stack, biomes)
}

# Assign names to each layer of the stack according to age (in millions of years)
names (biomes_stack) <- ma
plot (biomes_stack)

# Rasterise the reconstructed sediment shapefiles for each period
sed_ras <- rast()
for (i in 2:15) {
  shapefile <- read_sf(paste0("sediment_recons_", ma [i], "Ma.shp"))
  ras <- rasterize (x = shapefile, y = biomes_stack[[1]])
  sed_ras <- c (sed_ras, ras)
}


# 1)Calculate the areas occupied by each biome within the fossil sediment zones
biomes_areas <- data.frame (ma = ma [-1], tropical = NA, arid = NA, 
                            temperate = NA, cold = NA, polar = NA)
biomes_areas_percent <- data.frame (ma = ma [-1], tropical = NA, arid = NA, 
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
library(ggplot2)
library(tidyr)

# Convert the dataframe to long format where each row represents a combination of biome and its proportion in a given year
proportion_by_biome_long <- pivot_longer(proportion_by_biome, 
                                         cols = -ma, 
                                         names_to = "biome", 
                                         values_to = "proportion")

# Create stacked bar chart with proportions ("fill" position)
ggplot(proportion_by_biome_long, aes(x = ma, y = proportion, fill = biome)) +
  geom_col(position = "fill") +
  labs(x = "Million years (ma)", 
       y = "Proportion") +
  scale_y_continuous(labels = scales::percent) +  #Axis y as a percentage
  scale_x_continuous(breaks = c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69)) +
  theme_classic()
