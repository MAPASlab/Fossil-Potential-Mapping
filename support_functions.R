## METADATA ===============================================================
## Description: Support functions for fossils potential mapping
## 
## R version: 4.2.2 for Windows
## Date: 2024-10-15 19:00:34
## License: TODO MARTA
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##



#' Process Temperature and Precipitation Data for Geological Periods
#'
#' This function processes temperature and precipitation data for a given geological period, specified by `time_ID`. It reads in climate model data, applies necessary transformations (such as unit conversions), and returns a stack of rasters for each month.
#' 
#' @param time_ID Integer. Identifier for the geological period, valid values are c(0, 3, 11, 15, 20, 26, 31, 36, 40, 45, 52, 56, 61, 66, 69).
#' @param var_name Character. The variable name, either "1_5m_temp" for temperature or "precip" for precipitation.
#' @param dir Character. The directory containing the data files.
#'
#' @return A stack of raster objects, each representing the data for one month of the given period.
#' 
#' @details The function reads NetCDF files for each month, performs transformations depending on the variable (`var_name`), and applies a land mask to the resulting raster. Temperature is converted from Kelvin to Celsius, while precipitation is converted to mm/month.
#'
#' @examples
#' 
#' # Example usage
#' result <- monthly_stacks(15, "1_5m_temp", "./data/scotese/formatted_data")
#' 
#' @importFrom terra rast mask
#' @importFrom terra rotate
monthly_stacks <- function(time_ID, var_name, dir) {
  # setwd(file.path(dir, time_ID))
  land_mask <- rotate(rast(file.path(dir, time_ID, paste(time_ID, ".qrparm.mask.nc", sep = ""))))
  envar_stack <- NULL
  
  for (i in 1:length(month_ID)) {
    # i <- 1
    file_name <- file.path(dir, time_ID, paste(time_ID, "a.pdcl", month_ID[i], "_", var_name, ".nc", sep = ""))
    if (file.exists( file_name)) {
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


#' Köppen-Geiger Climate Classification Reclassification
#'
#' This function performs a Köppen-Geiger climate classification reclassification for given temperature (`Temp`) and precipitation (`Prec`) data.
#' The function can classify climate zones into either specific climate types (`class`) or broad climate categories (`broadclass`).
#' Fucntion authored by S. Galván, documented by O. Hagen.
#' 
#' @param Temp A three-dimensional array of temperature data (in degrees Celsius).
#' @param Prec A three-dimensional array of precipitation data (in mm/month).
#' @param type Character. Specifies the classification type, either "class" for specific classes (1-30) or "broadclass" for broad classes (1-5).
#' 
#' @return A two-dimensional array representing the climate classification for each cell in the spatial matrix.
#'
#' @details This function performs reclassification based on monthly temperature and precipitation data across different seasons.
#' It applies a variety of rules to determine the climatic category for each grid cell. The classification can be granular (`class`) or broader (`broadclass`).
#'
#' @examples
#'
#' # Example usage
#' temp_data <- array(runif(100, -10, 30), dim = c(10, 10, 12))
#' prec_data <- array(runif(100, 0, 500), dim = c(10, 10, 12))
#' result <- kg_reclass(temp_data, prec_data, type = "broadclass")
#'
#' @importFrom base apply
#' @importFrom terra rast
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
