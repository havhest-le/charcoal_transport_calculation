#####
# Script to test the OON distribution of the observed plume injection heights ##
#####

rm(list = ls(all= TRUE))

####
## Load packages ################################################################
#####

library(rgdal)
library(rgeos)
library(terra)
library(raster)
library(sf)
library(sp)
library(glue)
library(ggplot2)
library(tidyverse)
library(sp)
library(gdalUtils)
library(GISTools)  
library(dplyr)
library(MODISTools)
library(MODIS)
library(purrr)
library(plyr)
library(disk.frame)
library(stars)
library(Matrix)
library(data.table)
library(tidyselect)
library(dplyr)

suppressMessages(library(raster))
suppressMessages(library(rgdal))
suppressMessages(library(rgeos))
suppressMessages(library(plyr))
suppressMessages(library(sp))
suppressMessages(library(mapview))

## install.packages("devtools")
# devtools:::install_github("gearslaboratory/gdalUtils")
library(gdalUtils)

################################################################################
## Package MODIS ##
## website: https://cran.r-project.org/web/packages/MODIS/index.html
## website: https://github.com/MatMatt/MODIS
## install.packages("MODIS", repos="http://R-Forge.R-project.org")
## library(MODIS)

## Package h4h5tooles ##
## website: https://pjbartlein.github.io/REarthSysSci/hdf5_intro.html
## download: https://portal.hdfgroup.org/display/support/Download+h4h5tools
## install.packages("C:/Users/Vivi/Documents/R/win-library/4.1/h4h5tools-1.10.6-2.2.5-win10_64-vs15.zip")

#####
## Load data ###################################################################
#####

## Load hdf4 data, get the data sds information and change to raster-tiff 
## https://api.rpubs.com/shtien/408239 | function structure 

files_list <- list.files(path = "Data/PIH", pattern = "*hdf", 
                         all.files = T, full.names = T)

#####
## Get data information ########################################################
#####

read_hdf     <- function(file, n, extent) {
  sds        <- get_subdatasets(file)       # get data information
  data       <- tempfile(fileext = ".tif")  # the tiff file is saved as a temporary file
  gdal_translate(sds[n], data)              # select the PIH data set
  raster::brick(x = data)                   # the ouptut is a multi-layer raster file
}

## For loop for every day between 2001 and 2018 #
## Select the PIH data set, write a data.frame and and save as .csv file 

for(i in 2:length(files_list)){
  
  data <- files_list[i] ## select one file of one single day   
  
  ## Read the data structure with data sets
  sds  <- get_subdatasets(data) ## PIH layer is in position [8] of sds
  ## [8] "HDF4_EOS:EOS_GRID:Data/PIH/data:grid1km:Injection_Height"
  
  ## Select the dataset
  PIH  <- read_hdf(data, grep("grid1km:Injection_Height", sds))
  
  ## Create an other layer name
  names(PIH) <- paste0("PIH.", letters[1:nlayers(PIH)])
  
  ## Change crs and create a data frame
  pih_dataframe            <- as.data.frame(PIH, xy=TRUE)
  
  ## Define coordinates as lon and lat
  lat_lon                  <- sin_to_ll(x = pih_dataframe$x, y = pih_dataframe$y) 
  
  ## Merge PIH values with lon and lat coordinates
  pih_dataframe_WGS48      <- cbind(pih_dataframe, lat_lon) 
  pih_dataframe_WGS48[1:2] <- NULL
  
  ## Change columns names
  colnames(pih_dataframe_WGS48)[which(names(pih_dataframe_WGS48) == "longitude_ll")] <- "lon"
  colnames(pih_dataframe_WGS48)[which(names(pih_dataframe_WGS48) == "latitude_ll")]  <- "lat"
  
  ## Add the date of the file in extra column
  pih_dataframe_WGS48$date <- substr(data,12,21)
  
  ## Change positions of columns lat, lon and date
  pih_dataframe_WGS48_sub  <- pih_dataframe_WGS48      %>% relocate(lat,  .before = PIH.a) 
  pih_dataframe_WGS48_sub  <- pih_dataframe_WGS48_sub  %>% relocate(lon,  .before = PIH.a)
  pih_dataframe_WGS48_sub  <- pih_dataframe_WGS48_sub  %>% relocate(date, .before = PIH.a)
  
  ## Select the columns if the column name is starting with "PIH"
  PIH_select <- pih_dataframe_WGS48_sub %>% dplyr::select(starts_with("PIH"))
  col_num    <- ncol(PIH_select)-1 # defining the col numbers for remove rows with NA's

  ## Remove rows in data.frame which have NA's values in greater or equal columns numbers
  delete.na           <- function(data, n){ data[rowSums(is.na(data)) <= n,]}
  PIH_dataframe_final <- delete.na(data = pih_dataframe_WGS48_sub, n = col_num)

  ## Define the filename
  filename            <- substr(data,10,21)
  
  ## Print "{filename}_data.frame is empty" or save the csv-file in Results
  if(nrow(PIH_dataframe_final) == 0){
    print(glue("{filename}_data.frame is empty"))
  }else{
    write.csv(PIH_dataframe_final, file = glue("Results/PIH/PIH_csv/{filename}_PIH.csv"))
    }
}


#####
## Load csv files of PIH values ################################################
#####

## Load data and create a list of the .csv-files of PIH dataframes in R # 
PIH_data <- list.files(path = "Results/PIH/PIH_csv", pattern = ".csv", 
                       all.files = T, full.names = T)

## Selecte the different years (fire years, defined with MODIS FRP values)
PIH_select_years <- function(year, data){
                    Filter(function(e) grepl(year, e), data)}

## Classification of extent 1 and 2 of the study areas
PIH_1 <- PIH_select_years("1_2", PIH_data)
PIH_2 <- PIH_select_years("2_2", PIH_data)

## Select the different fire years
## 1. extent
PIH_1_2001 <- PIH_select_years("_2001", PIH_1)
PIH_1_2002 <- PIH_select_years("_2002", PIH_1)
PIH_1_2003 <- PIH_select_years("_2003", PIH_1)
PIH_1_2004 <- PIH_select_years("_2004", PIH_1)
PIH_1_2007 <- PIH_select_years("_2007", PIH_1)
PIH_1_2013 <- PIH_select_years("_2013", PIH_1)
PIH_1_2014 <- PIH_select_years("_2014", PIH_1)
PIH_1_2017 <- PIH_select_years("_2017", PIH_1)
PIH_1_2018 <- PIH_select_years("_2018", PIH_1)

## 2. extent
PIH_2_2001 <- PIH_select_years("_2001", PIH_2) # not first one
PIH_2_2003 <- PIH_select_years("_2003", PIH_2)
PIH_2_2004 <- PIH_select_years("_2004", PIH_2)
PIH_2_2007 <- PIH_select_years("_2007", PIH_2)
PIH_2_2013 <- PIH_select_years("_2013", PIH_2)
PIH_2_2014 <- PIH_select_years("_2014", PIH_2)

#####
## Merge the extent 1 and 2 of study area Lake Khamra ##########################
#####

Khamra_PIH_in_list_merge <- function(data_1, data_2, filename){ 
  
  list_1 <- list() ## create two empty lists
  list_2 <- list()
  
  for(o in PIH_1_2001){
    get_data    <- read.csv(o) ## open data.frames
    list_1[[o]] <- get_data    ## put data.frames for every single day in one list
    un <- unlist(list_1)
  } 
  
  ## matrix <- ldply(PIH_1_2001,readr::read_csv)
  csv_1    <- lapply(PIH_1_2001, read.csv)
  result_1 <- do.call(rbind.fill, csv_1)
  
  csv_2    <- lapply(PIH_2_2001, read.csv)
  result_2 <- do.call(rbind.fill, csv_2) 
  
  rbind_data    <- rbind.fill(result_1, result_2) 
  

  
  for(i in PIH_2_2001){
    get_data    <- read.csv(i) # Open data.frames
    list_2[[i]] <- get_data    # Put data.frames for every single day in one list
  }
  
  memory.limit(9999999999)
  rbind_data    <- rbind.fill(list_1, list_2) ## Merge list 1 and 2 to encompass 
                                              ## the entire study area
  rbind_data[1] <- NULL
  save(rbind_data, file = glue("Results/PIH/Khamra_PIH/{filename}_PIH.RData"))
  
  }

## Merge extent 1 and 2 for single years (Lake Khamra)
Khamra_2001 <- Khamra_PIH_in_list_merge(data_1 = PIH_1_2001, data_2 = PIH_2_2001, filename = "Khamra_2001")
Khamra_2003 <- Khamra_PIH_in_list_merge(data_1 = PIH_1_2003, data_2 = PIH_2_2003, filename = "Khamra_2003")
Khamra_2004 <- Khamra_PIH_in_list_merge(data_1 = PIH_1_2004, data_2 = PIH_2_2004, filename = "Khamra_2004")
Khamra_2007 <- Khamra_PIH_in_list_merge(data_1 = PIH_1_2007, data_2 = PIH_2_2007, filename = "Khamra_2007")
Khamra_2013 <- Khamra_PIH_in_list_merge(data_1 = PIH_1_2013, data_2 = PIH_2_2013, filename = "Khamra_2013")
Khamra_2014 <- Khamra_PIH_in_list_merge(data_1 = PIH_1_2014, data_2 = PIH_2_2014, filename = "Khamra_2014")

## Load merged PIH-RData for lake Khamra in R again
## Khamra fire years: 2001, 2003, 2004, 2007, 2013, 2014

Khamra_2001_PIH <- get(load("Results/PIH/Khamra_PIH/Khamra_2001_PIH.RData"))
Khamra_2003_PIH <- get(load("Results/PIH/Khamra_PIH/Khamra_2003_PIH.RData"))
Khamra_2004_PIH <- get(load("Results/PIH/Khamra_PIH/Khamra_2004_PIH.RData"))
Khamra_2007_PIH <- get(load("Results/PIH/Khamra_PIH/Khamra_2007_PIH.RData"))
Khamra_2013_PIH <- get(load("Results/PIH/Khamra_PIH/Khamra_2013_PIH.RData"))
Khamra_2014_PIH <- get(load("Results/PIH/Khamra_PIH/Khamra_2014_PIH.RData"))


## Merge all Khamra files
Khamra_all_PIH  <- rbind.fill(Khamra_2001_PIH, Khamra_2003_PIH, Khamra_2004_PIH, 
                              Khamra_2007_PIH, Khamra_2013_PIH, Khamra_2014_PIH)

## Load PIH-Data for lake Satagay and combine the single days for every year # 
Satagay_PIH_open <- function(data_1, filename){

                    list_1 <- list()             ## Create an empty list

                    for(p in data_1){
                    get_data    <- read.csv(p)   ## Open the data.frames
                    list_1[[p]] <- get_data  }   ## Put data.frames of single days in one list

                    ## Merge single days in a data.frame
                    rbind_data    <- rbind.fill(list_1, fill = TRUE) ## Merge rows of single days
                    rbind_data[1] <- NULL
                    save(rbind_data, file = glue("Results/PIH/Satagay_PIH/{filename}_PIH.RData"))

                    }

## Merge rows of single days for every year
Satagay_2002 <- Satagay_PIH_open(data_1 = PIH_1_2002, filename = "Satagay_2002")
Satagay_2013 <- Satagay_PIH_open(data_1 = PIH_1_2013, filename = "Satagay_2013")
Satagay_2014 <- Satagay_PIH_open(data_1 = PIH_1_2014, filename = "Satagay_2014")
Satagay_2017 <- Satagay_PIH_open(data_1 = PIH_1_2017, filename = "Satagay_2017")
Satagay_2018 <- Satagay_PIH_open(data_1 = PIH_1_2018, filename = "Satagay_2018")

## Load merged PIH-RData for lake Satagay in R again
## Satagay fire years: 2002, 2013, 2014, 2018
Satagay_2002_PIH <- get(load("Results/PIH/Satagay_PIH/Satagay_2002_PIH.RData"))
Satagay_2013_PIH <- get(load("Results/PIH/Satagay_PIH/Satagay_2013_PIH.RData"))
Satagay_2014_PIH <- get(load("Results/PIH/Satagay_PIH/Satagay_2014_PIH.RData"))
Satagay_2018_PIH <- get(load("Results/PIH/Satagay_PIH/Satagay_2018_PIH.RData"))

## Merge all Satagay files
Satagay_all_PIH <- rbind.fill(Satagay_2002_PIH, Satagay_2013_PIH, Satagay_2014_PIH,
                              Satagay_2018_PIH)

#####
## Crop PIH data files with the defined extents ################################
#####

## Khamra_buf_extent  # 111.1664, 114.8065, 59.08353, 60.89779  (xmin, xmax, ymin, ymax)
## Satagay_buf_extent # 116.0127, 120.009, 62.17579, 63.98115   (xmin, xmax, ymin, ymax)
 
crop_PIH_files <- function(data, lon_bigger, lon_smaller, lat_bigger, lat_smaller){
      data_sub <- data %>% filter(lon >= lon_bigger & lon <= lon_smaller &  
                                  lat >= lat_bigger & lat <= lat_smaller)}

## Crop Khamra extent
Khamra_all_PIH_sub  <- crop_PIH_files(data = Khamra_all_PIH,
                                      lon_bigger = 111.16, lon_smaller = 114.80,
                                      lat_bigger = 59.08 , lat_smaller = 60.89)

Khamra_all_PIH_sub$PIH.h <- as.numeric(Khamra_all_PIH_sub$PIH.h) ## set col "PIH.h" as.numeric
Khamra_all_PIH_sub$PIH.i <- as.numeric(Khamra_all_PIH_sub$PIH.i) ## set col "PIH.i" as.numeric

## Crop Satagay extent 
Satagay_all_PIH_sub <- crop_PIH_files(data = Satagay_all_PIH,
                                      lon_bigger = 116.01, lon_smaller = 120.00,
                                      lat_bigger = 62.17 , lat_smaller = 63.98)


Satagay_all_PIH_sub$PIH.h <- as.numeric(Satagay_all_PIH_sub$PIH.h) ## set col "PIH.h" as.numeric

#####
## Calculation #################################################################
#####

## Calculate the statistic values of PIH data for all orbit numbers #
library(matrixStats)

## Get the statistical information 
statistics_PIH          <- function(data, filename, col_n, col_n_select){
     col_pos            <- col_n
     col_pos_select     <- col_n_select
     data$range         <- matrixStats::rowRanges(as.matrix(data[,col_pos]),na.rm = T)
     data$max           <- matrixStats::rowMins(as.matrix(data[,col_pos]),na.rm = T)
     data$min           <- matrixStats::rowMaxs(as.matrix(data[,col_pos]),na.rm = T)
     data$mean          <- matrixStats::rowMeans2(as.matrix(data[,col_pos]), na.rm = T)
     data$median        <- matrixStats::rowMedians(as.matrix(data[,col_pos]),na.rm = T)
     data$median_select <- matrixStats::rowMedians(as.matrix(data[,col_pos_select]),na.rm = T)
     data$sd            <- matrixStats::rowSds(as.matrix(data[,col_pos]), na.rm = T)
     data               <- as.data.frame(data)
     save(data, file = glue("Results/Simulation/PIH/data/{filename}_PIH_statistics.rda"))
}

Khamra_PIH_stat  <- statistics_PIH(col_n = c(4,5,6,7,8,9,10,11,12),
                                   data = Khamra_all_PIH_sub, filename = "Khamra",
                                   col_n_select = c(6,7,8,9,10))

Satagay_PIH_stat <- statistics_PIH(col_n = c(4,5,6,7,8,9,10,11), 
                                   data = Satagay_all_PIH_sub, filename = "Satagay",
                                   col_n_select = c(6,7,8,9,10))

## Load PIH_stat files again in R 
Khamra_PIH_stat  <- get(load("Results/Simulation/PIH/data/Khamra_PIH_statistics.rda"))
Satagay_PIH_stat <- get(load("Results/Simulation/PIH/data/Satagay_PIH_statistics.rda")) 
Satagay_PIH_stat$median        <- na.omit(Satagay_PIH_stat$median)
Satagay_PIH_stat$median_select <- as.numeric(Satagay_PIH_stat$median_select)


## Define the time structure
Khamra_PIH_stat$date  <- as.POSIXct(x = Khamra_PIH_stat$date)
Satagay_PIH_stat$date <- as.POSIXct(x = Satagay_PIH_stat$date)

## Define the years
Khamra_PIH_stat$Year  <- as.numeric(format(Khamra_PIH_stat$date, "%Y"))
Satagay_PIH_stat$Year <- as.numeric(format(Satagay_PIH_stat$date, "%Y"))


## Filter relevant information of data frame
Khamra_PIH_stat  <- Khamra_PIH_stat %>% dplyr::select(lat, lon, date, Year, median)
Satagay_PIH_stat <- Satagay_PIH_stat %>% dplyr::select(lat, lon, date, Year, median)


Khamra_PIH_stat_2014     <- Khamra_PIH_stat %>% filter(Khamra_PIH_stat$Year == 2014)
Khamra_PIH_stat_2014     <- Khamra_PIH_stat_2014 %>% select(lat, lon, median)
Khamra_PIH_stat_2014$lat <- round(Khamra_PIH_stat_2014$lat, 2)
Khamra_PIH_stat_2014$lon <- round(Khamra_PIH_stat_2014$lon, 2)
colnames(Khamra_PIH_stat_2014)[which(names(Khamra_PIH_stat_2014) == "median")] <- "PIH"
Khamra_PIH_points_2014   <- st_as_sf(x = Khamra_PIH_stat_2014,
                            coords = c("lon", "lat"), # create coordinate geometry for points
                            crs = 4326)                          # change crs to WGS 48
save(Khamra_PIH_points_2014, file = "Data/PIH/Data_points/Khamra_PIH_points_2014.RData")

Khamra_PIH_stat$lat <- round(Khamra_PIH_stat$lat, 2)
Khamra_PIH_stat$lon <- round(Khamra_PIH_stat$lon, 2)
colnames(Khamra_PIH_stat)[which(names(Khamra_PIH_stat) == "median")]   <- "PIH"
colnames(Satagay_PIH_stat)[which(names(Satagay_PIH_stat) == "median")] <- "PIH"



## Convert data frame to points
Khamra_PIH_points  <- st_as_sf(x = Khamra_PIH_stat,
                               coords = c("lon", "lat"), ## Create coordinate geometry for points
                               crs = 4326)               ## Change crs to WGS 48

Satagay_PIH_points <- st_as_sf(x = Satagay_PIH_stat,
                               coords = c("lon", "lat"),
                               crs = 4326)

# # save(Khamra_PIH_points, file = "Results/Simulation/PIH/data/Khamra_PIH_points.RData")
# # save(Satagay_PIH_points, file = "Results/Simulation/PIH/data/Satagay_PIH_points.RData")
# 
# # Khamra_PIH_points  <- get(load("Results/Simulation/PIH/data/Khamra_PIH_points.RData"))
# # Satagay_PIH_points <- get(load("Results/Simulation/PIH/data/Satagay_PIH_points.RData"))


## Which atmosphere levels are important for calculation? 
PIH_min_Khamra  <- min(Khamra_PIH_stat$min)
PIH_max_Khamra  <- max(Khamra_PIH_stat$max)

PIH_min_Satagay <- min(Satagay_PIH_stat$min)
PIH_max_Satagay <- max(Satagay_PIH_stat$max)

## Calculate the variance of PIH values of the orbit overpass numbers 
Khamra_PIH_OBN_c    <- c(Khamra_PIH_stat$PIH.a, Khamra_PIH_stat$PIH.b, Khamra_PIH_stat$PIH.c,
                         Khamra_PIH_stat$PIH.d, Khamra_PIH_stat$PIH.e, Khamra_PIH_stat$PIH.f,
                         Khamra_PIH_stat$PIH.g, Khamra_PIH_stat$PIH.h, Khamra_PIH_stat$PIH.i)

Satagay_PIH_OBN_c   <- c(Satagay_PIH_stat$PIH.a, Satagay_PIH_stat$PIH.b, Satagay_PIH_stat$PIH.c,
                         Satagay_PIH_stat$PIH.d, Satagay_PIH_stat$PIH.e, Satagay_PIH_stat$PIH.f,
                         Satagay_PIH_stat$PIH.g, Satagay_PIH_stat$PIH.h)

var_PIH_OBN_Khamra  <- var(Khamra_PIH_OBN_c, na.rm = T)
var_PIH_OBN_Satagay <- var(Satagay_PIH_OBN_c, na.rm = T)

library(stats)
sd_PIH_OBN_Khamra   <- sd(Khamra_PIH_OBN_c, na.rm = T)  # 261.99
sd_PIH_OBN_Satagay  <- sd(Satagay_PIH_OBN_c, na.rm = T) # 359.44




#####
## Plot the dispersion of PIH values of single orbit numbers ###################
#####

col_PIH_orbit <- c("white", "azure2", "lightblue3", "lightblue3", "lightblue3", "lightblue3", "lightblue3", "white")
label_names   <- c(1:8)

## Lake Satagay with outliners
png(glue("Results/Simulation/PIH/Satagay/PIH_boxplot_orbit_numbers_Satagay.png"), width = 1000, height = 800)
  
boxplot(Satagay_PIH_stat[4:11], outline = T, na.action	= F,
                                      xlab = "Orbit overpass numbers (1 to 8)", 
                                      ylab = "Measured PIH [m]", col = col_PIH_orbit,
                                      ylim = c(0, 3000), cex.axis = 1.4, cex.lab = 1.7, names = label_names)
mytitle    = glue("Distribution of plume injection height values as a function\nof the number of overflights in the study area Lake Satagay")
legend("topright", legend=c("standard deviation = 359.44"), bty = "n", cex = 1.5)
mtext(line = 0.4,cex = 2.2, mytitle, font = 2)
# Save as png
dev.off()

## Lake Satagay without outliners
png(glue("Results/Simulation/PIH/Satagay/PIH_boxplot_orbit_numbers_Satagay_without_outliners.png"), width = 1000, height = 800)

boxplot(Satagay_PIH_stat[4:11], outline = F, na.action	= F, col = col_PIH_orbit,
        xlab = "Orbit overpass numbers (1 to 8)", 
        ylab = "Measured PIH [m]",
        ylim = c(0, 2000), cex.axis = 1.4, cex.lab = 1.7, names = label_names)
mytitle    = glue("Distribution of plume injection height values as a function\nof the number of overflights in the study area Lake Satagay")
legend("topright", legend=c("standard deviation = 359.44"), bty = "n", cex = 1.5)
mtext(line = 0.4,cex = 2.2, mytitle, font = 2)
# Save as png
dev.off()

summary(Satagay_PIH_stat[4:11])


## only orbit overpass numbers (c to g)
png(glue("Results/Simulation/PIH/Satagay/PIH_boxplot_orbit_numbers_Satagay_without_b.png"), width = 1000, height = 800)
boxplot(Satagay_PIH_stat[6:10], outline = F, na.action	= F, col = "lightblue3",
        xlab = "Orbit overpass numbers (a to h)", 
        ylab = "Measured PIH [m]",
        ylim = c(0, 2000), cex.axis = 1.4, cex.lab = 1.7)
mytitle    = glue("Distribution of plume injection height values as a function\nof the number of overflights in the study area Lake Satagay")
legend("topright", legend=c("standard deviation = 359.44"), bty = "n", cex = 1.5)
mtext(line = 0.4,cex = 2.2, mytitle, font = 2)
# Save as png
dev.off()



#####
## Multiple plot for orbit overpass numbers ####################################
#####

med_all_S_orbit <- round(median(Satagay_PIH_stat$median, na.rm = T),1)
med_sel_S_orbit <- round(median(Satagay_PIH_stat$median_select, na.rm = T),1)

col_PIH_orbit_2   <- c("skyblue2", "lightblue3")
label_names_2 <- c(glue("median = {med_all_S_orbit}"),glue("median = {med_sel_S_orbit}"))

png(glue("Results/Simulation/PIH/Satagay/PIH_boxplot_orbit_numbers_Satagay_overview.png"), width = 600, height = 800)
layout(matrix(c(1,1,2,2),nrow = 2, ncol = 2, byrow = TRUE))

boxplot(Satagay_PIH_stat[4:11], outline = F, na.action	= F, col = col_PIH_orbit,
        xlab = "OON (1 to 8)", 
        ylab = "Measured PIH [m]",
        ylim = c(0, 2000), cex.axis = 1.75, cex.lab = 1.75, names = label_names)

boxplot(Satagay_PIH_stat$median, Satagay_PIH_stat$median_select, outline = F, 
        names = label_names_2, col = col_PIH_orbit_2, cex.lab = 1.75, cex.axis = 1.75, ylim = c(0,2000),
        ylab = "Measured PIH [m]",  frame.plot = FALSE)

legend("topright", legend = c("PIH Distribution over all OON","PIH Distribution over OON 3 to 7", "OON 2"), 
       col=c("skyblue2", "lightblue3", "azure2"), pt.cex=3, pch=15, cex = 2)

dev.off()


PIH.b_2 <- na.omit(Satagay_PIH_stat$PIH.b)
length(PIH.b_2) ## 17954
PIH.c_3 <- na.omit(Satagay_PIH_stat$PIH.c)
length(PIH.c_3) ## 383101
PIH.d_4 <- na.omit(Satagay_PIH_stat$PIH.d)
length(PIH.d_4) ## 322364
PIH.e_5 <- na.omit(Satagay_PIH_stat$PIH.e)
length(PIH.e_5) ## 374703
PIH.f_5 <- na.omit(Satagay_PIH_stat$PIH.f)
length(PIH.f_5) ## 477469
PIH.g_6 <- na.omit(Satagay_PIH_stat$PIH.g)
length(PIH.g_6) ## 86145


#####
# Calculation the monthly median of PIH values #################################
#####

## Remove single PIH data of every orbit number 
Khamra_PIH_stat[4:12]  <- NULL
Satagay_PIH_stat[4:12] <- NULL

## Calculate the monthly median of PIH mean values for every single year #
stat_dataframe_PIH_mean <- function(data, filename){
              datalist  <- list()                                    ## create an empty list

                           for(y in 2000:2018){                      ## Define the time period
                                                                     ## y for year
                           monthList <- lapply(5:8, function(m) {    ## m for month (may to august)
                                        
                           cat(glue("\ryear {y} month {m}"))
                          
                           ## Select the time period of interest
                           subTab_PIH <- subset(data, as.numeric(format(data$date, "%Y")) %in% y & 
                                                      as.numeric(format(data$date, "%m")) %in% m ) 
                            
                           ## Calculate the statistic for every single month
                           median <- median(subTab_PIH$mean, nr.rm = TRUE)
                           min    <- min(subTab_PIH$mean, nr.rm = TRUE)
                           max    <- max(subTab_PIH$mean, nr.rm = TRUE)
                           sd     <- sd(subTab_PIH$mean)
                            
                           ## Create a data frame of statistic information of PIH 
                           df_PIH <- data.frame(median, min, max, sd , y, m)
                            
                           })
                          
                           ## Merge the data of single months  
                           stat_rbind_PIH <- rbind(monthList[[1]], monthList[[2]], 
                                                   monthList[[3]], monthList[[4]])
                           
                           ## Create a list with the data frames of single years
                           datalist[[y]]  <- stat_rbind_PIH
                            
                            }
                            
                           ## Merge all years and months in one data frame
                           big_data_PIH   <- do.call(rbind, datalist)
                            
                           ## save as .RData
                           save(big_data_PIH, file = glue("Results/PIH/PIH_mean_statistic_dataframe_{filename}.RData"))
                            
                           }


Khamra_PIH_dataframe_stat  <- stat_dataframe_PIH_mean(data = Khamra_PIH_stat, filename = "Khamra")
Satagay_PIH_dataframe_stat <- stat_dataframe_PIH_mean(data = Satagay_PIH_stat, filename = "Satagay")

## Load data frames of PIH mean stat 
Khamra_PIH_dataframe_stat  <- get(load("Results/PIH/PIH_mean_statistic_dataframe_Khamra.RData"))
Satagay_PIH_dataframe_stat <- get(load("Results/PIH/PIH_mean_statistic_dataframe_Satagay.RData"))


## Calculate the monthly median of PIH median values for every single year #
stat_dataframe_PIH_median <- function(data, filename){
                datalist  <- list()
                 
                             for(y in 2000:2018){                      # for the period 2000 to 2018
                          
                             monthList <- lapply(5:8, function(m) {    # for months 5 to 8
                            
                             cat(glue("\ryear {y} month {m}"))
                            
                             subTab_PIH <- subset(data, as.numeric(format(data$date, "%Y")) %in% y & 
                                                        as.numeric(format(data$date, "%m")) %in% m ) 
                            
                             # calculate the statistic for every single month
                             median <- median(subTab_PIH$median, nr.rm = TRUE)
                             min    <- min(subTab_PIH$median, nr.rm = TRUE)
                             max    <- max(subTab_PIH$median, nr.rm = TRUE)
                             sd     <- sd(subTab_PIH$median)
                            
                             # Create a data frame of statistic information of PIH 
                             df_PIH <- data.frame(median, min, max, sd , y, m)
                            
                             })
                      
                             # Merge the data of single months  
                             stat_rbind_PIH <- rbind(monthList[[1]], monthList[[2]], 
                                                     monthList[[3]], monthList[[4]])
                             
                             # Create a list
                             datalist[[y]]  <- stat_rbind_PIH
                            
                             }
  
                             # Merge all years and months in one data frame
                             big_data_PIH   <- do.call(rbind, datalist)
                            
                             # save as .RData
                             save(big_data_PIH, file = glue("Results/PIH/PIH_median_statistic_dataframe_{filename}.RData"))
                            
                             }

Khamra_PIH_dataframe_stat_median  <- stat_dataframe_PIH_median(data = Khamra_PIH_stat, filename = "Khamra")
Satagay_PIH_dataframe_stat_median <- stat_dataframe_PIH_median(data = Satagay_PIH_stat, filename = "Satagay")

## Load data frames of PIH median stat 
Khamra_PIH_dataframe_stat_median  <- get(load("Results/PIH/PIH_median_statistic_dataframe_Khamra.RData"))
Satagay_PIH_dataframe_stat_median <- get(load("Results/PIH/PIH_median_statistic_dataframe_Satagay.RData"))
