#####
# Script to calculate the plume injection height #
#####

rm(list = ls(all= TRUE))

# Load packages
library(rgdal)
library(rgeos)
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
library(hdf5r)
library(rhdf5)
library(MODIS)
library(purrr)
library(plyr)


suppressMessages(library(gdalUtils))
suppressMessages(library(raster))
suppressMessages(library(rgdal))
suppressMessages(library(rgeos))

suppressMessages(library(plyr))
suppressMessages(library(magrittr))
suppressMessages(library(sp))
suppressMessages(library(mapview))

# update.packages(ask = FALSE)

#################################################################################

## Package MODIS ##
# website: https://cran.r-project.org/web/packages/MODIS/index.html
# website: https://github.com/MatMatt/MODIS
# install.packages("MODIS", repos="http://R-Forge.R-project.org")
# library(MODIS)

## Package h4h5tooles ##
# website: https://pjbartlein.github.io/REarthSysSci/hdf5_intro.html
# download: https://portal.hdfgroup.org/display/support/Download+h4h5tools
# install.packages("C:/Users/Vivi/Documents/R/win-library/4.1/h4h5tools-1.10.6-2.2.5-win10_64-vs15.zip")


######################################################################################

# https://api.rpubs.com/shtien/408239 | function structure 


read_hdf     <- function(file, n, extent) {
  sds        <- get_subdatasets(file)
  data       <- tempfile(fileext = ".tif")        # the tiff file is saved as a temporary file
  gdal_translate(sds[n], data)
  raster::brick(x = data)           # the ouptut is a multi-layer raster file
}

files_list <- list.files(path = "Data/PIH", pattern = "*hdf", all.files = T, full.names = T)

##################################################################################

for(i in 2:length(files_list)){
  
  data <- files_list[i]
  
  # Read the data structure with datasets
  sds  <- get_subdatasets(data) # pih is in position [8]
  # [8] "HDF4_EOS:EOS_GRID:Data/PIH/data:grid1km:Injection_Height"
  
  # Select the necessary dataset
  PIH  <- read_hdf(data, grep("grid1km:Injection_Height", sds))
  
  # Create a different name for layer
  names(PIH) <- paste0("PIH.", letters[1:nlayers(PIH)])
  
  # Change crs and create a dataframe
  pih_dataframe            <- as.data.frame(PIH, xy=TRUE)
  lat_lon                  <- sin_to_ll(x = pih_dataframe$x, y = pih_dataframe$y)
  pih_dataframe_WGS48      <- cbind(pih_dataframe, lat_lon)
  pih_dataframe_WGS48[1:2] <- NULL
  colnames(pih_dataframe_WGS48)[which(names(pih_dataframe_WGS48) == "longitude_ll")] <- "lon"
  colnames(pih_dataframe_WGS48)[which(names(pih_dataframe_WGS48) == "latitude_ll")]  <- "lat"
  
  # Add filename
  pih_dataframe_WGS48$date <- substr(data,12,21)
  filename                 <- substr(data,10,21)
  
  # Save data as RData
  
  write.csv(pih_dataframe_WGS48, file = glue("E:/Bachelorarbeit_Wind/R-Scripte/SiberiaWind_R/Results/PIH/{filename}_PIH.csv"))
  
}


###################################################################################
# list the csv data in R
PIH_data <- list.files(path = "Results/PIH/", pattern = ".csv", all.files = T, full.names = T)

for(t in PIH_data){
  memory.limit(9999999999)
  read_data          <- read.csv(t)  # read csv data
  read_data[1]       <- NULL
  read_data_sub      <- read_data     %>% relocate(lat,  .before = PIH.a)
  read_data_sub      <- read_data_sub %>% relocate(lon,  .before = PIH.a)
  read_data_sub      <- read_data_sub %>% relocate(date, .before = PIH.a)
  filename           <- substr(t,13,24)
  save(read_data_sub, file = glue("Results/PIH/open_data/{filename}_PIH_open.RData"))
}


####################################################################################
file.remove("Results/PIH/open_data/2_2001-05-01_PIH_open.RData") 
PIH_RData <- list.files(path = "Results/PIH/open_data", pattern = ".RData", all.files = T, full.names = T)

# all_data_frames <- get(load(lapply((list.files(path = "Results/PIH/open_data", pattern = ".RData")))))


PIH_select_years <- function(year, data){
  Filter(function(e) grepl(year, e), data)
}

PIH_1 <- PIH_select_years("1_2", PIH_RData)
PIH_2 <- PIH_select_years("2_2", PIH_RData)

# Selecting fire years 
# 1. extent
PIH_1_2001 <- PIH_select_years("_2001", PIH_1)
PIH_1_2002 <- PIH_select_years("_2002", PIH_1)
PIH_1_2003 <- PIH_select_years("_2003", PIH_1)
PIH_1_2004 <- PIH_select_years("_2004", PIH_1)
PIH_1_2007 <- PIH_select_years("_2007", PIH_1)
PIH_1_2013 <- PIH_select_years("_2013", PIH_1)
PIH_1_2014 <- PIH_select_years("_2014", PIH_1)
PIH_1_2017 <- PIH_select_years("_2017", PIH_1)
PIH_1_2018 <- PIH_select_years("_2018", PIH_1)

# 2. extent
PIH_2_2001 <- PIH_select_years("_2001", PIH_2) # not first one
PIH_2_2003 <- PIH_select_years("_2003", PIH_2)
PIH_2_2004 <- PIH_select_years("_2004", PIH_2)
PIH_2_2007 <- PIH_select_years("_2007", PIH_2)
PIH_2_2013 <- PIH_select_years("_2013", PIH_2)
PIH_2_2014 <- PIH_select_years("_2014", PIH_2)

PIH_in_list <- function(data_1, data_2, days){
  
  list_1 <- list()
  list_2 <- list()
  
  for(o in data_1){
    memory.limit(9999999999)
    get_data    <- get(load(o))
    list_1[[o]] <- get_data
  }
  
  for(i in data_2){
    memory.limit(9999999999)
    get_data    <- get(load(i))
    list_2[[i]] <- get_data
  }
  
  memory.limit(9999999999)
  rbind.fill(list_1[days], list_2[days])
}


Khamra_2001 <- PIH_in_list(data_1 = PIH_1_2001, data_2 = PIH_2_2001, days = c(1:104))

Khamra_2003 <- PIH_in_list(data_1 = PIH_1_2003, data_2 = PIH_2_2003, days = c(1:123))

files <- files[which(str_detect(Khamra_2001,'PIH.R'))]


Khamra_2004 <- PIH_in_list(data_1 = PIH_1_2004, data_2 = PIH_2_2004, days = c(1:123))

Khamra_2007 <- PIH_in_list(data_1 = PIH_1_2007, data_2 = PIH_2_2007, days = c(1:123))
Khamra_2013 <- PIH_in_list(data_1 = PIH_1_2013, data_2 = PIH_2_2013, days = c(1:123))
Khamra_2014 <- PIH_in_list(data_1 = PIH_1_2014, data_2 = PIH_2_2014, days = c(1:123))







# list_2_2001$`Results/PIH/open_data/2_2001-05-01_PIH_open.RData` <- NULL


# library(data.table)
# try_data <- rbindlist(Map(merge, list_1_2001, list_2_2001, by = "date"))


# try_function <- function(x, y, days){
#                 rbind.fill(x[days], y[days])
#                 }
# 
# Khamra_2014 <- try_function(x= list_1_2014, y  = list_2_2014, days = c(1:123))
# 
# for(o in PIH_1_2001){
#     o <- PIH_1_2001[90]
#     i <- PIH_2_2001[91]
#     get_data_1          <- get(load(o))
#     get_data_2          <- get(load(i))
#     try <- rbind.fill(get_data_1, get_data_2)
#     
# }
# 
# 
# 
# for(o in PIH_2_2001){
#   memory.limit(9999999999)
#   get_data          <- get(load(o))
#   list_2_2001[[o]]  <- get_data
# }
# 
# 
# 
# 
# data_PIH_frame_2001_1 <- do.call(rbind.data.frame, list_1_2001)
# 
# view(data_2001_1_frame)
# 
# for(o in PIH_2_2001){
#   get_data         <- get(load(o))
#   list_2_2001[[o]] <- get_data
#   rbind(list_2_2001)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # 
# # get_data_PIH <- lapply(1:lenght(PIH_1_2001), function(list){
# #     get_data <- get(load(o))
# #     mylist   <- list
# #  mylist[[o]] <- get_data
# #  }))
# 
# PIH_2001_1 <- get_data_PIH(o = PIH_1_2001, list = list_1_2001)
# PIH_2001_2 <- get_data_PIH()
# 
# 
# get(load('C:/Users/isfar.Rdata'))
# 
# 
# 
# 
# 
# 
# 
# PIH_1_2001 <- PIH_select_years("1_2002", PIH_1)
# PIH_1_2001
# PIH_select_years <- function(year, data){
#   Filter(function(e) grepl(year, e), data)
# }
# PIH_1_2001 <- PIH_select_years("1_2001", PIH_1)
# PIH_1_2001
# PIH_select_years <- function(year, data){
#   Filter(function(e) grepl(year, e), data)
# }
# PIH_select_years <- function(year, data){
#   Filter(function(e) grepl(year, e), data)
# }
# PIH_1      <- PIH_select_years("1_2", PIH_data)
# PIH_2      <- PIH_select_years("2_2", PIH_data)
# # single years
# PIH_1_2001 <- PIH_select_years("_2002", PIH_1)
# PIH_1_2001
# # single years
# PIH_1_2001 <- PIH_select_years("_2001", PIH_1)
# PIH_1_2001
# # single years
# 




