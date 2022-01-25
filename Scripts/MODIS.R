rm(list = ls(all= TRUE))


library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(sf)
library(ncdf4)
library(rnaturalearth)
library(tidyverse)
library(data.table)
library(glue)


#####
# Input data
#####

map     <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
ext     <- extent(c(103.82, 180, 50.07, 80.56))

MODIS <- read.csv2("N:/bioing/data/Projekte/Fires/East Siberia/R_Data/MODIS_east_siberia.txt")

MODIS <- select(.data = MODIS,"FID", "LATITUDE","LONGITUDE", "ACQ_DATE", "FRP")
names(MODIS)[names(MODIS) == "ACQ_DATE"] <- "TIME"
MODIS$TIME <- as.POSIXct(x = MODIS$TIME, format = '%d.%m.%Y %H:%M:%S')
MODIS$FRP <- as.numeric(MODIS$FRP)
MODIS$Year <- as.numeric(format(MODIS$TIME, "%Y"))
MODIS_ext <- MODIS %>%
  filter(LONGITUDE >=103.82 & LONGITUDE <= 180 & LATITUDE >= 50.07  & LATITUDE <= 80.56 & FRP > 0)
head(MODIS_ext)

MODIS <- data.table(MODIS_ext)
save(MODIS, file = "Results/MODIS.RData")
