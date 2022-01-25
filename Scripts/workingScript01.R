rm(list = ls(all= TRUE))


library(sp)
library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(ncdf4)
library(rnaturalearth)
library(tidyverse)
library(glue)
library(mapview)
library(units)

####
# Input and clean data ####
####



Elgy <- read_sf("Data/Lakes/Elgygytgen/Lake_area.shp")
plot(Elgy$geometry, col  = "lightblue")
buffer <- st_transform(Elgy, crs = CRS("+proj=laea")) %>%
  st_buffer(400000) %>% st_transform(3857)
plot(buffer$geometry, add = T)

Poly_buffer <- st_buffer(x = Elgy, dist = rad)
plot(Poly_buffer, add = TRUE)

Khamra <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp")
plot(Khamra$geometry, col = "lightblue")

Ill <- read_sf("Data/Lakes/Illerney/Umriss_Polygon_Ilriney.shp")
plot(Ill$geometry, col = "lightblue")
