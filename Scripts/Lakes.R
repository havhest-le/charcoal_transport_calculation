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


#####
# Mapping 
#####

map     <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
ext     <- extent(c(103.82, 180, 50.07, 80.56))
mapCrop <- map %>% 
  st_intersection(st_as_sf(as(ext, "SpatialPolygons")) %>% 
                      st_set_crs(4326)) %>% st_geometry 

# Elgygytgen
Elgy <- read_sf("Data/Lakes/Elgygytgen/Lake_area.shp") %>% st_transform(4326)
proj <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Elgy))[,1]} +lat_0={st_coordinates(st_centroid(Elgy))[,2]}")
Elgy_buf <- st_transform(Elgy, crs = CRS(proj)) %>%
  st_buffer(100000) %>% st_transform(4326)
plot(Elgy$geometry, col = "blue")
# with 400 km buffer
plot(Elgy_buf$geometry, add = T)

# Khamra
Khamra <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp") %>% st_transform(4326)
proj <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Khamra))[,1]} +lat_0={st_coordinates(st_centroid(Khamra))[,2]}")
Kha_buf <- st_transform(Khamra, crs = CRS(proj)) %>%
    st_buffer(100000) %>% st_transform(4326)
plot(Khamra$geometry, col = "blue", add = T)
# with 400 km buffer
plot(Kha_buf$geometry, add = T)


# Illerney
Ill <- read_sf("Data/Lakes/Illerney/Umriss_Polygon_Ilriney.shp") %>% st_transform(4326)
proj <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Ill))[,1]} +lat_0={st_coordinates(st_centroid(Ill))[,2]}")
Ill_buf <- st_transform(Ill, crs = CRS(proj)) %>%
  st_buffer(100000) %>% st_transform(4326)
plot(Ill$geometry, col = "blue", add = T)
# with 400 km buffer
plot(Ill_buf$geometry, add = T)



#####
# Shapefile lakes
#####

lakes <- list(Elgy = list(lake = Elgy, buffer = Elgy_buf),
              Khamra = list(lake = Khamra, buffer = Kha_buf),
              Ill = list(lake = Ill, buffer = Ill_buf))

lakes_roh <-  list(Elgy = Elgy$geometry, Khamra =  Khamra$geometry, Ill =  Ill$geometry) 
lakes_buf <- list(Elgy = Elgy_buf$geometry, Khamra =  Kha_buf$geometry, Ill =  Ill_buf$geometry)

save(lakes, file = "Results/lakes.RData")
save(lakes_roh, file = "Results/lakes_roh.RData")
save(lakes_buf, file = "Results/lakes_buf.RData")



lakes_sf <- do.call("rbind", lapply(1:3, function(x) 
  data.frame(st_coordinates(lakes[[x]][[1]])[,1:2], lake = names(lakes)[x]))) %>%
  st_as_sf(coords = c("X", "Y"), dim = c("XY")) %>% st_set_crs(st_crs(lakes[[1]][[1]])) %>%
  group_by(lake) %>% summarise(do_union = F) %>% st_cast("POLYGON")

save(lakes_sf, file = "Data/lakesSHP/lakes_sf.shp")

