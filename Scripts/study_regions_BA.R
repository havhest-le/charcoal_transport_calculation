#####
## Script for visualization the study areas ####################################
#####

rm(list = ls(all= TRUE))

######
## Load the packages ###########################################################
######

library(sp)
library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(ncdf4)
library(rnaturalearth)
library(tidyverse)
library(mapview)
library(units)
library(lattice)
library(rasterVis)
library(ggplot2)
library(GISTools)  
library(tmap)
library(maps)
library(glue)

#####
## Load the data ###############################################################
#####

## Russia
russia            <- read_sf("Data/East_Siberia/gadm36_RUS_1.shp") %>% st_transform(4326)
russia_part       <- st_crop(russia, xmin = 100, xmax = 172, ymin = 77.10, ymax = 45) %>% st_transform(4326)

## Select Yakutia
yakutia           <- russia %>% filter(russia$NAME_1 == "Sakha")
save(yakutia, file = "Data/Jakutien/yakutia.shp")
yakutia_WGS       <- st_transform(yakutia, crs = st_crs(4326))

treeline_global   <- read_sf("Data/Vegetation/treeline_wgs84.shp")
treeline_siberia  <- st_intersection(treeline_global, yakutia)
veget             <- raster("Data/Borealer_Wald/boreal_forest_yakutia.tif")


## Create a data frame with coordinates
data_coord_lakes  <- data.frame(location = c("Khamra","Satagay"),
                                lon = c(112.98, 117.998), 
                                lat = c(59.99, 63.078)) 

yakustk_coord     <- data.frame(location = c("Yakusk"),
                                lon = c(129.44),
                                lat = c(62.20))

## Lake Khamra
Khamra      <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp") %>% st_transform(4326)
proj        <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Khamra))[,1]}
                     +lat_0={st_coordinates(st_centroid(Khamra))[,2]}")

## Create a 100 km radius around the lake as a center point 
Khamra_buf_100  <- st_transform(Khamra, crs = CRS(proj)) %>%
                   st_buffer(100000) %>% st_transform(4326)
Khamra_buf_50   <- st_transform(Khamra, crs = CRS(proj)) %>%
                   st_buffer(50000) %>% st_transform(4326)

## Lake Satagay 
Satagay     <- read_sf("Data/Satagay/satagay_poly_1.shp") %>% st_transform(4326)
proj        <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Satagay))[,1]}
                     +lat_0={st_coordinates(st_centroid(Satagay))[,2]}")

## Create a 100 km radius around the lake as a center point 
Satagay_buf_100 <- st_transform(Satagay, crs = CRS(proj)) %>%
                   st_buffer(100000) %>% st_transform(4326)

Satagay_buf_50  <- st_transform(Satagay, crs = CRS(proj)) %>%
                   st_buffer(50000) %>% st_transform(4326)


## Create study area extent as SpatialPolygons
Khamra_buf_extent  <- as(extent(st_bbox(Khamra_buf_100 %>% st_transform(4326) 
                                        %>% st_shift_longitude())[c(1,3,2,4)]), 
                                        "SpatialPolygons")

Satagay_buf_extent <- as(extent(st_bbox(Satagay_buf_100 %>% st_transform(4326) 
                                        %>% st_shift_longitude())[c(1,3,2,4)]), 
                                        "SpatialPolygons")


#####
# Visualization ################################################################
#####

## Vegetation reclassification
reclass_df                              <- c(1, 2, 1, 2, 4, 2)
reclass_m                               <- matrix(reclass_df, ncol = 3, byrow = TRUE)
veget_classified                        <- raster::reclassify(veget, reclass_m)
veget_classified[veget_classified == 0] <- NA
save(veget_classified, file = "Data/Borealer_Wald/veget_classified.tif")
colors                                  <- c("darkseagreen2", "darkgreen")
col                                     <- c("yellow1","red2")


png(glue("Results/study_area/study_area_yakutia_new.png"), width = 1000, height = 700)

# Make the window wider than taller
# windows(width = 4.5, height = 4)
# Save current graphical parameters
opar <- par(no.readonly = TRUE)
 
# Change the margins of the plot (the fourth is the right margin)
par(mar = c(5, 5, 4, 20))

plot(veget_classified,
     legend = FALSE,
     col = colors,
     box = FALSE)

plot(yakutia$geometry, lwd= 2, add =T)
plot(russia$geometry, add = T)
plot(treeline_siberia$geometry, add = T, col = "aquamarine4")
points(x = data_coord_lakes$lon, y = data_coord_lakes$lat, pch = c(19, 19), fill = col, col = col, cex = 3, lwd = 3)
points(x = yakustk_coord$lon, y = yakustk_coord$lat, pch = 4, col = "darkslategray", cex = 3, lwd = 2)

# Legend
legend("topright",
       inset = c(-0.26, 0),
       title = "     Legend",
       title.adj = 0.5,
       legend = c("   Lake Khamra", "   Lake Satagay"),
       pch = c(19, 19),
       col = col,
       #fill = col,
       #border = col,
       bty = "n",
       cex= 1.5,
       pt.cex= 2,
       xpd = TRUE)

legend("topright",
       inset = c(-0.19, 0.119),
       legend = c("   Yakutsk"),
       pch = 4,
       cex = 1.5, 
       col = "darkslategray", 
       bty = "n",
       pt.cex= 2,
       xpd = TRUE)

legend("topright",
       inset = c(-0.2, 0.16),
       legend = c("Treeline"),
       cex = 1.5, 
       col = "aquamarine4", 
       lwd = 2, 
       bty = "n",
       xpd = TRUE)

legend("topright",
       inset = c(-0.38, 0.24),
       title = "Forest type",
       title.adj = 0.5,
       legend = c("Summergreen needleleaf", "Everygreen needleleaf"),
       fill = colors,
       border = colors,
       bty = "n",
       cex= 1.5,
       xpd = TRUE)


# text(119.3, 64.2, "Satagay", cex = 1)
# text(115.7, 59.2, "Khamra", cex = 1)
# text(132.4, 62.2, "Yakutsk", cex = 1)

# North arrow and scale 
north.arrow(xb=106.5, yb=56.3, len = 0.4, lab = "N")
map.scale(x=146.8, y=56.5, ratio=F, relwidth=0.1) 

dev.off()



#####
## Visualization of each single lake with the 100 km buffer ####################
#####

# Lake Khamra 
plot(Khamra_buf_100$geometry, col = "cyan")
plot(Khamra$geometry, col = "darkblue", add = T)

# Lake Satagay
plot(Satagay_buf_100$geometry,col = "cyan")
plot(Satagay$geometry, col = "darkblue", add = T)


lakes <- list(Khamra  = list(lake = Khamra,  buffer = Khamra_buf_100),
              Satagay = list(lake = Satagay, buffer = Satagay_buf_100))

lakes_roh <-  list(Khamra = Khamra$geometry, Satagay =  Satagay$geometry) 
lakes_buf <-  list(Khamra = Khamra_buf_100$geometry, Satagay =  Satagay_buf_100$geometry)

save(lakes, file = "Results/lakes_Kham_Sat.RData")
save(lakes_roh, file = "Results/lakes_roh_Kham_Sat.RData")
save(lakes_buf, file = "Results/lakes_buf_Kham_Sat.RData")


lakes_sf <- do.call("rbind", lapply(1:2, function(x) 
            data.frame(st_coordinates(lakes[[x]][[1]])[,1:2], lake = names(lakes)[x]))) %>%
            st_as_sf(coords = c("X", "Y"), dim = c("XY")) %>% st_set_crs(st_crs(lakes[[1]][[1]])) %>%
            group_by(lake) %>% summarise(do_union = F) %>% st_cast("POLYGON")

st_write(lakes_sf, dsn = "Data/Lakes/lakes_sf.gbd",layer="lakes_sf", driver="ESRI Shapefile")



#####
## New study area map ##########################################################
#####
library(GISTools)
library(maps)


png(glue("Results/study_area/study_area_yakutia_new.png"), width = 1000, height = 700)

# Make the window wider than taller
# windows(width = 4.5, height = 4)

# Save current graphical parameters
opar <- par(no.readonly = TRUE)

# Change the margins of the plot (the fourth is the right margin)
par(mar = c(5, 5, 4, 10))

plot(veget_classified,
     legend = FALSE,
     col = colors,
     box = FALSE)

plot(yakutia$geometry, lwd= 2, add =T)
plot(russia$geometry, add = T)
plot(Khamra_buf_extent, col = alpha("violetred", 0.3), add = T)
plot(Satagay_buf_extent, col = alpha("violetred", 0.3), add = T)
points(x = data_coord_lakes$lon, y = data_coord_lakes$lat, pch = c(19, 19), fill = "maroon1", col = "maroon1", cex = 1, lwd = 2)
points(x = yakustk_coord$lon, y = yakustk_coord$lat, pch = 4, col = "darkslategray", cex = 2, lwd = 2)

# Legen
legend("topright",
       inset = c(0.0005, 0.05),
       legend = c("Summergreen needleleaf", "Everygreen needleleaf"),
       fill = colors,
       border = colors,
       bty = "n",
       cex= 1.4,
       xpd = TRUE)


legend("topright",
       inset = c(-0.035, 0.015),
       legend = "100 km radius                     ",
       fill = "violetred",
       border = "violetred",
       bty = "n",
       bg="transparent",
       cex= 1.4,
       xpd = TRUE)

text(117.4, 64.65, "Lake\n     Satagay", cex = 1.2)
text(112.7, 61.55, "Lake\n     Khamra", cex = 1.2)
text(132.5, 62.25, "Yakutsk", cex = 1.2)

# North arrow and scale 
GISTools::north.arrow(xb=102, yb=56.3, len = 0.4, lab = "N")
maps::map.scale(x=142, y=56.5, ratio=F, relwidth=0.1) 

dev.off()
