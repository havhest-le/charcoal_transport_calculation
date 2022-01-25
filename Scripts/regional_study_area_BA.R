###
# Script to map the regional area of lake Khamra and Satagay #
###

rm(list = ls(all= TRUE))

# load packages
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
library(plyr)
library(ggmap)
library(tmap)
library(plotly)

#####
# Load data ######################################################################
#####

# Khamra #########################################################################
Khamra                 <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp") %>% st_transform(4326)
proj                   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Khamra))[,1]}
                                +lat_0={st_coordinates(st_centroid(Khamra))[,2]}")
Khamra_buf_100         <- st_transform(Khamra, crs = CRS(proj)) %>%
                          st_buffer(100000) %>% st_transform(4326)
Khamra_data_frame      <- data.frame(lon = 59.99, lat =  112.98)

# Hydrology
Khamra_catchment       <- read_sf("Data/Khamra/khamra_catchment.shp") %>% st_transform(4326)
Khamra_streams         <- read_sf("Data/Khamra/streams.shp") %>% st_transform(4326)
Khamra_streams_new     <- Khamra_streams %>% filter(grid_code >= 5)

Khamra_depth           <- read_sf("Data/Khamra/KhamraWaterDepthAll.shp")
Khamra_cores           <- read_sf("Data/Khamra/Khamra_Cores.shp")
Khamra_bathemytrie     <- raster("Data/Khamra/Khamra_Bathymetry_test.tif")

# Digital elevation model
Khamra_dem             <- raster("Data/Khamra/fill_dem.tif")
Khamra_dem_proj        <- projectRaster(from = Khamra_dem, crs = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Satagay 2.0 ####################################################################
Satagay              <- read_sf("Data/Satagay/satagay_poly_1.shp") %>% st_transform(4326)
proj                 <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Satagay))[,1]}
                              +lat_0={st_coordinates(st_centroid(Satagay))[,2]}")
Satagay_buf_100      <- st_transform(Satagay, crs = CRS(proj)) %>%
                        st_buffer(100000) %>% st_transform(4326)

# Hydrology
Satagay_catchment    <- read_sf("Data/Satagay/satagay_catchment.shp") %>% st_transform(4326)
Satagay_streams      <- read_sf("Data/Satagay/satagay_streams.shp") %>% st_transform(4326)


# Digital elevation model
Satagay_dem_1       <- raster("Data/Satagay/TDM1_DEM__30_N63E118_DEM.tif")
Satagay_dem_2       <- raster("Data/Satagay/TDM1_DEM__30_N63E116_DEM.tif")
Satagay_dem_3       <- raster("Data/Satagay/TDM1_DEM__30_N62E118_DEM.tif")
Satagay_dem_4       <- raster("Data/Satagay/TDM1_DEM__30_N62E116_DEM.tif")


Satagay_dem_merge   <- merge(Satagay_dem_1, Satagay_dem_2, Satagay_dem_3, 
                             Satagay_dem_4,  tolerance = 0.5, overlap=TRUE, ext = NULL)

# Crop data with defined extent
Satagay_extent       <- extent(117.94, 118.13, 63.01, 63.10)
Satagay_dem_crop     <- crop(Satagay_dem_merge, Satagay_extent)
Satagay_extent_2     <- extent(117.94, 118.12, 63.01, 63.10)
Satagay_streams_crop <- st_crop(Satagay_streams, Satagay_extent_2)

#####
# Plotting ######################################################################
#####

# Satagay 2.0 ###################################################################

# DEM and catchment
plot(Satagay_dem_crop)
plot(Satagay_catchment$geometry, add = T)
plot(Satagay$geometry, add = T)
plot(Satagay_streams_crop$geometry, col = "dodgerblue3", add =T)

# Hillshade
slope_Sat     <- terrain(Satagay_dem_crop, opt='slope')
aspect_sat    <- terrain(Satagay_dem_crop, opt='aspect')
hillshade_sat <- hillShade(slope_Sat, aspect_sat, angle=45, direction=315)

# Plotting hillshade #
png(glue("Results/Satagay_hill.png"), width = 800, height = 700)
opar <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 10))

plot(hillshade_sat, col=grey.colors(100, start=0, end=1),legend=F)
plot(Satagay_catchment$geometry, lwd = 2, add = T)
plot(Satagay$geometry, col = "blue", add = TRUE)
plot(Satagay_streams_crop, col = "dodgerblue3", add = T)
title(main = "Hillshade of the catchment Lake Satagay 2.0",cex.main = 1.5)

legend("topright",
       inset = c(-0.26, 0),
       title = "Legend",
       title.adj = 0.5,
       legend = c("Lake Satagay", "Catchment"),
       pch = c(15, 0),
       col = c("blue", "black"),
       bty = "n",
       cex= 1.5,
       pt.cex= 2,
       xpd = TRUE)

north.arrow(xb=117.948, yb=63.09, len = 0.002, lab = "N")
map.scale(x=118.06, y=63.098, ratio=F, relwidth=0.2) 

dev.off()


# DEM and hillshade 
png(glue("Results/Satagay_dem_hill_2.png"), width = 800, height = 700)
opar <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 10))

plot(hillshade_sat,
     col = grey(1:100/100),
     legend = FALSE, xaxt = "n", yaxt = "n")

axis(1, at = c(117.95, 118.00, 118.05, 118.1))
axis(2, at = c(63.02, 63.04, 63.06, 63.08))

plot(Satagay_streams_crop, col = "deepskyblue1", lwd = 2.5, add = T)

plot(Satagay_dem_crop,
     main = "Lidar Digital Elevation Model (DEM)",
     add = TRUE, alpha = .5,col = grey(100:1/100),
     legend = F)

plot(Satagay_dem_crop, legend.only=TRUE, col = grey(100:1/100), legend.width=1, legend.shrink=0.9,
     smallplot=c(0.89,0.92, 0.1,0.7))
plot(Satagay_catchment$geometry, lwd = 2, add = T)
plot(Satagay$geometry, col = "blue", add = TRUE)
title(main = "Digital Elevation Model (DEM) overlayed on top \n of the hillshade of catchment Lake Satagay 2.0",cex.main = 1.5)

legend("topright",
       inset = c(-0.26, 0),
       title = "Legend",
       title.adj = 0.5,
       legend = c("Lake Satagay", "Catchment"),
       pch = c(15, 0),
       col = c("blue", "black"),
       bty = "n",
       cex= 1.5,
       pt.cex= 2,
       xpd = TRUE)

legend("topright",
       inset = c(-0.2, 0.19),
       legend = "DEM [m]",
       angle = 45,
       bty = "n",
       cex= 1.5,
       xpd = TRUE)

north.arrow(xb=118.1105, yb=63.09, len = 0.002, lab = "N")
map.scale(x=118.055, y=63.098, ratio=F, relwidth=0.2) 

dev.off()


# Create stack of dem and hillshade 
stack_dem_hill_Sat                <- stack(Satagay_dem_crop,hillshade_sat)
# Rasterstack to points
stack_dem_hill_points_Sat         <- rasterToPoints(stack_dem_hill)
# Write a data frame
stack_dem_hill_data_Sat           <- data.frame(stack_dem_hill_points)
colnames(stack_dem_hill_data_Sat) <- c("X","Y","DEM", "hillshade")

# function ploty() #
# DEM Satagay
p_Satagay <- plot_ly(x = stack_dem_hill_data_Sat$X, y = stack_dem_hill_data_Sat$Y, z = stack_dem_hill_data_Sat$DEM, type = "contour",
        colors = "Oranges") %>% colorbar(title = "DEM [in m]")  %>% 
        layout(title = 'DEM of area Lake Satagay 2.0', 
               plot_bgcolor = "#e5ecf6")
htmlwidgets::saveWidget(as_widget(p_Satagay), "Results/p_Satagay.html")
# conotur lines Lake Satagay
png(glue("Results/Satagay_contour.png"), width = 900, height = 800)
opar <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 10))

plot(hillshade_sat,
     col = grey(1:100/100),
     legend = FALSE, xaxt = "n", yaxt = "n")

axis(1, at = c(117.95, 118.00, 118.05, 118.1))
axis(2, at = c(63.02, 63.04, 63.06, 63.08))

plot(Satagay_streams_crop, col = "deepskyblue1", lwd = 2.5, add = T)
cols <- hcl.colors(1, "BrwnYl")
contour(Satagay_dem_crop,method = "flattest",
        col = "grey60", add = T) 
plot(Satagay_catchment$geometry, lwd = 2, add = T)
plot(Satagay$geometry, col = "blue", add = TRUE)
title(main = "Hillshade with contour lines of catchment Lake Satagay 2.0",cex.main = 1.5)

legend("topright",
       inset = c(-0.26, 0),
       title = "Legend",
       title.adj = 0.5,
       legend = c("Lake Satagay", "Catchment"),
       pch = c(15, 0),
       col = c("blue", "black"),
       bty = "n",
       cex= 1.5,
       pt.cex= 2,
       xpd = TRUE)

north.arrow(xb=118.1105, yb=63.09, len = 0.002, lab = "N")
map.scale(x=118.055, y=63.098, ratio=F, relwidth=0.2) 

dev.off()


# Khamra ##########################################################################

# DEM and catchment
Khamra_area_crop    <- extent(112.7, 113.08, 59.9, 60.03)
Khamra_dem_crop     <- crop(Khamra_dem_proj, Khamra_area_crop)
Khamra_area_crop_2  <- extent(112.7, 113.06, 59.9, 60.03)
Khamra_streams_crop <- st_crop(Khamra_streams_new,Khamra_area_crop_2) %>% st_transform(4326)

# Hillshade
slope              <- terrain(Khamra_dem_crop, opt='slope')
aspect             <- terrain(Khamra_dem_crop, opt='aspect')
hillshade_kahm     <- hillShade(slope, aspect, angle=45, direction=315)

# Plotting hillshade
plot(hillshade_kahm, col=grey.colors(100, start=0, end=1), legend=F)
plot(Khamra_catchment$geometry, add = T)
plot(Khamra$geometry, col = "blue", add = TRUE)
plot(Khamra_streams_crop$geometry,col = "blue", add = T)

# Create stack of dem and hillshade 
stack_dem_hill_Kham                <- stack(Khamra_dem_crop,hillshade_kahm)
# Rasterstack to points
stack_dem_hill_points_Kham         <- rasterToPoints(stack_dem_hill_Kham)
# Write a data frame
stack_dem_hill_data_Kham           <- data.frame(stack_dem_hill_points_Kham)
colnames(stack_dem_hill_data_Kham) <- c("X","Y","DEM", "hillshade")

# function ploty() #
# DEM Khamra
p_Khamra <- plot_ly(x = stack_dem_hill_data_Kham$X, y = stack_dem_hill_data_Kham$Y, z = stack_dem_hill_data_Kham$DEM, type = "contour",
            colors = "Oranges") %>% colorbar(title = "DEM [in m]")  %>% 
            layout(title = 'DEM of area Khamra', 
            plot_bgcolor = "#e5ecf6")
htmlwidgets::saveWidget(as_widget(p_Khamra), "Results/p_Khamra.html")

# DEM and hillshade 
png(glue("Results/Khamra_dem_hill_2.png"), width = 850, height = 600)
opar <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 10))

plot(hillshade_kahm,
     col = grey(1:100/100),
     legend = FALSE, xaxt = "n", yaxt = "n")

axis(1, at = c(112.75, 112.85, 112.95, 113.05))
axis(2, at = c(59.90, 59.95, 60, 60.05))

plot(Khamra_streams_crop$geometry,col = "deepskyblue1", lwd = 2.5, add = T)
plot(Khamra_dem_crop,
     main = "Lidar Digital Elevation Model (DEM)",
     add = TRUE, alpha = .5,col = grey(100:1/100),
     legend = F)

fun_color_range <- colorRampPalette(c("darkblue", "darkslategray1"))
my_colors <- fun_color_range(11)
fun_color_range_2 <- colorRampPalette(c("darkslategray1", "darkblue"))
my_colors_2 <- fun_color_range_2(11)

plot(Khamra_dem_crop, legend.only=TRUE, col = grey(100:1/100), legend.width=1, legend.shrink=0.75,
     smallplot=c(0.87,0.9, 0.5,0.75))
plot(Khamra_dem_crop, legend.only=TRUE, legend.width=1, legend.shrink=0.75,
     smallplot=c(0.87,0.9, 0.1,0.45), col = my_colors_2)

plot(Khamra_bathemytrie, add = T, col = my_colors, legend = F)
plot(Khamra_catchment$geometry, lwd = 2, add = T)
title(main = "Digital Elevation Model (DEM) overlayed on top \n of the hillshade of catchment Lake Khamra",
      cex.main = 1.5)


legend("topright",
       inset = c(-0.17, 0),
       title = "Legend",
       title.adj = 0.5,
       legend = "Catchment",
       pch = 0,
       col = "black",
       bty = "n",
       cex= 1.5,
       pt.cex= 2,
       xpd = TRUE)

legend("topright",
       inset = c(-0.13, 0.135),
       legend = "DEM [m]",
       angle = 45,
       bty = "n",
       cex= 1.2,
       xpd = TRUE)

legend("topright",
       inset = c(-0.17, 0.52),
       legend = "Bathymetry [m]",
       angle = 45,
       bty = "n",
       cex= 1.2,
       xpd = TRUE)

north.arrow(xb=112.71, yb=60.013, len = 0.004, lab = "N")
map.scale(x=112.725, y=60.025, ratio=F, relwidth=0.1) 

dev.off()

###############################################################################
