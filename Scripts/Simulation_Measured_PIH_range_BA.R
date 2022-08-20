#####
## Script to simulate potential and actual modern source areas of charcoal using 
## observed plume injection heights between the range of Q[0.2] and Q[0.8] #####
#####


rm(list = ls(all= TRUE))

#####
## Load packages ###############################################################
#####

library(tidyverse)
library(sf)
library(terra)
library(stars)
library(dplyr)
library(glue)
library(base)
library(raster)

## Plotting
library(maps)
library(viridis)
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(ggnewscale)
library(MASS)
library(survival)
library(fitdistrplus)
library(geosphere)

#####
## Load and clean data #########################################################
#####

## Lakes
data_coord_lakes        <- data.frame(location = c("Khamra","Satagay"),
                                      lon = c(112.98, 117.998), 
                                      lat = c(59.99, 63.078))
## Khamra
Khamra                  <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp") %>% st_transform(4326)
proj_K                  <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Khamra))[,1]}
                                            +lat_0={st_coordinates(st_centroid(Khamra))[,2]}")
Khamra_catchment        <- read_sf("Data/Khamra/khamra_catchment.shp") %>% st_transform(4326)
Khamra_buf_100          <- st_transform(Khamra, crs = st_crs(proj_K))  %>%  st_buffer(100000) %>% st_transform(4326) 
Khamra_buf_extent       <- as(extent(st_bbox(Khamra_buf_100 %>% 
                                               st_transform(4326)  %>% 
                                               st_shift_longitude())[c(1,3,2,4)]), "SpatialPolygons")
crs(Khamra_buf_extent)  <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

## Satagay 2.0
Satagay                 <- read_sf("Data/Satagay/satagay_poly_1.shp") %>% st_transform(4326) 
Satagay$InPoly_FID      <- NULL
Satagay$MaxSimpTol      <- NULL
Satagay$MinSimpTol      <- NULL
Satagay$SimPgnFlag      <- NULL
proj_S                  <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Satagay))[,1]}
                                            +lat_0={st_coordinates(st_centroid(Satagay))[,2]}")
Satagay_buf_100         <- st_transform(Satagay, crs = st_crs(proj_S)) %>%
  st_buffer(100000) %>% st_transform(4326)
Satagay_catchment       <- read_sf("Data/Satagay/satagay_catchment.shp") %>% st_transform(4326)
Satagay_buf_extent      <- as(extent(st_bbox(Satagay_buf_100 %>% 
                                               st_transform(4326) %>% 
                                               st_shift_longitude())[c(1,3,2,4)]), "SpatialPolygons")
crs(Satagay_buf_extent) <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

PIH_range_Khamra  <- get(load("Results/Simulation/data/PIH_range_Khamra.rda"))
PIH_range_Satagay <- get(load("Results/Simulation/data/PIH_range_Satagay.rda"))

#####
## Simulation for lake Khamra ###### ###########################################
#####
# 
# ## Wind data
wndTab_Khamra <- get(load("Results/Simulation/data/wndTab_Khamra.rda"))

## FRP and PIH data
# outMODIS_PIH_Kahmra_geom         <- get(load("Results/Simulation/data/outMODIS_PIH_Kahmra_geom.rda"))
# outMODIS_PIH_Kahmra_geom_na      <- na.omit(outMODIS_PIH_Kahmra_geom)
# outMODIS_PIH_Kahmra_geom_na$TIME <- as.Date(outMODIS_PIH_Kahmra_geom_na$TIME)

# #####
# ## Creating a tibble with FRP, PIH and wind ####################################
# #####
# 
# tibble_frp_pih_wind_Khamra <- do.call("rbind", lapply(unique(outMODIS_PIH_Kahmra_geom_na$TIME), function(d) {
# 
#   cat(sprintf('\rDate %s of %s  ', which(d == unique(outMODIS_PIH_Kahmra_geom_na %>%
#                                                        pull(TIME) %>% as.Date())), length(unique(outMODIS_PIH_Kahmra_geom_na %>% pull(TIME) %>% as.Date()))))
# 
#   wndTab_Khamra$date <- as.Date(wndTab_Khamra$date)
#   tmp_wind    <- wndTab_Khamra %>% group_by(date) %>% filter(date == d) %>%
#     dplyr::select("date", "lon", "lat", "wind_spd", "wind_dir")
#   tmp_wind_2  <- st_as_sf(tmp_wind, coords = c("lon", "lat")) %>% st_set_crs(4326) # define coordinates
#   tmp_frp_pih <- outMODIS_PIH_Kahmra_geom_na %>% group_by(TIME) %>% filter(TIME == d)
#   tmp_frp_pih <- st_as_sf(tmp_frp_pih, crs = st_crs(4326))
# 
#   # Calculating the distance between tmp2 and x
#   distM <- st_distance(tmp_frp_pih, tmp_wind_2, by_element = F)
# 
#   # Selecting the minimum distance between PIH pixel to FRP pixel
#   ind <- apply(distM, 1, function(y) which.min(y))
# 
#   # Creating a tibble of coordinates, FRP, PIH and dist
#   dist <- sapply(1:length(ind), function(z) distM[z, ind[z]]/1000)
# 
#   df <- tmp_frp_pih %>% mutate(w_spd = tmp_wind_2$wind_spd[ind],
#                                w_dir = tmp_wind_2$wind_dir[ind],
#                                dist = as.numeric(dist))
# 
# }))
# save(tibble_frp_pih_wind_Khamra, file = "Results/Simulation/data/tibble_frp_pih_wind_Khamra.rda")
tibble_frp_pih_wind_Khamra  <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Khamra.rda"))

## Convert MW to cal/s
tibble_frp_pih_wind_Khamra$FRP_cal <- tibble_frp_pih_wind_Khamra$FRP*238845.8966275

#####
## PIH values | range between Q[0.2] and Q[0.8] ################################
#####

#####
## Calculation of the travel distance of spherical particle ####################
#####

## Cropping tibble with 100 km radius box around the lake
tibble_frp_pih_wind_Khamra <- st_crop(tibble_frp_pih_wind_Khamra, Khamra_buf_extent)

## Plotting empirical_density for measured PIH
png(glue("Results/Simulation/PIH/Khamra/range/Empirical_density_PIH_range_Khamra.png"), width = 1000, height = 700)
plotdist(as.numeric(PIH_range_Khamra), histo = TRUE, demp = TRUE, )
dev.off()

## Gamma distribution of measured PIH
PIH_distr_Khamra_range <- MASS::fitdistr(as.numeric(PIH_range_Khamra), densfun = "gamma")

#####
## Function of the potential modern origin of charcoal #########################
#####

origin  <- function(spd, dir, pih, diam, lon_orig, lat_orig){
  vt <- (((0.5-0.00127)*981*diam^2)/(18*(1.8*10^-4)))/100 # cm/s in m/s
  t  <- pih/(vt)
  delta_x <- t*spd
  geosphere::destPoint(c(lon_orig, lat_orig),
                       (dir+180)%%360, delta_x)
}

## Creating a charcoal diameter range from 150 to 500 micrometer based on the
## lake sediment samples analyzed by Glückler et al.
charc_diam_cm <- seq(0.015, 0.05, length = 100)

## Calculation of the travel distance ##########################################

sumCrds_Khamra_PIH_range <- do.call("rbind", lapply(unique(wndTab_Khamra$date), function(d) {

  cat(sprintf('\rDate %s of %s  ', which(d == unique(wndTab_Khamra %>% pull(date) %>% as.Date())),
              length(unique(wndTab_Khamra %>% pull(date) %>% as.Date()))))

  tmp <- wndTab_Khamra %>% filter(date == d) %>% ungroup() %>% st_as_sf(coords = c("lon", "lat")) %>%
         st_set_crs(4326) %>% dplyr::select(wind_spd, wind_dir)

  tmp <- tmp[apply(st_distance(tmp, data_coord_lakes       %>%
                                 st_as_sf(coords = c("lon", "lat"))                %>%
                                 st_set_crs(4326)), 2, function(x) which.min(x)),] %>%
                                 st_drop_geometry()

  t<- do.call("rbind", lapply(1:20, function(r) {
    tibble(lake = 1, as.data.frame(origin(spd  = tmp$wind_spd[1],
                                          dir  = tmp$wind_dir[1],
                                          pih  = rgamma(1, PIH_distr_Khamra_range$estimate[1], PIH_distr_Khamra_range$estimate[2]),
                                          diam = charc_diam_cm[sample(1:length(charc_diam_cm),1)],
                                          lon_orig = data_coord_lakes$lon[1],
                                          lat_orig = data_coord_lakes$lat[1]))) %>% setNames(c("lake", "lon", "lat"))
  })) %>% dplyr::mutate(date = d) %>% dplyr::select(date, lake, lon, lat)
}))
save(sumCrds_Khamra_PIH_range, file = "Results/Simulation/data/sumCrds_Khamra_PIH_range.rda")
sumCrds_Khamra_PIH_range <- get(load("Results/Simulation/data/sumCrds_Khamra_PIH_range.rda"))

#####
## Calculating the distance between the coordinates, put it into tibble and 
## convert to points
#####

dist_data_Khamra_PIH_range <- sumCrds_Khamra_PIH_range %>% filter(!is.na(lon) & !is.na(lat)) %>% st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>%
  group_split(lake) %>% lapply(function(x) {
    x %>% mutate(dist = as.numeric(st_distance(geometry, (Khamra_buf_100 %>% st_centroid() %>% st_transform(4326))[unique(lake),])))
  }) %>% do.call("rbind", .)

save(dist_data_Khamra_PIH_range, file = "Results/Simulation/data/dist_data_Khamra_PIH_range.rda")
dist_data_Khamra_PIH_range <- get(load("Results/Simulation/data/dist_data_Khamra_PIH_range.rda"))

#####
## Plotting the travel distance of spherical particle in form of histogram #####
##### 

## 100 breaks
png(glue("Results/Simulation/PIH/Khamra/range/dist/travel_distance_Khamra_PIH_4000_100_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Khamra_PIH_range$dist[dist_data_Khamra_PIH_range$lake==1 & dist_data_Khamra_PIH_range$dist<400000]/1000, breaks= 100, 
     xlim = c(0,20), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks
png(glue("Results/Simulation/PIH/Khamra/range/dist/travel_distance_Khamra_PIH_2000_200_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Khamra_PIH_range$dist[dist_data_Khamra_PIH_range$lake==1 & dist_data_Khamra_PIH_range$dist<400000]/1000, breaks= 200, 
     xlim = c(0,20), ylim = c(0, 3000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 100 breaks | 40 km
png(glue("Results/Simulation/PIH/Khamra/range/dist/travel_distance_Khamra_PIH_4000_100_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Khamra_PIH_range$dist[dist_data_Khamra_PIH_range$lake==1 & dist_data_Khamra_PIH_range$dist<400000]/1000, breaks= 100, 
     xlim = c(0,40), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks | 40 km
png(glue("Results/Simulation/PIH/Khamra/range/dist/travel_distance_Khamra_PIH_2000_200_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Khamra_PIH_range$dist[dist_data_Khamra_PIH_range$lake==1 & dist_data_Khamra_PIH_range$dist<400000]/1000, breaks= 200, 
     xlim = c(0,40), ylim = c(0, 3000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

#####
## Plotting a boxplot of calculated PIH ########################################
#####

n_dist_K_range <- length(dist_data_Khamra_PIH_range$dist) # 14760
dist_data_Khamra_PIH_range$dist <- dist_data_Khamra_PIH_range$dist/1000
summary(dist_data_Khamra_PIH_range$dist) ## in km 
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##   0.00118  0.56387   1.15756   2.13798   2.48253  46.51162 
median_dist_Khamra_PIH_range <- summary(dist_data_Khamra_PIH_range$dist)[3]

## Boxplot
## 20 km
png(glue("Results/Simulation/PIH/Khamra/range/dist/travel_distance_Khamra_PIH_range_boxplot_20km.png"), width = 900, height = 800)
boxplot(x = dist_data_Khamra_PIH_range$dist, outline = F, ylim = c(0,20),
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra", 
        ylab = "Δx [km]", xlab = glue("n = {n_dist_K_range}"), cex.main=1.7, cex.lab=1.5, cex.axis=1.5, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 12 km
png(glue("Results/Simulation/PIH/Khamra/range/dist/travel_distance_Khamra_PIH_range_boxplot_12km.png"), width = 900, height = 800)
boxplot(x = dist_data_Khamra_PIH_range$dist, outline = F, ylim = c(0,12),
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra", 
        ylab = "Δx [km]", xlab = glue("n = {n_dist_K_range}"), cex.main=1.7, cex.lab=1.5, cex.axis=1.5, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 8 km 
png(glue("Results/Simulation/PIH/Khamra/range/dist/travel_distance_Khamra_PIH_range_boxplot_8km.png"), width = 900, height = 800)
boxplot(x = dist_data_Khamra_PIH_range$dist, outline = F, ylim = c(0,8),
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra", 
        ylab = "Δx [km]", xlab = glue("n = {n_dist_K_range}"), cex.main=1.7, cex.lab=1.5, cex.axis=1.5, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## Cropping distance tibble with 100 km radius box around the lake
dist_data_crop_Khamra_PIH_range  <- st_crop(dist_data_Khamra_PIH_range, Khamra_buf_extent)

## Selecting a FRP-range
MODIS_Khamra_PIH_range  <- tibble_frp_pih_wind_Khamra %>% filter(FRP_cal > 0 & FRP_cal <= 71653768.9882488) 
## 300 MW = 71653768.9882488 cal/s

#####
## Plotting the backward simulation of potential source areas of charcoal particles
#####

## FRP (cal/s)
png(glue("Results/Simulation/PIH/Khamra/range/psource/Khamra_psource_PIH_range_cal.png"), width = 1000, height = 900)
ggplot(Khamra_buf_extent) +
  geom_sf(data = dist_data_crop_Khamra_PIH_range, aes(colour  = "Locations of\np-source areas\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Khamra_PIH_range, aes(color = FRP_cal), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(title = "Backward simulation of potential source areas\nof charcoal particles in the study area Lake Khamra",
       color = "FRP [cal/s]") +
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.5) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA, alpha = 0.3))) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[1]], y = lat[[1]],
                           shape = location[[1]]), data = data_coord_lakes, colour = "turquoise1",
             size = 3, stroke = 1.5) +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 20, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()


## FRP (MW)
MODIS_Khamra_PIH_range  <- tibble_frp_pih_wind_Khamra %>% filter(FRP > 0 & FRP <= 300) 

png(glue("Results/Simulation/PIH/Khamra/range/psource/Khamra_psource_PIH_range_MW.png"), width = 1000, height = 900)
ggplot(Khamra_buf_extent) +
  geom_sf(data = dist_data_crop_Khamra_PIH_range, aes(colour  = "Locations of\np-source areas\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Khamra_PIH_range, aes(color = FRP), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(title = "Backward simulation of potential source areas\nof charcoal particles in the study area Lake Khamra",
       color = "FRP [MW]") +
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.5) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA, alpha = 0.3))) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[1]], y = lat[[1]],
                           shape = location[[1]]), data = data_coord_lakes, colour = "turquoise1",
             size = 3, stroke = 1.5) +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 20, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()

#####
## Converting to raster and plotting p-source areas ############################
#####

## Creating an empty raster 
r0_Khamra_PIH_range <- raster::raster(raster::extent(as(data_coord_lakes[1,2:3] %>% 
                                                       st_as_sf(coords = c("lon", "lat")) %>% st_buffer(0.5)    %>%   # 0.5 degress around the lake
                                                       st_bbox() %>% st_as_sfc(), "Spatial")), nrow = 50, ncol = 50) 

## Put the data into the empty raster
raster_Khamra_PIH_range <- rasterize(dist_data_crop_Khamra_PIH_range %>% st_coordinates(), r0_Khamra_PIH_range, fun='count')
crs(raster_Khamra_PIH_range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

## Converting the raster to data.frame
df_r_K_PIH_range     <- as.data.frame(raster_Khamra_PIH_range, xy = TRUE)
df_r_K_PIH_range_na  <- na.omit(df_r_K_PIH_range)

## Rasterize the fire radiative power data and creating MULTIPOLYGON
FRP_Khamra_PIH_range <- rasterize(tibble_frp_pih_wind_Khamra %>% filter(!is.na(FRP)) %>% st_coordinates(),
                               raster_Khamra_PIH_range, fun = 'count') %>% raster::rasterToPolygons() %>%
  st_as_sf() %>% filter(layer > 0) %>% st_union() %>% 
  st_set_crs(4326) %>% st_transform(4326)

## Normalizing the p-source	probabilities 
df_r_K_PIH_range_na$normalized <- (df_r_K_PIH_range_na$layer-min(df_r_K_PIH_range_na$layer))/(max(df_r_K_PIH_range_na$layer)-min(df_r_K_PIH_range_na$layer))

## Set crs 4326 for Khamra catchment
Khamra_catchment <- Khamra_catchment %>% st_set_crs(4326) %>% st_transform(crs=4326)

## Plotting
png(glue("Results/Simulation/PIH/Khamra/range/psource/Khamra_psource_raster_PIH_range.png"), width = 1500, height = 1000)
ggplot(Khamra_buf_extent) +
  geom_tile(data = df_r_K_PIH_range_na , aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\np-source area\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])")+
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.2) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill") +
  geom_sf(data = FRP_Khamra_PIH_range, mapping = aes(fill = "Burned area"), alpha = 0.2, colour = NA) +
  scale_fill_manual(values = c("Burned area" = "firebrick2"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill")+
  geom_sf(data = Khamra, aes(fill = "Lake Khamra"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake Khamra" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75))
dev.off()

#####
## Masking p-source raster data with modern FRP data ###########################
######
FRP_sp_Khamra_PIH_range <- as(FRP_Khamra_PIH_range, "Spatial")  ## MULTIPOLYGON to SpatialPolygons
r_mask_Khamra_PIH_range <- raster::mask(raster_Khamra_PIH_range, FRP_sp_Khamra_PIH_range)

df_r_K_mask_PIH_range   <- as.data.frame(r_mask_Khamra_PIH_range, xy = TRUE)
df_r_K_mask_PIH_range   <- na.omit(df_r_K_mask_PIH_range)

df_r_K_mask_PIH_range$normalized <- (df_r_K_mask_PIH_range$layer-min(df_r_K_mask_PIH_range$layer))/(max(df_r_K_mask_PIH_range$layer)-min(df_r_K_mask_PIH_range$layer))

## Plotting
png(glue("Results/Simulation/PIH/Khamra/range/psource/Khamra_psource_raster_extract_PIH_range.png"), width = 1500, height = 1000)
ggplot(Khamra_buf_extent) +
  geom_tile(data = df_r_K_mask_PIH_range, aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\nactual source area\n(measured PIH |\nRange between\n Q[0.2] and Q[0.8])")+
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.2) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  
  new_scale("fill")+
  geom_sf(data = Khamra, aes(fill = "Lake Khamra")) +
  scale_fill_manual(values = c("Lake Khamra" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75))
dev.off()


#####
## Simulation for lake Satagay ##################################################
#####

## Wind data
wndTab_Satagay <- get(load("Results/Simulation/data/wndTab_Satagay.rda"))
# 
# ## FRP and PIH data
# outMODIS_PIH_Satagay_geom         <- get(load("Results/Simulation/data/outMODIS_PIH_Satagay_geom.rda"))
# outMODIS_PIH_Satagay_geom_na      <- na.omit(outMODIS_PIH_Satagay_geom)
# outMODIS_PIH_Satagay_geom_na$TIME <- as.Date(outMODIS_PIH_Satagay_geom_na$TIME)
# 
# #####
# ## Creating a tibble with FRP, PIH and wind ####################################
# #####
# 
# tibble_frp_pih_wind_Satagay <- do.call("rbind", lapply(unique(outMODIS_PIH_Satagay_geom_na$TIME), function(d) { 
#   
#   cat(sprintf('\rDate %s of %s  ', which(d == unique(outMODIS_PIH_Satagay_geom_na %>% 
#                                                        pull(TIME) %>% as.Date())), length(unique(outMODIS_PIH_Satagay_geom_na %>% pull(TIME) %>% as.Date()))))
#   
#   wndTab_Satagay$date <- as.Date(wndTab_Satagay$date)
#   tmp_wind    <- wndTab_Satagay %>% group_by(date) %>% filter(date == d) %>%
#     dplyr::select("date", "lon", "lat", "wind_spd", "wind_dir")
#   tmp_wind_2  <- st_as_sf(tmp_wind, coords = c("lon", "lat")) %>% st_set_crs(4326) # define coordinates
#   tmp_frp_pih <- outMODIS_PIH_Satagay_geom_na %>% group_by(TIME) %>% filter(TIME == d) 
#   tmp_frp_pih <- st_as_sf(tmp_frp_pih, crs = st_crs(4326))
#   
#   # Calculating the distance between tmp2 and x 
#   distM <- st_distance(tmp_frp_pih, tmp_wind_2, by_element = F)
#   
#   # Selecting the minimum distance between PIH pixel to FRP pixel
#   ind <- apply(distM, 1, function(y) which.min(y))
#   
#   # Creating a tibble of coordinates, FRP, PIH and dist
#   dist <- sapply(1:length(ind), function(z) distM[z, ind[z]]/1000)
#   
#   df <- tmp_frp_pih %>% mutate(w_spd = tmp_wind_2$wind_spd[ind], 
#                                w_dir = tmp_wind_2$wind_dir[ind],
#                                dist = as.numeric(dist))
#   
# }))
# save(tibble_frp_pih_wind_Satagay, file = "Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda")
tibble_frp_pih_wind_Satagay  <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))

## Convert MW to cal/s
tibble_frp_pih_wind_Satagay$FRP_cal <- tibble_frp_pih_wind_Satagay$FRP*238845.8966275

#####
## Calculation of the travel distance of spherical particle ####################
#####

## Cropping tibble with 100 km radius box around the lake
tibble_frp_pih_wind_Satagay <- st_crop(tibble_frp_pih_wind_Satagay, Satagay_buf_extent)

## Plotting empirical_density for measured PIH
png(glue("Results/Simulation/PIH/Satagay/range/Empirical_density_PIH_range_Satagay.png"), width = 1000, height = 700)
plotdist(as.numeric(PIH_range_Satagay), histo = TRUE, demp = TRUE, )
dev.off()

## Gamma distribution of measured PIH
PIH_distr_Satagay_range <- MASS::fitdistr(as.numeric(PIH_range_Satagay), densfun = "gamma")
save(PIH_distr_Satagay_range, file = "Results/Simulation/data/PIH_distr_Satagay_range.rda")
PIH_distr_Satagay_range <- get(load("Results/Simulation/data/PIH_distr_Satagay_range.rda"))

#####
## Function of the potential modern origin of charcoal #########################
#####

origin  <- function(spd, dir, pih, diam, lon_orig, lat_orig){
  vt <- (((0.5-0.00127)*981*diam^2)/(18*(1.8*10^-4)))/100 # cm/s in m/s
  t  <- pih/(vt)
  delta_x <- t*spd
  geosphere::destPoint(c(lon_orig, lat_orig),
                       (dir+180)%%360, delta_x)
}

## Creating a charcoal diameter range from 150 to 500 micrometer based on the
## lake sediment samples analyzed by Glückler et al.
charc_diam_cm <- seq(0.015, 0.05, length = 100)

## Calculation of the travel distance ##########################################

sumCrds_Satagay_PIH_range <- do.call("rbind", lapply(unique(wndTab_Satagay$date), function(d) {

  cat(sprintf('\rDate %s of %s  ', which(d == unique(wndTab_Satagay %>% pull(date) %>% as.Date())),
              length(unique(wndTab_Satagay %>% pull(date) %>% as.Date()))))

  tmp <- wndTab_Satagay %>% filter(date == d) %>% ungroup() %>% st_as_sf(coords = c("lon", "lat")) %>%
    st_set_crs(4326) %>% dplyr::select(wind_spd, wind_dir)

  tmp <- tmp[apply(st_distance(tmp, data_coord_lakes       %>%
                                 st_as_sf(coords = c("lon", "lat"))                %>%
                                 st_set_crs(4326)), 2, function(x) which.min(x)),] %>% st_drop_geometry()

  t <- do.call("rbind", lapply(1:20, function(r) {
    tibble(lake = 2, as.data.frame(origin(spd  = tmp$wind_spd[2],
                                          dir  = tmp$wind_dir[2],
                                          pih  = rgamma(1, PIH_distr_Satagay_range$estimate[1], PIH_distr_Satagay_range$estimate[2]),
                                          diam = charc_diam_cm[sample(1:length(charc_diam_cm),1)],
                                          lon_orig = data_coord_lakes$lon[2],
                                          lat_orig = data_coord_lakes$lat[2]))) %>% setNames(c("lake", "lon", "lat"))
  })) %>% dplyr::mutate(date = d) %>% dplyr::select(date, lake, lon, lat)
}))
save(sumCrds_Satagay_PIH_range, file = "Results/Simulation/data/sumCrds_Satagay_PIH_range.rda")
sumCrds_Satagay_PIH_range <- get(load("Results/Simulation/data/sumCrds_Satagay_PIH_range.rda"))

#####
## Calculating the distance between the coordinates, put it into tibble and
## convert to points
#####

dist_points_Satagay_PIH_range  <- sumCrds_Satagay_PIH_range %>% filter(!is.na(lon) & !is.na(lat)) %>%
                                  st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326)

dist_data_Satagay_PIH_range <- do.call("rbind", lapply(unique(dist_points_Satagay_PIH_range$date), function(x) {

  cat(sprintf('\rDate %s of %s  ', which(x == unique(dist_points_Satagay_PIH_range %>% pull(date) %>% as.Date())),
  length(unique(dist_points_Satagay_PIH_range %>% pull(date) %>% as.Date()))))

  tmp <- dist_points_Satagay_PIH_range %>% filter(date == x)

  tmp %>% dplyr::mutate(dist = as.numeric(st_distance(geometry, (Satagay_buf_100 %>%
                 st_centroid() %>% st_transform(4326)))))

}))
save(dist_data_Satagay_PIH_range, file = "Results/Simulation/data/dist_data_Satagay_PIH_range.rda")
dist_data_Satagay_PIH_range <- get(load("Results/Simulation/data/dist_data_Satagay_PIH_range.rda"))
dist_data_Satagay_PIH_range <- na.omit(dist_data_Satagay_PIH_range)

#####
## Plotting the travel distance of spherical particle in form of histogram #####
##### 

## 100 breaks
png(glue("Results/Simulation/PIH/Satagay/range/dist/travel_distance_Satagay_PIH_4000_100_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Satagay_PIH_range$dist[dist_data_Satagay_PIH_range$lake==2 & dist_data_Satagay_PIH_range$dist<400000]/1000, breaks= 100, 
     xlim = c(0,20), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks
png(glue("Results/Simulation/PIH/Satagay/range/dist/travel_distance_Satagay_PIH_2000_200_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Satagay_PIH_range$dist[dist_data_Satagay_PIH_range$lake==2 & dist_data_Satagay_PIH_range$dist<400000]/1000, breaks= 200, 
     xlim = c(0,20), ylim = c(0, 3000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 100 breaks | 40 km
png(glue("Results/Simulation/PIH/Satagay/range/dist/travel_distance_Satagay_PIH_4000_100_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Satagay_PIH_range$dist[dist_data_Satagay_PIH_range$lake==2 & dist_data_Satagay_PIH_range$dist<400000]/1000, breaks= 100, 
     xlim = c(0,40), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks | 40 km
png(glue("Results/Simulation/PIH/Satagay/range/dist/travel_distance_Satagay_PIH_2000_200_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Satagay_PIH_range$dist[dist_data_Satagay_PIH_range$lake==2 & dist_data_Satagay_PIH_range$dist<400000]/1000, breaks= 200, 
     xlim = c(0,40), ylim = c(0, 3000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

#####
## Plotting a boxplot of measured PIH ##########################################
#####

n_dist_S_range <- length(dist_data_Satagay_PIH_range$dist) #9840
dist_data_Satagay_PIH_range$dist <- dist_data_Satagay_PIH_range$dist/1000
summary(dist_data_Satagay_PIH_range$dist) ## in km 
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  0.02977  0.89088     1.66699  2.78027   3.25907  66.26952 
median_dist_Satagay_PIH_range <- summary(dist_data_Satagay_PIH_range$dist)[3]

## Boxplot

## 12 km
png(glue("Results/Simulation/PIH/Satagay/range/dist/travel_distance_Satagay_PIH_range_boxplot_12km.png"), width = 900, height = 800)
boxplot(x = dist_data_Satagay_PIH_range$dist, outline = F, ylim = c(0,12),
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay", 
        ylab = "Δx [km]", xlab = glue("n = {n_dist_S_range}"), cex.main=1.7, cex.lab=1.5, cex.axis=1.5, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## 8 km 
png(glue("Results/Simulation/PIH/Satagay/range/dist/travel_distance_Satagay_PIH_range_boxplot_8km.png"), width = 900, height = 800)
boxplot(x = dist_data_Satagay_PIH_range$dist, outline = F, ylim = c(0,8),
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay", 
        ylab = "Δx [km]", xlab = glue("n = {n_dist_S_range}"), cex.main=1.7, cex.lab=1.5, cex.axis=1.5, col = "lightblue3")
legend("topright", legend=c("Measured PIH | Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)
dev.off()

## Boxplot dist Khamra and Satagay parallel ####################################

png(glue("Results/Simulation/PIH/boxplot_travel_distances_K_S_range.png"), width = 1200, height = 800)
layout(matrix(c(1,2),nrow = 1, ncol = 2, byrow = TRUE))
boxplot(x = dist_data_Khamra_PIH_range$dist, outline = F, 
        xlab =  "Study area Lake Khamra", 
        ylab = "Δx [km]", cex.main=1.5, cex.lab=1.5, cex.axis=1.5, 
        col = "lightblue3", ylim = c(0,12))
boxplot(x = dist_data_Satagay_PIH_range$dist, outline = F, 
        ylim = c(0,12),
        xlab = "Study area Lake Satagay", cex.main=1.2, cex.lab=1.5, cex.axis=1.5, col = "lightblue3", ylim = c(0,8))
mtext("Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density) into the atmosphere\nwith measured PIH [m] | Range between Q[0.2] and Q[0.8]", 
      line = -3.2 ,outer = TRUE, cex = 1.6, font = 2, adj = 0.6)
dev.off()

png(glue("Results/Simulation/PIH/boxplot_travel_distances_K_S_range_8km.png"), width = 1200, height = 800)
layout(matrix(c(1,2),nrow = 1, ncol = 2, byrow = TRUE))
boxplot(x = dist_data_Khamra_PIH_range$dist, outline = F, 
        xlab =  "Study area Lake Khamra", 
        ylab = "Δx [km]", cex.main=1.5, cex.lab=1.5, cex.axis=1.5, 
        col = "lightblue3", ylim = c(0,8))
boxplot(x = dist_data_Satagay_PIH_range$dist, outline = F, 
        ylim = c(0,8),
        xlab = "Study area Lake Satagay", cex.main=1.2, cex.lab=1.5, cex.axis=1.5, col = "lightblue3", ylim = c(0,8))
mtext("Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density) into the atmosphere\nwith measured PIH [m] | Range between Q[0.2] and Q[0.8]", 
      line = -3.2 ,outer = TRUE, cex = 1.6, font = 2, adj = 0.6)
dev.off()


## Multiple plot for Khamra and Satagay ########################################
dist_data_Satagay_PIH_range <- get(load("Results/Simulation/data/dist_data_Satagay_PIH_range.rda"))
dist_data_Satagay_PIH_range <- na.omit(dist_data_Satagay_PIH_range)
dist_data_Khamra_PIH_range <- get(load("Results/Simulation/data/dist_data_Khamra_PIH_range.rda"))

hist_dist_S_range <- hist(dist_data_Satagay_PIH_range$dist[dist_data_Satagay_PIH_range$lake== 2 & dist_data_Satagay_PIH_range$dist<400000]/1000,  breaks= 100, 
                       xlim = c(0,20), ylim = c(0, 4000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")

hist_dist_K_range <- hist(dist_data_Khamra_PIH_range$dist[dist_data_Khamra_PIH_range$lake==1 & dist_data_Khamra_PIH_range$dist<400000]/1000, breaks = 100,
                       xlim = c(0,20), ylim = c(0, 4000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
legend("topright", legend=c("Measured PIH |  Range between Q[0.2] and Q[0.8]"), bty = "n", cex = 1.5)


png(glue("Results/Simulation/PIH/boxplot_hist_travel_distances_K_S_range.png"), width = 1000, height = 800)
layout(matrix(c(1,2,3,3),nrow = 2, ncol = 2, byrow = TRUE))
boxplot(x = dist_data_Khamra_PIH_range$dist/1000, outline = F, 
        ylab = "Δx [km]", cex.main=1.5, cex.lab=1.8, cex.axis=1.8, 
        col = rgb(1,0,1,1/4), ylim = c(0,8))
boxplot(x = dist_data_Satagay_PIH_range$dist/1000, outline = F, 
        ylim = c(0,8),cex.main=1.2, cex.lab=1.8, cex.axis=1.8,
        col = rgb(0,0,1,1/4))

plot(hist_dist_S_range, col = rgb(0,1,1, 1/4), main = NULL, xlab = "Δx [km]", cex.lab = 1.7, xlim = c(0, 20), ylim = c(0,4000))
plot(hist_dist_K_range, col = rgb(1,0,1,1/4), add = T, main = NULL, xlab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5,  xlim = c(0, 20), ylim = c(0,4000))

legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col=c(rgb(1,0,1,1/4), rgb(0,0,1,1/4)), pt.cex=3, pch=15, cex = 1.8)

mtext("Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density) into the atmosphere\nMeasured PIH [m] | Range between Q[0.2] and Q[0.8]",
      line = -3.5 ,outer = TRUE, cex = 1.6, font = 2, adj = 0.6)
dev.off()


## Selecting all distance values above the median 
# dist_above_median_Satagay_PIH_20 <- dist_data_Satagay_PIH_20 %>% group_by(dist) %>% filter(dist >= median_dist_Satagay_PIH_20)
# length(dist_above_median_Satagay_PIH_20$dist) #7380

## Cropping distance tibble with 100 km radius box around the lake
dist_data_crop_Satagay_PIH_range  <- st_crop(dist_data_Satagay_PIH_range, Satagay_buf_extent)

## Selecting a FRP-range
MODIS_Satagay_PIH_range  <- tibble_frp_pih_wind_Satagay %>% filter(FRP_cal > 0 & FRP_cal <= 71653768.9882488) 
## 300 MW = 71653768.9882488 cal/s

#####
## Plotting the backward simulation of potential source areas of charcoal particles
#####

## FRP (cal/s)
png(glue("Results/Simulation/PIH/Satagay/range/psource/Satagay_psource_PIH_range_cal.png"), width = 1000, height = 900)
ggplot(Satagay_buf_extent) +
  geom_sf(data = dist_data_crop_Satagay_PIH_range, aes(colour  = "Locations of\np-source area\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Satagay_PIH_range, aes(color = FRP_cal), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(title = "Backward simulation of potential source areas\nof charcoal particles in the study area Lake Satagay",
       color = "FRP [cal/s]") +
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.8) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA, alpha = 0.3))) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[2]], y = lat[[2]],
                           shape = location[[2]]), data = data_coord_lakes, colour = "turquoise1",
             size = 2, stroke = 1.5) +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 20, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()


## FRP (MW)
MODIS_Satagay_PIH_range  <- tibble_frp_pih_wind_Satagay %>% filter(FRP > 0 & FRP <= 300) 

png(glue("Results/Simulation/PIH/Satagay/range/psource/Satagay_psource_PIH_range_MW.png"), width = 1000, height = 900)
ggplot(Satagay_buf_extent) +
  geom_sf(data = dist_data_crop_Satagay_PIH_range, aes(colour  = "Locations of\np-source area\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source area\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Satagay_PIH_range, aes(color = FRP), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(title = "Backward simulation of potential source areas\nof charcoal particles in the study area Lake Satagay",
       color = "FRP [MW]") +
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.8) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA, alpha = 0.3))) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[2]], y = lat[[2]],
                           shape = location[[2]]), data = data_coord_lakes, colour = "turquoise1",
             size = 2, stroke = 1.5) +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 20, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()

#####
## Converting to raster and plotting p-source areas ############################
#####

## Creating an empty raster 
r0_Satagay_PIH_range <- raster::raster(raster::extent(as(data_coord_lakes[2,2:3] %>% 
                                                        st_as_sf(coords = c("lon", "lat")) %>% st_buffer(0.5)    %>%   # 0.5 degress around the lake
                                                        st_bbox() %>% st_as_sfc(), "Spatial")), nrow = 50, ncol = 50) 

## Put the data into the empty raster
raster_Satagay_PIH_range <- rasterize(dist_data_crop_Satagay_PIH_range %>% st_coordinates(), r0_Satagay_PIH_range, fun='count')
crs(raster_Satagay_PIH_range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

## Converting the raster to data.frame
df_r_S_PIH_range      <- as.data.frame(raster_Satagay_PIH_range, xy = TRUE)
df_r_S_PIH_range_na   <- na.omit(df_r_S_PIH_range)

## Rasterize the fire radiative power data and creating MULTIPOLYGON
FRP_Satagay_PIH_range <- rasterize(tibble_frp_pih_wind_Satagay %>% filter(!is.na(FRP)) %>% st_coordinates(),
                         raster_Satagay_PIH_range, fun = 'count') %>% raster::rasterToPolygons() %>%
                         st_as_sf() %>% filter(layer > 0) %>% st_union() %>% 
                         st_set_crs(4326) %>% st_transform(4326)

## Normalizing the p-source	probabilities 
df_r_S_PIH_range_na$normalized <- (df_r_S_PIH_range_na$layer-min(df_r_S_PIH_range_na$layer))/(max(df_r_S_PIH_range_na$layer)-min(df_r_S_PIH_range_na$layer))

## Set crs 4326 for Khamra catchment
Satagay_catchment <- Satagay_catchment %>% st_set_crs(4326) %>% st_transform(crs=4326)

## Plotting
png(glue("Results/Simulation/PIH/Satagay/range/psource/Satagay_psource_raster_PIH_range.png"), width = 1500, height = 1000)
ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_PIH_range_na , aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\np-source area\n(measured PIH |\n Range between\n Q[0.2] and Q[0.8])")+
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.4) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill") +
  geom_sf(data = FRP_Satagay_PIH_range, mapping = aes(fill = "Burned area"), alpha = 0.2, colour = NA) +
  scale_fill_manual(values = c("Burned area" = "firebrick2"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake Satagay"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake Satagay" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75))
dev.off()

#####
## Masking p-source raster data with modern FRP data ###########################
######
FRP_sp_Satagay_PIH_range <- as(FRP_Satagay_PIH_range, "Spatial")  ## MULTIPOLYGON to SpatialPolygons
r_mask_Satagay_PIH_range <- raster::mask(raster_Satagay_PIH_range, FRP_sp_Satagay_PIH_range)

df_r_S_mask_PIH_range   <- as.data.frame(r_mask_Satagay_PIH_range, xy = TRUE)
df_r_S_mask_PIH_range   <- na.omit(df_r_S_mask_PIH_range)

df_r_S_mask_PIH_range$normalized <- (df_r_S_mask_PIH_range$layer-min(df_r_S_mask_PIH_range$layer))/(max(df_r_S_mask_PIH_range$layer)-min(df_r_S_mask_PIH_range))

## Plotting
png(glue("Results/Simulation/PIH/Satagay/range/psource/Satagay_psource_raster_extract_PIH_range.png"), width = 1500, height = 1000)
ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_mask_PIH_range, aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\nactual source area\n(measured PIH |\nRange between\n Q[0.2] and Q[0.8])")+
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.2) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake Satagay")) +
  scale_fill_manual(values = c("Lake Satagay" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75))
dev.off()
