#####
## Script to simulate potential and actual modern source areas of charcoal #####
## using modeled plume injection heights #######################################
#####

## The term "H" is used for the calculated plume injection height based on the 
## calculation of Vachula & Richter (2018). In the figures, the term has been 
## changed to "calculated PIH" to allow better comparison with the 
## satellite-measured PIH from MODIS 

## rm(list = ls(all= TRUE))

#####
## Load packages ###############################################################
#####

library(tidyverse)
library(sf)
library(terra)
library(stars)
library(dplyr)
library(glue)

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
library(graphics)

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
Khamra_buf_100          <- st_transform(Khamra, crs = CRS(proj_K)) %>% 
                           st_buffer(100000) %>% st_transform(4326) 
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
Satagay_buf_100         <- st_transform(Satagay, crs = CRS(proj_S)) %>%
                           st_buffer(100000) %>% st_transform(4326)
Satagay_catchment       <- read_sf("Data/Satagay/satagay_catchment.shp") %>% st_transform(4326)
Satagay_buf_extent      <- as(extent(st_bbox(Satagay_buf_100 %>% 
                                               st_transform(4326) %>% 
                                               st_shift_longitude())[c(1,3,2,4)]), "SpatialPolygons")
crs(Satagay_buf_extent) <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")



#####
## Simulation for Lake Khamra ##################################################
#####

## Fire radiative power [FRP] ##################################################
## Fire years for Lake Khamra: 2001, 2003, 2004, 2007, 2013, 2014

sf_MODIS_Khamra_box      <- get(load("Results/FRP/FRP_Khamra_points_box.RData"))
sf_MODIS_Khamra_box$date <- sf_MODIS_Khamra_box %>% pull(TIME) %>% as.Date()
sf_MODIS_Khamra_box$TIME <- NULL
sf_MODIS_Khamra_box      <- as.tibble(sf_MODIS_Khamra_box)

## Filter the time period 2001-2018 and the classified fire years
sf_MODIS_Khamra_box_fy   <- sf_MODIS_Khamra_box %>% group_by(date)                  %>% 
                            mutate(Year  = as.numeric(format(date, "%Y")),
                                   Month = as.numeric(format(date, "%m")))          %>% 
                            filter(Year %in% c(2001, 2003, 2004, 2007, 2013, 2014)) %>% 
                            filter(Month >= 5 & Month <= 8)


## Wind data ###################################################################
wndTabSpdDir  <- get(load("Data/ERA5/wndTabSpdDir.rda"))
wndTab_Khamra <- wndTabSpdDir %>% group_by(date)  %>% 
                 mutate(Year  = as.numeric(format(date, "%Y")),
                        Month = as.numeric(format(date, "%m"))) %>% 
                 filter(Year %in% c(2001, 2003, 2004, 2007, 2013, 2014))
save(wndTab_Khamra, file = "Results/Simulation/data/wndTab_Khamra.rda")
wndTab_Khamra <- get(load("Results/Simulation/data/wndTab_Khamra.rda"))


#####
## Create a tibble with FRP and wind values ####################################
#####

tibble_wind_frp_Khamra <- do.call("rbind", lapply(unique(sf_MODIS_Khamra_box_fy$date), function(d) { 
  
cat(sprintf('\rDate %s of %s  ', which(d == unique(sf_MODIS_Khamra_box_fy %>% 
pull(date) %>% as.Date())), length(unique(sf_MODIS_Khamra_box_fy %>% pull(date) %>% as.Date()))))
  
wndTab_Khamra$date <- as.Date(wndTab_Khamra$date)
       tmp_wind    <- wndTab_Khamra %>% group_by(date) %>% filter(date == d) %>%
                      dplyr::select("date", "lon", "lat", "wind_spd", "wind_dir")
       ## Wind variables were converted from data type st to sf
       tmp_wind_2  <- st_as_sf(tmp_wind, coords = c("lon", "lat")) %>% st_set_crs(4326) # define coordinates
       tmp_frp_pih <- sf_MODIS_Khamra_box_fy %>% group_by(date) %>% filter(date == d) 
       tmp_frp_pih <- st_as_sf(tmp_frp_pih, crs = st_crs(4326))

  ## Calculate the distance between FRP and wind data points
  distM <- st_distance(tmp_frp_pih, tmp_wind_2, by_element = F)
  
  ## Determine the minimum distance between these points
  ind <- apply(distM, 1, function(y) which.min(y))
  
  ## Create a tibble of geographical coordinates, FRP, PIH and distances
  dist <- sapply(1:length(ind), function(z) distM[z, ind[z]]/1000)
  
  df <- tmp_frp_pih %>% mutate(w_spd = tmp_wind_2$wind_spd[ind], 
                               w_dir = tmp_wind_2$wind_dir[ind],
                               dist = as.numeric(dist))
  
}))
save(tibble_wind_frp_Khamra, file = "Results/Simulation/data/tibble_wind_frp_Khamra.rda")
tibble_wind_frp_Khamra     <- get(load("Results/Simulation/data/tibble_wind_frp_Khamra.rda"))
tibble_wind_frp_Khamra     <- tibble_wind_frp_Khamra %>% mutate(H = NA)

## Convert MW to cal/s
tibble_wind_frp_Khamra$FRP <- tibble_wind_frp_Khamra$FRP*238845.8966275

#####
## Calculation of the plume injection height, as a function of fire intensities
## and wind speed based on the simulation model of Vachula & Richter (2018) 
#####

H_function_Khamra <- do.call("rbind", lapply(unique(tibble_wind_frp_Khamra$date), function(j) { 
  
                     cat(sprintf('\rDate %s of %s  ', which(j == unique(tibble_wind_frp_Khamra %>% pull(date) %>% as.Date())),
                     length(unique(tibble_wind_frp_Khamra %>% pull(date) %>% as.Date()))))
        
              tmp <- tibble_wind_frp_Khamra %>% group_by(date) %>% filter(date == j) 
             
             for(r in 1:length(tmp$FRP)) {
             calc <- if((tmp$FRP[r]) < (1.4*10^6)){
                      (0.01*tmp$FRP[r])^0.75/tmp$w_spd[r]}
                      else{(0.085*tmp$FRP[r])^0.6/tmp$w_spd[r]}
                      tmp$H[r] <- calc
                      }
                      tmp
    
}))
save(H_function_Khamra, file = "Results/Simulation/data/H_function_Khamra.rda")
H_function_Khamra <- get(load("Results/Simulation/data/H_function_Khamra.rda"))
H_function_Khamra$FRP[H_function_Khamra$FRP == 0] <- NA
H_function_Khamra$H[H_function_Khamra$H == 0]     <- NA
H_function_Khamra                                 <- na.omit(H_function_Khamra)


#####
## Visualization of the calculated plume injection heights #####################
#####


## Histogramm
n_K <- length(H_function_Khamra$H) # [1] 4816

png(glue("Results/Simulation/H/Khamra/H/calculated_PIH_VR_H_Khamra.png"), width = 1000, height = 700)
hist(H_function_Khamra$H, xlim = c(0, 10000), ylim = c(0,1000), 
     breaks= 100, xlab = "Calculated PIH [m]", 
     main = "Calculated plume injection heights (n = 4816) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "grey40") # in m 
dev.off()

## Summary
summary(H_function_Khamra$H)
## Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 63.71   399.07   691.85  1207.22  1315.94 26725.81 
## Boxplot
png(glue("Results/Simulation/H/Khamra/H/calculated_PIH_VR_H_Khamra_boxplot.png"), width = 1000, height = 800)
boxplot(H_function_Khamra$H, ylab = "Calculated PIH [m]",outline=FALSE,
     main = "Calculated plume injection heights (n = 4816) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, xlab = glue("n = {n_K}"),ylim = c(0,3000),  col = "grey40") # in m 
dev.off()



#####
## Visualization the link between FRP, calculated PIH and wind speed ###########
#####

## with unit cal/s
png(glue("Results/Simulation/H/Khamra/linkage/Calculated_H_Multiple_2000_Khamra_71653768.9882488.png"), width = 1000, height = 900)
ggplot(data = H_function_Khamra,
       aes(x = FRP, y = H, color = w_spd))+
       geom_point() +
       scale_color_gradientn(colours = rainbow(10), limit = c(1,20))+ 
                             #breaks = c(2,4,6,8,10,12,14,16)) +
       labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Khamra",
            color = "spd [m/s]")+
       xlab("FRP [cal/s]") +
       ylab("Calculated PIH [m]")+
       xlim(c(0,71653768.9882488))+
       ylim(c(0, 2000))+
       theme_minimal() +
         theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
               axis.title    = element_text(size = 20, hjust = 0.5),
               legend.title  = element_text(size = 20, vjust = 0.5),
               legend.text   = element_text(size = 18, vjust = 0.75),
               axis.text     = element_text(size = 16, vjust = 0.75))
dev.off() 

## with ylim(c(0, 1000))
png(glue("Results/Simulation/H/Khamra/linkage/Calculated_H_Multiple_10000_Khamra_71653768.9882488.png"), width = 1000, height = 900)
ggplot(data = H_function_Khamra,
       aes(x = FRP, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10), limit = c(1,20))+
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Khamra",
       color = "spd [m/s]")+
  xlab("FRP [cal/s]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 71653768.9882488))+
  ylim(c(0, 1000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()       

## with unit MW
H_function_Khamra$FRP_MW    <- H_function_Khamra$FRP/238845.8966275
png(glue("Results/Simulation/H/Khamra/linkage/Calculated_H_Multiple_2000_Khamra_MW_300.png"), width = 1000, height = 900)
ggplot(data = H_function_Khamra,
       aes(x = FRP_MW, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10), limit = c(1,20))+
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Khamra",
       color = "spd [m/s]")+
  xlab("FRP [MW]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 300))+
  ylim(c(0, 2000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off() 

## with ylim(c(0, 1000))
png(glue("Results/Simulation/H/Khamra/linkage/Calculated_H_Multiple_1000_Khamra_MW_3000.png"), width = 1000, height = 900)
ggplot(data = H_function_Khamra,
       aes(x = FRP_MW, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10), limit = c(1,20))+
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Khamra",
       color = "spd [m/s]")+
  xlab("FRP [MW]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 300))+
  ylim(c(0, 1000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()       




#####
## Calculation of the horizontal travel distance of spherical particles within the air
#####

## Crop the tibble (H_function) with 100 km radius box around the lake
Khamra_data <- st_crop(H_function_Khamra, Khamra_buf_extent)
Khamra_data$H[Khamra_data$H==0] <- NA
Khamra_data <- na.omit(Khamra_data)

## Visualization of the empirical_density of calculated PIH's
png(glue("Results/Simulation/H/Khamra/H/Empirical_density_VR_H_Khamra.png"), width = 1000, height = 700)
plotdist(Khamra_data$H, histo = TRUE, demp = TRUE, )
dev.off()

## Gamma distribution of calculated PIH's
H_distr_Khamra <- MASS::fitdistr(Khamra_data$H, densfun = "gamma")



#####
## Function for the calculation of the potential modern origin of charcoal #####
#####

origin  <- function(spd, dir, pih, diam, lon_orig, lat_orig){
     vt <- (((0.5-0.00127)*981*diam^2)/(18*(1.8*10^-4)))/100 # cm/s in m/s
     t  <- pih/(vt)
           delta_x <- t*spd # maximum charcoal dispersion distances
           # Determination of the distance between lake as center and detected
           # potential source areas
           geosphere::destPoint(c(lon_orig, lat_orig),
                               (dir+180)%%360, delta_x)
}


## In this study, it was assumed a spherical particle shape with a density of 0.5 g/cm3
## (like Clark, 1988; Gilgen et al., 2018) and a diameter range of 150 to 500 μm 
## (0.015 to 0.050 cm) based on the sampled lake sediment records (Glückler et al. 2021)

charc_diam_cm <- seq(0.015, 0.05, length = 100) 

#####
## Calculation of the travel distance ##########################################
#####

sumCrds_Khamra_VR <- do.call("rbind", lapply(unique(wndTab_Khamra$date), function(d) { 
  
cat(sprintf('\rDate %s of %s  ', which(d == unique(wndTab_Khamra %>% pull(date) %>% as.Date())),
length(unique(wndTab_Khamra %>% pull(date) %>% as.Date()))))
  
tmp <- wndTab_Khamra %>% filter(date == d) %>% ungroup() %>% st_as_sf(coords = c("lon", "lat")) %>%
       st_set_crs(4326) %>% dplyr::select(wind_spd, wind_dir)
              
tmp <- tmp[apply(st_distance(tmp, data_coord_lakes       %>% 
       st_as_sf(coords = c("lon", "lat"))                %>% 
       st_set_crs(4326)), 2, function(x) which.min(x)),] %>%
       st_drop_geometry()
              
do.call("rbind", lapply(1:20, function(r) {
                 tibble(lake = 1,
                 as.data.frame(origin(spd  = tmp$wind_spd[1],
                                      dir  = tmp$wind_dir[1],
                                      pih  = rgamma(1, H_distr_Khamra$estimate[1], H_distr_Khamra$estimate[2]),
                                      diam = charc_diam_cm[sample(1:length(charc_diam_cm),1)],
                                      lon_orig = data_coord_lakes$lon[1],
                                      lat_orig = data_coord_lakes$lat[1]))) %>% setNames(c("lake", "lon", "lat"))
})) %>% dplyr::mutate(date = d) %>% dplyr::select(date, lake, lon, lat)
}))
save(sumCrds_Khamra_VR, file = "Results/Simulation/data/sumCrds_Khamra_VR.rda") 
sumCrds_Khamra_VR <- get(load("Results/Simulation/data/sumCrds_Khamra_VR.rda"))


#####
## Calculation of the distance between detecteed potential charcoal source areas
## and the sampled lake to determine the actual charcoal dispersion distances
#####

dist_data_Khamra <- sumCrds_Khamra_VR %>% filter(!is.na(lon) & !is.na(lat)) %>% st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>%
                    group_split(lake) %>% lapply(function(x) {
                    x %>% mutate(dist = as.numeric(st_distance(geometry, (Khamra_buf_100 %>% st_centroid() %>% st_transform(4326))[unique(lake),])))
                    }) %>% do.call("rbind", .)
save(dist_data_Khamra, file = "Results/Simulation/data/dist_data_Khamra.rda")
dist_data_Khamra <- get(load("Results/Simulation/data/dist_data_Khamra.rda"))


#####
## Visualization of the horizontal travel distances of spherical particles in the air in form of histogram
##### 

## 100 breaks | 20 km 
png(glue("Results/Simulation/H/Khamra/dist/travel_distance_Khamra_VR_H_100_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Khamra$dist[dist_data_Khamra$lake==1 & dist_data_Khamra$dist<400000]/1000, breaks= 100, 
     xlim = c(0,20), ylim = c(0, 6000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey40")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()

## 100 breaks | 40 km
png(glue("Results/Simulation/H/Khamra/dist/travel_distance_Khamra_VR_H_100_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Khamra$dist[dist_data_Khamra$lake==1 & dist_data_Khamra$dist<400000]/1000, breaks= 100, 
     xlim = c(0,40), ylim = c(0, 6000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey40")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks | 20 km
png(glue("Results/Simulation/H/Khamra/dist/travel_distance_Khamra_VR_H_200_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Khamra$dist[dist_data_Khamra$lake==1 & dist_data_Khamra$dist<400000]/1000, breaks= 200, 
     xlim = c(0,20), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1,col = "grey40")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks | 40 km
png(glue("Results/Simulation/H/Khamra/dist/travel_distance_Khamra_VR_H_200_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Khamra$dist[dist_data_Khamra$lake==1 & dist_data_Khamra$dist<400000]/1000, breaks= 200, 
     xlim = c(0,40), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey40")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()


#####
## Boxplot of calculated PIH's #################################################
#####

n_dist_K <- length(dist_data_Khamra$dist)
dist_data_Khamra$dist <- dist_data_Khamra$dist/1000
summary(dist_data_Khamra$dist) ## in km 
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  0.01686   1.18493   3.15735   6.97611   7.73558 192.08712 
median_dist_Khamra_VR <- summary(dist_data_Khamra$dist)[3]
sd(dist_data_Khamra$dist)


## 20 km 
png(glue("Results/Simulation/H/Khamra/dist/travel_distance_Khamra_VR_H_boxplot.png"), width = 900, height = 800)
boxplot(x = dist_data_Khamra$dist, outline = F, ylim= c(0,20),
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Khamra", 
        ylab = "Δx [km]", xlab = glue("n = {n_dist_K}"), cex.main=1.7, cex.lab=1.5, cex.axis=1.5, col = "grey40")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()


## Select all travel distance values above the average median 
dist_above_median_Khamra <- dist_data_Khamra %>% group_by(dist) %>% filter(dist >= median_dist_Khamra_VR)
length(dist_above_median_Khamra$dist) #7380

## Crop the tibble with 100 km radius box around the lake
dist_data_crop_Khamra     <- st_crop(dist_data_Khamra, Khamra_buf_extent)
# dist_data_crop_Khamra_am  <- st_crop(dist_above_median_Khamra, Khamra_buf_extent)

## Select the FRP-range
MODIS_Khamra_VR        <- H_function_Khamra %>% filter(FRP > 0 & FRP <= 71653768.9882488) 
## 300 MW = 71653768.9882488 cal/s

#####
## Visualization of the backward simulation of potential source areas of charcoal particles
#####

## FRP (cal/s)
png(glue("Results/Simulation/H/Khamra/psource/Khamra_psource_VR_H_cal.png"), width = 1000, height = 900)
ggplot(Khamra_buf_extent) +
 geom_sf(data = dist_data_crop_Khamra, aes(colour  = "Locations of\np-source areas\n(calculated PIH)"), show.legend = "point") +
 scale_colour_manual(values = c("Locations of\np-source areas\n(calculated PIH)" = "darkslategrey"), name = NULL,
                     guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                              shape = 16))) + # do not plot the color in legend
 new_scale("colour") +
 geom_sf(data = MODIS_Khamra_VR, aes(color = FRP), size = 2,  alpha = 0.5, na.rm = T) +
 scale_color_continuous_sequential(palette = "Heat") +
 labs(color = "FRP [cal/s]") +
 new_scale("fill") +
 geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.4) +
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
 theme(legend.title  = element_text(size = 22, vjust = 0.5),
       legend.text   = element_text(size = 22, vjust = 0.75),
       axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()


## FRP (MW)
MODIS_Khamra_VR$FRP_MW <- MODIS_Khamra_VR$FRP/238845.8966275

png(glue("Results/Simulation/H/Khamra/psource/Khamra_psource_VR_H_MW.png"), width = 1000, height = 900)
ggplot(Khamra_buf_extent) +
  geom_sf(data = dist_data_crop_Khamra, aes(colour  = "Locations of\np-source areas\n(calculated PIH)"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(calculated PIH)" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Khamra_VR, aes(color = FRP_MW), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(color = "FRP [MW]") +
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.4) +
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
  theme(legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()

## All distances above the median
# png(glue("Results/Simulation/p_source/Khamra_psource_VR_H_cal_am.png"), width = 1000, height = 900)
# ggplot(Khamra_buf_extent) +
#   geom_sf(data = dist_data_crop_Khamra_am, aes(colour  = "p-source\n(with calculated PIH)"), show.legend = "point") +
#   scale_colour_manual(values = c("p-source\n(with calculated PIH)" = "darkslategrey"), name = NULL,
#                       guide = guide_legend(override.aes = list(linetype = c("blank"),
#                                                                shape = 16))) + # do not plot the color in legend
#   new_scale("colour") +
#   geom_sf(data = MODIS_Khamra_VR, aes(color = FRP), size = 2,  alpha = 0.5, na.rm = T) +
#   scale_color_continuous_sequential(palette = "Heat") +
#   labs(title = "Backward simulation of potential source areas\nof charcoal particles in the study area lake Khamra",
#        color = "FRP [cal/s]") +
#   new_scale("fill") +
#   geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.3) +
#   scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
#                     guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA, alpha = 0.3))) +
#   theme_minimal() +
#   new_scale("colour") +
#   geom_point(mapping = aes(x = lon[[1]], y = lat[[1]],
#                            shape = location[[1]]), data = data_coord_lakes, colour = "turquoise1",
#              size = 3, stroke = 1.5) +
#   xlab("") +
#   ylab("") +
#   labs(shape = "Lake") +
#   theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
#         legend.title  = element_text(size = 20, vjust = 0.5),
#         legend.text   = element_text(size = 20, vjust = 0.75),
#         axis.text     = element_text(size = 12, vjust = 0.75))
# dev.off()

#####
## Convert to raster and plot p-source areas again ##############################
#####

## Create an empty raster 
r0_Khamra <- raster::raster(raster::extent(as(data_coord_lakes[1,2:3] %>% 
             st_as_sf(coords = c("lon", "lat")) %>% st_buffer(0.5)    %>%   # 0.5 degress around the lake
             st_bbox() %>% st_as_sfc(), "Spatial")), nrow = 50, ncol = 50) 

## Put the data into the empty raster
raster_Khamra      <- rasterize(dist_data_Khamra %>% st_coordinates(), r0_Khamra, fun='count')
crs(raster_Khamra) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

## Convert the raster to data.frame
df_r_K_VR     <- as.data.frame(raster_Khamra, xy = TRUE)
df_r_K_na_VR  <- na.omit(df_r_K_VR)

## Rasterize the fire radiative power data and create a MULTIPOLYGON
FRP_Khamra_VR <- rasterize(H_function_Khamra %>% filter(!is.na(FRP)) %>% st_coordinates(),
                           raster_Khamra, fun = 'count') %>% raster::rasterToPolygons() %>%
                           st_as_sf() %>% filter(layer > 0) %>% st_union()    %>% 
                           st_set_crs(4326) %>% st_transform(4326)

## Normalize the p-source	probabilities 
df_r_K_na_VR$normalized <- (df_r_K_na_VR$layer-min(df_r_K_na_VR$layer))/(max(df_r_K_na_VR$layer)-min(df_r_K_na_VR$layer))

## Set crs 4326 for Khamra catchment
Khamra_catchment <- Khamra_catchment %>% st_set_crs(4326) %>% st_transform(crs=4326)

#####
## Visualization ################################################################
#####

png(glue("Results/Simulation/H/Khamra/psource/Khamra_psource_raster_VR_H.png"), width = 1500, height = 1000)
ggplot(Khamra_buf_extent) +
 geom_tile(data = df_r_K_na_VR , aes(x = x, y = y, fill = normalized)) +
 scale_fill_gradientn(colours = rev(viridis::mako(99))) +
 labs(fill = "Probability of\np-source areas\n(calculated PIH)")+
 new_scale("fill") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.2) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill") +
  geom_sf(data = FRP_Khamra_VR, mapping = aes(fill = "Burned area"), alpha = 0.2, colour = NA) +
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
  theme(legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()

#####
## Mask the p-source raster data with modern FRP data ###########################
#####

FRP_sp_Khamra      <- as(FRP_Khamra_VR, "Spatial")  ## MULTIPOLYGON to SpatialPolygons
r_mask_Khamra      <- raster::mask(raster_Khamra, FRP_sp_Khamra)

df_r_K_VR_mask     <- as.data.frame(r_mask_Khamra, xy = TRUE)
df_r_K_VR_mask     <- na.omit(df_r_K_VR_mask)

## Normalization
df_r_K_VR_mask$normalized <- (df_r_K_VR_mask$layer-min(df_r_K_VR_mask$layer))/(max(df_r_K_VR_mask$layer)-min(df_r_K_VR_mask$layer))

#####
## Visualization ################################################################
#####

png(glue("Results/Simulation/H/Khamra/psource/Khamra_psource_raster_VR_H_extract.png"), width = 1500, height = 1000)
ggplot(Khamra_buf_extent) +
  geom_tile(data = df_r_K_VR_mask, aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\nactual source areas\n(calculated PIH)")+
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
  theme(legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()




#####
## Simulation for Lake Satagay ##################################################
#####

## Fire radiative power ########################################################
## Fire years for Lake Satagay: 2002, 2013, 2014, 2018

## Lake Satagay
sf_MODIS_Satagay_box      <- get(load("Results/FRP/FRP_Satagay_points_box.RData"))
sf_MODIS_Satagay_box$date <- sf_MODIS_Satagay_box %>% pull(TIME) %>% as.Date()
sf_MODIS_Satagay_box$TIME <- NULL
sf_MODIS_Satagay_box$Year <- NULL
sf_MODIS_Satagay_box      <- as.tibble(sf_MODIS_Satagay_box)

## Filter the time period 2001 to 2018 and the classified fire years
sf_MODIS_Satagay_box_fy <- sf_MODIS_Satagay_box %>% group_by(date)        %>% 
                           mutate(Year  = as.numeric(format(date, "%Y")),
                                  Month = as.numeric(format(date, "%m"))) %>% 
                           filter(Year %in% c(2002, 2013, 2014, 2018))    %>% 
                           filter(Month >= 5 & Month <= 8)

## Wind data ###################################################################
wndTabSpdDir   <- get(load("Data/ERA5/wndTabSpdDir.rda"))
wndTab_Satagay <- wndTabSpdDir %>% group_by(date)  %>% 
                  mutate(Year  = as.numeric(format(date, "%Y")),
                         Month = as.numeric(format(date, "%m"))) %>% 
                  filter(Year %in% c(2002, 2013, 2014, 2018))
save(wndTab_Satagay, file = "Results/Simulation/data/wndTab_Satagay.rda")
wndTab_Satagay <- get(load("Results/Simulation/data/wndTab_Satagay.rda"))


#####
## Create a tibble with FRP and wind values ####################################
#####


tibble_wind_frp_Satagay <- do.call("rbind", lapply(unique(sf_MODIS_Satagay_box_fy$date), function(d) { 
  
cat(sprintf('\rDate %s of %s  ', which(d == unique(sf_MODIS_Satagay_box_fy %>% 
pull(date) %>% as.Date())), length(unique(sf_MODIS_Satagay_box_fy %>% pull(date) %>% as.Date()))))
  
  wndTab_Satagay$date <- as.Date(wndTab_Satagay$date)
  tmp_wind    <- wndTab_Satagay %>% group_by(date) %>% filter(date == d) %>%
                 dplyr::select("date", "lon", "lat", "wind_spd", "wind_dir")
  ## Wind variables were converted from data type st to sf
  tmp_wind_2  <- st_as_sf(tmp_wind, coords = c("lon", "lat")) %>% st_set_crs(4326) # define coordinates
  tmp_frp_pih <- sf_MODIS_Satagay_box_fy %>% group_by(date) %>% filter(date == d) 
  tmp_frp_pih <- st_as_sf(tmp_frp_pih, crs = st_crs(4326))
  
  # Calculate the distance between FRP and wind data points
  distM <- st_distance(tmp_frp_pih, tmp_wind_2, by_element = F)
  
  # Determine the minimum distance between these points
  ind <- apply(distM, 1, function(y) which.min(y))
  
  # Create a tibble of geographical coordinates, FRP, PIH and distances
  dist <- sapply(1:length(ind), function(z) distM[z, ind[z]]/1000)
  
  df <- tmp_frp_pih %>% mutate(w_spd = tmp_wind_2$wind_spd[ind], 
                               w_dir = tmp_wind_2$wind_dir[ind],
                               dist = as.numeric(dist))
  
}))
save(tibble_wind_frp_Satagay, file = "Results/Simulation/data/tibble_wind_frp_Satagay.rda")
tibble_wind_frp_Satagay     <- get(load("Results/Simulation/data/tibble_wind_frp_Satagay.rda"))
tibble_wind_frp_Satagay     <- tibble_wind_frp_Satagay %>% mutate(H = NA)

## Convert MW to cal/s
tibble_wind_frp_Satagay$FRP <- tibble_wind_frp_Satagay$FRP*238845.8966275

#####
## Calculation of the plume injection height, as a function of fire intensities
## and wind speed based on the simulation model of Vachula & Richter (2018) 
#####

H_function_Satagay <- do.call("rbind", lapply(unique(tibble_wind_frp_Satagay$date), function(j) { 
  
cat(sprintf('\rDate %s of %s  ', which(j == unique(tibble_wind_frp_Satagay %>% pull(date) %>% as.Date())),
length(unique(tibble_wind_frp_Satagay %>% pull(date) %>% as.Date()))))
  
  tmp <- tibble_wind_frp_Satagay %>% group_by(date) %>% filter(date == j) 
  
  for(r in 1:length(tmp$FRP)) {
    calc <- if((tmp$FRP[r]) < (1.4*10^6)){
      (0.01*tmp$FRP[r])^0.75/tmp$w_spd[r]}
    else{(0.085*tmp$FRP[r])^0.6/tmp$w_spd[r]}
    tmp$H[r] <- calc
  }
  tmp
  
}))
save(H_function_Satagay, file = "Results/Simulation/data/H_function_Satagay.rda")
H_function_Satagay <- get(load("Results/Simulation/data/H_function_Satagay.rda"))
H_function_Satagay$FRP[H_function_Satagay$FRP == 0] <- NA
H_function_Satagay$H[H_function_Satagay$H == 0]     <- NA
H_function_Satagay                                  <- na.omit(H_function_Satagay)


#####
## Visualization of the calculated plume injection heights #####################
#####

n_S <- length(H_function_Satagay$H) # [1] 20758 ## number of data points


## Histogramm
png(glue("Results/Simulation/H/Satagay/H/calculated_PIH_VR_H_Satagay.png"), width = 1000, height = 700)
hist(H_function_Satagay$H, xlim = c(0, 10000), 
     breaks= 100, xlab = "Calculated PIH [m]", 
     main = "Calculated plume injection height for the study area Lake Satagay", cex.main = 2,
     cex.lab = 1.5, col = "grey80") # in m 
dev.off()

## Summary
summary(H_function_Satagay$H)
## Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 61.08   316.08   516.58   838.23   914.48 20137.35 
## Boxplot
png(glue("Results/Simulation/H/Satagay/H/calculated_PIH_VR_H_Satagay_boxplot.png"), width = 1000, height = 800)
boxplot(H_function_Satagay$H, ylab = "Calculated PIH [m]",outline=FALSE, 
        main = "Calculated plume injection heights for the study area Lake Satagay", cex.main = 2,
        cex.lab = 1.5, xlab = glue("n = {n_S}"), ylim = c(0,3000), col = "grey80") # in m 
dev.off()

#####
## Boxplot of PIH's for Khamra and Satagay parallel ############################
#####

n_H_K <- length(H_function_Khamra$H)
n_H_S <- length(H_function_Satagay$H)

# png(glue("Results/Simulation/H/boxplot_H_Khamra_Satagay.png"), width = 1200, height = 800)
# layout(matrix(c(1,2),nrow = 1, ncol = 2, byrow = TRUE))
# boxplot(H_function_Khamra$H, ylab = "Calculated PIH [m]",outline=FALSE,
#         xlab = "Study area Lake Khamra", main = NULL,
#         cex.lab = 1.5, xlab = glue("n = {n_H_K}"), ylim = c(0,3000)) # in m 
# boxplot(H_function_Satagay$H, outline=FALSE,
#         xlab = "Study area Lake Satagay", main = NULL,
#         cex.lab = 1.5, ylim = c(0,3000),  xlab = glue("n = {n_H_S}"), ylim = c(0,3000)) # in m 
# mtext("Calculated plume injection heights", line = -1.5,outer = TRUE, cex = 1.8, font = 2)
# dev.off()

H_df <- merge(data.frame(H_function_Khamra, row.names=NULL), data.frame(H_function_Satagay, row.names=NULL), 
            by = 0, all = TRUE)[-1]

median(H_function_Khamra$H)  ##  691.8474 m
median(H_function_Satagay$H) ##  516.5799 m

label_names_H     <- c(glue("n = {n_H_K}"), glue("n = {n_H_S}"))
colors_H_boxplot <- c("grey40", "grey80")

png(glue("Results/Simulation/H/boxplot_H_Khamra_Satagay.png"), width = 800, height = 800)
boxplot(H_df$H.x, H_df$H.y, ylim = c(0,3000),  outline=FALSE, names = label_names_H, col = colors_H_boxplot, 
        main = NULL, cex.main = 2, ylab = "Calculated PIH [m]",
        cex.lab = 1.75, cex.axis = 1.75)
legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col=c("grey40", "grey80"), pt.cex=3, pch=15, cex = 2)
dev.off()


## Multiple plot for calculated PIH's ###########################################

hist_calc_PIH_K <- hist(H_function_Khamra$H, xlim = c(0, 6000), ylim = c(0,4000), 
                        breaks= 200, xlab = "Calculated PIH [m]", main = NULL, cex.axis = 1.5,
                        cex.lab = 1.7, col = "grey40") # in m 
hist_calc_PIH_S <- hist(H_function_Satagay$H, xlim = c(0, 6000),  ylim = c(0,4000),cex.axis = 1.5,
                        breaks= 200, xlab = "Calculated PIH [m]", main = NULL,
                        cex.lab = 1.7, col = "grey80") # in m 


png(glue("Results/Simulation/H/box_hist_calculated_PIH.png"), width = 600, height = 800)
layout(matrix(c(1,1,2,2),nrow = 2, ncol = 2, byrow = TRUE))
boxplot(H_df$H.x, H_df$H.y, ylim = c(0,4000),  outline=FALSE, names = label_names_H, col = colors_H_boxplot, 
        main = NULL, cex.main = 2, ylab = "Calculated PIH [m]",
        cex.lab = 1.75, cex.axis = 1.75)
plot(hist_calc_PIH_S, col = "grey80",  main = NULL, xlab = "Calculated PIH [m]", 
     cex.lab = 1.75, cex.axis = 1.75, xlim = c(0,4000), ylim = c(0,4000))
plot(hist_calc_PIH_K, col = "grey40", add = T,main = NULL, xlab = "Calculated PIH [m]",
     cex.lab = 1.75, cex.axis = 1.75,  xlim = c(0,4000), ylim = c(0,4000))

legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col=c("grey40", "grey80"), pt.cex=3, pch=15, cex = 2)

dev.off()



#####
## Visualization the link between FRP, wind speed and calculated PIH ########
#####

## with unit cal/s
png(glue("Results/Simulation/H/Satagay/linkage/Calculated_H_Multiple_2000_Satagay_71653768.9882488.png"), width = 1000, height = 900)
ggplot(data = H_function_Satagay,
       aes(x = FRP, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10),limit = c(1,20))+
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Satagay",
       color = "spd [m/s]")+
  xlab("FRP [cal/s]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 71653768.9882488))+
  ylim(c(0, 2000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off() 

## with ylim(c(0, 1000))
png(glue("Results/Simulation/H/Satagay/linkage/Calculated_H_Multiple_1000_Satagay_71653768.9882488.png"), width = 1000, height = 900)
ggplot(data = H_function_Satagay,
       aes(x = FRP, y = H, color = w_spd),limit = c(1,20))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10))+#, 
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Satagay",
       color = "spd [m/s]")+
  xlab("FRP [cal/s]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 71653768.9882488))+
  ylim(c(0, 1000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()       

## with unit MW
H_function_Satagay$FRP_MW    <- H_function_Satagay$FRP/238845.8966275
png(glue("Results/Simulation/H/Satagay/linkage/Calculated_H_Multiple_2000_Satagay_MW_300.png"), width = 1000, height = 900)
ggplot(data = H_function_Satagay,
       aes(x = FRP_MW, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10),limit = c(1,20))+ 
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Satagay",
       color = "spd [m/s]")+
  xlab("FRP [MW]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 300))+
  ylim(c(0, 2000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off() 

## with ylim(c(0, 1000))
png(glue("Results/Simulation/H/Satagay/linkage/Calculated_H_Multiple_1000_Satagay_MW_300.png"), width = 1000, height = 900)
ggplot(data = H_function_Satagay,
       aes(x = FRP_MW, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10),limit = c(1,20))+
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "The link between fire radiative power and wind speed,\nas a function of plume injection height for the study area Lake Satagay",
       color = "spd [m/s]")+
  xlab("FRP [MW]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 300))+
  ylim(c(0, 1000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()       



#####
## Calculation of the horizontal travel distance of spherical particles within the air
#####

## Crop the tibble (H_function) with 100 km radius box around the lakeSatagay_data <- st_crop(H_function_Satagay, Satagay_buf_extent)
Satagay_data$H[Satagay_data$H==0] <- NA
Satagay_data <- na.omit(Satagay_data)

## Visualization of the empirical_density of calculated PIH's
png(glue("Results/Simulation/H/Satagay/H/Empirical_density_VR_H_Satagay.png"), width = 1000, height = 700)
plotdist(Satagay_data$H, histo = TRUE, demp = TRUE, )
dev.off()

## Gamma distribution of calculated PIH's
H_distr_Satagay <- MASS::fitdistr(Satagay_data$H, densfun = "gamma")



#####
## Function for the calculation of the potential modern origin of charcoal #####
#####

origin  <- function(spd, dir, pih, diam, lon_orig, lat_orig){
  vt <- (((0.5-0.00127)*981*diam^2)/(18*(1.8*10^-4)))/100 # cm/s in m/s
  t  <- pih/(vt)
  delta_x <- t*spd                                        
  geosphere::destPoint(c(lon_orig, lat_orig), 
                       (dir+180)%%360, delta_x)
}


## In this study, it was assumed a spherical particle shape with a density of 0.5 g/cm3
## (like Clark, 1988; Gilgen et al., 2018) and a diameter range of 150 to 500 μm 
## (0.015 to 0.050 cm) based on the sampled lake sediment records (Glückler et al. 2021)
charc_diam_cm <- seq(0.015, 0.05, length = 100) 



#####
## Calculation of the travel distance ##########################################
#####


sumCrds_Satagay_VR <- do.call("rbind", lapply(unique(wndTab_Satagay$date), function(d) { 
  
cat(sprintf('\rDate %s of %s  ', which(d == unique(wndTab_Satagay %>% pull(date) %>% as.Date())),
length(unique(wndTab_Satagay %>% pull(date) %>% as.Date()))))
  
tmp <- wndTab_Satagay %>% filter(date == d) %>% ungroup() %>% st_as_sf(coords = c("lon", "lat")) %>%
       st_set_crs(4326) %>% dplyr::select(wind_spd, wind_dir)
  
tmp <- tmp[apply(st_distance(tmp, data_coord_lakes       %>% 
       st_as_sf(coords = c("lon", "lat"))                %>% 
       st_set_crs(4326)), 2, function(x) which.min(x)),] %>%
       st_drop_geometry()
  
do.call("rbind", lapply(1:20, function(r) {
        tibble(lake = 2,
        as.data.frame(origin(spd  = tmp$wind_spd[2],
                             dir  = tmp$wind_dir[2],
                             pih  = rgamma(1, H_distr_Satagay$estimate[1], H_distr_Satagay$estimate[2]),
                             diam = charc_diam_cm[sample(1:length(charc_diam_cm),1)],
                             lon_orig = data_coord_lakes$lon[2],
                             lat_orig = data_coord_lakes$lat[2]))) %>% setNames(c("lake", "lon", "lat"))
})) %>% dplyr::mutate(date = d) %>% dplyr::select(date, lake, lon, lat)
}))
save(sumCrds_Satagay_VR, file = "Results/Simulation/data/sumCrds_Satagay_VR.rda") 
sumCrds_Satagay_VR <- get(load("Results/Simulation/data/sumCrds_Satagay_VR.rda"))



#####
## Calculation of the distance between detecteed potential charcoal source areas
## and the sampled lake to determine the actual charcoal dispersion distances
#####

p_Satagay <- sumCrds_Satagay_VR %>% filter(!is.na(lon) & !is.na(lat)) %>% 
             st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326)

dist_data_Satagay <- do.call("rbind", lapply(unique(p_Satagay$date), function(x) { 
  
cat(sprintf('\rDate %s of %s  ', which(x == unique(p_Satagay %>% pull(date) %>% as.Date())),
length(unique(p_Satagay %>% pull(date) %>% as.Date()))))
  
              tmp <- p_Satagay %>% filter(date == x)
  
              tmp %>% dplyr::mutate(dist = as.numeric(st_distance(geometry, (Satagay_buf_100 %>% 
              st_centroid() %>% st_transform(4326)))))
            
}))
save(dist_data_Satagay, file = "Results/Simulation/data/dist_data_Satagay.rda")
dist_data_Satagay <- get(load("Results/Simulation/data/dist_data_Satagay.rda"))



#####
## Visualization of the horizontal travel distances of spherical particles in the air in form of histogram
##### 


## 100 breaks | 20 km
png(glue("Results/Simulation/H/Satagay/dist/travel_distance_Satagay_VR_H_100_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Satagay$dist[dist_data_Satagay$lake==2 & dist_data_Satagay$dist<400000]/1000, breaks= 100, 
     xlim = c(0,20), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey80")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()

## 100 breaks | 40 km
png(glue("Results/Simulation/H/Satagay/dist/travel_distance_Satagay_VR_H_100_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Satagay$dist[dist_data_Satagay$lake==2 & dist_data_Satagay$dist<400000]/1000, breaks= 100, 
     xlim = c(0,40), ylim = c(0, 4000), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey80")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks | 20 km
png(glue("Results/Simulation/H/Satagay/dist/travel_distance_Satagay_VR_H_200_breaks_20km.png"), width = 1000, height = 700)
hist(dist_data_Satagay$dist[dist_data_Satagay$lake==2 & dist_data_Satagay$dist<400000]/1000, breaks= 200, 
     xlim = c(0,20), ylim = c(0, 2500), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey80")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()

## 200 breaks | 40 km
png(glue("Results/Simulation/H/Satagay/dist/travel_distance_Satagay_VR_H_200_breaks_40km.png"), width = 1000, height = 700)
hist(dist_data_Satagay$dist[dist_data_Satagay$lake==2 & dist_data_Satagay$dist<400000]/1000, breaks= 200, 
     xlim = c(0,40), ylim = c(0, 2500), 
     main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay",
     xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey80")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()


#####
## Boxplot of calculated PIH's #################################################
#####


n_S_dist <- length(dist_data_Satagay$dist) # 9840
summary(dist_data_Satagay$dist/1000) ## in km 
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  0.01127   1.21230   2.80081   5.74208   6.52604 151.50041  
sd(dist_data_Satagay$dist/1000)

## Boxplot
png(glue("Results/Simulation/H/Satagay/dist/travel_distance_Satagay_VR_H_boxplot.png"), width = 900, height = 800)
boxplot(x = dist_data_Satagay$dist/1000, outline = F, ylim = c(0,20),
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere for the study area Lake Satagay", 
        ylab = "Δx [km]", xlab = glue("n = {n_S_dist}"), cex.main=1.7, cex.lab=1.5, cex.axis=1.5,col = "grey80")
legend("topright", legend=c("Calculated PIH"), bty = "n", cex = 1.5)
dev.off()


#####
## Boxplot dist Khamra and Satagay parallel ####################################
#####

dist_data_Khamra_H  <- get(load("Results/Simulation/data/dist_data_Khamra.rda"))
dist_data_Satagay_H <- get(load("Results/Simulation/data/dist_data_Satagay.rda"))
dist_data_Satagay_H <- na.omit(dist_data_Satagay)

n_dist_K <- length(dist_data_Khamra_H$dist)
n_dist_S <- length(dist_data_Satagay_H$dist)

dist_df <- merge(data.frame(dist_data_Khamra_H, row.names=NULL), data.frame(dist_data_Satagay_H, row.names=NULL), 
              by = 0, all = TRUE)[-1]

label_names_dist     <- c(glue("n = {n_dist_K}"), glue("n = {n_dist_S}"))
colors_dist_boxplot  <- c("grey40", "grey80")

png(glue("Results/Simulation/H/boxplot_H_Khamra_Satagay.png"), width = 800, height = 800)
boxplot(dist_df$dist.x/1000, dist_df$dist.y/1000, ylim = c(0,20),  outline=FALSE, names = label_names_dist, col = colors_dist_boxplot, 
        main = "Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density)\ninto the atmosphere with calculated PIH [m]", 
        cex.main = 1.7, ylab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5)
legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col= c("grey40", "grey80"), pt.cex=3, pch=15, cex = 1.5)
dev.off()


#####
## Multiple plot for Khamra and Satagay ########################################
#####

hist_dist_S_H <- hist(dist_data_Satagay_H$dist[dist_data_Satagay_H$lake== 2 & dist_data_Satagay_H$dist<400000]/1000, breaks =100,
                       xlim = c(0,40), ylim = c(0, 6000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey80")

hist_dist_K_H <- hist(dist_data_Khamra_H$dist[dist_data_Khamra_H$lake==1 & dist_data_Khamra_H$dist<400000]/1000,breaks = 100,
                       xlim = c(0,40), ylim = c(0, 6000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "grey40")


png(glue("Results/Simulation/H/travel_distance_H_K_S_multipe.png"), width = 600, height = 800)
layout(matrix(c(1,1,2,2),nrow = 2, ncol = 2, byrow = TRUE))
boxplot(dist_df$dist.x/1000, dist_df$dist.y/1000, ylim = c(0,40),  outline=FALSE, names = label_names_dist, col = colors_dist_boxplot, 
        main = NULL, ylab = "Δx [km]", cex.lab = 1.75, cex.axis = 1.75)

plot(hist_dist_K_H, col = "grey40", main = NULL, xlab = "Δx [km]", cex.lab = 1.75, cex.axis = 1.75,  xlim = c(0, 40), ylim = c(0,6000))
plot(hist_dist_S_H, col = "grey80", add = T, main = NULL, xlab = "Δx [km]", cex.lab = 1.75, xlim = c(0, 40), ylim = c(0,6000))

legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col= c("grey40", "grey80"), pt.cex=3, pch=15, cex = 2)
dev.off()


png(glue("Results/Simulation/H/travel_distance_H_K_S_multipe_20.png"), width = 600, height = 800)
layout(matrix(c(1,1,2,2),nrow = 2, ncol = 2, byrow = TRUE))
boxplot(dist_df$dist.x/1000, dist_df$dist.y/1000, ylim = c(0,20),  outline=FALSE, col = colors_dist_boxplot, 
        main = NULL, ylab = "Δx [km]", cex.lab = 1.75, cex.axis = 1.75)

plot(hist_dist_K_H, col = "grey40", main = NULL, xlab = "Δx [km]", cex.lab = 1.75, cex.axis = 1.75,  xlim = c(0, 40), ylim = c(0,6000))
plot(hist_dist_S_H, col = "grey80", add = T, main = NULL, xlab = "Δx [km]", cex.lab = 1.75, xlim = c(0, 40), ylim = c(0,6000))

legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col= c("grey40", "grey80"), pt.cex=3, pch=15, cex = 2)
dev.off()


## Crop the tibble with 100 km radius box around the lake
dist_data_crop_Satagay <- st_crop(dist_data_Satagay, Satagay_buf_extent)


## Select the FRP-range
MODIS_Satagay_VR       <- H_function_Satagay %>% filter(FRP > 0 & FRP <= 71653768.9882488) 
## 300 MW = 71653768.9882488 cal/s



#####
## Visualization of the backward simulation of potential source areas of charcoal particles
#####

## FRP (cal/s)
png(glue("Results/Simulation/H/Satagay/psource/Satagay_psource_VR_H_cal.png"), width = 1000, height = 900)
ggplot(Satagay_buf_extent) +
  geom_sf(data = dist_data_crop_Satagay, aes(colour  = "Locations of\np-source areas\n(calculated PIH)"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(calculated PIH)" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Satagay_VR, aes(color = FRP), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(color = "FRP [cal/s]") +
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
  theme(legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()


## FRP (MW)
MODIS_Satagay_VR$FRP_MW <- MODIS_Satagay_VR$FRP/238845.8966275

png(glue("Results/Simulation/H/Satagay/psource/Satagay_psource_VR_H_MW.png"), width = 1000, height = 900)
ggplot(Satagay_buf_extent) +
  geom_sf(data = dist_data_crop_Satagay, aes(colour  = "Locations of\np-source areas\n(calculated PIH)"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(calculated PIH)" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Satagay_VR, aes(color = FRP_MW), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(color = "FRP [MW]") +
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
  theme(legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()


#####
## Convert to raster and plot p-source areas again ##############################
#####


## Create an empty raster 
r0_Satagay <- raster::raster(raster::extent(as(data_coord_lakes[2,2:3] %>% 
                                                st_as_sf(coords = c("lon", "lat")) %>% st_buffer(0.5)    %>%   # 0.5 degress around the lake
                                                st_bbox() %>% st_as_sfc(), "Spatial")), nrow = 50, ncol = 50) 

## Put the data into the empty raster
raster_Satagay     <- rasterize(dist_data_Satagay %>% st_coordinates(), r0_Satagay, fun='count')
crs(raster_Satagay) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

## Convert the raster to data.frame
df_r_S_VR     <- as.data.frame(raster_Satagay, xy = TRUE)
df_r_S_na_VR  <- na.omit(df_r_S_VR)

## Rasterize the fire radiative power data and create a MULTIPOLYGON
FRP_Satagay_VR <- rasterize(H_function_Satagay %>% filter(!is.na(FRP)) %>% st_coordinates(),
                            raster_Satagay, fun = 'count') %>% raster::rasterToPolygons() %>%
                            st_as_sf() %>% filter(layer > 0) %>% st_union() %>% 
                            st_set_crs(4326) %>% st_transform(4326)

## Normalize the p-source	probabilities  
df_r_S_na_VR$normalized <- (df_r_S_na_VR$layer-min(df_r_S_na_VR$layer))/(max(df_r_S_na_VR$layer)-min(df_r_S_na_VR$layer))

## Set crs 4326 for Khamra catchment
Satagay_catchment <- Satagay_catchment %>% st_set_crs(4326) %>% st_transform(crs=4326)



#####
## Visualization ################################################################
#####

png(glue("Results/Simulation/H/Satagay/psource/Satagay_psource_raster_VR_H.png"), width = 1500, height = 1000)
ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_na_VR , aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\np-source areas\n(calculated PIH)")+
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.4) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill") +
  geom_sf(data = FRP_Satagay_VR, mapping = aes(fill = "Burned area"), alpha = 0.2, colour = NA) +
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
  theme(legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()


#####
## Mask the p-source raster data with modern FRP data ###########################
#####

FRP_sp_Satagay     <- as(FRP_Satagay_VR, "Spatial")  ## MULTIPOLYGON to SpatialPolygons
r_mask_Satagay     <- raster::mask(raster_Satagay, FRP_sp_Satagay)

df_r_S_VR_mask     <- as.data.frame(r_mask_Satagay, xy = TRUE)
df_r_S_VR_mask     <- na.omit(df_r_S_VR_mask)
FRP_Khamra_VR_mask <- rasterize(H_function_Satagay %>% filter(!is.na(FRP)) %>% st_coordinates(),
                                r_mask_Satagay, fun = 'count') %>% raster::rasterToPolygons() %>%
                                st_as_sf() %>% filter(layer > 0) %>% st_union()    %>% 
                                st_set_crs(4326) %>% st_transform(4326)

## Normalization
df_r_S_VR_mask$normalized <- (df_r_S_VR_mask$layer-min(df_r_S_VR_mask$layer))/(max(df_r_S_VR_mask$layer)-min(df_r_S_VR_mask$layer))



#####
## Visualization ################################################################
#####

png(glue("Results/Simulation/H/Satagay/psource/Satagay_psource_raster_VR_H_extract.png"), width = 1500, height = 1000)
ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_VR_mask, aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\nactual source areas\n(calculated PIH)")+
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
  theme(legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))
dev.off()


#####
## Mapping calculated PIH [m] for Lake Khamra and Satagay ######################
#####

H_function_Khamra_2500  <- H_function_Khamra %>% filter(H <= 2000)
H_function_Satagay_2500 <- H_function_Satagay %>% filter(H <= 2000)

H_points_Khamra  <- st_as_sf(H_function_Khamra_2500)
H_points_Satagay <- st_as_sf(H_function_Satagay_2500)

png(glue("Results/Simulation/H/Khamra/H/study_area_H_Khamra_2000.png"), width = 1000, height = 800)
ggplot(Khamra_buf_extent) +
  geom_sf(data = H_points_Khamra, aes(color = H), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limit = c(1,2000))+
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.4) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Khamra, aes(fill = "Lake Khamra"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake Khamra" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  xlab("") +
  ylab("") +
  labs(title = "Calculated plume injection heights for study area Lake Khamra",
       color = "Calculated PIH [m]", shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 16, vjust = 0.5),
        legend.text   = element_text(size = 16, vjust = 0.75))
dev.off()


png(glue("Results/Simulation/H/Satagay/H/study_area_H_Satagay_2000.png"), width = 1000, height = 800)
ggplot(Satagay_buf_extent) +
  geom_sf(data = H_points_Satagay, aes(color = H), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limit = c(1,2000)) +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.8) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake Satagay"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake Satagay" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  xlab("") +
  ylab("") +
  labs(title = "Calculated plume injection heights for study area Lake Satagay",
       color = "Calculated PIH [m]", shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 16, vjust = 0.5),
        legend.text   = element_text(size = 16, vjust = 0.75))
dev.off()
