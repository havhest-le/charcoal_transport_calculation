#####
## Script to select the data set "Injection_Height" with daily values #
## from MAIAC MODIS data product, calculate the different quantiles ###
## and visualize the plume injection height distribution ##############
#####

## rm(list = ls(all= TRUE))

#####
## Load packages ###############################################################
#####
library(base)
library(dplyr)
library(tidyverse)
library(sf)
library(sp)
library(stars)
library(stats)
library(terra)
library(glue)
library(tibble)

## optional
library(raster)

## Visualization
library(ggnewscale)
library(ggplot2)
library(graphics)

library(maps)
library(viridis)
library(RColorBrewer)
library(colorspace)


#####
## Load and clean data #########################################################
#####

## Study area
buffer_lakes_200_UTM <- data.frame(location = c("Khamra","Satagay"),
                                   lon = c(112.98, 117.998),
                                   lat = c(59.99, 63.078)) %>%
                                   st_as_sf(coords = c("lon", "lat")) %>%
                                   st_set_crs(4326) %>%
                                   st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") %>%
                                   st_buffer(200000)

buffer_lakes_200     <- data.frame(location = c("Khamra","Satagay"),
                                   lon = c(112.98, 117.998),
                                   lat = c(59.99, 63.078)) %>%
                                   st_as_sf(coords = c("lon", "lat")) %>%
                                   st_set_crs(4326) %>%
                                   st_buffer(200000)

## Lake Khamra
Khamra                  <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp") %>% st_transform(4326)
proj_K                  <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Khamra))[,1]}
                                            +lat_0={st_coordinates(st_centroid(Khamra))[,2]}")
Khamra_catchment        <- read_sf("Data/Khamra/khamra_catchment.shp") %>% st_transform(4326)
Khamra_buf_100          <- st_transform(Khamra, crs = CRS(proj_K)) %>%
                           st_buffer(100000) %>% st_transform(4326)
Khamra_buf_extent       <- as(extent(st_bbox(Khamra_buf_100 %>%
                           st_transform(4326)               %>%
                           st_shift_longitude())[c(1,3,2,4)]), "SpatialPolygons")
crs(Khamra_buf_extent)  <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

## Lake Satagay 2.0
Satagay                 <- read_sf("Data/Satagay/satagay_poly_1.shp") %>% st_transform(4326)
proj_S                  <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Satagay))[,1]}
                                      +lat_0={st_coordinates(st_centroid(Satagay))[,2]}")
Satagay_catchment       <- read_sf("Data/Satagay/satagay_catchment.shp") %>% st_transform(4326)
Satagay_buf_100         <- st_transform(Satagay, crs = CRS(proj_S)) %>%
                           st_buffer(100000) %>% st_transform(4326)
Satagay_buf_extent      <- as(extent(st_bbox(Satagay_buf_100 %>%
                           st_transform(4326)                %>%
                           st_shift_longitude())[c(1,3,2,4)]), "SpatialPolygons")
crs(Satagay_buf_extent) <- crs("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")


## Fire radiative power [FRP] ##################################################
sf_MODIS_FRP_Khamra_fy  <- get(load("Results/FRP/FRP_Khamra_points_box.RData")) %>%
                           filter(Year == "2001" | Year == "2003" |
                                  Year == "2004" | Year == "2007" |
                                  Year == "2013" | Year == "2014" )

sf_MODIS_FRP_Satagay_fy <- get(load("Results/FRP/FRP_Satagay_points_box.RData")) %>%
                           filter(Year == "2002" | Year == "2013" |
                                  Year == "2014" | Year == "2018" )


## Plume injection height [PIH] ################################################
PIH_files <- tibble(path = list.files(path = "Data/PIH/hdf4/", pattern = "hdf", full.names = T)) %>%
             mutate(extent = substr(sapply(strsplit(path, "/"), function(x) x[length(x)]), 1, 1),
                    date   = as.POSIXct(substr(sapply(strsplit(path, "/"), function(x) x[length(x)]), 3, 12))) %>%
                    arrange(date)


#####
## Get the PIH data for each lon lat coordinate, where FRP is located ##########
#####

## To establish the spatial relationship between detected active fires and observed PIH′s, the
## PIH values were assigned to the FRP points. A 450 m radius was created around the FRP
## points and the nearest neighbor method was used to assign the nearest PIH values to the FRP points.

## Study area Lake Khamra 
sf_MODIS_PIH_FRP_Khamra <- sf_MODIS_FRP_Khamra_fy %>% mutate(PIH = NA)

## Create an empty tibble
outMODIS_PIH_Kahmra_geom <- tibble(index = NA, TIME = NA, FRP = NA, `2.5%` = NA, `20%` = NA, `50%` = NA,  `80%` = NA, `97.5%` = NA)

for(date in unique(sf_MODIS_FRP_Khamra_fy %>% pull(TIME) %>% as.Date())) {

cat(sprintf('\rDate %s of %s  ', which(date == unique(sf_MODIS_FRP_Khamra_fy %>% pull(TIME) %>% as.Date())),
           length(unique(sf_MODIS_FRP_Khamra_fy %>% pull(TIME) %>% as.Date()))))

tmp_modis <- sf_MODIS_PIH_FRP_Khamra %>% rownames_to_column(var = "index") %>% filter(as.Date(TIME) == date)

if(length(PIH_files %>% filter(as.Date(date)==as.Date(tmp_modis$TIME[1])) %>% pull(path)) > 0) {

      ras <- do.call("st_mosaic", lapply(PIH_files %>% filter(as.Date(date) == as.Date(tmp_modis$TIME[1])) %>% pull(path),
             function(x) terra::rast(x, subds=8) %>% st_as_stars()))
             ## get the data layer "Injection height" from the layer position 8

     out <- tryCatch(st_extract(ras, tmp_modis %>% st_transform(st_crs(ras)) 
                                %>% st_buffer(45000)),    ## Create a 450 m radius around the FRP points
                                error = function(e) NULL)

rasExtr <- tryCatch(ras[tmp_modis %>% st_transform(st_crs(ras)) %>% st_buffer(4500) %>% st_union(),,,], error = function(e) NULL)

           if(!is.null(rasExtr)) {
           outMODIS_PIH_Kahmra_geom <- outMODIS_PIH_Kahmra_geom %>% bind_rows(
           tmp_modis %>% dplyr::select(index, TIME, FRP, Year)  %>% bind_cols(do.call("rbind", lapply(1:nrow(tmp_modis), function(x) {
           rasExtr <- tryCatch(ras[tmp_modis[x,] %>% st_transform(st_crs(ras)) 
                                                 %>% st_buffer(50000) ## Create a 500 m radius around the FRP points
                                                 %>% st_union(),,,], error = function(e) NULL)
           if(is.null(rasExtr)) quantile(c(NA, NA), probs = c(0.025, 0.2, 0.5, 0.8, 0.975), na.rm = T) else {
           ## Calculate the quantiles Q[0.025], Q[0.2], Q[0.5], Q[0.8], and Q[0.975] of PIH
           quantile(c(rasExtr[[1]][]), probs = c(0.025, 0.2, 0.5, 0.8, 0.975), na.rm = T) 
          }
          }))) %>% filter(!apply(., 1, function(y) all(is.na(y[-c(1:3)])))) ## Put the values into the created empty tibble
          )
          }
           if(!is.null(out)) sf_MODIS_PIH_FRP_Khamra$PIH[as.numeric(tmp_modis %>% pull(index))] <- apply(out[[1]], 1, max, na.rm = T)
          }
}
save(outMODIS_PIH_Kahmra_geom, file = "Results/Simulation/data/outMODIS_PIH_Kahmra_geom.rda")
outMODIS_PIH_Kahmra_geom    <- get(load("Results/Simulation/data/outMODIS_PIH_Kahmra_geom.rda"))
outMODIS_PIH_Kahmra_geom_na <- na.omit(outMODIS_PIH_Kahmra_geom)



#####
## Create a range between Q[0.2] and Q[0.8] of PIH values ######################
#####

## To compare the transport results of PIH quantile ranges with single
## quantiles, an additional matrix of self-generated PIH value sequences between Q[0.2] and
## Q[0.8] with a length of 1000 was created

length(outMODIS_PIH_Kahmra_geom_na$`20%`) ## 3372

## Empty matrix with a length of 1000
mat_Khamra = matrix(nrow = 3372, ncol = 1000)

## Fill the matrix with values
for(i in c(1:3372)){
  seq <- seq(from = outMODIS_PIH_Kahmra_geom_na$`20%`[i], to = outMODIS_PIH_Kahmra_geom_na$`80%`[i], length.out = 1000)
  mat_Khamra[i,] <- seq
  mat_Khamra
}

save(mat_Khamra, file = "Results/Simulation/data/PIH_range_Khamra.rda")
PIH_range_Khamra <- get(load("Results/Simulation/data/PIH_range_Khamra.rda"))




## Study area Lake Satagay 
sf_MODIS_PIH_FRP_Satagay  <- sf_MODIS_FRP_Satagay_fy %>% mutate(PIH = NA)

outMODIS_PIH_Satagay_geom <- tibble(index = NA, TIME = NA, FRP = NA, `2.5%` = NA, `20%` = NA, `50%` = NA,  `80%` = NA, `97.5%` = NA)

for(date in unique(sf_MODIS_FRP_Satagay_fy %>% pull(TIME) %>% as.Date())) {

  cat(sprintf('\rDate %s of %s  ', which(date == unique(sf_MODIS_FRP_Satagay_fy %>% pull(TIME) %>% as.Date())),
              length(unique(sf_MODIS_FRP_Satagay_fy %>% pull(TIME) %>% as.Date()))))

tmp_modis <- sf_MODIS_PIH_FRP_Satagay %>% rownames_to_column(var = "index") %>% filter(as.Date(TIME) == date)

             if(length(PIH_files %>% filter(as.Date(date)==as.Date(tmp_modis$TIME[1])) %>% pull(path)) > 0) {


      ras <- do.call("st_mosaic", lapply(PIH_files %>% filter(as.Date(date) == as.Date(tmp_modis$TIME[1])) %>% pull(path),
                     function(x) terra::rast(x, subds=8) %>% st_as_stars()))


      out <- tryCatch(st_extract(ras, tmp_modis %>% st_transform(st_crs(ras)) %>% st_buffer(45000)), 
             error = function(e) NULL)

  rasExtr <- tryCatch(ras[tmp_modis %>% st_transform(st_crs(ras)) %>% st_buffer(4500) %>% st_union(),,,], error = function(e) NULL)

             if(!is.null(rasExtr)) {
             outMODIS_PIH_Satagay_geom <- outMODIS_PIH_Satagay_geom %>% bind_rows(
             tmp_modis %>% dplyr::select(index, TIME, FRP, Year)    %>% bind_cols(do.call("rbind", lapply(1:nrow(tmp_modis), function(x) {
             rasExtr <- tryCatch(ras[tmp_modis[x,] %>% st_transform(st_crs(ras)) 
                                                   %>% st_buffer(50000)
                                                   %>% st_union(),,,], error = function(e) NULL)
             if(is.null(rasExtr)) quantile(c(NA, NA), probs = c(0.025, 0.2, 0.5, 0.8, 0.975), na.rm = T) else {
             quantile(c(rasExtr[[1]][]), probs = c(0.025, 0.2, 0.5, 0.8, 0.975), na.rm = T)
             }
             }))) %>% filter(!apply(., 1, function(y) all(is.na(y[-c(1:3)]))))
             )
             }
             if(!is.null(out)) sf_MODIS_PIH_FRP_Satagay$PIH[as.numeric(tmp_modis %>% pull(index))] <- apply(out[[1]], 1, max, na.rm = T)
             }
}

save(outMODIS_PIH_Satagay_geom, file = "Results/Simulation/data/outMODIS_PIH_Satagay_geom.rda")
outMODIS_PIH_Satagay_geom    <- get(load("Results/Simulation/data/outMODIS_PIH_Satagay_geom.rda"))
outMODIS_PIH_Satagay_geom_na <- na.omit(outMODIS_PIH_Satagay_geom)


## Create a range between Q[0.2] and Q[0.8] of PIH values ########################
length(outMODIS_PIH_Satagay_geom_na$`20%`) ## 16956

## Empty matrix
mat_Satagay = matrix(nrow = 16956, ncol = 1000)

## Fill matrix with values
for(i in c(1:16956)){
  seq <- seq(from = outMODIS_PIH_Satagay_geom_na$`20%`[i], to = outMODIS_PIH_Satagay_geom_na$`80%`[i], length.out = 1000)
  mat_Satagay[i,] <- seq
  mat_Satagay
}

save(mat_Satagay, file = "Results/Simulation/data/PIH_range_Satagay.rda")
PIH_range_Satagay <- get(load("Results/Simulation/data/PIH_range_Satagay.rda"))



#####
## Distribution of the PIH values for each individual quantile #################
#####

summary(outMODIS_PIH_Kahmra_geom_na$`80%`)  ## max 1802.7 m
summary(outMODIS_PIH_Satagay_geom_na$`80%`) ## max 2139.437 m

summary(outMODIS_PIH_Kahmra_geom_na$`97.5%`)  ## max 2210.21 m
summary(outMODIS_PIH_Satagay_geom_na$`97.5%`) ## max 2422.072 m 

hist(outMODIS_PIH_Kahmra_geom_na$`50%`)    ## max 1389.80 m
summary(outMODIS_PIH_Kahmra_geom_na$`50%`)

hist(outMODIS_PIH_Satagay_geom_na$`50%`)
summary(outMODIS_PIH_Satagay_geom_na$`50%`) ## max 1940.993 m 



## Visualization ###############################################################


## Study area Lake Khamra
length(outMODIS_PIH_Kahmra_geom_na$`50%`) # [1]  3372

png(glue("Results/Simulation/PIH/Khamra/2.5/Measured_PIH_Khamra_025.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Kahmra_geom_na$`2.5%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.025]",
     main = "Measured plume injection heights (n = 3372) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Khamra/20/Measured_PIH_Khamra_20.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Kahmra_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.2]",
     main = "Measured plume injection heights (n = 3372) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Khamra/50/Measured_PIH_Khamra_50.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Kahmra_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.5]",
     main = "Measured plume injection heights (n = 3372) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Khamra/80/Measured_PIH_Khamra_80.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Kahmra_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.8]",
     main = "Measured plume injection heights (n = 3372) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Khamra/97.5/Measured_PIH_Khamra_975.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Kahmra_geom_na$`97.5%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.975]",
     main = "Measured plume injection heights (n = 3372) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Khamra/range/Measured_PIH_Khamra_range_800000.png"), width = 1000, height = 700)
hist(PIH_range_Khamra, xlim = c(0, 2000), ylim = c(0,800000),
     breaks= 100, xlab = "Measured PIH [m] | Range between Q[0.2] and Q[0.8]",
     main = "Measured plume injection heights (n = 3372) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Khamra/range/Measured_PIH_Khamra_range_200000.png"), width = 1000, height = 700)
hist(PIH_range_Khamra, xlim = c(0, 2000), ylim = c(0,200000),
     breaks= 100, xlab = "Measured PIH [m] | Range between Q[0.2] and Q[0.8]",
     main = "Measured plume injection heights (n = 3372) for the study area Lake Khamra", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()



## Multiple plot ###############################################################
png(glue("Results/Simulation/PIH/Khamra/Measured_PIH_Khamra.png"), width = 1000, height = 700)
layout(matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE))

hist(outMODIS_PIH_Kahmra_geom_na$`2.5%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.025]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`97.5%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.975]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.2]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.8]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,1200),
     breaks= 100, xlab = "Measured PIH [m]",  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.975]", line = -1)

mtext("Measured plume injection heights (n = 3372) for the study area Lake Khamra", side = 3, line = - 2, outer = TRUE, cex = 1.5)
dev.off()

 
## Multiple plot ################################################################
png(glue("Results/Simulation/PIH/Khamra/Measured_PIH_Khamra_300.png"), width = 1000, height = 700)
layout(matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE))

hist(outMODIS_PIH_Kahmra_geom_na$`2.5%`, xlim = c(0, 2000), ylim = c(0,300),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.025]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`97.5%`, xlim = c(0, 2000), ylim = c(0,300),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.975]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,300),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.2]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,300),
     breaks= 100, xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.8]", line = -1)

hist(outMODIS_PIH_Kahmra_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,300),
     breaks= 100, xlab = "Measured PIH [m]",  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.975]", line = -1)

mtext("Measured plume injection heights (n = 3372) for the study area Lake Khamra", side = 3, line = - 2, outer = TRUE, cex = 1.5)

dev.off()


################################################################################

## Study area Lake Satagay
length(outMODIS_PIH_Satagay_geom_na$`50%`) # [1] 16956

png(glue("Results/Simulation/PIH/Satagay/2.5/Measured_PIH_Satagay_025.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Satagay_geom_na$`2.5%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.025]",
     main = "Measured plume injection heights (n = 16956) for the study area Lake Satagay", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Satagay/20/Measured_PIH_Satagay_20.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Satagay_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.2]",
     main = "Measured plume injection heights (n = 16956) for the study area Lake Satagay", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Satagay/50/Measured_PIH_Satagay_50.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Satagay_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.5]",
     main = "Measured plume injection heights (n = 16956) for the study area Lake Satagay", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Satagay/80/Measured_PIH_Satagay_80.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Satagay_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.8]",
     main = "Measured plume injection heights (n = 16956) for the study area Lake Satagay", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Satagay/97.5/Measured_PIH_Satagay_975.png"), width = 1000, height = 700)
hist(outMODIS_PIH_Satagay_geom_na$`97.5%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = "Measured PIH [m] | Q[0.975]",
     main = "Measured plume injection heights (n = 16956) for the study area Lake Satagay", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()

png(glue("Results/Simulation/PIH/Satagay/range/Measured_PIH_Satagay_range_800000.png"), width = 1000, height = 700)
hist(PIH_range_Satagay, xlim = c(0, 2000), ylim = c(0,800000),
     breaks= 100, xlab = "Measured PIH [m] | Range between Q[0.2] and Q[0.8]",
     main = "Measured plume injection heights (n = 16956) for the study area Lake Satagay", cex.main = 2,
     cex.lab = 1.5, col = "lightblue3") # in m
dev.off()


## Multiple plot ###############################################################
png(glue("Results/Simulation/PIH/Satagay/Measured_PIH_Satagay_7000.png"), width = 1000, height = 700)
layout(matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE))

hist(outMODIS_PIH_Satagay_geom_na$`2.5%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.025]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`97.5%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.975]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.2]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.8]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,7000),
     breaks= 100, xlab = "Measured PIH [m]",  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.5]", line = -1)

mtext("Measured plume injection heights (n = 16956) for the study area Lake Satagay", side = 3, line = - 2, outer = TRUE, cex = 1.5)

dev.off()

## ylim = c(0,2500) ############################################################
png(glue("Results/Simulation/PIH/Satagay/Measured_PIH_Satagay_2500.png"), width = 1000, height = 700)
layout(matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE))

hist(outMODIS_PIH_Satagay_geom_na$`2.5%`, xlim = c(0, 2000), ylim = c(0,2500),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.025]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`97.5%`, xlim = c(0, 2000), ylim = c(0,2500),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.975]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,2500),
     breaks= 100, xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.2]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,2500),
     breaks= 100, xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.8]", line = -1)

hist(outMODIS_PIH_Satagay_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,2500),
     breaks= 100, xlab = "Measured PIH [m]",  main = NULL,cex.lab = 1.5, col = "lightblue3")
title(main = "Q[0.5]", line = -1)

mtext("Measured plume injection heights (n = 16956) for the study area Lake Satagay", side = 3, line = - 2, outer = TRUE, cex = 1.5)

dev.off()





#####
## Overlapped hist() ###########################################################
#####

## Satagay
h_02_S    <- hist(outMODIS_PIH_Satagay_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,8000),
             xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
h_05_S    <- hist(outMODIS_PIH_Satagay_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,8000),
              xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
h_08_S    <- hist(outMODIS_PIH_Satagay_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,8000),
             xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
h_range_S <- hist(PIH_range_Satagay, xlim = c(0, 2000), ylim = c(0,4000000),
             xlab = "Measured PIH [m] | Range between Q[0.2] and Q[0.8]",
             main = NULL,cex.lab = 1.5, col = "lightblue3")

## Khamra
h_02_K <- hist(outMODIS_PIH_Kahmra_geom_na$`20%`, xlim = c(0, 2000), ylim = c(0,2000),
               xlab = NULL, main = NULL, cex.lab = 1.5, col = "lightblue3")
h_05_K <- hist(outMODIS_PIH_Kahmra_geom_na$`50%`, xlim = c(0, 2000), ylim = c(0,2000),
               xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
h_08_K <- hist(outMODIS_PIH_Kahmra_geom_na$`80%`, xlim = c(0, 2000), ylim = c(0,2000),
               xlab = NULL,  main = NULL,cex.lab = 1.5, col = "lightblue3")
h_range_K <- hist(PIH_range_Khamra, xlim = c(0, 2000), ylim = c(0,4000000),
                  xlab = "Measured PIH [m] | Range between Q[0.2] and Q[0.8]",
                  main = NULL,cex.lab = 1.5, col = "lightblue3")


#####
##  Boxplots of PIH values for Lake Khamra and Satagay #########################
#####

## Study area Lake Khamra
n_PIH_K <- length(outMODIS_PIH_Kahmra_geom_na$`20%`)

label_names_PIH_K     <- c(glue("n = {n_PIH_K}"), glue("n = {n_PIH_K}"), glue("n = {n_PIH_K}"))
colors_PIH_boxplot  <- c(rgb(0,0,1,1/4), col = rgb(0,1,1,1/4), rgb(1,0,1,1/4))

## ylim = c(0,300)
png(glue("Results/Simulation/PIH/Khamra/boxplots_histo_PIH_K.png"), width = 600, height = 800)
layout(matrix(c(1,1,2,2),nrow = 2, ncol = 2, byrow = TRUE))

boxplot(outMODIS_PIH_Kahmra_geom_na$`20%`, outMODIS_PIH_Kahmra_geom_na$`50%`, outMODIS_PIH_Kahmra_geom_na$`80%`,
        col = colors_PIH_boxplot, outline = FALSE, 
        ylim =c(0,2000), cex.axis = 1.75, 
        ylab = "Measured PIH [m]", cex.lab = 1.75, names = label_names_PIH_K)

plot(h_02_K, col = rgb(0,0,1,1/4), main = NULL, xlab = "Measured PIH [m]", 
     cex.lab = 1.75, cex.axis = 1.75,  xlim = c(0, 2000), ylim = c(0,2000))
plot(h_05_K, col = rgb(0,1,1,1/4), add = T, main = NULL, xlab = "Measured PIH [m]", cex.lab = 1.75)
plot(h_08_K, col = rgb(1,0,1,1/4), add = T, main = NULL, xlab = "Measured PIH [m]", cex.lab = 1.75)

legend("topright", legend = c("Q[0.2]","Q[0.5]", "Q[0.8]"), 
       col= colors_PIH_boxplot, pt.cex=3, pch=15, cex = 2)
dev.off()

################################################################################

## Study area Lake Satagay
n_PIH_S           <- length(outMODIS_PIH_Satagay_geom_na$`20%`)
label_names_PIH_S <- c(glue("n = {n_PIH_S}"), glue("n = {n_PIH_S}"), glue("n = {n_PIH_S}"))


png(glue("Results/Simulation/PIH/Satagay/boxplots_histo_PIH_S.png"), width = 600, height = 800)
layout(matrix(c(1,1,2,2),nrow = 2, ncol = 2, byrow = TRUE))

boxplot(outMODIS_PIH_Satagay_geom_na$`20%`, outMODIS_PIH_Satagay_geom_na$`50%`, outMODIS_PIH_Satagay_geom_na$`80%`,
        col = colors_PIH_boxplot, outline = FALSE, 
        ylim =c(0,2000), cex.axis = 1.75, 
        ylab = "Measured PIH [m]", cex.lab = 1.7, names = label_names_PIH_S)

plot(h_02_S, col = rgb(0,0,1,1/4), main = NULL, xlab = "Measured PIH [m]", 
     cex.lab = 1.75, cex.axis = 1.75,  xlim = c(0, 2000), ylim = c(0,8000))
plot(h_05_S, col = rgb(0,1,1,1/4), add = T, main = NULL, xlab = "Measured PIH [m]", cex.lab = 1.75)
plot(h_08_S, col = rgb(1,0,1,1/4), add = T, main = NULL, xlab = "Measured PIH [m]", cex.lab = 1.75)

legend("topright", legend = c("Q[0.2]","Q[0.5]", "Q[0.8]"), 
       col= colors_PIH_boxplot, pt.cex=3, pch=15, cex = 2)
dev.off()



## Range between Q[0.2] and Q[0.8] ##############################################
png(glue("Results/Simulation/PIH/histogram_PIH_ranges.png"), width = 1000, height = 800)
plot(h_range_K, col = rgb(1,0,1,1/4), main = NULL, xlab = "Measured PIH [m]", 
     cex.lab = 1.6, cex.axis = 1.5,  xlim = c(0, 2000), ylim = c(0,4000000))
plot(h_range_S, col = rgb(0,1,1,1/4), add = T, main = NULL, cex.lab = 1.6)

legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col=c(rgb(0,0,1,1/4), rgb(0,1,1,1/4)), pt.cex=3, pch=15, cex = 1.5)
dev.off()


## Histogram and boxplot #######################################################

png(glue("Results/Simulation/PIH/histogram_boxplot_PIH_ranges.png"), width = 600, height = 800)
layout(matrix(c(1,2,3,3),nrow = 2, ncol = 2, byrow = TRUE))

boxplot(PIH_range_Khamra[1:3372000], outline = FALSE, col =  rgb(0,0,1,1/4), 
        ylim = c(0,2000), ylab = "Measured PIH [m]",
        cex.lab = 1.75, cex.axis = 1.75)
        
boxplot(PIH_range_Satagay[1:16956000], outline = FALSE, col =  rgb(0,1,1,1/4), 
        ylim = c(0,2000), 
        cex.lab = 1.75, cex.axis = 1.75)      

plot(h_range_K, col = rgb(1,0,1,1/4), main = NULL, xlab = "Measured PIH [m]", 
     cex.lab = 1.75, cex.axis = 1.75,  xlim = c(0, 2000), ylim = c(0,4000000))
plot(h_range_S, col = rgb(0,1,1,1/4), add = T, main = NULL, cex.lab = 1.75)

legend("topright", legend = c("Study area Lake Khamra","Study area Lake Satagay"), 
       col=c(rgb(0,0,1,1/4), rgb(0,1,1,1/4)), pt.cex=3, pch=15, cex = 2)
dev.off()


#####
## Spatial visualization of the PIH points within the study regions ############
#####

PIH_points_Khamra  <- st_as_sf(outMODIS_PIH_Kahmra_geom_na)
save(PIH_points_Khamra, file = "Results/Simulation/PIH/data/PIH_points_Khamra.rda")
PIH_points_Satagay <- st_as_sf(outMODIS_PIH_Satagay_geom_na)
save(PIH_points_Satagay, file = "Results/Simulation/PIH/data/PIH_points_Satagay.rda")


## PIH values | 20 % ###########################################################

png(glue("Results/Simulation/PIH/Khamra/20/Map/study_area_PIH_Khamra_fire_years_20.png"), width = 1000, height = 800)
ggplot(Khamra_buf_extent) +
  geom_sf(data = PIH_points_Khamra, aes(color = `20%`), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000))+
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
  labs(title = "Measured plume injection heights for study area Lake Khamra",
       color = "Measured PIH [m] | Q[0.2]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 18, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75), 
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()


png(glue("Results/Simulation/PIH/Satagay/20/Map/study_area_PIH_Satagay_fire_years_20.png"), width = 1000, height = 800)
ggplot(Satagay_buf_extent) +
  geom_sf(data = PIH_points_Satagay, aes(color = `20%`), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000))+
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
  labs(title = "Measured plume injection heights for study area Lake Satagay",
       color = "Measured PIH [m] | Q[0.2]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 18, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()


## PIH values | 50 % ###########################################################


png(glue("Results/Simulation/PIH/Khamra/50/Map/study_area_PIH_Khamra_fire_years_50.png"), width = 1000, height = 800)
ggplot(Khamra_buf_extent) +
  geom_sf(data = PIH_points_Khamra, aes(color = `50%`), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000))+
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
  labs(title = "Measured plume injection heights for study area Lake Khamra",
       color = "Measured PIH [m] | Q[0.5]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 18, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()

png(glue("Results/Simulation/PIH/Satagay/50/Map/study_area_PIH_Satagay_fire_years_50.png"), width = 1000, height = 800)
ggplot(Satagay_buf_extent) +
  geom_sf(data = PIH_points_Satagay, aes(color = `50%`), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000))+
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
  labs(title = "Measured plume injection heights for study area Lake Satagay",
       color = "Measured PIH [m] | Q[0.5]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 18, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()

## PIH values | 80 % ###########################################################

png(glue("Results/Simulation/PIH/Khamra/80/Map/study_area_PIH_Khamra_fire_years_80.png"), width = 1000, height = 800)
ggplot(Khamra_buf_extent) +
  geom_sf(data = PIH_points_Khamra, aes(color = `80%`), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000))+
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
  labs(title = "Measured plume injection heights for study area Lake Khamra",
       color = "Measured PIH [m] | Q[0.8]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 18, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()


png(glue("Results/Simulation/PIH/Satagay/80/Map/study_area_PIH_Satagay_fire_years_80.png"), width = 1000, height = 800)
ggplot(Satagay_buf_extent) +
  geom_sf(data = PIH_points_Satagay, aes(color = `80%`), size = 2, na.rm = T, show.legend = T) +
  scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000))+
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.4) + 
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake Satagay"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake Satagay" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  xlab("") +
  ylab("") +
  labs(title = "Measured plume injection heights for study area Lake Satagay",
       color = "Measured PIH [m] | Q[0.8]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 18, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 12, vjust = 0.75))
dev.off()




