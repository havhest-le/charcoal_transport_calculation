#####
## Script for wind field calculation over Yakutia, Russia ######################
## Based on the ERA5 reanalysis wind data, the modern wind patterns were calculated
#####

rm(list = ls(all= TRUE))

#####
## Load packages ###############################################################
#####

library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(sf)
library(parallel)
library(ncdf4)
library(rnaturalearth)
library(tidyverse)
library(data.table)
library(glue)
library(circular)
library(ggpubr)
library(gridExtra)
library(ggsn)
library(ggnewscale)
library(tidyr)


#####
## Load and clear data #########################################################
#####

## Maps
russia   <- read_sf("Data/East_Siberia/gadm36_RUS_1.shp") %>% st_transform(4326)
yakutia  <- russia %>% filter(russia$NAME_1 == "Sakha") %>% 
            st_shift_longitude() %>% st_geometry() %>% st_transform(4326)

## Lakes 
lakes    <- read_sf("Data/Lakes/lakes_sf.gbd/lakes_sf.shp")

## Lake Khamra
Khamra   <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp") %>% st_transform(4326)
proj_K   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Khamra))[,1]}
                             +lat_0={st_coordinates(st_centroid(Khamra))[,2]}")

## Lake Satagay
Satagay  <- read_sf("Data/Satagay/satagay_poly_1.shp") %>% st_transform(4326)
proj_S   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Satagay))[,1]}
                             +lat_0={st_coordinates(st_centroid(Satagay))[,2]}")

## Lake coordinates
data_coord_lakes  <- data.frame(location = c("Khamra","Satagay"),
                                lon = c(112.98, 117.998), 
                                lat = c(59.99, 63.078)) 


## Study regions
buffer_lakes_200_UTM <- data.frame(location = c("Khamra","Satagay"),
                                   lon = c(112.98, 117.998), 
                                   lat = c(59.99, 63.078)) %>% 
                        st_as_sf(coords = c("lon", "lat")) %>%
                        st_set_crs(4326) %>% 
                        st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") %>%
                        st_buffer(200000)

buffer_lakes_200 <- data.frame(location = c("Khamra","Satagay"),
                                    lon = c(112.98, 117.998), 
                                    lat = c(59.99, 63.078)) %>% 
                    st_as_sf(coords = c("lon", "lat")) %>%
                    st_set_crs(4326) %>% 
                    st_buffer(200000)


#####
## Definition of study regions #################################################
#####

## 100 km radius around the lake as center point
Khamra_buf_100  <- st_transform(Khamra, crs = CRS(proj_K)) %>% 
                   st_buffer(100000) %>% st_transform(4326) 

Satagay_buf_100 <- st_transform(Satagay, crs = CRS(proj_S)) %>%
                   st_buffer(100000) %>% st_transform(4326)

## Define the extent of study areas
Khamra_buf_extent  <- as(extent(st_bbox(Khamra_buf_100 %>% st_transform(4326) 
                                                       %>% st_shift_longitude())[c(1,3,2,4)]), 
                                                           "SpatialPolygons")

Satagay_buf_extent <- as(extent(st_bbox(Satagay_buf_100 %>% st_transform(4326) 
                                                        %>% st_shift_longitude())[c(1,3,2,4)]), 
                                                            "SpatialPolygons")

##  200 km radius around the lake as center point
buffer_lakes_200 <- data.frame(location = c("Khamra","Satagay"),
                               lon = c(112.98, 117.998), 
                               lat = c(59.99, 63.078)) %>% 
                    st_as_sf(coords = c("lon", "lat")) %>%
                    st_set_crs(4326) %>% 
                    st_buffer(200000)

## Convert the study area extent into a SpatialPolygon | 200 km 
Khamra_study_area  <- as(extent(st_bbox(buffer_lakes_200$geometry[1]  %>% st_transform(4326) 
                                        %>% st_shift_longitude())[c(1,3,2,4)]), 
                                        "SpatialPolygons")

Satagay_study_area <- as(extent(st_bbox(buffer_lakes_200$geometry[2]  %>% st_transform(4326) 
                                        %>% st_shift_longitude())[c(1,3,2,4)]), 
                                        "SpatialPolygons")
 
#####
## ERA 5 data product ##########################################################
#####

files   <- list.files(path = "Z:/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", 
                      pattern = "*.nc", all.files = T, full.names = T)

## Open the data and define the time structure
fls_Tab <- do.call("rbind", lapply(files, function(x) {
     nf <- nc_open(x)
    tms <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, origin = "1900-01-01")
    out <- data.frame(path = x, date = tms, id = 1:length(tms))
           nc_close(nf)
           out 
           }))

#####
## Calculation of the monthly average median of the u and v wind vectors #######
#####

monthly_median_wind_v_u_vector <- function(extent, lake_name) {

for(y in 2000:2018) {                  # Define the time period
                                       # y for year
monthList <- lapply(5:8, function(m) { # m for month (5:8, may to august)

             cat(glue("\ryear {y} month {m} for lake {lake_name}"))
  
   ## Select the period of the study
   subTab <- subset(fls_Tab, as.numeric(format(date, "%Y")) %in% y & 
                             as.numeric(format(date, "%m")) %in% m)
    
   ## Create a list of u and v wind vectors for each individual pressure level for the study area
   levelList <- lapply(1:7, function(level) { # level for pressure level 
           u <- raster::crop(brick(as.character(unique(subTab$path)), varname=  "u", level = level), extent)
           v <- raster::crop(brick(as.character(unique(subTab$path)), varname = "v", level = level), extent)
                list(u, v)}) ## 1:7 levels = 1000, 950, 900, 850, 800, 750 and 700 hPa
    
   ## Create a list of raster bricks of the the daily median values of u and v wind vectors for each individual pressure level 
   rasterList <- lapply(1:nlayers(levelList[[1]][[1]]), function(dts) { # dts for dates
       levTmp <- lapply(1:7, function(level) {                          # level for pressure level 
                 list(levelList[[level]][[1]][[dts]],                   # [[1]] for u wind vector 
                      levelList[[level]][[2]][[dts]])                   # [[2]] for v wind vector
                 })
       
        uDate <- calc(do.call("brick", lapply(X = levTmp, FUN = function(l) l[[1]])), median, na.rm = T) # l for pressure level
        vDate <- calc(do.call("brick", lapply(X = levTmp, FUN = function(l) l[[2]])), median, na.rm = T)
                 list(uDate, vDate)
                 })
  
   ## Create a list of raster bricks of the the monthly median values of u and v wind vectors
   uMonth <- calc(do.call("brick", lapply(X = rasterList, FUN = function(k) k[[1]])), median, na.rm = T) # k for pressure level
   vMonth <- calc(do.call("brick", lapply(X = rasterList, FUN = function(k) k[[2]])), median, na.rm = T)
             list(uMonth, vMonth)
             })
  
   ## Put the raster bricks of each individual year into a list
   yearBrick <- list(do.call("brick", lapply(monthList, function(j) j[[1]])),                            # j for pressure level
                     do.call("brick", lapply(monthList, function(j) j[[2]])))
  
   ## Save the list as RData file
   if(lake_name == "Khamra"){
     save(yearBrick, file = glue("Results/yearBricks/Khamra/{y}_brick_{lake_name}.RData"))
   }else{
     save(yearBrick, file = glue("Results/yearBricks/Satagay/{y}_brick_{lake_name}.RData"))
  
   }
   }
   }

## Run the For Loop for the study region of Lake Khamra and Satagay 
Khamra_year_brick  <- monthly_median_wind_v_u_vector(extent = Khamra_study_area,  lake_name = "Khamra")
Satagay_year_brick <- monthly_median_wind_v_u_vector(extent = Satagay_study_area, lake_name = "Satagay")

#####
## Change the data information #################################################
#####

## Load data again
files_Khamra     <- list.files(path = "Results/yearBricks/Khamra/", 
                               pattern = ".RData", all.files = T, full.names = F)

files_Satagay    <- list.files(path = "Results/yearBricks/Satagay/", 
                               pattern = ".RData", all.files = T, full.names = F)

## Split the single years out of the file name and set as.numeric
years_K          <- as.numeric(sapply(strsplit(files_Khamra, "_"),  function(x) x[[1]]))
years_S          <- as.numeric(sapply(strsplit(files_Satagay, "_"), function(x) x[[1]]))

## Select the period of the study 
fls_list_Khamra  <- lapply(which(years_K %in% c(2000:2018)), function(x) {
                    get(load(glue("Results/yearBricks/Khamra/{years_K[x]}_brick_Khamra.RData")))
                    })

fls_list_Satagay <- lapply(which(years_S %in% c(2000:2018)), function(x) {
                    get(load(glue("Results/yearBricks/Satagay/{years_S[x]}_brick_Satagay.RData")))
                    })

## Save the edited data as RData file
save(fls_list_Khamra,  file = "Results/wind_values/file_lists/list_files_K.RData")
save(fls_list_Satagay, file = "Results/wind_values/file_lists/list_files_S.RData")

fls_list_Khamra  <- get(load("Results/wind_values/file_lists/list_files_K.RData"))
fls_list_Satagay <- get(load("Results/wind_values/file_lists/list_files_S.RData"))




#####
## Calculation of wind variables: wind speed, direction and variance ###########
#####

## For 100 km radius around the lake ###########################################

wind_calculation <- function(lake_brick, extent, lake_number, years, lake, 
                             years_2, lon_min, lon_max, lat_min, lat_max){
  
  # lake_brick  = fls_list_Khamra[c(2,4,5,8,14,15)]
  # extent = Khamra_buf_extent
  # lake_number = 1
  # years = "fire_years"
  # lake = "Khamra"
  # years_2 = "Fire years"
  # lon_min = 111.0
  # lon_max = 114.8
  # lat_min = 59.125
  # lat_max = 60.9

  # lake_brick  = fls_list_Satagay[c(3,14,15,19)]
  # extent = Satagay_buf_extent
  # lake_number = 2
  # years = "fire_years"
  # lake = "Satagay"
  # years_2 = "Fire years"
  # lon_min     = 116.0127
  # lon_max = 120.009
  # lat_min = 62.0
  # lat_max = 64

  ## Create u and v wind vector raster bricks 
  uBrick    <- raster::crop(brick(lapply(lake_brick, function(x) x[[1]])), extent) # x for single year
  vBrick    <- raster::crop(brick(lapply(lake_brick, function(x) x[[2]])), extent)
  
  ## Calculation of wind speed (spd)
  spdBrick   <- sqrt((uBrick^2)+(vBrick^2))
  
  ## Calculation of wind direction (dir)
  dirBrick   <-  raster::crop(brick(lapply(1:length(lake_brick), function(x) {
                 brick(lapply(1:4, function(y) {  # 1:4 for month 5 to 8 
                  
                  ## Create a raster brick with u and v wind vector values for every month and year          
                  u_v   <- brick(lake_brick[[x]][[1]][[y]],       # x for year
                                 lake_brick[[x]][[2]][[y]])       # y for month
                  dir   <- atan2(u_v[[1]], u_v[[2]])*(180/pi)     # [[1]] for u wind vector                         
                  dir[] <- ifelse(dir[]< 0, 360 + dir[], dir[])   # [[2]] for v wind vector                      
                  dir                                                                 
                }))
                })), extent)

  ## Create a list with the monthly median values of wind spd and dir
  medSpd  <- raster::calc(x = spdBrick, fun = median)
  ## Khamra fire years   2.781827 max | 1.263348 min
  ## Satagay fire years  2.248793 max | 1.855766 min
  
  medDir  <- raster::calc(x = dirBrick, fun = median)
  
  ## Create a data frame with the wind spd values and their coordinates
  data_map           <- as(medSpd, "SpatialPixelsDataFrame")
  data_spd           <- as.data.frame(data_map)
  colnames(data_spd) <- c("value", "x", "y")
  
  ## Create wind arrows 
  aggrR_dir_spd <- aggregate(brick(medSpd, medDir), 2) # nrow 7 and ncell 15
  
  ## Crop the raster brick with the study area extent
  aggrR_dir_spd <- crop(aggrR_dir_spd, extent)    
  
  ## Define the coordinates of the wind arrows
  r_coord       <- coordinates(aggrR_dir_spd)
  
  ## Create a data frame with the calculated start and end points of wind arrows 
  ## Each wind arrow is 5 m/h long
  arrow         <- data.frame(r_coord, 
                              geosphere::destPoint(r_coord, 
                                                   aggrR_dir_spd[[2]], 
                                                   aggrR_dir_spd[[1]]*60*60*5))
                                                   ## *60*60*5 means 5 meter per hours
  
  
  #####
  ## Visualization of the monthly average median of wind spd and dir ###########
  #####
  
  plotDirSpeedMap <- ggplot(data = extent) +  
                     geom_tile(data = data_spd, aes(x = x, y = y, fill = value), alpha = 0.8) +
                     scale_fill_gradientn(colours = rev(viridis::plasma(99)),
                                          breaks=c(1, 2, 3, 4), 
                                          labels=c( "1", "2", "3", "4"),
                                          limits=c(0.5,4))+
                     geom_point(mapping = aes(x = lon[[lake_number]], y = lat[[lake_number]], 
                                              shape = location[[lake_number]]), data = data_coord_lakes, 
                                              colour  = "turquoise1",size = 3.5, stroke = 1.5, show.legend = FALSE) +
                     new_scale("colour") +
                     geom_segment(data  = subset(arrow, lon>0), aes(x = x, xend = lon, y = y, 
                                   yend = lat, colour = "5 m/h"),
                                  arrow = arrow(length = unit(0.1, "cm"))) +
                     scale_colour_manual(values = c("5 m/h" = "black"), name = NULL)+
                     theme_minimal() +
                     labs(subtitle = "The average wind direction and speed", 
                              fill = "spd [m/s]")+
                     xlab("") +
                     ylab("") +
                     xlim(c(lon_min, lon_max)) +
                     ylim(c(lat_min, lat_max)) +
                     theme(plot.subtitle = element_text(size = 26, hjust = 0.5, vjust = -3),
                           legend.title  = element_text(size = 18, vjust = 1),
                           legend.text   = element_text(size = 16, vjust = 0.75),
                           axis.text     = element_text(size = 14))
                  
  ## Save the plot as .png file
  png(glue("Results/wind_calculation/100_km/wind_speed/windspeed_{years}_{lake}_100_km.png"), width = 1000, height = 800)
  plot(plotDirSpeedMap)
  memory.limit(size = 9999999999)
  dev.off()
  
  
  #####
  ## Calculation of the variance of wind dir ###################################
  #####
  
  var_dir <- calc(x = dirBrick, function(x) {
     crc  <- circular(x, type = "angles", units = "degrees") # Create a circle 
     out  <- quantile(crc, probs = c(0.2,0.8))               # Define the quantiles
             if(diff(out) < 0) {
             diff(out) + 360
             } else diff(out)
             })
  
  ## Create a data frame of the wind dir variation values and their coordinates 
  map_var                 <- as(var_dir, "SpatialPixelsDataFrame")
  map_var_data            <- as.data.frame(map_var)
  colnames(map_var_data)  <- c("value", "x", "y")
  
  ## Normalization #############################################################
  ## The wind var of wind dir was normalized and thus has a range of values between 0 and 1
  map_var_data$normalized <- (map_var_data$value-min(map_var_data$value))/(max(map_var_data$value)-min(map_var_data$value))
  
  
  #####
  ## Visualization of the calculated var of wind dir ###########################
  #####
  
  plotVarMap <- ggplot(extent) +
                geom_tile(data = map_var_data, aes(x = x, y = y, fill = normalized), alpha = 0.8) +
                scale_fill_gradientn(colours = rev(viridis::viridis(99)),breaks=c(0,0.5,1)) +
                geom_point(mapping = aes(x = lon[[lake_number]], y = lat[[lake_number]], 
                             shape = location[[lake_number]]), data = data_coord_lakes, 
                           colour  = "turquoise1",size = 3.5, stroke = 1.5) +
                scale_color_manual(name = "Lake", values = c(1))+         
                theme_minimal() +
                labs(subtitle = "The variance of wind direction", 
                     fill = "wind variation") +
                labs(shape = "Lake") +
                xlab("") +
                ylab("") +
                xlim(c(lon_min, lon_max)) +
                ylim(c(lat_min, lat_max)) +
                theme(plot.subtitle = element_text(size = 26, hjust = 0.5, vjust = -3),
                      legend.title  = element_text(size = 18, vjust = 1),
                      legend.text   = element_text(size = 16, vjust = 0.75), 
                      axis.text     = element_text(size = 14))
  
  
  ## Save the plot as .png file
  png(glue("Results/wind_calculation/100_km/wind_variance/windvariance_{years}_{lake}_100_km.png"), width = 1000, height = 800)
  plot(plotVarMap)
  memory.limit(size = 9999999999)
  dev.off()
  
  
  figure_1   <- ggarrange(plotDirSpeedMap, plotVarMap, nrow = 1, ncol = 2, widths = c(1,1),heights = c(1,1),
                          common.legend = F, legend = "bottom")
  
  ## Save the plot as .png file
  png(glue("Results/wind_calculation/100_km/wind_spd_dir_{years}_{lake}_100_km.png"), width = 1200, height = 600)
  plot(figure_1)
  memory.limit(size = 9999999999)
  dev.off() 
  
  
  #####
  ## Wind roses ################################################################
  #####
  
  ## Extract the raster bricks of wind spd and dir from the study areas
  extr_dir <- raster::extract(dirBrick, data_coord_lakes[lake_number, 2:3])
  extr_spd <- raster::extract(spdBrick, data_coord_lakes[lake_number, 2:3])
  
  ## Create a data frame for selected study areas, date, wind spd and dir values  
  extr_Tab <- data.frame(loc  = rep(data_coord_lakes$location[lake_number], ncol(extr_dir)),
                         time = 1:length(lake_brick),
                         each = nrow(extr_dir),
                         dir  = c(extr_dir), spd = c(extr_spd)) 
  
  ## Visualization
  plts <- lapply(data_coord_lakes$location[lake_number], function(i) {
          c_sub     <- subset(extr_Tab, extr_Tab$loc == i) # i for study location
          breaks    <- seq(0, 360, 10)                     # create the wind rose structure 
          bins      <- cut(c_sub$dir, breaks)
          data_bins <- c_sub %>%
                       mutate(bins = cut(c_sub$dir, breaks, labels = seq(5, 360, 10)), 
                              bins = as.numeric(as.character(bins))) %>%
                       group_by(bins) %>%
                       summarise(count = n(), spd = median(spd))
          
          
ggplot(data_bins) +
geom_bar(aes(bins, count, fill = spd), stat = "identity", show.legend = F) +
scale_fill_gradientn(colours = "grey20", 
                     breaks = round(seq(3, 16, length = 4), 0), 
                     limits  = c(0,16)) +
coord_polar(start = 0) + 
theme_minimal() +
labs(subtitle = "Monthly median values from month 5 to 8") +
ylab("Count") +
theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = -2)) +
scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 315, by = 45),
                   labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) +
scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4))+
theme(plot.subtitle = element_text(size = 20, hjust = 2.5, vjust = -130),
        axis.text.x = element_text(size = 18, face = 'bold'), 
         axis.title = element_text(size = 20, face = 'bold'),
          axis.text = element_text(size = 18))
})

## Save the plot as .png file
png(glue("Results/wind_calculation/100_km/wind_roses/wind_rose_{years}_{lake}_100_km.png"), width = 1000, height = 800)
plot(plts[[1]])
memory.limit(size = 9999999999)
dev.off()
      
  
  #####
  ## Create an png-file of all figures #########################################
  #####
  
  png(glue("Results/wind_calculation/100_km/Overview/Windplot_{years}_{lake}_100_km.png"), width = 1200, height = 1200)
  g <- grid.arrange(figure_1, plts[[1]], nrow = 2,layout_matrix = rbind(c(1,1), c(2,2)))
  print(annotate_figure(g, top = text_grob(glue("Wind calculation of study area Lake {lake} | {years_2}"), 
                                           vjust = 0.8, hjust = 0.5, face = "bold", size = 30)))
  memory.limit(size = 9999999999)
  dev.off()
  
  
  ## THE END ##
}



#####
## Run the function for each study area and time period of interest #############
#####

Khamra_2000_2009_100_km   <- wind_calculation(lake_brick  = fls_list_Khamra[1:10] , extent = Khamra_buf_extent,
                                              lake_number = 1, years = "2000_2009_100_km", lake = "Khamra", years_2 = "2000 to 2009",
                                              lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2010_2018_100_km   <- wind_calculation(lake_brick  = fls_list_Khamra[11:19] , extent = Khamra_buf_extent,
                                              lake_number = 1, years = "2010_2018_100_km", lake = "Khamra",years_2 = "2010 to 2018",
                                              lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2000_2018_100_km   <- wind_calculation(lake_brick  = fls_list_Khamra[1:19] , extent = Khamra_buf_extent,
                                              lake_number = 1, years = "2000_2018_100_km", lake = "Khamra",years_2 = "2000 to 2018",
                                              lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)


Satagay_2000_2009_100_km  <- wind_calculation(lake_brick  = fls_list_Satagay[1:10] , extent = Satagay_buf_extent,
                                              lake_number = 2, years = "2000_2009_100_km", lake = "Satagay", years_2 = "2000 to 2009",
                                              lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

Satagay_2010_2018_100_km  <- wind_calculation(lake_brick  = fls_list_Satagay[11:19] , extent = Satagay_buf_extent,
                                              lake_number = 2, years = "2010_2018_100_km", lake = "Satagay",years_2 = "2010 to 2018",
                                              lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

Satagay_2000_2018_100_km  <- wind_calculation(lake_brick  = fls_list_Satagay[1:19] , extent = Satagay_buf_extent,
                                              lake_number = 2, years = "2000_2018_100_km", lake = "Satagay",years_2 = "2000 to 2018",
                                              lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

## Fire years for Lake Khamra: 2001, 2003, 2004, 2007, 2013, 2014
Khamra_fire_years   <- wind_calculation(lake_brick  = fls_list_Khamra[c(2,4,5,8,14,15)], extent = Khamra_buf_extent,
                                        lake_number = 1, years = "fire_years", lake = "Khamra", years_2 = "Fire years",
                                        lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2001   <- wind_calculation(lake_brick  = fls_list_Khamra[2], extent = Khamra_buf_extent,
                                  lake_number = 1, years = "2001", lake = "Khamra", years_2 = "2001",
                                  lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2003   <- wind_calculation(lake_brick  = fls_list_Khamra[4], extent = Khamra_buf_extent,
                                  lake_number = 1, years = "2003", lake = "Khamra", years_2 = "2003",
                                  lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2004   <- wind_calculation(lake_brick  = fls_list_Khamra[5], extent = Khamra_buf_extent,
                                  lake_number = 1, years = "2004", lake = "Khamra", years_2 = "2004",
                                  lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2007   <- wind_calculation(lake_brick  = fls_list_Khamra[8], extent = Khamra_buf_extent,
                                  lake_number = 1, years = "2007", lake = "Khamra", years_2 = "2007",
                                  lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2013   <- wind_calculation(lake_brick  = fls_list_Khamra[14], extent = Khamra_buf_extent,
                                  lake_number = 1, years = "2013", lake = "Khamra", years_2 = "2013",
                                  lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)

Khamra_2014   <- wind_calculation(lake_brick  = fls_list_Khamra[15], extent = Khamra_buf_extent,
                                  lake_number = 1, years = "2014", lake = "Khamra", years_2 = "2014",
                                  lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)



## Fire years of lake Satagay: 2002, 2013, 2014, 2018
Satagay_fire_years  <- wind_calculation(lake_brick  = fls_list_Satagay[c(3,14,15,19)] , extent = Satagay_buf_extent,
                                        lake_number = 2, years = "fire_years", lake = "Satagay",years_2 = "Fire years",
                                        lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

Satagay_2002  <- wind_calculation(lake_brick  = fls_list_Satagay[3] , extent = Satagay_buf_extent,
                                  lake_number = 2, years = "2002", lake = "Satagay",years_2 = "2002",
                                  lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

Satagay_2013  <- wind_calculation(lake_brick  = fls_list_Satagay[14] , extent = Satagay_buf_extent,
                                  lake_number = 2, years = "2013", lake = "Satagay",years_2 = "2013",
                                  lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

Satagay_2014  <- wind_calculation(lake_brick  = fls_list_Satagay[15] , extent = Satagay_buf_extent,
                                  lake_number = 2, years = "2014", lake = "Satagay",years_2 = "2014",
                                  lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

Satagay_2018  <- wind_calculation(lake_brick  = fls_list_Satagay[19] , extent = Satagay_buf_extent,
                                  lake_number = 2, years = "2018", lake = "Satagay",years_2 = "2018",
                                  lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)




#####
## Take the daily wind values ##################################################
## Daily values were used for the following charcoal dispersal calculation #####
#####

## Load the u and v wind vector data by ERA 5 ##################################

wind_files <- tibble(path = list.files("Z:/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", full.names = T)) %>%
              mutate(year  = as.numeric(sapply(strsplit(path, "/"), function(x) substring(x[7], 11,14))),
                     month = as.numeric(sapply(strsplit(path, "/"), function(x) substring(x[7], 16,17)))) %>%
              filter(year%in%c(2001:2018) & month%in%c(5:8))


wndTab <- do.call("rbind", lapply(wind_files$path, function(k) {
  
suppressMessages({do.call("rbind", lapply(1:7, function(level) { ## level = pressure level
     ## 1:7 levels = 1000, 950, 900, 850, 800, 750 and 700 hPa
    
     ## u wind components 
     u <- mask(crop(brick(k, varname =  "u", level = level), as(buffer_lakes_200 %>% st_combine(), "Spatial")),
               as(buffer_lakes_200 %>% st_combine(), "Spatial"))
               # A radius of 200 km around the center of the lake was chosen, which was later reduced
    
    uTab <- u[] %>% as_tibble() %>% bind_cols(coordinates(u) %>% as_tibble() %>% rename(lon = x, lat = y)) %>%
            pivot_longer(starts_with("X")) %>% filter(!is.na(value)) %>%
            mutate(date = as.POSIXct(name, format = "X%Y.%m.%d.%H.%M.%S", tz = "GMT")) %>%
            group_by(date = as.Date(date), lon, lat) %>%
            summarise(wnd = median(value)) %>% mutate(type = "u", level = level)
          
    ## v wind components 
    v <- mask(crop(brick(k, varname =  "v", level = level), as(buffer_lakes_200 %>% st_combine(), "Spatial")),
              as(buffer_lakes_200 %>% st_combine(), "Spatial"))
    
    vTab <- v[] %>% as_tibble() %>% bind_cols(coordinates(u) %>% as_tibble() %>% rename(lon = x, lat = y)) %>%
            pivot_longer(starts_with("X")) %>% filter(!is.na(value)) %>%
            mutate(date = as.POSIXct(name, format = "X%Y.%m.%d.%H.%M.%S", tz = "GMT")) %>%
            group_by(date = as.Date(date), lon, lat) %>%
            summarise(wnd = median(value)) %>% mutate(type = "v", level = level)
    
    ## Bind the rows of the u and v wind components together
    bind_rows(uTab, vTab)})) %>% pivot_wider(id_cols = c(date, lon, lat), names_from = c(type, level), values_from = wnd) })
  
}))

## Save the u and v wind vector components as .RData
save(wndTab, file = "Data/ERA5/windTab.RData") 
wndTab <- get(load("Z:/bioing/user/slisovsk/windTab.rda")) 
## The function did not run on my laptop because the RAM is only 8 GB.


#####
## Calculation of the wind variables: wind spd and dir #########################
#####

wndTabSpdDir <- wndTab %>%
                mutate(med_u    = pmap_dbl(list(u_1, u_2, u_3, u_4, u_5, u_6, u_7), ~median(c(...))), ## calculate the median 
                       med_v    = pmap_dbl(list(v_1, v_2, v_3, v_4, v_5, v_6, v_7), ~median(c(...))),
                       ## wind speed
                       wind_spd = sqrt(med_u^2) + sqrt(med_v^2),
                       ## wind direction
                       wind_dir = ifelse(atan2(med_u, med_v)*(180/pi) < 0, 360 +  atan2(med_u, med_v)*(180/pi),atan2(med_u, med_v)*(180/pi)))

## Save the data frame as .RData
save(wndTabSpdDir, file = "Data/ERA5/wndTabSpdDir.rda")
wndTabSpdDir <- get(load("Data/ERA5/wndTabSpdDir.rda"))



