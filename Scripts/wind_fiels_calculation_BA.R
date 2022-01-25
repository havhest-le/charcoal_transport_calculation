#####
# Script to identify the yearly median of wind fiels over Yakutia, Russia #
#####

rm(list = ls(all= TRUE))

# Load packages 
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

#####
# Load and clear data #
#####

# Map # 
russia            <- read_sf("Data/East_Siberia/gadm36_RUS_1.shp") %>% st_transform(4326)
yakutia           <- russia %>% filter(russia$NAME_1 == "Sakha") %>% st_shift_longitude() %>% st_geometry() %>% st_transform(4326)

# Lakes # 
lakes             <- read_sf("Data/Lakes/lakes_sf.gbd/lakes_sf.shp")
Khamra            <- read_sf("Data/Lakes/Khamra/Khamra_polygon.shp") %>% st_transform(4326)
proj_K            <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Khamra))[,1]}
                           +lat_0={st_coordinates(st_centroid(Khamra))[,2]}")
Satagay           <- read_sf("Data/Satagay/satagay_poly_1.shp") %>% st_transform(4326)
proj_S            <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(Satagay))[,1]}
                           +lat_0={st_coordinates(st_centroid(Satagay))[,2]}")

# Lakes Puffer 100 km # 
Khamra_buf_100     <- st_transform(Khamra, crs = CRS(proj_K)) %>%
                      st_buffer(100000) %>% st_transform(4326) 
Satagay_buf_100    <- st_transform(Satagay, crs = CRS(proj_S)) %>%
                      st_buffer(100000) %>% st_transform(4326)


# Lakes Puffer 50 km # 
Khamra_buf_50      <- st_transform(Khamra, crs = CRS(proj_K)) %>%
                      st_buffer(50000) %>% st_transform(4326) 
Satagay_buf_50     <- st_transform(Satagay, crs = CRS(proj_S)) %>%
                      st_buffer(50000) %>% st_transform(4326)

# Create extent for Lakes Puffer # 
Khamra_buf_extent  <- as(extent(st_bbox(Khamra_buf_100 %>% st_transform(4326) 
                                       %>% st_shift_longitude())[c(1,3,2,4)]), 
                                       "SpatialPolygons")
Satagay_buf_extent <- as(extent(st_bbox(Satagay_buf_100 %>% st_transform(4326) 
                                       %>% st_shift_longitude())[c(1,3,2,4)]), 
                                       "SpatialPolygons")
Khamra_buf_50_ext  <- as(extent(st_bbox(Khamra_buf_50 %>% st_transform(4326) 
                                        %>% st_shift_longitude())[c(1,3,2,4)]), 
                                        "SpatialPolygons")
Satagay_buf_50_ext <- as(extent(st_bbox(Satagay_buf_50 %>% st_transform(4326) 
                                        %>% st_shift_longitude())[c(1,3,2,4)]), 
                                        "SpatialPolygons")

# Wind files #
files   <- list.files(path = "Z:/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", 
                                pattern = "*.nc", all.files = T, full.names = T)

fls_Tab <- do.call("rbind", lapply(files, function(x) {
     nf <- nc_open(x)
    tms <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, origin = "1900-01-01")
    out <- data.frame(path = x, date = tms, id = 1:length(tms))
           nc_close(nf)
           out 
}))



#################################################################################

# Script from Simeon with Lake Khamra and 100 km Puffer # 

for(y in 2000:2018) {

  monthList <- lapply(5:8, function(m) {
    
    cat(glue("\ryear {y} month {m}"))
    
    subTab <- subset(fls_Tab, as.numeric(format(date, "%Y")) %in% y &
                                         as.numeric(format(date, "%m")) %in% m)
    
    
    # Creating a list for every level #
    levelList <- lapply(1:7, function(level) {
            u <- raster::crop(brick(as.character(unique(subTab$path)), varname=  "u", level = level), Khamra_buf_extent)
            v <- raster::crop(brick(as.character(unique(subTab$path)), varname = "v", level = level), Khamra_buf_extent)
                 list(u, v)})
    
    rasterList <- lapply(1:nlayers(levelList[[1]][[1]]), function(dts) {
      levTmp <- lapply(1:7, function(level) {
        list(levelList[[level]][[1]][[dts]],
             levelList[[level]][[2]][[dts]])
      })
      uDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[1]])), median, na.rm = T)
      vDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[2]])), median, na.rm = T)
      list(uDate, vDate)
    })
    
    uMonth <- calc(do.call("brick", lapply(rasterList, function(y) y[[1]])), median, na.rm = T)
    vMonth <- calc(do.call("brick", lapply(rasterList, function(y) y[[2]])), median, na.rm = T)
    list(uMonth, vMonth)
  })
  
  yearBrick_Khamra <- list(do.call("brick", lapply(monthList, function(y) y[[1]])),
                      do.call("brick", lapply(monthList, function(y) y[[2]])))
  
  save(yearBrick_Khamra, file = glue("Results/yearBricks/Khamra/{y}_brick_Khamra.RData"))
}


###################################################################################
for(y in 2000:2020){
  
  monthList <- lapply(5:8, function(m) {
    
    cat(glue("\ryear {y} month {m} for lake Satagay"))
    
    subTab <- subset(fls_Tab, as.numeric(format(date, "%Y")) %in% y &
                       as.numeric(format(date, "%m")) %in% m)
    
    # Creating a list for every level #
    levelList <- lapply(1:7, function(level) {
      u <- raster::crop(brick(as.character(unique(subTab$path)), varname=  "u", level = level), Satagay_buf_extent)
      v <- raster::crop(brick(as.character(unique(subTab$path)), varname = "v", level = level), Satagay_buf_extent)
      list(u, v)})
    
    rasterList <- lapply(1:nlayers(levelList[[1]][[1]]), function(dts) {
        levTmp <- lapply(1:7, function(level) {
        list(levelList[[level]][[1]][[dts]],
             levelList[[level]][[2]][[dts]])
      })
      uDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[1]])), median, na.rm = T)
      vDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[2]])), median, na.rm = T)
      list(uDate, vDate)
    })
    
    uMonth <- calc(do.call("brick", lapply(rasterList, function(y) y[[1]])), median, na.rm = T)
    vMonth <- calc(do.call("brick", lapply(rasterList, function(y) y[[2]])), median, na.rm = T)
    list(uMonth, vMonth)
  })
  
  yearBrick_Satagay <- list(do.call("brick", lapply(monthList, function(y) y[[1]])),
                       do.call("brick", lapply(monthList, function(y) y[[2]])))
  
  save(yearBrick_Satagay, file = glue("Results/yearBricks/Satagay/{y}_brick_Satagay.RData"))
}

####################################################################################

#####
# Calculation of wind values # 
#####

files_Khamra     <- list.files(path = "Results/yearBricks/Khamra/", 
                               pattern = "*.RData", all.files = T, full.names = F)
files_Satagay    <- list.files(path = "Results/yearBricks/Satagay/", 
                               pattern = "*.RData", all.files = T, full.names = F)

years_K          <- as.numeric(sapply(strsplit(files_Khamra, "_"), function(x) x[[1]]))
years_S          <- as.numeric(sapply(strsplit(files_Satagay, "_"), function(x) x[[1]]))

fls_list_Khamra  <- lapply(which(years_K%in%c(2000:2020)), function(x) {
  load(glue("Results/yearBricks/Khamra/{years_K[x]}_brick_Khamra.RData"))
  yearBrick_Khamra
})

fls_list_Satagay <- lapply(which(years_K%in%c(2000:2020)), function(x) {
  load(glue("Results/yearBricks/Satagay/{years_K[x]}_brick_Satagay.RData"))
  yearBrick_Satagay
})


newproj         <- "+proj=longlat +datum=WGS84 +no_defs"
extent_area_K   <- extent(c(111.125, 114.875, 59.125, 60.875))
extent_area_S   <- extent(c(116.125, 120.125, 62.125, 63.875))


###################################################################################

wind_calculation <- function(lake_brick, extent, lake_number, years, lake, 
                             years_2, lon_min, lon_max, lat_min, lat_max){

 
  # Creating u and v wind vector bricks 
  uBrick    <- brick(lapply(lake_brick, function(x) x[[1]]))
  vBrick    <- brick(lapply(lake_brick, function(x) x[[2]]))
  
  # Calculation wind speed # 
  spdBrick    <- sqrt((uBrick^2)+(vBrick^2))
  
  # Calculation wind direction #
  dirBrick   <- brick(lapply(1:length(lake_brick), function(x) {
                brick(lapply(1:4, function(y) {
       rb    <- brick(lake_brick[[x]][[1]][[y]], lake_brick[[x]][[2]][[y]])
       dir   <- atan2(rb[[1]], rb[[2]])*(180/pi)
       dir[] <- ifelse(dir[]< 0, 360 + dir[], dir[])
       dir
       }))
  }))
  
  # Creating wind list with median values of direction and speed
  medSpd  <- calc(spdBrick, median)
  medDir  <- calc(dirBrick, median)
  
  # Creating data frame for values and coordinates
  data_map           <- as(medSpd, "SpatialPixelsDataFrame")
  data_spd           <- as.data.frame(data_map)
  colnames(data_spd) <- c("value", "x", "y")
  
  # Creating arrows
  aggrR_dir_spd    <- aggregate(brick(medSpd, medDir), 2) # nrow 7 and ncell 15
  aggrR_dir_spd    <- crop(aggrR_dir_spd, extent)
  r_coord          <- coordinates(aggrR_dir_spd) 
  arrow            <- data.frame(r_coord, geosphere::destPoint(r_coord, 
                                                               aggrR_dir_spd[[2]], 
                                                               aggrR_dir_spd[[1]]*60*60*5)) 

  #####
  # Plotting wind speed and direction #
  #####

  plotDirSpeedMap <- ggplot(data = extent) +  
                     geom_tile(data = data_spd, aes(x = x, y = y, fill = value), alpha = 0.8) +
                     scale_fill_gradientn(colours = rev(viridis::plasma(99)),
                                          breaks  = round(seq(min(medSpd[]), max(medSpd[]), length = 4), 1)) +
                     geom_point(mapping = aes(x = lon[[lake_number]], y = lat[[lake_number]], 
                                          shape = location[[lake_number]]), data = data_coord_lakes, colour = "black",
                                          fill  = "white", size = 3, stroke = 1.5, show.legend = F) +
                     geom_segment(data  = subset(arrow, lon>0), aes(x = x, xend = lon, y = y, yend = lat),
                                  arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
                     scale_shape_manual(name = "Lake", values = c(1))+
                     theme_minimal() +
                     labs(subtitle = "The average wind direction and speed", fill = "Wind speed [m/s]")+
                     xlab("") +
                     ylab("") +
                     xlim(c(lon_min, lon_max)) +
                     ylim(c(lat_min, lat_max))+
                     theme(plot.subtitle = element_text(size = 20, hjust = 0.5, vjust = -3),
                           legend.title  = element_text(size = 12, vjust = 1),
                           legend.text   = element_text(size = 10, vjust = 0.75))
  
  # save as png #
  
  png(glue("Results/windspeed_{years}_{lake}.png"), width = 800, height = 700)
  plot(plotDirSpeedMap)
  memory.limit(size = 9999999999)
  dev.off()
  
  
  #####
  # Calculation of the variance of wind direction and speed #
  #####
  
  var_dir <- calc(x = dirBrick, function(x) {
     crc  <- circular(x, type = "angles", units = "degrees") # create circular with values on circle
     out  <- quantile(crc, probs = c(0.2,0.8))
     if(diff(out)<0) {
     diff(out)+360
     } else diff(out)
     })
  
  # wind variance in data frame #
  map_var                <- as(var_dir, "SpatialPixelsDataFrame")
  map_var_data           <- as.data.frame(map_var)
  colnames(map_var_data) <- c("value", "x", "y")
  
  #####
  # Plotting the wind variance #
  #####
  plotVarMap <- ggplot(extent) +
                geom_tile(data = map_var_data, aes(x = x, y = y, fill = value), alpha = 0.8) +
                scale_fill_gradientn(colours = rev(viridis::viridis(99)),
                                     breaks  = round(seq(min(map_var_data[, 1]), max(map_var_data[, 1]), 
                                     length  = 6), 0)) +
                geom_point(mapping = aes(x = lon[[lake_number]], y = lat[[lake_number]], 
                             shape = location[[lake_number]]), data = data_coord_lakes, colour = "black",
                             fill  = "white", size = 3, stroke = 1.5, show.legend = T) +
                scale_shape_manual(name = "Lake", values = c(1))+
                theme_minimal() +
                labs(subtitle = "The variance of wind direction", fill = "Interquantile\n[Q0.8-Q0.2]") +
                xlab("") +
                ylab("") +
                xlim(c(lon_min, lon_max)) +
                ylim(c(lat_min, lat_max))+
                theme(plot.subtitle = element_text(size = 20, hjust = 0.5, vjust = -3),
                      legend.title  = element_text(size = 12, vjust = 0.5),
                      legend.text   = element_text(size = 9, vjust = 0.75))
              
              
  # save as png #
  png(glue("Results/windvariance_{years}_{lake}.png"), width = 800, height = 700)
  plot(plotVarMap)
  memory.limit(size = 9999999999)
  dev.off()
  
  
  figure_1   <- ggarrange(plotDirSpeedMap, plotVarMap, nrow = 1, ncol = 2, widths = c(1,1),heights = c(1,1),
                          common.legend = F, legend = "bottom")
            
  # save as png #
  png(glue("Results/wind_spd_dir_{years}_{lake}.png"), width = 1200, height = 600)
  plot(figure_1)
  memory.limit(size = 9999999999)
  dev.off() 
  
  
  #####
  # Plotting the windrose #
  #####
  
  extr_dir <- raster::extract(dirBrick, data_coord_lakes[lake_number, 2:3])
  extr_spd <- raster::extract(spdBrick, data_coord_lakes[lake_number, 2:3])
  extr_Tab <- data.frame(loc  = rep(data_coord_lakes$location[lake_number], ncol(extr_dir)),
                         time = 1:length(lake_brick),
                         each = nrow(extr_dir),
                         dir  = c(extr_dir), spd = c(extr_spd)) 


  plts <- lapply(data_coord_lakes$location[lake_number], function(i) {
          c_sub     <- subset(extr_Tab, extr_Tab$loc == i) # bei i wird auf data_coord zugriffen!!! 
          breaks    <- seq(0, 360, 10)
          bins      <- cut(c_sub$dir, breaks)
          data_bins <- c_sub %>%
                       mutate(bins = cut(c_sub$dir, breaks, labels = seq(5, 360, 10)), 
                              bins = as.numeric(as.character(bins))) %>%
                              group_by(bins) %>%
                              summarise(count = n(), spd = median(spd))
          ggplot(data_bins) +
          geom_bar(aes(bins, count , fill = spd), stat = "identity", show.legend = F) +
          scale_fill_gradientn(colours = "grey20", breaks = round(seq(3, 16, length = 4), 0), 
                               limits  = c(0,16)) +
          coord_polar(start = 0) + 
          theme_minimal() +
          ylab("Count") +
          theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = -2)) +
          scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 315, by = 45),
                                 labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) +
          theme(axis.text.x = element_text(size = 14, face = 'bold'))
          })   
      
  
  #####
  # Creating the entire png #
  #####
  
  png(glue("Results/Windplot_{years}_{lake}.png"), width = 1200, height = 1200)
  g <- grid.arrange(figure_1, plts[[1]], nrow = 2,layout_matrix = rbind(c(1,1), c(2,2)))
  print(annotate_figure(g, top = text_grob(glue("Wind calculation of lake {lake} ({years_2})"), 
                                           vjust = 0.8, hjust = 0.5, face = "bold", size = 31)))
  memory.limit(size = 9999999999)
  dev.off()
  

  # THE END #
  
  }


Khamra_2000_2009_100_km   <- wind_calculation(lake_brick  = fls_list_Khamra[1:10] , extent = Khamra_buf_extent, 
                                              lake_number = 1, years = "2000_2009_100_km", lake = "Khamra", years_2 = "2000 to 2009 | 100 km", 
                                              lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)
Khamra_2010_2018_100_km   <- wind_calculation(lake_brick  = fls_list_Khamra[11:19] , extent = Khamra_buf_extent, 
                                              lake_number = 1, years = "2010_2018_100_km", lake = "Khamra",years_2 = "2010 to 2018 | 100 km", 
                                              lon_min     = 111.0, lon_max = 114.8, lat_min = 59.125, lat_max = 60.9)
Satagay_2000_2009_100_km  <- wind_calculation(lake_brick  = fls_list_Satagay[1:10] , extent = Satagay_buf_extent,
                                              lake_number = 2, years = "2000_2009_100_km", lake = "Satagay", years_2 = "2000 to 2009 | 100 km",
                                              lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)
Satagay_2010_2018_100_km  <- wind_calculation(lake_brick  = fls_list_Satagay[11:19] , extent = Satagay_buf_extent, 
                                              lake_number = 2, years = "2010_2018_100_km", lake = "Satagay",years_2 = "2010 to 2018 | 100 km",
                                              lon_min     = 116.0127, lon_max = 120.009, lat_min = 62.0, lat_max = 64)

