rm(list = ls(all= TRUE))

library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(sp)
library(sf)
library(ncdf4)
library(rnaturalearth)
library(tidyverse)
library(glue)
library(dplyr)
library(circular)
library(scales)
library(viridis)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(plotrix)
library(glue)
library(cowplot)
library(ggthemes)
library(ggplot2)



map     <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
ext     <- extent(c(103.82, 180, 50.07, 80.56))

####
# Input data ####
####

files <- list.files(path = "N:/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", pattern = "*.nc", all.files = T, full.names = T)

fls_Tab <- do.call("rbind", lapply(files, function(x) {
  nf <- nc_open(x)
  tms <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, origin = "1900-01-01")
  out <- data.frame(path = x, date = tms, id = 1:length(tms))
  nc_close(nf)
  out
}))


# data_coord <- data.frame(location = c("Khamra","Satagay 2.0","Malaya Chabyda","Illirney","Rauchagytgyn", "Elgygytgyn"),
#                          lon = c(112.98, 117.99, 129.41,168.20,168.71, 172.90), lat = c(59.99, 63.08,61.92,67.21, 67.80, 67.30))

data_coord <- data.frame(location = c("Khamra","Illirney", "Elgygytgyn"), lon = c(112.98,168.20, 172.90), lat = c(59.99,67.21, 67.30))

map <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
ext <- as(extent(c(103.82, 180, 50.07, 80.56)), "SpatialPolygons")

mapCrop <- map %>% 
  st_intersection(st_as_sf(as(ext, "SpatialPolygons")) %>% 
                    st_set_crs(4326)) %>% st_geometry 

for(y in 1:2){
  
  # cat(glue("\r{y}"))
 
  if(y==1) years <- 1980:2000 else years <- 2001:2018
  
  subTab <- subset(fls_Tab, as.numeric(format(date, "%Y")) %in% years &
                     as.numeric(format(date, "%m")) %in% c(6:8)) 

  # Creating a list for every date with wind direction and speed
  rasterList <- lapply(unique(subTab$path), function(x) {
    
    levelList <- lapply(1:7, function(level) {
      u <- raster::crop(brick(as.character(x), varname=  "u", level = level), ext)
      v <- raster::crop(brick(as.character(x), varname = "v", level = level), ext)
      list(u, v) })
    
    lapply(1:nlayers(levelList[[1]][[1]]), function(dts) {
      levTmp <- lapply(1:7, function(level) {
        list(levelList[[level]][[1]][[dts]],
             levelList[[level]][[2]][[dts]])})
    
    uDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[1]])), median, na.rm = T)
    vDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[2]])), median, na.rm = T)
    list(uDate, vDate)
    })
  })
      
      
    u <- raster::crop(brick(x, varname=  "u", level = 1), ext)
    v <- raster::crop(brick(x, varname = "v", level = 1), ext)
    
    dir <- atan2(u, v)*(180/pi)
    dir[] <- ifelse(dir[]<0, 360+dir[], dir[])
    spd <- sqrt((u^2)+(v^2))
    
    list <- list(dir, spd)
  })
  
  # Creating a brick for wind direction and speed 
  dirBrick <- brick(lapply(rasterList, function(x) x[[1]]))
  spdBrick <- brick(lapply(rasterList, function(x) x[[2]]))
  medSpd <- calc(spdBrick, median)
  medDir <- calc(dirBrick, median)
  
  
  data_map <- as(medSpd, "SpatialPixelsDataFrame")
  data_spd <- as.data.frame(data_map)
  colnames(data_spd) <- c("value", "x", "y")
  
  aggrR   <- aggregate(brick(medSpd, medDir), 20)
  r_coord <- coordinates(aggrR)
  r_coord[,1] <- ifelse(r_coord[,1]>180,NA, r_coord[,1])
  arrow <- data.frame(r_coord, geosphere::destPoint(r_coord, aggrR[[2]], aggrR[[1]]*60*60*5))
  
  
  plotDirSpeedMap <- ggplot() +  
    geom_tile(data = data_spd, aes(x = x, y = y, fill = value), alpha = 0.8) + 
    geom_sf(data = mapCrop)+
    scale_fill_gradientn(colours = rev(viridis::plasma(99)),
                         breaks = round(seq(min(medSpd[]), max(medSpd[]), length = 5), 0))+
    geom_point(mapping = aes(x = lon, y = lat, shape = location), data = data_coord, colour = "black",
               fill = "white", size = 3, stroke = 1.5)+
    geom_segment(data = subset(arrow, lon>0), aes(x = x, xend = lon, y = y, yend = lat),
                 arrow = arrow(length = unit(0.1, "cm")), colour = "black")+
    scale_shape_manual(name = "Lakes", values = c(0:2))+
    theme_minimal() +
    labs(subtitle = "The average wind direction and speed", fill = "Wind speed\n[m/s]")+
    xlab("") +
    ylab("") +
    ylim(c(50, 80))+
    theme(plot.subtitle = element_text(size = 20, hjust = 0.5, vjust = -3),
          legend.title = element_text(size = 12, vjust = 1),
          legend.text = element_text(size = 10, vjust = 0.75))
  
  
  ####
  # Calculating the variance of wind direction and speed ####
  ####
  
  var_dir <- calc(x = dirBrick, function(x) {
    crc <- circular(x, type = "angles", units = "degrees") # circular erstellen mit Datenpunkten auf Kreis
    out <- quantile(crc, probs = c(0.2,0.8)) #Quantile festlegen
    if(diff(out)<0) { #damit nicht Minuswerte
      diff(out)+360
    } else diff(out)
  })
  
  map_var <- as(var_dir, "SpatialPixelsDataFrame")
  map_var_data <- as.data.frame(map_var)
  colnames(map_var_data) <- c("value", "x", "y")
  
  plotVarMap <- ggplot() +
    geom_tile(data = map_var_data, aes(x = x, y = y, fill = value), alpha = 0.8)+
    geom_sf(data = mapCrop)+
    scale_fill_gradientn(colours = rev(viridis::viridis(99)),
                         breaks = round(seq(min(map_var_data[, 1]), max(map_var_data[, 1]), length = 6), 0))+
    geom_point(mapping = aes(x = lon, y = lat, shape = location), data = data_coord, colour = "black",
               fill = "white", size = 3, stroke = 1.5, show.legend = F)+
    scale_shape_manual(name = "Lakes", values = c(0:2))+
    theme_minimal() +
    labs(subtitle = "The variance of wind direction", fill = "Interquantile\n[Q0.8-Q0.2]")+
    xlab("") +
    ylab("") +
    ylim(c(50, 80))+
    theme(plot.subtitle = element_text(size = 20, hjust = 0.5, vjust = -3),
          legend.title = element_text(size = 12, vjust = 0.5),
          legend.text = element_text(size = 8, vjust = 0.75))
  
  
  
  figure_1 <- ggarrange(plotDirSpeedMap, plotVarMap, nrow = 1, ncol = 2, widths = c(1,1),heights = c(1,1),
                        common.legend = F, legend = "bottom")
  
  
  ####
  # Plotting windrose ####
  ####
  extr_dir <- raster::extract(dirBrick, data_coord[, 2:3])
  extr_spd <- raster::extract(spdBrick, data_coord[, 2:3])
  extr_Tab <- data.frame(loc = rep(data_coord$location, ncol(extr_dir)),
                         time = rep(subTab$date, each = nrow(extr_dir)), 
                         dir = c(extr_dir), spd = c(extr_spd)) 
  
  head(extr_Tab)
  
  plts <- lapply(data_coord$location, function(i) {
    c_sub <- subset(extr_Tab, extr_Tab$loc == i) # bei i wird auf data_coord zugriffen!!! 
    breaks <- seq(0, 360, 10)
    bins <- cut(c_sub$dir, breaks)
    data_bins <- c_sub %>%
      mutate(bins = cut(c_sub$dir, breaks, labels = seq(5, 360, 10)), 
             bins = as.numeric(as.character(bins))) %>%
      group_by(bins) %>%
      summarise(count = n(), spd = median(spd))
    ggplot(data_bins)+
      geom_bar(aes(bins, count , fill = spd), stat = "identity", show.legend = F)+
      scale_fill_gradientn(colours = "grey20", breaks = round(seq(3, 16, length = 4), 0), 
                           limits = c(0,16))+
      coord_polar(start = 0) + theme_minimal()+
      ylab("Count")+
      theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = -2))+
      scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 315, by = 45),
                         labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))+
      theme(axis.text.x = element_text(size = 14, face = 'bold'))
  })
  
  
  figure_2 <- ggarrange(plts[[1]],plts[[2]],plts[[3]],nrow = 1, ncol = 3,
                        labels = c("Khamra","Illirney","Elgygytgyn"), vjust = 10,
                        font.label = list(size = 14, face = "italic",common.legend = TRUE, legend = "bottom"))
  
  
  png(glue("Results/WindPlot{y}.png"), width = 1200, height = 1200)
  g <- grid.arrange(figure_1, figure_2, nrow = 2, layout_matrix = rbind(c(1,1), c(2,2)))
  print(annotate_figure(g, top = text_grob(glue("East Siberia"), vjust = 0.8, hjust = 0.5, face = "bold", size = 30)))
  memory.limit(size = 9999999999)
  dev.off()
  
} ### end


