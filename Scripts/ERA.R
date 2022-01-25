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
library(ggnewscale)
library(RCurl)

data_coord <- data.frame(location = c("Khamra","Illirney","Rauchagytgyn", "Elgygytgyn"),
                         lon = c(112.98, 167.57, 168.71, 172.15), lat = c(59.99, 67.15, 67.80, 67.58))
map     <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
ext     <- extent(c(103.82, 180, 50.07, 80.56))
load("Results/lakes.RData")
lakes
load("Results/lakes_roh.RData")
lakes_roh
load("Results/buffer.RData")
buffer


####
# Input data ####
####

files <- list.files(path = "Z:/data/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", pattern = "*.nc", all.files = T, full.names = T)

fls_Tab <- do.call("rbind", lapply(files, function(x) {
  nf <- nc_open(x)
  tms <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, origin = "1900-01-01")
  out <- data.frame(path = x, date = tms, id = 1:length(tms))
  nc_close(nf)
  out
}))


for(y in 1979:1983) {
  
  # y <- 1982
  cat(glue("\rWir befinden uns im Jahre {y} nach Christus. Ganz Gallien ist nicht mehr von den R?mern besetzt."))
  
  subTab <- subset(fls_Tab, as.numeric(format(date, "%Y")) %in% y &
                     as.numeric(format(date, "%m")) %in% c(6)) #### !!!!! only one month for test
   
  # Creating a list for every level
  rasterList <- lapply(unique(subTab$path), function(x) {
    
    # x <- unique(subTab$path)[[1]]
    
    levelList <- lapply(1:7, function(level) {
      # level = 1
      u <- raster::crop(brick(x, varname=  "u", level = level), ext)
      v <- raster::crop(brick(x, varname = "v", level = level), ext)
      
      dir <- atan2(u, v)*(180/pi)
      dir[] <- ifelse(dir[]<0, 360+dir[], dir[])
      spd <- sqrt((u^2)+(v^2))
      
      list(dir, spd)
    })
    

    lapply(1:nlayers(levelList[[1]][[1]]), function(dts) {
      
      levTmp <- lapply(1:7, function(level) {
        list(levelList[[level]][[1]][[dts]],
             levelList[[level]][[2]][[dts]])
      })
     dirDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[1]])), median, na.rm = T)
     spdDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[2]])), median, na.rm = T)
     
     list(dirDate, spdDate)
    })
    
  })
  
  # Creating a brick for wind direction and speed 
  dirBrick <- brick(lapply(rasterList, function(x) brick(lapply(x, function(y) y[[1]]))))
  spdBrick <- brick(lapply(rasterList, function(x) brick(lapply(x, function(y) y[[2]]))))
  
  crdsTab <- do.call("rbind", lapply(1:nlayers(dirBrick), function(z) {
    
    # 1 get layers
    spd <- spdBrick[[z]]
    dir <- dirBrick[[z]]
    
    # 2 get wind arrows
    windLines <- cbind(coordinates(spd), geosphere::destPoint(coordinates(spd), dir[], spd[]*60*60*7))
    head(windLines)
    plot(windLines)
    
    
    # 3 get st_lines
    sf_lines <- st_sfc(lapply(1:nrow(windLines), function(l) st_linestring(matrix(windLines[l,], ncol = 2, byrow = TRUE))), crs = 4326) %>%
                  st_geometry()
    head(sf_lines)
    plot(sf_lines)
  
    
    # (4) make buffer around each line
    
    # 5 merge lakes
    list <- list(buffer, lakes_roh)
  
    # 6 st_intersect    
    # st_intersects doesn't work
    
    for(i in list[[1]]){
      inters <- st_intersection(i, sf_lines)
      plot(i, col = "grey")
       for(j in list[[2]]) {
        plot(j, col = "lightblue")
        plot(i, add = TRUE)
        plot(inters, add = TRUE)
      } }

  
    # 7 data.table(lonOrigin, latOrigin, date, lakee
    # how can I select lon and lat of LINESTRING in "inters"?
        tab <- data.frame(lon = inters[1], lat = inters[2], lake = list[[2]])
    
    
  }))
    
  
  crdsTab
  
  
}

### rasterize output





data_map <- as(medSpd, "SpatialPixelsDataFrame")
data_spd <- as.data.frame(data_map)
colnames(data_spd) <- c("value", "x", "y")

data_brick   <- brick(medSpd, medDir)
brick_coord <- coordinates(data_brick)
brick_coord[,1] <- ifelse(brick_coord[,1]>180,NA, brick_coord[,1])
arrow <- data.frame(brick_coord, geosphere::destPoint(brick_coord, data_brick[[2]], data_brick[[1]]*60*60*7))

ggplot() +  
    geom_tile(data = data_spd, aes(x = x, y = y, fill = value), alpha = 0.8) + 
    geom_sf(data = lakes$Ill$lake, mapping = aes(colours = "white", size = 10), show.legend = F)+
    scale_fill_gradientn(colours = rev(viridis::plasma(99)),
                       breaks = round(seq(min(medSpd[]), max(medSpd[]), length = 5), 0))+
    geom_segment(data = subset(arrow, lon>0), aes(x = x, xend = lon, y = y, yend = lat),
                 arrow = arrow(length = unit(0.1, "cm")), colour = "black")+
    theme_minimal() +
    labs(subtitle = "The average wind direction and speed", fill = "Wind speed\n[m/s]")+
    xlab("") +
    ylab("") +
    theme(plot.subtitle = element_text(size = 20, hjust = 0.5, vjust = -3),
          legend.title = element_text(size = 12, vjust = 1),
          legend.text = element_text(size = 8, vjust = 0.75))


rows <- split(arrow, seq(nrow(arrow)))
lines <- lapply(rows, function(row) { 
  lmat <- matrix(unlist(row[1:4]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)})
lines <- st_sfc(lines)
lines_sf <- st_sf('geometry' = lines)
lines_CRS <- st_set_crs(lines_sf, "+proj=longlat +datum=WGS84") %>% st_geometry()
head(lines_CRS)

# Lake Khamra
Khamra_lake <- st_set_crs(lakes$Khamra$lake, "+proj=longlat +datum=WGS84") %>% st_geometry()
Khamra_buffer <- st_set_crs(lakes$Khamra$buffer, "+proj=longlat +datum=WGS84") %>% st_geometry()
intersect_Kham <- st_intersection(x = Khamra_lake, y = lines_CRS)
jpeg("Results/Khamra_dir.jpeg")
plot(Khamra_lake)
plot(intersect_Kham, add = T)
dev.off()

# Lake Elgy
Ely_lake <- Ely_lake <- st_set_crs(lakes$Elgy$lake, "+proj=longlat +datum=WGS84") %>% st_geometry()
intersect_Ely <- st_intersection(x = Ely_lake, y = lines_CRS)
jpeg("Results/Ely_dir.jpeg")
plot(Ely_lake)
plot(intersect_Ely, add = T)
dev.off()

# Lake Ill
Ill_lake <- st_set_crs(lakes$Ill$lake, "+proj=longlat +datum=WGS84") %>% st_geometry()
intersect_Ill <- st_intersection(x = Ill_lake, y = lines_CRS)
jpeg("Results/Ill_dir.jpeg")
plot(Ill_lake)
plot(intersect_Ill, add = T)
dev.off()



# Creating Wingroses
extr_dir <- raster::extract(dirBrick, data_coord[, 2:3])
extr_spd <- raster::extract(spdBrick, data_coord[, 2:3])
extr_Tab <- data.frame(loc = rep(data_coord$location, ncol(extr_dir)), # f?r jedes Datum einmal 1:4 durchlaufen lassen
                       time = rep(subTab$date, each = nrow(extr_dir)), 
                       dir = c(extr_dir), spd = c(extr_spd)) 

head(extr_Tab)
i = "Illirney" 
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

print(plts)

figure_2 <- ggarrange(plts[[6]],plts[[5]],plts[[3]],plts[[1]],plts[[2]],plts[[4]],
                      labels = c("Satagay 2.0","Malaya Chabyda","Rauchagytgyn","Khamra","Illirney","Elgygytgyn"), 
                      vjust = 2.25, font.label = list(size = 14, face = "italic",common.legend = TRUE, legend = "bottom"))


png(glue("C:/Users/vreichel/Documents/GitHub/SiberianWindSOMs/WindPlot{y}.png"), width = 1200, height = 1200)
g <- grid.arrange(figure_1, figure_2, nrow = 3,
                  ncol = 3, layout_matrix = rbind(c(1,1), c(2,2)))
print(annotate_figure(g, top = text_grob(glue("Wind direction and speed in East Siberia {y}"),vjust = 0.8, hjust = 0.5, face = "bold", size = 30)))

dev.off()

} ### end


# ####
# # Calculating the variance of wind direction and speed ####
# ####
# 
# var_dir <- calc(x = dirBrick, function(x) {
#   crc <- circular(x, type = "angles", units = "degrees") # circular erstellen mit Datenpunkten auf Kreis
#   out <- quantile(crc, probs = c(0.2,0.8)) #Quantile festlegen
#   if(diff(out)<0) { #damit nicht Minuswerte
#     diff(out)+360
#   } else diff(out)
# })
# 
# map_var <- as(var_dir, "SpatialPixelsDataFrame")
# map_var_data <- as.data.frame(map_var)
# colnames(map_var_data) <- c("value", "x", "y")
# 
# plotVarMap <- ggplot() +
#   geom_tile(data = map_var_data, aes(x = x, y = y, fill = value), alpha = 0.8)+
#   geom_sf(data = mapCrop)+
#   scale_fill_gradientn(colours = rev(viridis::viridis(99)),
#                        breaks = round(seq(min(map_var_data[, 1]), max(map_var_data[, 1]), length = 6), 0))+
#   geom_point(mapping = aes(x = LONGITUDE,y = LATITUDE, alpha= FRP, size = FRP),
#              data = subFire, colour = "white", show.legend = F) +
#   scale_alpha(range = c(1/100,1/1000))+
#   geom_point(mapping = aes(x = lon, y = lat, shape = location), data = data_coord,
#              size = 4, stroke = 2, show.legend = F)+
#   scale_shape_manual(values = c(0:5))+
#   theme_minimal() +
#   labs(subtitle = "The variance of wind direction", fill = "Interquantile\n[Q0.8-Q0.2]")+
#   xlab("") +
#   ylab("") +
#   ylim(c(50, 80))+
#   theme(plot.subtitle = element_text(size = 20, hjust = 0.5, vjust = -3),
#         legend.title = element_text(size = 12, vjust = 0.5),
#         legend.text = element_text(size = 8, vjust = 0.75))
# 
# 
# 
# 
# 
# figure_1 <- ggarrange(plotDirSpeedMap, plotVarMap, nrow = 1, ncol = 2, widths = c(1,1),heights = c(1,1),
#                       common.legend = F, legend = "bottom")
#         

####
# Plotting windrose ####
####