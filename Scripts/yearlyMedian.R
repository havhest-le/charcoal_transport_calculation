#####
# Script for wind fields # 
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


#
map0     <- read_sf("~/Google Drive/GeoDat/NaturalEarth/10m_cultural/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp") %>%
              st_shift_longitude()

map  <- map0[map0$adm1_code %in% c("RUS-2612", "RUS-2615", "RUS-2321"),] %>% st_geometry()
proj <- glue("+proj=moll +lon_0={median(st_bbox(map)[c(1,3)])}")

mapS    <- map  %>% st_transform(crs = proj)
map_all <- map0 %>% st_transform(crs = proj) %>% st_geometry() %>% st_buffer(0) %>% st_intersection(st_as_sfc(as(extent(st_bbox(mapS)[c(1,3,2,4)] + c(-250000, 250000, -250000, 250000)), "SpatialPolygons")) 
                                                                                                    %>% st_set_crs(proj))

lakes   <- read_sf("Data/lakesSHP/lakes_sf.shp") 
# plot(lakes, add = T)


### files
files <- list.files(path = "/Volumes/potsdam/data/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", pattern = "*.nc", all.files = T, full.names = T)

fls_Tab <- do.call("rbind", lapply(files, function(x) {
  nf <- nc_open(x)
  tms <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, origin = "1900-01-01")
  out <- data.frame(path = x, date = tms, id = 1:length(tms))
  nc_close(nf)
  out
}))

bfr <- sib %>% st_buffer(250000) %>% st_union()
ext <- as(extent(st_bbox(bfr %>% st_transform(4326) %>% st_shift_longitude())[c(1,3,2,4)]), "SpatialPolygons")

for(y in 1980:2020) {
  
  cat(glue("\ryear {y}  "))
  
  monthList <- mclapply(5:8, function(mth) {
    
    subTab <- as.character(unique(subset(fls_Tab, as.numeric(format(date, "%Y")) %in% y &
                                                  as.numeric(format(date, "%m")) %in% m)$path))
    
    # Creating a list for every level
    levelList <- mclapply(1:7, function(level) {
      u <- raster::crop(brick(subTab, varname=  "u", level = level), ext)
      v <- raster::crop(brick(subTab, varname=  "v", level = level), ext)
      list(u, v)}, mc.cores = 3)
    
    
    rasterList <- mclapply(1:nlayers(levelList[[1]][[1]]), function(dts) {
      levTmp <- lapply(1:7, function(level) {
        list(levelList[[level]][[1]][[dts]],
             levelList[[level]][[2]][[dts]])
      })
      
      uDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[1]])), median, na.rm = T)
      vDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[2]])), median, na.rm = T)
      list(uDate, vDate)
      
    }, mc.cores = 3)
    
    
    uMonth <- calc(do.call("brick", lapply(rasterList, function(y) y[[1]])), median, na.rm = T)
    vMonth <- calc(do.call("brick", lapply(rasterList, function(y) y[[2]])), median, na.rm = T)
    
    list(uMonth, vMonth)
    
  
  }, mc.silent = 2)

  yearBrick <- list(do.call("brick", lapply(monthList, function(y) y[[1]])),
                    do.call("brick", lapply(monthList, function(y) y[[2]])))

  save(yearBrick, file = glue("Results/yearBricks/{y}_brick.RData"))
}



#### Plots
ext <- sib %>% st_buffer(250000) %>% st_union() %>% st_transform(4326) %>% st_shift_longitude()


## 1980-2003
opar <- par(mar = c(1,1,1,1), mfrow = c(2,1), bty = "n")

fls <- list.files("Results/yearBricks/")
  years <- as.numeric(sapply(strsplit(fls, "_"), function(x) x[[1]]))

allList <- lapply(which(years%in%c(1980:2003)), function(x) {
  load(glue("Results/yearBricks/{years[x]}_brick.RData"))
  yearBrick
})

uRaster <- crop(mask(projectRaster(rotate(calc(brick(lapply(allList, function(x) x[[1]])), median, na.rm = T)), crs = CRS(proj)), as(bfr, "Spatial")), as(bfr, "Spatial"))
vRaster <- crop(mask(projectRaster(rotate(calc(brick(lapply(allList, function(x) x[[2]])), median, na.rm = T)), crs = CRS(proj)), as(bfr, "Spatial")), as(bfr, "Spatial"))
spd  <- sqrt((uRaster^2)+(vRaster^2))

aggrDir <- projectRaster(rotate(brick(lapply(1:length(allList), function(x) {
  brick(lapply(1:4, function(y) {
    rb <- brick(allList[[x]][[1]][[y]], allList[[x]][[2]][[y]])
    dir <- atan2(rb[[1]], rb[[2]])*(180/pi)
    dir[] <- ifelse(dir[]<0, 360+dir[], dir[])
    dir
  }))
}))), crs = CRS(proj))

coord <- st_sample(st_as_sf(as(extent(spd), "SpatialPolygons")) %>% st_set_crs(proj), size = 250, type = "regular")
coord <- coord[!is.na(raster::extract(spd, st_coordinates(coord))),]

extr <- raster::extract(aggrDir, as(coord %>% st_buffer(100000), "Spatial"))
extrSpd <- raster::extract(spd, as(coord %>% st_buffer(100000), "Spatial"), fun = median, na.rm = T)

## PLOT
plot(map_all, col = "grey90", border = "grey90")
plot(mapS, col = "grey80", border = "grey50", add = T)

invisible(lapply(1:length(coord), function(d) {

# plot(coord[d,] %>% st_buffer(100000), add = T, col = "orange")
  
tt   <- circular(c(extr[[d]]), units = "degrees")
med  <- as.numeric(median.circular(tt[!is.na(tt)])) %% 360
qun  <- as.numeric(quantile(tt[!is.na(tt)], probs = c(0.2, 0.8))) %% 360

orig     <- st_coordinates(coord[d,] %>% st_transform(4326)) 
dest_med <- project(geosphere::destPoint(matrix(orig, ncol = 2), med, 24*60*60*extrSpd[d]), proj)
# arrows(st_coordinates(coord[d,])[,1], st_coordinates(coord[d,])[,2], dest_med[,1], dest_med[,2], length = 0.05)

dest_q <- project(geosphere::destPoint(matrix(orig, ncol = 2), qun, 1200000), proj)
# arrows(st_coordinates(coord[d,])[,1], st_coordinates(coord[d,])[,2], dest_q[,1], dest_q[,2], length = 0)

dest_s <- project(geosphere::destPoint(matrix(orig, ncol = 2), med, 1200000), proj)
pie <- coord[d,] %>% st_buffer(85000) %>% st_intersection(
          st_polygon(list(matrix(c(st_coordinates(coord[d,]), dest_q[1,], dest_s, dest_q[2,], st_coordinates(coord[d,])), ncol = 2, byrow = T))) %>% st_sfc(crs = proj))
  

plot(pie, add = T, col = adjustcolor("cornflowerblue", alpha.f = 0.4), border = adjustcolor("cornflowerblue", alpha.f = 0.7))
arrows(st_coordinates(coord[d,])[,1], st_coordinates(coord[d,])[,2], dest_med[,1], dest_med[,2], length = 0.05, lwd = 1.5)

}))



allList <- lapply(which(years%in%c(2000:2018)), function(x) {
  load(glue("Results/yearBricks/{years[x]}_brick.RData"))
  yearBrick
})

uRaster <- crop(mask(projectRaster(rotate(calc(brick(lapply(allList, function(x) x[[1]])), median, na.rm = T)), crs = CRS(proj)), as(bfr, "Spatial")), as(bfr, "Spatial"))
vRaster <- crop(mask(projectRaster(rotate(calc(brick(lapply(allList, function(x) x[[2]])), median, na.rm = T)), crs = CRS(proj)), as(bfr, "Spatial")), as(bfr, "Spatial"))
spd  <- sqrt((uRaster^2)+(vRaster^2))

aggrDir <- projectRaster(rotate(brick(lapply(1:length(allList), function(x) {
  brick(lapply(1:4, function(y) {
    rb <- brick(allList[[x]][[1]][[y]], allList[[x]][[2]][[y]])
    dir <- atan2(rb[[1]], rb[[2]])*(180/pi)
    dir[] <- ifelse(dir[]<0, 360+dir[], dir[])
    dir
  }))
}))), crs = CRS(proj))

coord <- st_sample(st_as_sf(as(extent(spd), "SpatialPolygons")) %>% st_set_crs(proj), size = 250, type = "regular")
coord <- coord[!is.na(raster::extract(spd, st_coordinates(coord))),]

extr <- raster::extract(aggrDir, as(coord %>% st_buffer(100000), "Spatial"))
extrSpd <- raster::extract(spd, as(coord %>% st_buffer(100000), "Spatial"), fun = median, na.rm = T)

## PLOT
plot(map_all, col = "grey90", border = "grey90")
plot(mapS, col = "grey80", border = "grey50", add = T)

invisible(lapply(1:length(coord), function(d) {
  
  # plot(coord[d,] %>% st_buffer(100000), add = T, col = "orange")
  
  tt   <- circular(c(extr[[d]]), units = "degrees")
  med  <- as.numeric(median.circular(tt[!is.na(tt)])) %% 360
  qun  <- as.numeric(quantile(tt[!is.na(tt)], probs = c(0.2, 0.8))) %% 360
  
  orig     <- st_coordinates(coord[d,] %>% st_transform(4326)) 
  dest_med <- project(geosphere::destPoint(matrix(orig, ncol = 2), med, 24*60*60*extrSpd[d]), proj)
  # arrows(st_coordinates(coord[d,])[,1], st_coordinates(coord[d,])[,2], dest_med[,1], dest_med[,2], length = 0.05)
  
  dest_q <- project(geosphere::destPoint(matrix(orig, ncol = 2), qun, 1200000), proj)
  # arrows(st_coordinates(coord[d,])[,1], st_coordinates(coord[d,])[,2], dest_q[,1], dest_q[,2], length = 0)
  
  dest_s <- project(geosphere::destPoint(matrix(orig, ncol = 2), med, 1200000), proj)
  pie <- coord[d,] %>% st_buffer(85000) %>% st_intersection(
    st_polygon(list(matrix(c(st_coordinates(coord[d,]), dest_q[1,], dest_s, dest_q[2,], st_coordinates(coord[d,])), ncol = 2, byrow = T))) %>% st_sfc(crs = proj) %>% st_buffer(0))
  
  
  plot(pie, add = T, col = adjustcolor("cornflowerblue", alpha.f = 0.4), border = adjustcolor("cornflowerblue", alpha.f = 0.7))
  arrows(st_coordinates(coord[d,])[,1], st_coordinates(coord[d,])[,2], dest_med[,1], dest_med[,2], length = 0.05, lwd = 1.5)
  
}))

plot(st_centroid(lakes %>% st_transform(proj)), add = T)

par(opar)
