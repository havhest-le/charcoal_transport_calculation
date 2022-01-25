rm(list = ls(all= TRUE))

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

map0     <- read_sf("~/Google Drive/GeoDat/NaturalEarth/10m_cultural/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp") %>%
  st_shift_longitude()

map  <- map0[map0$adm1_code %in% c("RUS-2612", "RUS-2615", "RUS-2321"),] %>% st_geometry()
proj <- glue("+proj=moll +lon_0={median(st_bbox(map)[c(1,3)])}")

mapS    <- map  %>% st_transform(crs = proj)
map_all <- map0 %>% st_transform(crs = proj) %>% st_geometry() %>% st_buffer(0) %>% st_intersection(st_as_sfc(as(extent(st_bbox(mapS)[c(1,3,2,4)] + c(-250000, 250000, -250000, 250000)), "SpatialPolygons")) 
                                                                                                    %>% st_set_crs(proj))

lakes   <- read_sf("Data/lakesSHP/lakes_sf.shp") 

ext <- as(extent(st_bbox(sib %>% st_buffer(250000) %>% st_union() %>% st_transform(4326) %>% st_shift_longitude())[c(1,3,2,4)]), "SpatialPolygons")


### files
files <- list.files(path = "/Volumes/potsdam/data/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", pattern = "*.nc", all.files = T, full.names = T)

fls_Tab <- do.call("rbind", lapply(files, function(x) {
  nf <- nc_open(x)
  tms <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, origin = "1900-01-01")
  out <- data.frame(path = x, date = tms, id = 1:length(tms))
  nc_close(nf)
  out
}))



mod <- read.table("/Volumes/potsdam/data/bioing/data/Projekte/Fires/East Siberia/R_Data/MODIS_east_siberia.txt", sep = ";", header = T, dec = ",")
mod$datetime <- as.POSIXct(mod$ACQ_DATE, format = '%d.%m.%Y %H:%M:%S')
mod$LONGITUDE[mod$LONGITUDE<0] <- 180 + (180+mod$LONGITUDE[mod$LONGITUDE<0])

mod_sf <- st_as_sf(mod[,c("datetime", "LATITUDE", "LONGITUDE", "FRP")], coords = c("LONGITUDE", "LATITUDE"), crs = 4326) 
ind    <- mod_sf %>% st_intersects(sib %>% st_buffer(250000) %>% st_union() %>% st_transform(4326) %>% st_shift_longitude(), sparse = FALSE)


dates <- unique(mod_sf$datetime[ind])
proj   <- glue("+proj=laea +lon_0={mean(extent(ext)[1:2])} +lat_0={mean(extent(ext)[3:4])}")
lks    <- lakes %>% st_transform(proj)

# plot(mod_sf$geometry[ind,])
# plot(sib %>% st_buffer(250000) %>% st_union() %>% st_transform(4326) %>% st_shift_longitude(), add = T)

for(d in 1:length(dates)) {
  
  cat(glue("\rdate {d} of {length(dates)}"))
  
  era <- which.min(abs(dates[d]-fls_Tab$date))
 
  levelList <- parallel::mclapply(1:7, function(level) {
    u <- raster::crop(brick(as.character(fls_Tab$path[era]), varname=  "u", level = level)[[fls_Tab$id[era]]], ext)
    v <- raster::crop(brick(as.character(fls_Tab$path[era]), varname=  "v", level = level)[[fls_Tab$id[era]]], ext)
    list(u, v) }, mc.cores = 7)
  
  windBrick <- brick(calc(do.call("brick", lapply(levelList, function(y) y[[1]])), median, na.rm = T),
                  calc(do.call("brick", lapply(levelList, function(y) y[[2]])), median, na.rm = T))


  
  if(sum(ind & mod_sf$datetime==dates[d])>1) {
    
      pts    <- st_as_sf(data.frame(st_coordinates(mod_sf)[ind & mod_sf$datetime==dates[d],]), coords = c("X", "Y"), crs = 4326) %>% st_transform(proj) %>% st_coordinates()
      wBrick <- projectRaster(windBrick, crs = proj)
      
      for(i in 1:12) {
          if(i==1) {
            extrWnd <- raster::extract(wBrick, pts)
            windTab <- data.table(tm = i, id = 1:nrow(pts), lon = pts[,1], lat = pts[,2], 
                                  u = extrWnd[,1], v = extrWnd[,2])
          } else {
            new <- with(windTab[windTab$tm==(i-1) & !is.na(windTab$v),], 
                        data.table(tm  = i, id = id, lon = lon+u*60*60, lat = lat+v*60*60))
            new[, c('u', 'v') := data.table(raster::extract(wBrick, new[,c("lon", "lat")]))]
            windTab <- rbind(windTab, new)      
          }}
        
      windTab <- subset(windTab, !is.na(id))
      windTab <- windTab[windTab$id%in%which(sapply(unique(windTab$id), function(p) sum(windTab$id==p))>2) & !is.na(windTab$id),]
        
      tracks_sf <- sfheaders::sf_linestring(
          obj = windTab[!is.na(windTab$lon),c("lon", "lat", "id")][order(windTab$id),],
          x = "lon",
          y = "lat",
          linestring_id = "id"
      ) %>% st_set_crs(proj)
    
      ### intersect
      out <- data.frame(cbind(project(as.matrix(data.table(st_coordinates(tracks_sf))[, head(.SD,1), by = "L1"][,c("X", "Y")]), proj, inv = T), 
                              dates[d],
                              do.call("cbind", lapply(seq(0, 5000, 1000), function(d) {
                              t(st_intersects(lks %>% st_buffer(d), tracks_sf, sparse = F))
                              }))))
      names(out) <- c("lon", "lat", "date", paste0(lakes$lake, "_", rep(seq(0, 5000, 1000), each = 3)))
      
      if(!file.exists("Results/crdsTab.csv")) {
        write.table(out, "Results/crdsTab.csv", sep = ",", row.names = F)
      } else write.table(out, "Results/crdsTab.csv", sep = ",", row.names = F, col.names = F, append = TRUE)
  }
}




trackMap <- read.csv("Results/crdsTab.csv", header = TRUE)
  trackMap$date <- as.POSIXct(trackMap$date, origin = "1970-01-01", tz = "GMT")
  
sf_tracks <- st_as_sf(trackMap, coords = c("lon", "lat"), crs = 4326)

# pdf("lakesProbs_ill.pdf", height = 10, width = 10)
# opar <- par(mfrow = c(1,1), mar = c(1,1,1,1))

### Illarney
library(oce)
library(ceramic)



projIll   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(lakes[3,]))[,1]} +lat_0={st_coordinates(st_centroid(lakes[3,]))[,2]}")
ill       <- lakes[3,] %>% st_transform(CRS(projIll))
illBuffer <- ill %>% st_buffer(150*1000) %>% st_geometry()

sub <- raster::extent(as(illBuffer, "Spatial"))
im  <- cc_location(as(illBuffer %>% st_buffer(1500), "Spatial"), type =  "mapbox.satellite", zoom = 7)

illMap <- map %>% st_transform(CRS(projIll)) %>% st_buffer(0) %>% st_intersection(illBuffer)
illDat <- sf_tracks[,c(1, seq(4, 19, 3))] %>% st_transform(CRS(projIll)) %>% st_intersection(illBuffer)

plotRGB(mask(projectRaster(im, crs = CRS("+proj=laea +lon_0=168.32417546563 +lat_0=67.3583097714494 +ellps=WGS84")), as(illBuffer, "Spatial")), alpha = 160)
plot(illBuffer, add = T, lwd = 2, border = "grey40")
plot(ill, add = T, col = "blue", border = "grey30")
points(st_coordinates(illDat), pch = 16, cex = 1.6, col = adjustcolor("grey20", alpha.f = 0.5))
points(st_coordinates(illDat)[apply(as.data.frame(illDat)[,-c(1, ncol(illDat))], 1, sum)>0,], pch = 21, cex = 1.6, col = "white", bg = adjustcolor("orange", alpha.f = 0.68))



projIll   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(lakes[1,]))[,1]} +lat_0={st_coordinates(st_centroid(lakes[1,]))[,2]}")
ill       <- lakes[1,] %>% st_transform(CRS(projIll))
illBuffer <- ill %>% st_buffer(150*1000) %>% st_geometry()

sub <- raster::extent(as(illBuffer, "Spatial"))
im  <- cc_location(as(illBuffer %>% st_buffer(1500), "Spatial"), type =  "mapbox.satellite", zoom = 7)

illDat <- sf_tracks[,c(1, seq(4, 19, 3))] %>% st_transform(CRS(projIll)) %>% st_intersection(illBuffer)

plotRGB(mask(projectRaster(im, crs = CRS(projIll), method = "ngb"), as(illBuffer, "Spatial")), alpha = 160)
plot(illBuffer, add = T, lwd = 2, border = "grey40")
plot(ill, add = T, col = "blue", border = "grey30")
points(st_coordinates(illDat), pch = 16, cex = 1.6, col = adjustcolor("grey20", alpha.f = 0.5))
points(st_coordinates(illDat)[apply(as.data.frame(illDat)[,-c(1, ncol(illDat))], 1, sum)>0,], pch = 21, cex = 1.6, col = "white", bg = adjustcolor("orange", alpha.f = 0.68))





projIll   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(lakes[2,]))[,1]} +lat_0={st_coordinates(st_centroid(lakes[2,]))[,2]}")
ill       <- lakes[2,] %>% st_transform(CRS(projIll))
illBuffer <- ill %>% st_buffer(150*1000) %>% st_geometry()

sub <- raster::extent(as(illBuffer, "Spatial"))
im  <- cc_location(as(illBuffer %>% st_buffer(1500), "Spatial"), type =  "mapbox.satellite", zoom = 7)

illDat <- sf_tracks[,c(1, seq(4, 19, 3))] %>% st_transform(CRS(projIll)) %>% st_intersection(illBuffer)

plotRGB(mask(projectRaster(im, crs = CRS(projIll), method = "ngb"), as(illBuffer, "Spatial")), alpha = 210)
plot(illBuffer, add = T, lwd = 2, border = "grey40")
plot(ill, add = T, col = "blue", border = "grey30")

points(st_coordinates(illDat), pch = 16, cex = 1.6, col = adjustcolor("grey20", alpha.f = 0.5))
points(st_coordinates(illDat)[apply(as.data.frame(illDat)[,-c(1, ncol(illDat))], 1, sum)>0,], pch = 21, cex = 1.6, col = "white", bg = adjustcolor("orange", alpha.f = 0.68))












rast1 <- rasterize(st_coordinates(illDat)[apply(as.data.frame(illDat)[,-c(1, ncol(illDat))], 1, sum)>0,], r0, fun = "count")
rast2 <- rasterize(st_coordinates(illDat), r0, fun = "count")

cexF  <- approxfun(range(rast2[], na.rm = T), c(0.5, 6))


plot(illMap)
plot(ill, add = T, col = "darkblue", border = "darkblue", lwd = 2)
plot(illBuffer, add =T)
# points(coordinates(rast1), pch = 21, cex = cexF(rast2[]), bg = adjustcolor("firebrick", alpha.f = 0.5), col = adjustcolor("firebrick", alpha.f = 0.8))

firePie <- cbind(coordinates(rast1), ifelse(is.na(rast1[]), 0, rast1[]), ifelse(is.na(rast2[]), 0, rast2[]))
points(coordinates(rast1), pch = 16, cex = cexF(rast2[]), col = adjustcolor(ifelse(firePie[,3]>0, "firebrick", "grey30"), alpha.f = 0.8))
par(opar)
dev.off()

pdf("lakesProbs_elgy.pdf", height = 10, width = 10)
opar <- par(mfrow = c(1,1), mar = c(1,1,1,1))
### Elgy
projIll   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(lakes[1,]))[,1]} +lat_0={st_coordinates(st_centroid(lakes[1,]))[,2]}")
ill       <- lakes[1,] %>% st_transform(CRS(projIll))
illBuffer <- ill %>% st_buffer(500*1000) %>% st_geometry()

illMap <- map %>% st_transform(CRS(projIll)) %>% st_buffer(0) %>% st_intersection(illBuffer)
illDat <- sf_tracks[,c(1, seq(2, 17, 3))] %>% st_transform(CRS(projIll)) %>% st_intersection(illBuffer)


r0 <- raster(extent(as(illBuffer, "Spatial")), crs = CRS(projIll), res = 30000)
# plot(rasterToPolygons(r0), add = T)

rast1 <- rasterize(st_coordinates(illDat)[apply(as.data.frame(illDat)[,-c(1, ncol(illDat))], 1, sum)>0,], r0, fun = "count")
rast2 <- rasterize(st_coordinates(illDat), r0, fun = "count")

cexF  <- approxfun(range(rast2[], na.rm = T), c(0.5, 6))


plot(illMap)
plot(ill, add = T, col = "darkblue", border = "darkblue", lwd = 2)
plot(illBuffer, add =T)
# points(coordinates(rast1), pch = 21, cex = cexF(rast2[]), bg = adjustcolor("firebrick", alpha.f = 0.5), col = adjustcolor("firebrick", alpha.f = 0.8))

firePie <- cbind(coordinates(rast1), ifelse(is.na(rast1[]), 0, rast1[]), ifelse(is.na(rast2[]), 0, rast2[]))
points(coordinates(rast1), pch = 16, cex = cexF(rast2[]), col = adjustcolor(ifelse(firePie[,3]>0, "firebrick", "grey30"), alpha.f = 0.8))
par(opar)
dev.off()


pdf("lakesProbs_kham.pdf", height = 10, width = 10)
opar <- par(mfrow = c(1,1), mar = c(1,1,1,1))
### Khamra
projIll   <- glue("+proj=laea +lon_0={st_coordinates(st_centroid(lakes[2,]))[,1]} +lat_0={st_coordinates(st_centroid(lakes[2,]))[,2]}")
ill       <- lakes[2,] %>% st_transform(CRS(projIll))
illBuffer <- ill %>% st_buffer(500*1000) %>% st_geometry()

illMap <- map %>% st_transform(CRS(projIll)) %>% st_buffer(0) %>% st_intersection(illBuffer)
illDat <- sf_tracks[,c(1, seq(3, 18, 3))] %>% st_transform(CRS(projIll)) %>% st_intersection(illBuffer)


r0 <- raster(extent(as(illBuffer, "Spatial")), crs = CRS(projIll), res = 30000)
# plot(rasterToPolygons(r0), add = T)

rast1 <- rasterize(st_coordinates(illDat)[apply(as.data.frame(illDat)[,-c(1, ncol(illDat))], 1, sum)>0,], r0, fun = "count")
rast2 <- rasterize(st_coordinates(illDat), r0, fun = "count")

cexF  <- approxfun(range(rast2[], na.rm = T), c(0.5, 6))


plot(illMap)
plot(ill, add = T, col = "darkblue", border = "darkblue", lwd = 2)
plot(illBuffer, add =T)
# points(coordinates(rast1), pch = 21, cex = cexF(rast2[]), bg = adjustcolor("firebrick", alpha.f = 0.5), col = adjustcolor("firebrick", alpha.f = 0.8))

firePie <- cbind(coordinates(rast1), ifelse(is.na(rast1[]), 0, rast1[]), ifelse(is.na(rast2[]), 0, rast2[]))
points(coordinates(rast1), pch = 16, cex = cexF(rast2[]), col = adjustcolor(ifelse(firePie[,3]>0, "firebrick", "grey30"), alpha.f = 0.8))
par(opar)
dev.off()
