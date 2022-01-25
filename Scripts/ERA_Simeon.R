rm(list = ls(all= TRUE))

library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(sf)
library(ncdf4)
library(rnaturalearth)
library(tidyverse)
library(data.table)
library(glue)


data_coord <- data.frame(location = c("Khamra","Illirney","Rauchagytgyn", "Elgygytgyn"),
                         lon = c(112.98, 167.57, 168.71, 172.15), lat = c(59.99, 67.15, 67.80, 67.58))
map     <- read_sf("Data/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
ext     <- extent(c(95, 200, 50, 80))

## I created a simple sf polygon collection with the lakes
## We will create the buffer later
lakes   <- read_sf("Data/lakesSHP/lakes_sf.shp") 

####
# Input data ####
####

# for Linux4 Server
files <- list.files(path = "/Volumes/potsdam-1/data/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", pattern = "*.nc", all.files = T, full.names = T)

# for vivis pc
# files <- list.files(path = "Z:/data/bioing/data/Data_Reanalyse/ERA5/Wind_u_v_pressureLevels", pattern = "*.nc", all.files = T, full.names = T)

fls_Tab <- do.call("rbind", lapply(files, function(x) {
  nf <- nc_open(x)
  tms <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, origin = "1900-01-01")
  out <- data.frame(path = x, date = tms, id = 1:length(tms))
  nc_close(nf)
  out
}))


for(y in 1979:2020) {
  
  for(m in 5:8) {
  
  cat(glue("\ryear {y} month {m}"))
    
  subTab <- subset(fls_Tab, as.numeric(format(date, "%Y")) %in% y &
                     as.numeric(format(date, "%m")) %in% m)
   
  # Creating a list for every level
  levelList <- parallel::mclapply(1:7, function(level) {
      u <- raster::crop(brick(as.character(unique(subTab$path)), varname=  "u", level = level), ext)
      v <- raster::crop(brick(as.character(unique(subTab$path)), varname = "v", level = level), ext)
      list(u, v) }, mc.cores = 7)
    

   rasterList <- parallel::mclapply(1:nlayers(levelList[[1]][[1]]), function(dts) {
      levTmp <- lapply(1:7, function(level) {
        list(levelList[[level]][[1]][[dts]],
             levelList[[level]][[2]][[dts]])
      })
     uDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[1]])), median, na.rm = T)
     vDate <- calc(do.call("brick", lapply(levTmp, function(y) y[[2]])), median, na.rm = T)
     list(uDate, vDate)
   }, mc.cores = 4)
    

  # Creating a brick for wind direction and speed 
  uBrick <- rotate(brick(lapply(rasterList, function(x) x[[1]])))
  vBrick <- rotate(brick(lapply(rasterList, function(x) x[[2]])))
  
  proj   <- glue("+proj=laea +lon_0={mean(ext[1:2])} +lat_0={mean(ext[3:4])}")
  crds   <- data.frame(coordinates(uBrick))
  pts    <- st_as_sf(crds, coords = c("x", "y"), crs = 4326) %>% st_transform(proj)
  lks    <- lakes %>% st_transform(proj)
  lkb    <- lks %>% st_buffer(700000)
  ind    <- apply(st_intersects(pts, lkb, sparse = F), 1, any) & apply(st_intersects(pts, map %>% st_transform(proj), sparse = F), 1, any)
  
  trackPts <- st_coordinates(pts)[ind,]
  endPts   <- crds[ind,]
  
  # plot(map %>% st_transform(proj) %>% st_buffer(0) %>% st_crop(extent(c(apply(trackPts, 2, function(x) c(min(x), max(x)))))), border = NA)
  # plot(mask(projectRaster(uBrick[[1]], crs = CRS(proj)), as(map %>% st_transform(proj), "Spatial")), add = T, legend = F)
  # plot(lakes %>% st_transform(proj), add = T)
  # points(trackPts, pch = 16, cex = 0.2)
  # plot(map %>% st_transform(proj), add = T)
  
  crdsTab <- do.call("rbind", parallel::mclapply(1:nlayers(uBrick), function(z) {

    wBrick <- projectRaster(brick(uBrick[[z]], vBrick[[z]]), crs = CRS(proj))
    
    for(i in 1:12) {
      if(i==1) {
        extrWnd <- raster::extract(wBrick, trackPts)
        windTab <- data.table(tm = i, id = 1:nrow(trackPts), lon = trackPts[,1], lat = trackPts[,2], 
                              u = extrWnd[,1], v = extrWnd[,2])
      } else {
        new <- with(windTab[windTab$tm==(i-1) & !is.na(windTab$v),], 
                    data.table(tm  = i, id = id, lon = lon+u*60*60, lat = lat+v*60*60))
        new[, c('u', 'v') := data.table(raster::extract(wBrick, new[,c("lon", "lat")]))]
        windTab <- rbind(windTab, new)      
      }   
    }
    
    windTab <- windTab[windTab$id%in%which(sapply(unique(windTab$id), function(p) sum(windTab$id==p))>2),]
    
    tracks_sf <- sfheaders::sf_linestring(
      obj = windTab[!is.na(windTab$lon),c("lon", "lat", "id")][order(windTab$id),],
      x = "lon",
      y = "lat",
      linestring_id = "id"
    ) %>% st_set_crs(proj)

    
    ### intersect
    out <- do.call("rbind", lapply(seq(0, 5000, 500), function(d) {
      inters <- st_intersects(lks %>% st_buffer(d), tracks_sf, sparse = F)
      # plot(lks[1,] %>% st_buffer(550000) %>% st_geometry())
      # points(trackPts, cex = 0.2, pch = 16)
      # plot(tracks_sf[inters[1,],], add = T)
      # plot(lks[1,] %>% st_buffer(d), add = T)
      # points(outTmp[,c("X", "Y")], pch = 21)
      do.call("rbind", lapply(1:nrow(inters), function(x) {
        if(sum(inters[x,])>0) {
        outTmp <- project(as.matrix(data.table(st_coordinates(tracks_sf[inters[x,],]))[, head(.SD,1), by = "L1"][,c("X", "Y")]), proj, inv = T)
        if(nrow(outTmp)>0) data.frame(lake = lakes$lake[x], date = subTab$date[z], lon = outTmp[,1], lat = outTmp[,2], buffer = d)
        } else NULL
      }))
    }))
    
    out[!duplicated(glue("{out[,1]}_{out[,2]}_{out[,3]}")),]
    
  }, mc.cores = 4))
  
  if(!file.exists("/Volumes/potsdam-1/data/bioing/user/slisovsk/crdsTab.csv")) {
    write.table(crdsTab, "/Volumes/potsdam-1/data/bioing/user/slisovsk/crdsTab.csv", sep = ",")
  } else write.table(crdsTab, "/Volumes/potsdam-1/data/bioing/user/slisovsk/crdsTab.csv", sep = ",", col.names = F, append = TRUE)
  
  }
  
}

### output plot
# plot(sqrt(uBrick[[1]]^2 + vBrick[[1]]^2))
# points(crdsTab[,c("lon", "lat")], pch = 16, cex = 0.4)
# plot(lakes$geometry, lwd = 13, add = T, col = "red")


## test
# r0 <- raster(extent(as(lakes[1,] %>% st_transform(proj) %>% st_buffer(310000) %>% st_transform(4326), "Spatial")),
#              res = .5, crs = 4326)
# 
# rLake <- rasterize(crdsTab[,c("lon", "lat")], r0, fun = 'count')
# plot(rLake)
# plot(lakes[1,], add = T)