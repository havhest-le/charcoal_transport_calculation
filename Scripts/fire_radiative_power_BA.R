#####
## Script for determination of geographical locations, the burn time and the 
## number of active fire data points as well as their thermal energy output based on the
## Thermal Anomalies/ Fire locations Collection 6.1 of MODIS
#####

rm(list = ls(all= TRUE))

#####
# Load packages ################################################################
#####

library(utils)
library(dplyr)
library(base)
library(tidyverse)
library(sf)
library(ggplot2)
library(glue)
library(graphics)
library(ggpubr)
library(stats)

#####
## Load and clean data #########################################################
#####

## Study area 

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


## Lake coordinates
data_coord_lakes  <- data.frame(location = c("Khamra","Satagay"),
                                lon = c(112.98, 117.998), 
                                lat = c(59.99, 63.078)) 


#####
## Fire radiative power [FRP] ##################################################
#####

MODIS_txt <- read.csv2("Data/MODIS_FRP/MODIS_east_siberia.txt")       %>%
             mutate(TIME =  as.POSIXct(ACQ_DATE, format = '%d.%m.%Y %H:%M:%S'),
                    FRP = as.numeric(FRP), YEAR = as.numeric(format(TIME, "%Y")), 
                    MONTH = as.numeric(format(TIME, "%m")))           %>% 
                    filter(YEAR <= 2018 & MONTH >= 5 & MONTH <= 8)    %>%
                    dplyr::select("LATITUDE","LONGITUDE", "TIME", "FRP") 


#####
## Crop FRP data with defined extent of the study area #########################
#####

## Convcert data type st to sf with spatial reference | box ####################
## Change unit MW/km² to MW
sf_MODIS_100 <- st_as_sf(MODIS_txt, coords = c("LONGITUDE", "LATITUDE")) %>% 
                st_set_crs(4326)                                         %>%
                mutate(lake = apply(st_intersects(., buffer_lakes_100    %>% 
                ## 100 km radius around the lake as center point
                st_transform(4326), sparse = F), 1, function(x) 
                ifelse(any(x), which(x), NA))) %>% filter(!is.na(lake)) 
                ## lake 1 == Khamra

save(sf_MODIS_100, file = "Results/FRP/sf_MODIS_100.RData")
get(load("Results/FRP/sf_MODIS_100.RData"))


## Convcert data type st to sf with spatial reference | circle #################

## Khamra 
sf_MODIS_K_circ      <- st_as_sf(MODIS_txt, coords = c("LONGITUDE", "LATITUDE")) %>% 
                        st_set_crs(4326)                                         %>%
                        mutate(lake = apply(st_intersects(., Khamra_buf_100      %>%
                        st_transform(4326), sparse = F), 1, function(x) 
                        ifelse(any(x), which(x), NA))) %>% filter(!is.na(lake)) 

sf_MODIS_K_circ$Year <- as.numeric(format(sf_MODIS_K_circ$TIME, "%Y"))
save(sf_MODIS_K_circ, file = "Results/FRP/FRP_Khamra_points_circ.RData")
get(load("Results/FRP/FRP_Khamra_points_circ.RData"))

## Satagay
sf_MODIS_S_circ      <- st_as_sf(MODIS_txt, coords = c("LONGITUDE", "LATITUDE")) %>% 
                        st_set_crs(4326)                                         %>%
                        mutate(lake = apply(st_intersects(., Satagay_buf_100     %>%
                        st_transform(4326), sparse = F), 1, function(x) 
                        ifelse(any(x), which(x), NA))) %>% filter(!is.na(lake)) 

sf_MODIS_S_circ$Year <- as.numeric(format(sf_MODIS_S_circ$TIME, "%Y"))
save(sf_MODIS_S_circ, file = "Results/FRP/FRP_Satagay_points_circ.RData")



## Crop the data with the study area extent ####################################
MODIS_for_region <- function(lon_bigger, lon_smaller, lat_bigger, lat_smaller){
                    MODIS_txt %>% 
                    filter (LONGITUDE >= lon_bigger & LONGITUDE <= lon_smaller &  
                            LATITUDE  >= lat_bigger & LATITUDE  <= lat_smaller)} 


## Khamra
Khamra_MODIS_FRP_box   <- MODIS_for_region(lon_bigger  = 111.1664, 
                                           lon_smaller = 114.8065,
                                           lat_bigger  = 59.08353, 
                                           lat_smaller = 60.89779)

## Satagay
Satagay_MODIS_FRP_box  <- MODIS_for_region(lon_bigger  = 116.0127, 
                                           lon_smaller = 120.009,
                                           lat_bigger  = 62.17579, 
                                           lat_smaller = 63.98115)


## Convcert data type st to sf with spatial reference (extent as a box)
sf_MODIS_Khamra_box       <- st_as_sf(Khamra_MODIS_FRP_box, 
                             coords = c("LONGITUDE", "LATITUDE")) %>%
                             st_set_crs(4326)

sf_MODIS_Khamra_box$Year  <- as.numeric(format(sf_MODIS_Khamra_box$TIME, "%Y"))
save(sf_MODIS_Khamra_box,  file = "Results/FRP/FRP_Khamra_points_box.RData")
get(load("Results/FRP/FRP_Khamra_points_box.RData"))

sf_MODIS_Satagay_box      <- st_as_sf(Satagay_MODIS_FRP_box, 
                             coords = c("LONGITUDE", "LATITUDE")) %>%
                             st_set_crs(4326)
sf_MODIS_Satagay_box$Year <- as.numeric(format(sf_MODIS_Satagay_box$TIME, "%Y"))
save(sf_MODIS_Satagay_box, file = "Results/FRP/FRP_Satagay_points_box.RData")
get(load("Results/FRP/FRP_Satagay_points_box.RData"))

#####
## Load saved FRP data again ###################################################
#####
 
# MODIS_data           <- get(load("Results/FRP/FRP_sub.RData")) ## entire MODIS data
# MODIS_data$Year      <- as.numeric(format(MODIS_data$TIME, "%Y"))

## extent as a circle
sf_MODIS_K_circ      <- get(load("Results/FRP/FRP_Khamra_points_circ.RData"))
sf_MODIS_S_circ      <- get(load("Results/FRP/FRP_Satagay_points_circ.RData"))

## extent as a box
sf_MODIS_Khamra_box  <- get(load("Results/FRP/FRP_Khamra_points_box.RData")) %>% 
                        filter(FRP > 0)
length(sf_MODIS_Khamra_box$FRP) ## 6673 fires

sf_MODIS_Satagay_box <- get(load("Results/FRP/FRP_Satagay_points_box.RData")) %>% 
                        filter(FRP > 0)
length(sf_MODIS_Satagay_box$FRP) ## 23927 fires


## Create a list of FRP data of Khamra and Satagay (extent as a box) 
points_lakes_box     <- list(sf_MODIS_Khamra_box, sf_MODIS_Satagay_box)


#####
## Statistical calculation of FRP ##############################################
#####

## Summary 
statstic_FRP_Khamra  <- summary(sf_MODIS_Khamra_box$FRP)
statstic_FRP_Satagay <- summary(sf_MODIS_Satagay_box$FRP)

## Select the median of FRP from the summary and round the median on the first position
statstic_FRP_Khamra[4]  <- round(statstic_FRP_Khamra[4], 1) 
statstic_FRP_Satagay[4] <- round(statstic_FRP_Satagay[4], 1)

## Create a list with the monthly median of FRP 
FRP_boxplot_list        <- list(sf_MODIS_Khamra_box$FRP, sf_MODIS_Satagay_box$FRP)
names(FRP_boxplot_list) <- c(paste("                   Khamra\n                    n = " , 
                                   length(sf_MODIS_Khamra_box$FRP) , sep=""), 
                             
                             paste("                  Satagay\n                     n = " , 
                                   length(sf_MODIS_Satagay_box$FRP) , sep=""))


#####
# Classification of fire years and non-fire years ##############################
#####

## Average median across all years between 2001-2018
med_FRP_Khamra  <- statstic_FRP_Khamra[3]   # 19.4
med_FRP_Satagay <- statstic_FRP_Satagay[3]  # 33.1

## How many days are above the average of FRP? #################################
fire_years_Khamra_days   <- sf_MODIS_Khamra_box %>%
                            filter(FRP >= med_FRP_Khamra) %>%
                            group_by(TIME) %>%
                            summarise("med_FRP" = median(FRP, na.rm = TRUE))

fire_years_Satagay_days  <- sf_MODIS_Satagay_box %>%
                            filter(FRP >= med_FRP_Satagay) %>%
                            group_by(TIME) %>%
                            summarise("med_FRP" = median(FRP, na.rm = TRUE))

length(fire_years_Khamra_days$TIME)  # 159 days
length(fire_years_Satagay_days$TIME) # 196 days


## How many days are under the 1st and above the 3rd quantile? #################
first_qua_FRP_Khamra    <- statstic_FRP_Khamra[2]  # [2] position of 1st quantile
first_qua_FRP_Satagay   <- statstic_FRP_Satagay[2]
third_qua_FRP_Khamra    <- statstic_FRP_Khamra[5]  # [5] position of 3rd quantile
third_qua_FRP_Satagay   <- statstic_FRP_Satagay[5]

## Under the 1st quantile 
under_first_qua_Khamra_days   <- sf_MODIS_Khamra_box %>%
                                 filter(FRP <= first_qua_FRP_Khamra) %>%
                                 group_by(TIME) %>%
                                 summarise("med_FRP" = median(FRP, na.rm = TRUE))

under_first_qua_Satagay_days  <- sf_MODIS_Satagay_box %>%
                                 filter(FRP <= first_qua_FRP_Satagay) %>%
                                 group_by(TIME) %>%
                                 summarise("med_FRP" = median(FRP, na.rm = TRUE))

length(under_first_qua_Khamra_days$TIME)  # 252 days
length(under_first_qua_Satagay_days$TIME) # 313 days

## Above the 3rd quantile 
above_third_qua_Khamra_days   <- sf_MODIS_Khamra_box %>%
                                 filter(FRP >= third_qua_FRP_Khamra) %>%
                                 group_by(TIME) %>%
                                 summarise("med_FRP" = median(FRP, na.rm = TRUE))

above_third_qua_Satagay_days  <- sf_MODIS_Satagay_box %>%
                                 filter(FRP >= third_qua_FRP_Satagay) %>%
                                 group_by(TIME) %>%
                                 summarise("med_FRP" = median(FRP, na.rm = TRUE))

length(above_third_qua_Khamra_days$TIME)  # 113 days
length(above_third_qua_Satagay_days$TIME) # 129 days



## Classification of fire years ################################################
## Calculation of the monthly median of FRP for each individual year
## Years with FRP values greater than or equal to the median of FRP over 
## the entire time period 2001 to 2018 were defined as fire years.

fire_years_Khamra_med   <- sf_MODIS_Khamra_box                              %>%
                           group_by(Year)                                   %>% 
                           summarise("med_FRP" = median(FRP, na.rm = TRUE)) %>%
                           filter(med_FRP >= med_FRP_Khamra)

fire_years_Satagay_med  <- sf_MODIS_Satagay_box                             %>%
                           group_by(Year)                                   %>% 
                           summarise("med_FRP" = median(FRP, na.rm = TRUE)) %>%
                           filter(med_FRP >= med_FRP_Satagay)

fire_years_Khamra_med  
## Fire years of study area Lake Khamra = 2001, 2003, 2004, 2007, 2013, 2014
fire_years_Satagay_med 
## Fire years of study area Lake Satagay = 2002, 2013, 2014, 2018


## Select the fire years 
fire_years_Khamra  <- sf_MODIS_Khamra_box                      %>%
                      group_by(Year)                           %>%
                      filter(Year == "2001" | Year == "2003" |  
                             Year == "2004" | Year == "2007" | 
                             Year == "2013" | Year == "2014" )
save(fire_years_Khamra, file = "Results/Simulation/data/fire_years_FRP_Khamra.rda")
fire_years_Khamra <- get(load("Results/Simulation/data/fire_years_FRP_Khamra.rda"))


fire_years_Satagay <- sf_MODIS_Satagay_box                     %>%
                      group_by(Year)                           %>%
                      filter(Year == "2002" | Year == "2013" | 
                             Year == "2014" | Year == "2018" )
save(fire_years_Satagay, file = "Results/Simulation/data/fire_years_FRP_Satagay.rda")
fire_years_Satagay <- get(load("Results/Simulation/data/fire_years_FRP_Satagay.rda"))


#####
## Visualization of the monthly median of FRP with boxplot() ###################
#####

# number_fires_Khamra <- aggregate(FRP~Year,FUN=length,data=sf_MODIS_Khamra_box)
# number_fire

boxplot_FRP <- function(lakes, fire_years_data, vjust, lake_name, median) {

  # lakes = 1
  # fire_years_data = fire_years_Khamra
  # vjust = -0.2
  # lake_name = "Khamra"
  # median = med_FRP_Khamra

  data            <- points_lakes_box[[lakes]]
  data$Year       <- as.factor(data$Year)
  median_FRP      <- data %>%
                     group_by(Year) %>% dplyr::summarise("med_FRP" = median(FRP, na.rm = TRUE))

  fire_years      <- fire_years_data
  fire_years$Year <- as.factor(fire_years$Year)
  

  cols            <- c("Fire years" = "firebrick3")
  
  cols_2          <- c("Non-fire years" = "white")
  
  boxplot         <-   ggplot(data = data, mapping = aes(x = Year, y = FRP)) +
                       stat_boxplot(geom = 'errorbar') +
                       geom_boxplot(outlier.size = 0.1, colour = "black", fatten = 0.1, show.legend = T) +
                       geom_text(data     = median_FRP, aes(y = med_FRP, label = round(med_FRP, 2)),
                                 position = position_dodge(width = 1),
                                 vjust    = vjust, size = 5, check_overlap = TRUE) + 
    
                       scale_colour_manual(values = cols_2, labels = glue("FRP median < {median} MW")) +
                       labs(colour = "Non-fire years") +
                       geom_boxplot(data = fire_years,,mapping = aes(x = Year, y = FRP, fill = "Fire years"),
                                    alpha= 0.4, fatten = 0.0001,
                                    show.legend = T) +
                       labs(fill = "Fire years", colour = "Non-fire years") +
                       scale_fill_manual(values = cols, labels = glue("FRP median >= {median} MW")) +
                       ylim(c(0,100)) + 
                       theme_minimal() +
                       scale_alpha_manual(values= c(1,0.1)) +
                       xlab("Year") +
                       ylab("FRP [MW]") +
                       coord_cartesian(clip = "off") +
                       theme(legend.title = element_text(size= 22),
                             legend.text = element_text(size= 20), 
                             legend.key.size = unit(3, 'cm'),
                             axis.title = element_text(size = 24),
                             axis.text = element_text(size = 18))
  
  png(glue("Results/FRP/Statistic/boxplot_median_FRP_{lake_name}_new.png"), width = 1300, height = 900)
  plot(boxplot)
  dev.off()
}


## Run the function for both study areas
Khamra_boxplot_FRP  <- boxplot_FRP(lakes = 1, fire_years_data = fire_years_Khamra,
                                   vjust = -0.2, lake_name = "Khamra", 
                                   median = med_FRP_Khamra)
Satagay_boxplot_FRP <- boxplot_FRP(lakes = 2, fire_years_data = fire_years_Satagay, 
                                   vjust = -0.5, lake_name = "Satagay", 
                                   median = med_FRP_Satagay)

#####
## Boxplot over all years ######################################################
#####


## Time period 2001 to 2018 ####################################################

summary(sf_MODIS_Khamra_box$FRP)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 2.10   10.90   19.40   49.61   42.00 2998.60
summary(sf_MODIS_Satagay_box$FRP) 
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 2.50   17.30   33.10   79.10   72.65 5112.50 


## Study area Lake Khamra
png(glue("Results/FRP/Statistic/boxplot_FRP_median_Khamra.png"), width = 800, height = 700)
boxplot(sf_MODIS_Khamra_box$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4))
title(main = "Fire radiative power of study area Lake Khamra | 2001 to 2018", 
      xlab = glue("n  =  {length(sf_MODIS_Khamra_box$FRP)}"), 
      ylab = "FRP [MW]", cex.main = 1.8, cex.lab = 1.5)
legend("topright", 
       legend=c("Min.         2.10","1st Qu.    10.90", "Median    19.40", "Mean       49.61", "3rd Qu.    42.00", "Max.        2998.60"), 
       title="Summary", cex = 1.5)
dev.off()

## Study area Lake Satagay
png(glue("Results/FRP/Statistic/boxplot_FRP_median_Satagay.png"), width = 800, height = 700)
boxplot(sf_MODIS_Satagay_box$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4))
title(main = "Fire radiative power of study area Lake Satagay | 2001 to 2018", 
      xlab = glue("n  =  {length(sf_MODIS_Satagay_box$FRP)}"),
      ylab = "FRP [MW]", cex.main = 1.8, cex.lab = 1.5)
legend("topright", 
       legend=c("Min.         2.50","1st Qu.    17.30", "Median    33.10", "Mean       79.10", "3rd Qu.    72.65", "Max.        5112.50"), 
       title="Summary", cex = 1.5)
dev.off()

## Multiple plot of study area Lake Khamra and Satagay
png(glue("Results/FRP/Statistic/boxplot_FRP_median_parallel.png"), width = 1200, height = 800)
layout(matrix(c(1,2),nrow = 1, ncol = 2, byrow = TRUE))
boxplot(sf_MODIS_Khamra_box$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4), ylim = c(0,150))
title(main = "Study area Lake Khamra", 
      xlab = glue("n  =  {length(sf_MODIS_Khamra_box$FRP)}"), 
      ylab = "FRP [MW]", cex.main = 1.4, cex.lab = 1.5)
boxplot(sf_MODIS_Satagay_box$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4))
title(main = "Study area Lake Satagay", 
      xlab = glue("n  =  {length(sf_MODIS_Satagay_box$FRP)}"),
      cex.main = 1.4, cex.lab = 1.5)
dev.off()


## Fire years ##################################################################

summary(fire_years_Khamra$FRP)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 2.70   12.50   22.60   57.70   50.23 2998.60
summary(fire_years_Satagay$FRP) 
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 2.50   17.70   33.90   78.89   74.55 4435.70 

## Study area Lake Khamra
png(glue("Results/FRP/Statistic/boxplot_FRP_median_Khamra_fire_years.png"), width = 800, height = 700)
boxplot(fire_years_Khamra$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4))
title(main = "Fire radiative power of study area Lake Khamra | Fire years", 
      xlab = glue("n  =  {length(fire_years_Khamra$FRP)}"), 
      ylab = "FRP [MW]", cex.main = 1.8, cex.lab = 1.5)
legend("topright", 
       legend=c("Min.         2.70","1st Qu.    12.50", "Median    22.60", "Mean       57.70", "3rd Qu.    50.23", "Max.        2998.60"), 
       title="Summary", cex = 1.5)
dev.off()

## Study area Lake Satagay
png(glue("Results/FRP/Statistic/boxplot_FRP_median_Satagay_fire_years.png"), width = 800, height = 700)
boxplot(fire_years_Satagay$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4))
title(main = "Fire radiative power of study area Lake Satagay | Fire years", 
      xlab = glue("n  =  {length(fire_years_Satagay$FRP)}"),
      ylab = "FRP [MW]", cex.main = 1.8, cex.lab = 1.5)
legend("topright", 
       legend=c("Min.         2.50","1st Qu.    17.70", "Median    33.90", "Mean       78.89", "3rd Qu.    74.55", "Max.        4435.70"), 
       title="Summary", cex = 1.5)
dev.off()

## Multiple plot of study area Lake Khamra and Satagay
png(glue("Results/FRP/Statistic/boxplot_FRP_median_parallel_fire_years.png"), width = 1200, height = 800)
layout(matrix(c(1,2),nrow = 1, ncol = 2, byrow = TRUE))
boxplot(fire_years_Khamra$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4), ylim = c(0,150))
title(main = "Study area Lake Khamra", 
      xlab = glue("n  =  {length(fire_years_Khamra$FRP)}"), 
      ylab = "FRP [MW]", cex.main = 1.4, cex.lab = 1.5)
boxplot(fire_years_Satagay$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4))
title(main = "Study area Lake Satagay", 
      xlab = glue("n  =  {length(fire_years_Satagay$FRP)}"),
      cex.main = 1.4, cex.lab = 1.5)
mtext("Fire radiative power | Fire years", line = -1.5,outer = TRUE, cex = 1.8, font = 2)
dev.off()


#####
## Histogram of FRP ############################################################
#####

## Multiple plot of study area Lake Khamra and Satagay | Fire years
png(glue("Results/FRP/Statistic/hist_FRP_median_parallel_fire_years.png"), width = 1200, height = 800)
layout(matrix(c(1,2),nrow = 1, ncol = 2, byrow = TRUE))
hist(fire_years_Khamra$FRP, col = rgb(1, 0, 0, alpha = 0.4), xlim = c(0,600), 
     breaks = 200, xlab = "FRP [MW]", main = NULL, cex.lab = 1.5)
title(main = "Study area Lake Khamra", cex.main = 1.4)
hist(fire_years_Satagay$FRP, col = rgb(1, 0, 0, alpha = 0.4), xlim = c(0,600), 
     breaks = 200,  xlab = "FRP [MW]", main = NULL, cex.lab = 1.5)
title(main = "Study area Lake Satagay", cex.main = 1.4)
mtext("Fire radiative power | Fire years", line = -1.5,outer = TRUE, cex = 1.8, font = 2)
dev.off()


png(glue("Results/FRP/Statistic/hist_FRP_median_parallel_fire_years_300.png"), width = 1200, height = 800)
layout(matrix(c(1,2),nrow = 1, ncol = 2, byrow = TRUE))
hist(fire_years_Khamra$FRP, col = rgb(1, 0, 0, alpha = 0.4), xlim = c(0,300), 
     breaks = 200, xlab = "FRP [MW]", main = NULL, cex.lab = 1.5)
title(main = "Study area Lake Khamra", cex.main = 1.4)
hist(fire_years_Satagay$FRP, col = rgb(1, 0, 0, alpha = 0.4), xlim = c(0,300), 
     breaks = 200,  xlab = "FRP [MW]", main = NULL, cex.lab = 1.5)
title(main = "Study area Lake Satagay", cex.main = 1.4)
mtext("Fire radiative power | Fire years", line = -1.5,outer = TRUE, cex = 1.8, font = 2)
dev.off()


#####
## Visualization | One layout ##################################################
#####

png(glue("Results/FRP/Statistic/entire_boxplot.png"), width = 1000, height = 1000)
par(mfrow=c(2,2))

boxplot(sf_MODIS_Khamra_box$FRP, outline = FALSE, col = "lightgrey", ylim = c(0,200))
title(main = "Study area Lake Khamra", 
      xlab = glue("n  =  {length(sf_MODIS_Khamra_box$FRP)}"), 
      ylab = "FRP [MW]", cex.main = 1.6, cex.lab =2)
boxplot(sf_MODIS_Satagay_box$FRP, outline = FALSE, col = "lightgrey", ylim = c(0,200))
title(main = "Study area Lake Satagay", 
      xlab = glue("n  =  {length(sf_MODIS_Satagay_box$FRP)}"),
      cex.main = 1.6, cex.lab =2)
mtext("2001 to 2018", line = -2,outer = TRUE, cex = 2, font = 2)


boxplot(fire_years_Khamra$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4), ylim = c(0,200))
title(xlab = glue("n  =  {length(fire_years_Khamra$FRP)}"), 
      ylab = "FRP [MW]", cex.main = 1.6, cex.lab =2)
boxplot(fire_years_Satagay$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4), ylim = c(0,200))
title(xlab = glue("n  =  {length(fire_years_Satagay$FRP)}"),
      cex.main = 1.6, cex.lab = 2)
mtext("Fire years", line = -43,outer = TRUE, cex = 2, font = 2)
dev.off()


#####
## Comparison between the FRP values of both study areas in form of a ##########
## multiple plot (histogram, boxplot) ##########################################
#####

## 2001 to 2018 ################################################################
hist_FRP_K_2001 <- hist(sf_MODIS_Khamra_box$FRP,  breaks= 200, 
                       xlim = c(0,300), ylim = c(0, 8000), main = NULL,
                       xlab = "FRP [MW]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightgrey")

hist_FRP_S_2001 <- hist(sf_MODIS_Satagay_box$FRP,  breaks= 200, 
                       xlim = c(0,300), ylim = c(0, 8000), main = NULL,
                       xlab = "FRP [MW]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightgrey")


## Fire years ##################################################################
hist_FRP_K_fy   <- hist(fire_years_Khamra$FRP,  breaks= 200, 
                       xlim = c(0,300), ylim = c(0, 8000), main = NULL,
                       xlab = "FRP [MW]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = rgb(1, 0, 0, alpha = 0.4))

hist_FRP_S_fy  <- hist(fire_years_Satagay$FRP,  breaks= 200, 
                       xlim = c(0,300), ylim = c(0, 8000), main = NULL,
                       xlab = "FRP [MW]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = rgb(1, 0, 0, alpha = 0.4))


n_FRP_K    <- length(sf_MODIS_Khamra_box$FRP)
median(sf_MODIS_Khamra_box$FRP) # 19.4
n_FRP_S    <- length(sf_MODIS_Satagay_box$FRP)
median(sf_MODIS_Satagay_box$FRP) # 33.1
n_FRP_K_fy <- length(fire_years_Khamra$FRP)
median(fire_years_Khamra$FRP) # 22.6
n_FRP_S_fy <- length(fire_years_Satagay$FRP)
median(fire_years_Satagay$FRP) # 33.9


FRP_df_K    <- merge(data.frame(sf_MODIS_Khamra_box, row.names=NULL), data.frame(fire_years_Khamra, row.names=NULL),
                   by = 0, all = TRUE)[-1]
FRP_df_S  <- merge(data.frame(sf_MODIS_Satagay_box, row.names=NULL), data.frame(fire_years_Satagay, row.names=NULL),
                   by = 0, all = TRUE)[-1]

label_names_FRP_K    <- c(glue("n = {n_FRP_K}"), glue("n = {n_FRP_K_fy}"))
label_names_FRP_S    <- c(glue("n = {n_FRP_S}"), glue("n = {n_FRP_S_fy}"))

colors_FRP_boxplot <- c("whitesmoke", rgb(1, 0, 0, alpha = 0.4))


################################################################################

png(glue("Results/FRP/Statistic/entire_boxplot_histo_FRP.png"), width = 900, height = 800)
layout(matrix(c(1,2,3,4),nrow = 2, ncol = 2, byrow = TRUE))

boxplot(FRP_df_K$FRP.x, FRP_df_K$FRP.y, ylim = c(0,300), outline = FALSE, 
        names = label_names_FRP_K, col = colors_FRP_boxplot, ylab = "FRP [MW]",
        cex.axis = 1.75, cex.lab = 1.75,  main = "Study area Lake Khamra", cex.main = 2.5)  
boxplot(FRP_df_S$FRP.x, FRP_df_S$FRP.y, ylim = c(0,300), outline = FALSE, 
        names = label_names_FRP_S, col = colors_FRP_boxplot,
        cex.axis = 1.75, cex.lab = 1.75, main = "Study area Lake Satagay", cex.main = 2.5)  

plot(hist_FRP_K_2001, col = "lightgrey", xlim= c(0,300), ylim = c(0, 4000), 
     xlab= "FRP [MW]", main = NULL, cex.lab = 1.75, cex.axis = 1.75)
plot(hist_FRP_K_fy, add = T, col=  rgb(1, 0, 0, alpha = 0.4), xlim= c(0,300), 
     ylim = c(0, 8000),main = NULL,  cex.axis = 1.75)

plot(hist_FRP_S_2001,col = "lightgrey", xlim= c(0,300), ylim = c(0, 8000), main = NULL,
     xlab= "FRP [MW]", ylab= NULL,cex.lab = 1.75,  cex.axis = 1.75)
plot(hist_FRP_S_fy,  add= T, col = rgb(1, 0, 0, alpha = 0.4), xlim= c(0,300), main = NULL,
     ylim = c(0, 8000),cex.axis = 1.75)

legend("topright", legend = c("2001 to 2018","Fire years"), 
       col=c("lightgrey", rgb(1, 0, 0, alpha = 0.4)), pt.cex=3, pch=15, cex = 2)
dev.off()

################################################################################

png(glue("Results/FRP/Statistic/entire_boxplot_histo_FRP_200_400.png"), width = 900, height = 800)
layout(matrix(c(1,2,3,4),nrow = 2, ncol = 2, byrow = TRUE))

boxplot(FRP_df_K$FRP.x, FRP_df_K$FRP.y, ylim = c(0,200), outline = FALSE, 
        names = label_names_FRP_K, col = colors_FRP_boxplot, ylab = "FRP [MW]",
        cex.axis = 1.75, cex.lab = 1.75,  main = "Study area Lake Khamra", cex.main = 2.5)  
boxplot(FRP_df_S$FRP.x, FRP_df_S$FRP.y, ylim = c(0,200), outline = FALSE, 
        names = label_names_FRP_S, col = colors_FRP_boxplot,
        cex.axis = 1.75, cex.lab = 1.75, main = "Study area Lake Satagay", cex.main = 2.5)  

plot(hist_FRP_K_2001, col = "lightgrey", xlim= c(0,400), ylim = c(0, 4000), 
     xlab= "FRP [MW]", main = NULL, cex.lab = 1.75, cex.axis = 1.75)
plot(hist_FRP_K_fy, add = T, col=  rgb(1, 0, 0, alpha = 0.4), xlim= c(0,400), 
     ylim = c(0, 8000),main = NULL,  cex.axis = 1.75)

plot(hist_FRP_S_2001,col = "lightgrey", xlim= c(0,400), ylim = c(0, 8000), main = NULL,
     xlab= "FRP [MW]", ylab= NULL,cex.lab = 1.75,  cex.axis = 1.75)
plot(hist_FRP_S_fy,  add= T, col = rgb(1, 0, 0, alpha = 0.4), xlim= c(0,400), main = NULL,
     ylim = c(0, 8000),cex.axis = 1.75)

legend("topright", legend = c("2001 to 2018","Fire years"), 
       col=c("lightgrey", rgb(1, 0, 0, alpha = 0.4)), pt.cex=3, pch=15, cex = 2)
dev.off()



## Multiple plot for study area Lake Khamra
# png(glue("Results/FRP/Statistic/entire_boxplot_histo_Khamra.png"), width = 1000, height = 900)
# layout(matrix(c(1,2,3,3),nrow = 2, ncol = 2, byrow = TRUE))
# 
# boxplot(sf_MODIS_Khamra_box$FRP, outline = FALSE, col = "lightgrey", ylim = c(0,200))
# title(main = NULL,
#       xlab = glue("n  =  {length(sf_MODIS_Khamra_box$FRP)}"), 
#       ylab = "FRP [MW]", cex.main = 1.6, cex.lab =1.5)
# 
# boxplot(fire_years_Khamra$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4), ylim = c(0,200))
# title(xlab = glue("n  =  {length(fire_years_Khamra$FRP)}"), 
#       cex.main = 1.6, cex.lab =1.5)
# 
# plot(hist_FRP_K_2001, col = "lightgrey", xlim= c(0,300), ylim = c(0, 8000), xlab= "FRP [MW]", main = NULL, cex.main = 2,cex.lab = 1.5)
# plot(hist_FRP_K_fy, add = T, col=  rgb(1, 0, 0, alpha = 0.4), xlim= c(0,300), ylim = c(0, 8000),main = NULL)
# 
# legend("topright", legend = c("2001 to 2018","Fire years"), 
#        col=c("lightgrey", rgb(1, 0, 0, alpha = 0.4)), pt.cex=3, pch=15, cex = 1.8)
# 
# mtext("Fire radiative power for study area Lake Khamra",
#       line = -2.2 ,outer = TRUE, cex = 1.6, font = 2, adj = 0.6)
# 
# dev.off()

## Study area Lake Khamra | 4000
# png(glue("Results/FRP/Statistic/entire_boxplot_histo_Khamra_4000.png"), width = 1000, height = 900)
# layout(matrix(c(1,2,3,3),nrow = 2, ncol = 2, byrow = TRUE))
# 
# boxplot(sf_MODIS_Khamra_box$FRP, outline = FALSE, col = "lightgrey", ylim = c(0,200))
# title(main = NULL,
#       xlab = glue("n  =  {length(sf_MODIS_Khamra_box$FRP)}"), 
#       ylab = "FRP [MW]", cex.main = 1.6, cex.lab =1.5)
# 
# boxplot(fire_years_Khamra$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4), ylim = c(0,200))
# title(xlab = glue("n  =  {length(fire_years_Khamra$FRP)}"), 
#       cex.main = 1.6, cex.lab =1.5)
# 
# plot(hist_FRP_K_2001, col = "lightgrey", xlim= c(0,300), ylim = c(0, 4000), xlab= "FRP [MW]", main = NULL, cex.main = 2,cex.lab = 1.5)
# plot(hist_FRP_K_fy, add = T, col=  rgb(1, 0, 0, alpha = 0.4), xlim= c(0,300), ylim = c(0, 4000),main = NULL)
# 
# legend("topright", legend = c("2001 to 2018","Fire years"), 
#        col=c("lightgrey", rgb(1, 0, 0, alpha = 0.4)), pt.cex=3, pch=15, cex = 1.8)
# 
# mtext("Fire radiative power for study area Lake Khamra",
#       line = -2.2 ,outer = TRUE, cex = 1.6, font = 2, adj = 0.6)
# 
# dev.off()

## Multiple plot for study area Lake Satagay
# png(glue("Results/FRP/Statistic/entire_boxplot_histo_Satagay.png"), width = 1000, height = 900)
# layout(matrix(c(1,2,3,3),nrow = 2, ncol = 2, byrow = TRUE))
# 
# boxplot(sf_MODIS_Satagay_box$FRP, outline = FALSE, col = "lightgrey", ylim = c(0,200))
# title(main = NULL,
#       xlab = glue("n  =  {length(sf_MODIS_Satagay_box$FRP)}"),
#       cex.main = 1.6, cex.lab =1.5)
# 
# boxplot(fire_years_Satagay$FRP, outline = FALSE, col = rgb(1, 0, 0, alpha = 0.4), ylim = c(0,200))
# title(xlab = glue("n  =  {length(fire_years_Satagay$FRP)}"),
#       cex.main = 1.6, cex.lab = 1.5)
# 
# plot(hist_FRP_S_2001,col = "lightgrey", xlim= c(0,300), ylim = c(0, 8000),main = NULL,xlab= "FRP [MW]", ylab= NULL,  cex.main = 2,cex.lab = 1.5)
# plot(hist_FRP_S_fy,  add= T, col = rgb(1, 0, 0, alpha = 0.4), xlim= c(0,300), ylim = c(0, 8000),main = NULL)
# 
# legend("topright", legend = c("2001 to 2018","Fire years"), 
#        col=c("lightgrey", rgb(1, 0, 0, alpha = 0.4)), pt.cex=3, pch=15, cex = 1.8)
# 
# mtext("Fire radiative power for study area Lake Satagay",
#       line = -2.2 ,outer = TRUE, cex = 1.6, font = 2, adj = 0.6)
# 
# dev.off()




#####
## Spatial visualization of the FRP points within the study regions ############
#####

##Points with FRP values less than or equal to 300 MW were plotted 
## FRP limited by the value range from 1-300 MW according to Rogers et al. (2015)
MODIS_Khamra_FRP_300_2001_2018  <- sf_MODIS_Khamra_box  %>% filter(FRP > 0 & FRP <= 300)
MODIS_Satagay_FRP_300_2001_2018 <- sf_MODIS_Satagay_box %>% filter(FRP > 0 & FRP <= 300)


## 2001 to 2018 ################################################################
library(colorspace)


## Study area Lake Khamra
png(glue("Results/FRP/Maps/study_area_FRP_Khamra.png"), width = 1000, height = 800)
ggplot(Khamra_buf_extent) +
  geom_sf(data = MODIS_Khamra_FRP_300_2001_2018, aes(color = FRP), size = 2, na.rm = T, show.legend = T) +
  scale_color_continuous_sequential(palette = "Heat") +
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
  labs(title = "Fire radiative power for study area Lake Khamra | 2001 to 2018",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 16, vjust = 0.5),
        legend.text   = element_text(size = 16, vjust = 0.75))
dev.off()


## Study area Lake Satagay
png(glue("Results/FRP/Maps/study_area_FRP_Satagay.png"), width = 1000, height = 800)
ggplot(Satagay_buf_extent) +
  geom_sf(data = MODIS_Satagay_FRP_300_2001_2018, aes(color = FRP), size = 2, na.rm = T, show.legend = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.6) + 
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake Satagay"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake Satagay" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  xlab("") +
  ylab("") +
  labs(title = "Fire radiative power for study area Lake Satagay | 2001 to 2018",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 16, vjust = 0.5),
        legend.text   = element_text(size = 16, vjust = 0.75))
dev.off()



## Multiple Plot of FRP points into the study areas ############################
plot_FRP_1 <- 
  ggplot(Khamra_buf_extent) +
  geom_sf(data = MODIS_Khamra_FRP_300_2001_2018, aes(color = FRP), size = 2, na.rm = T, show.legend = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.4) + 
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Khamra, aes(fill = "Lake"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  xlab("") +
  ylab("") +
  labs(title = "Study area Lake Khamra",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
        legend.title     = element_text(size = 24, vjust = 1),
        legend.text      = element_text(size = 24),
        axis.text        = element_text(size = 18),
        legend.key.width = unit(1, 'cm'))

plot_FRP_2 <- 
  ggplot(Satagay_buf_extent) +
  geom_sf(data = MODIS_Satagay_FRP_300_2001_2018, aes(color = FRP), size = 2, na.rm = T, show.legend = F) +
  scale_color_continuous_sequential(palette = "Heat") +
  geom_sf(data = Satagay_catchment, fill = "mediumorchid1", alpha = 0.6) + 
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Satagay, fill = "turquoise1", alpha = 0.9) +
  xlab("") +
  ylab("") +
  labs(title = "Study area Lake Satagay",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
        legend.title     = element_text(size = 24, vjust = 1),
        legend.text      = element_text(size = 24),
        axis.text        = element_text(size = 18),
        legend.key.width = unit(1, 'cm'))


figure_1   <- ggarrange(plot_FRP_1, plot_FRP_2, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")


png(glue("Results/FRP/Maps/FRP_Map.png"), width = 1400, height = 800)
plot(figure_1)
dev.off()



## Fire years ###################################################################

MODIS_Khamra_FRP_300  <- fire_years_Khamra  %>% filter(FRP > 0 & FRP <= 300)
save(MODIS_Khamra_FRP_300, file = "Results/Simulation/data/MODIS_Khamra_FRP_300.rda")
MODIS_Khamra_FRP_300 <- get(load("Results/Simulation/data/MODIS_Khamra_FRP_300.rda"))

MODIS_Satagay_FRP_300 <- fire_years_Satagay %>% filter(FRP > 0 & FRP <= 300)
save(MODIS_Satagay_FRP_300, file = "Results/Simulation/data/MODIS_Satagay_FRP_300.rda")
MODIS_Satagay_FRP_300 <- get(load("Results/Simulation/data/MODIS_Satagay_FRP_300.rda"))

## Graphic-packages
library(maps)
library(colorspace)
library(ggnewscale)

## Study Lake Khamra
png(glue("Results/FRP/Maps/study_area_FRP_Khamra_fire_years.png"), width = 1000, height = 800)
ggplot(Khamra_buf_extent) +
  geom_sf(data = MODIS_Khamra_FRP_300, aes(color = FRP), size = 2, na.rm = T, show.legend = T) +
  scale_color_continuous_sequential(palette = "Heat") +
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
  labs(title = "Fire radiative power for study area Lake Khamra | Fire years",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 16, vjust = 0.5),
        legend.text   = element_text(size = 16, vjust = 0.75))
dev.off()


## Study area Lake Satagay
png(glue("Results/FRP/Maps/study_area_FRP_Satagay_fire_years.png"), width = 1000, height = 800)
ggplot(Satagay_buf_extent) +
  geom_sf(data = MODIS_Satagay_FRP_300, aes(color = FRP), size = 2, na.rm = T, show.legend = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.6) + 
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake Satagay"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake Satagay" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  xlab("") +
  ylab("") +
  labs(title = "Fire radiative power for study area Lake Satagay | Fire years",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title    = element_text(size = 24, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 16, vjust = 0.5),
        legend.text   = element_text(size = 16, vjust = 0.75))
dev.off()



## Multiple plot ###############################################################
FRP_plot_1 <- ggplot(Khamra_buf_extent) +
  geom_sf(data = MODIS_Khamra_FRP_300, aes(color = FRP), size = 2, na.rm = T, show.legend = F) +
  scale_color_continuous_sequential(palette = "Heat") +
  geom_sf(data = Khamra_catchment, fill = "mediumorchid1", alpha = 0.4)+
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Khamra, fill = "turquoise1", alpha = 0.9) +
  xlab("") +
  ylab("") +
  labs(title = "Study area Lake Khamra",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 17.5, vjust = 0.9),
        legend.text   = element_text(size = 17.5, vjust = 0.75),
        axis.text     = element_text(size = 14))

FRP_plot_2 <- ggplot(Satagay_buf_extent) +
  geom_sf(data = MODIS_Satagay_FRP_300, aes(color = FRP), size = 2, na.rm = T, show.legend = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.6) + 
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  xlab("") +
  ylab("") +
  labs(title = "Study area Lake Satagay",
       color = "FRP [MW]", shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 17.5, vjust = 0.9),
        legend.text   = element_text(size = 17, vjust = 0.75),
        axis.text     = element_text(size = 14))


figure_1   <- ggarrange(FRP_plot_1, FRP_plot_2, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")

png(glue("Results/FRP/Maps/FRP_fire_years_map.png"), width = 1400, height = 800)
plot(figure_1)
dev.off()



#####
## Calculation of the maximum of FRP ###########################################
#####



max(fire_years_Khamra$FRP)
max(fire_years_Satagay$FRP)

length(fire_years_Khamra$FRP > 300)
length(fire_years_Satagay$FRP > 300)

FRP_max_Khamra  <- max(sf_MODIS_Khamra_box$FRP)
FRP_max_Satagay <- max(sf_MODIS_Satagay_box$FRP)

## The maximum of FRP for each individual year
Khamra_FRP_max  <- sf_MODIS_Khamra_box %>% group_by(Year) %>%
                   summarise("max_FRP" = max(FRP, na.rm = TRUE))
Satagay_FRP_max <- sf_MODIS_Satagay_box %>% group_by(Year) %>%
                   summarise("max_FRP" = max(FRP, na.rm = TRUE))


## Filter the years with high max of FRP for classified fire years
fire_years_max_Khamra  <- Khamra_FRP_max %>%
                          group_by(Year) %>%
                          filter(Year == "2001" | Year == "2003" |
                                 Year == "2004" | Year == "2007" |
                                 Year == "2013" | Year == "2014" )

fire_years_max_Satagay <- Satagay_FRP_max %>%
                          group_by(Year)  %>%
                          filter(Year == "2002" | Year == "2013" |
                                 Year == "2014" | Year == "2018" )



#####
## Visualization of the years with highest max of FRP ##########################
#####

max_FRP_over_years <- function(lake_data, FRP_max_years_data, lake_name){
  
      max_FRP_plot <- ggplot() +
                      geom_bar(data = lake_data, mapping = aes(x = Year, y = max_FRP),
                               stat = 'identity', width = 0.6, colour = "black",
                               fill = "white", show.legend = T) +
                      geom_bar(data = FRP_max_years_data, mapping = aes(x = Year, y = max_FRP, fill = "Fire years"),
                               stat = 'identity', width = 0.6, alpha = 0.4, show.legend = T) +
                     scale_fill_manual(values = c("Fire years" = "firebrick3"), name = NULL) +
                     geom_text(data     = lake_data, aes(x = Year, y = max_FRP, label = round(max_FRP, 2)),
                               position = position_dodge(width = 1),
                               vjust    = -0.3, size = 4) +
                      scale_x_continuous(breaks = seq(2001,2018, by = 1)) +
                      theme_minimal()   +
                      labs(title = glue("The maximum values of FRP for study area Lake {lake_name}")) +
                      xlab("Year") +
                      ylab("Maximum of FRP [MW]") +
                      theme(plot.title = element_text(size = 24, hjust = 0.5, vjust = -1),
                            axis.text.x = element_text(size = 12,  hjust = .5),
                            axis.text.y = element_text(size = 18, hjust = 1),
                            axis.title.x = element_text(size = 20, hjust = .5, vjust = -0.5),
                            axis.title.y = element_text(size = 20, hjust = .5, vjust = 0.5),
                            axis.line = element_line(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            panel.background = element_blank(),legend.text = element_text(size = 20),
                                        text = element_text(size = 28)
                            )
                    png(glue("Results/FRP/Statistic/FRP_max_{lake_name}.png"), width = 1000, height = 800)
                    plot(max_FRP_plot)
                    dev.off()
}

## Run the function for study areas
Khamra_max_FRP_plot <- max_FRP_over_years(lake_data = Khamra_FRP_max,
                                          FRP_max_years_data = fire_years_max_Khamra,
                                          lake_name = "Khamra")
Satagy_max_FRP_plot <- max_FRP_over_years(lake_data = Satagay_FRP_max,
                                          FRP_max_years_data = fire_years_max_Satagay,
                                          lake_name = "Satagay")


#####
## Further calculations ########################################################
#####

## Sort the entire data to FRP
## Study area Lake Khamra
MODIS_Khamra_data       <- sf_MODIS_Khamra_box
MODIS_Khamra_data$Year  <- as.factor(MODIS_Khamra_data$Year)
Khamra_sort_FRP         <- MODIS_Khamra_data %>%
                           arrange(desc(FRP))
head(Khamra_sort_FRP)

## Study area Lake Satagay
MODIS_Satagay_data      <- sf_MODIS_Satagay_box
MODIS_Satagay_data$Year <- as.factor(MODIS_Satagay_data$Year)
Satagay_sort_FRP        <- MODIS_Satagay_data %>%
                           arrange(desc(FRP))
head(Khamra_sort_FRP)


## Calculate the modern active fire locations ##################################
leng_Khamra  <- length(sf_MODIS_Khamra_box$FRP)
leng_Satagay <- length(sf_MODIS_Satagay_box$FRP)

data_leng <- data.frame(name   = c("Khamra", "Satagay"),
                       values = c(leng_Khamra, leng_Satagay))


## The number of modern active fire locations
png(glue("Results/FRP/Statistic/boxplot_sum_FRP_lakes.png"), width = 1000, height = 800)
ggplot() +
  geom_bar(data = data_leng, mapping = aes(x = name, y = values),
           stat = 'identity',  width = 0.6, colour = "black", fill = "lightcoral") +
  theme_minimal()+
  labs(title = "The number of modern active fire locations")+
  xlab("Study area of Lake") +
  ylab("Number of FRP points") +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title = element_text(size = 20, hjust = 0.5),
        axis.text  = element_text(size = 18, hjust = 0.5))
dev.off()

