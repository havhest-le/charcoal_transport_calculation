###
# Script to calculate the burned areas #
###

rm(list = ls(all= TRUE))

# load packages
library(utils)
library(dplyr)
library(base)
library(tidyverse)
library(sf)
library(ggplot2)
library(glue)
library(graphics)

#####
# Load data ######################################################################
#####

# Wildfires #

#####
# MODIS data #
#####

MODIS_txt        <- read.csv2("Data/MODIS/MODIS_east_siberia.txt")
MODIS_txt        <- dplyr::select(.data = MODIS_txt, "FID", "LATITUDE","LONGITUDE", "ACQ_DATE", "FRP")
names(MODIS_txt)[names(MODIS_txt) == "ACQ_DATE"] <- "TIME"
MODIS_txt$TIME   <- as.POSIXct(x = MODIS_txt$TIME, format = '%d.%m.%Y %H:%M:%S')
MODIS_txt$FRP    <- as.numeric(MODIS_txt$FRP)
MODIS_txt$Year   <- as.numeric(format(MODIS_txt$TIME, "%Y"))


MODIS_for_region <- function(lon_bigger, lon_smaller, lat_bigger, lat_smaller){
                    MODIS_txt %>% 
                    filter (LONGITUDE >= lon_bigger & LONGITUDE <= lon_smaller &  
                            LATITUDE  >= lat_bigger & LATITUDE  <= lat_smaller & 
                            FRP > 0 & Year <= 2018)}

# Khamra MODIS # 
Khamra_MODIS          <- MODIS_for_region(lon_bigger = 111.16, 
                                          lon_smaller = 114.80,
                                          lat_bigger = 59.08, 
                                          lat_smaller = 60.9)
str(Khamra_MODIS)


Khamra_MODIS_point    <- st_as_sf(x = Khamra_MODIS, 
                                  coords = c("LONGITUDE", "LATITUDE"),
                                  crs = 4326)
str(Khamra_MODIS_point)

save(Khamra_MODIS_point, file = "Results/Khamra_MODIS_point.RData")

# Satagay MODIS #
Satagay_MODIS          <- MODIS_for_region(lon_bigger = 116.01, 
                                           lon_smaller = 120.00,
                                           lat_bigger = 62.17, 
                                           lat_smaller = 63.98)
Satagay_MODIS_point    <- st_as_sf(x = Satagay_MODIS, 
                                   coords = c("LONGITUDE", "LATITUDE"),
                                   crs = 4326)

save(Satagay_MODIS_point, file = "Results/Satagay_MODIS_point.RData")

points_lakes           <- list(Khamra_MODIS_point, Satagay_MODIS_point)

#####################################################################################

#####
# Calculation of data #
#####

# Median of FRP and identification of fire years #
median_FRP_Khamra    <- Khamra_MODIS_point %>% 
                        group_by(Year)     %>% 
                        summarise("med_FRP" = median(FRP, na.rm = TRUE))
median_FRP_Satagay   <- Satagay_MODIS_point %>% 
                        group_by(Year)      %>% 
                        summarise("med_FRP" = median(FRP, na.rm = TRUE))

statstic_FRP_Khamra  <- summary(Khamra_MODIS_point$FRP)
statstic_FRP_Satagay <- summary(Satagay_MODIS_point$FRP)

statstic_FRP_Khamra[4]  <- round(statstic_FRP_Khamra[4], 1)
statstic_FRP_Satagay[4] <- round(statstic_FRP_Satagay[4], 1)


FRP_boxplot_list        <- list(Khamra_MODIS_point$FRP, Satagay_MODIS_point$FRP)
names(FRP_boxplot_list) <- c(paste("                   Khamra\n                    n = " , length(Khamra_MODIS_point$FRP) , sep=""), 
                             paste("                  Satagay\n                     n = " , length(Satagay_MODIS_point$FRP) , sep=""))

 

#####
# Fire years #####################################################################
#####

# Filtering all values over the median #
med_FRP_Khamra      <- statstic_FRP_Khamra[3]
med_FRP_Satagay     <- statstic_FRP_Satagay[3]

# How many days are above the averrage of FRP? #
fire_years_Khamra_days  <- Khamra_MODIS_point %>%
                           filter(FRP >= med_FRP_Khamra) %>%
                           group_by(TIME) %>%
                           summarise("med_FRP" = median(FRP, na.rm = TRUE))
length(fire_years_Khamra_days$TIME) # 174

fire_years_Satagay_days  <- Satagay_MODIS_point %>%
                            filter(FRP >= med_FRP_Satagay) %>%
                            group_by(TIME) %>%
                            summarise("med_FRP" = median(FRP, na.rm = TRUE))
length(fire_years_Satagay_days$TIME) # 275

# How many days are under the 1st and above the 3rd quartile? #
first_qua_FRP_Khamra         <- statstic_FRP_Khamra[2]
first_qua_FRP_Satagay        <- statstic_FRP_Satagay[2]
third_qua_FRP_Khamra         <- statstic_FRP_Khamra[5]
third_qua_FRP_Satagay        <- statstic_FRP_Satagay[5]

under_first_qua_Khamra_days  <- Khamra_MODIS_point %>%
                                filter(FRP <= first_qua_FRP_Khamra) %>%
                                group_by(TIME) %>%
                                summarise("med_FRP" = median(FRP, na.rm = TRUE))
length(under_first_qua_Khamra_days$TIME) # 271

under_first_qua_Satagay_days  <- Satagay_MODIS_point %>%
                                 filter(FRP <= first_qua_FRP_Satagay) %>%
                                 group_by(TIME) %>%
                                 summarise("med_FRP" = median(FRP, na.rm = TRUE))
length(under_first_qua_Satagay_days$TIME) # 450

above_third_qua_Khamra_days  <- Khamra_MODIS_point %>%
                                filter(FRP >= third_qua_FRP_Khamra) %>%
                                group_by(TIME) %>%
                                summarise("med_FRP" = median(FRP, na.rm = TRUE))
length(above_third_qua_Khamra_days$TIME) # 122

above_third_qua_Satagay_days  <- Satagay_MODIS_point %>%
                                 filter(FRP >= third_qua_FRP_Satagay) %>%
                                 group_by(TIME) %>%
                                 summarise("med_FRP" = median(FRP, na.rm = TRUE))
length(above_third_qua_Satagay_days$TIME) # 171


# Identification of fire years # 
fire_years_Khamra_med   <- Khamra_MODIS_point                               %>%
                           group_by(Year)                                   %>% 
                           summarise("med_FRP" = median(FRP, na.rm = TRUE)) %>%
                           filter(med_FRP >= med_FRP_Khamra)
fire_years_Khamra_med

fire_years_Satagay_med  <- Satagay_MODIS_point                              %>%
                           group_by(Year)                                   %>% 
                           summarise("med_FRP" = median(FRP, na.rm = TRUE)) %>%
                           filter(med_FRP >= med_FRP_Satagay)
fire_years_Satagay_med


# Identification of fire years #
fire_years_Khamra  <- Khamra_MODIS_point                               %>%
                      group_by(Year)                                   %>%
                      filter(Year == "2001" | Year == "2003" | 
                             Year == "2004" | Year == "2007" | 
                             Year == "2013" | Year == "2014" )

fire_years_Satagay <- Satagay_MODIS_point                              %>%
                      group_by(Year)                                   %>%
                      filter(Year == "2002" | Year == "2013" | 
                             Year == "2014" | Year == "2017" | 
                             Year == "2018" )


# PLotting the boxplot with median of FRP #

boxplot_FRP <- function(lakes, fire_years_data, vjust, lake_name) {
  
               data            <- points_lakes[[lakes]]
               data$Year       <- as.factor(data$Year)
               median_FRP      <- data %>%
                                  group_by(Year) %>%
                                  summarise("med_FRP" = median(FRP, na.rm = TRUE))
                                
               fire_years      <- fire_years_data
               fire_years$Year <- as.factor(fire_years$Year)
               
               cols            <- c("with median\nabove the averrage\nmedian over all years" = "grey80")
               
               boxplot         <-   ggplot(data = data, mapping = aes(x = Year, y = FRP)) +
                                    stat_boxplot(geom = 'errorbar') +
                                    geom_boxplot(outlier.size = 0.1, colour = "black", fatten = 0.1) +
                                    geom_text(data     = median_FRP, aes(y = med_FRP, label = round(med_FRP, 2)),
                                              position = position_dodge(width = 1),
                                              vjust    = vjust, size = 3) +
                                    geom_boxplot(mapping = aes(x = Year, y = FRP),
                                                 data    = fire_years, fill = "grey80", alpha= 0.4, fatten = 0.0001,
                                                 show.legend = T) + 
                                    scale_fill_manual(name = "fire years ", values = cols) + 
                                    ylim(c(0,100))  +
                                    theme_minimal() +
                                    scale_alpha_manual(values= c(1,0.1)) +
                                    labs(title = glue("Fire Radiative Power of lake {lake_name} 2001-2018"),
                                         subtitle = "FRP with median values and outliners over the 3rd quartile") +
                                    xlab("Year") +
                                    ylab("FRP [MW / km²]") +
                                    theme(plot.title = element_text(size = 24, hjust = 0.5, vjust = -0.5),
                                          plot.subtitle = element_text(size = 14, hjust = 0.5, vjust = -2),
                                          legend.text = element_text(size= 12),
                                          legend.title = element_text(size= 14,face = "bold"))
                                  

               # Creating boxplot #
               png(glue("Results/boxplot_median_FRP_{lake_name}.png"), width = 900, height = 600)
               plot(boxplot)
               dev.off()
               }
 
  

Khamra_boxplot_FRP  <- boxplot_FRP(lakes = 1, fire_years_data = fire_years_Khamra,
                                   vjust = -0.2, lake_name = "Khamra")
Satagay_boxplot_FRP <- boxplot_FRP(lakes = 2, fire_years_data = fire_years_Satagay, 
                                   vjust = -0.5, lake_name = "Satagay")


# Calculation the max of fire intensiv [FRP] #
FRP_max_Khamra  <- max(Khamra_MODIS_point$FRP)  # 2998.6 in year 2014
FRP_max_Satagay <- max(Satagay_MODIS_point$FRP) # 5112.5 in year 2017

Khamra_FRP_max  <- Khamra_MODIS_point %>% group_by(Year) %>%
                   summarise("max_FRP" = max(FRP, na.rm = TRUE))
Satagay_FRP_max <- Satagay_MODIS_point %>% group_by(Year) %>%
                   summarise("max_FRP" = max(FRP, na.rm = TRUE))

fire_years_max_Khamra  <- Khamra_FRP_max %>%
                          group_by(Year) %>%
                          filter(Year == "2001" | Year == "2003" | 
                                 Year == "2006" | Year == "2007" | 
                                 Year == "2014" | Year == "2016" )
fire_years_max_Satagay <- Satagay_FRP_max %>%
                          group_by(Year) %>%
                          filter(Year == "2002" | Year == "2013" |
                                 Year == "2014" | Year == "2015" | 
                                 Year == "2017" | Year == "2018" )


# Plotting
max_FRP_over_years <- function(lake_data, FRP_max_years_data, lake_name){

max_FRP_plot <- ggplot() +
                geom_bar(data = lake_data, mapping = aes(x = Year, y = max_FRP), 
                         stat = 'identity', width = 0.6, colour = "black", 
                         fill = "white", show.legend = T) +
                geom_bar(data = FRP_max_years_data, mapping = aes(x = Year, y = max_FRP),
                         stat = 'identity', width = 0.6, fill = "firebrick3", alpha = 0.4, show.legend = T) +
                geom_text(data     = lake_data, aes(x = Year, y = max_FRP, label = round(max_FRP, 2)),
                          position = position_dodge(width = 1),
                          vjust    = -0.3, size = 4) +
                scale_x_continuous(breaks = seq(2000,2018, by = 1)) +
                
                theme_minimal()   +
                labs(title = glue("The maximum of FRP of the 100 km puffer around lake {lake_name}")) +
                xlab("Year") +
                ylab("max of FRP") +
                theme(plot.title = element_text(size = 22, hjust = 0.5, vjust = -1),
                      axis.text.x = element_text(size = 10,  hjust = .5),
                      axis.text.y = element_text(size = 10, hjust = 1),
                      axis.title.x = element_text(size = 15, hjust = .5, vjust = -0.5),
                      axis.title.y = element_text(size = 15, hjust = .5, vjust = 0.5),
                      axis.line = element_line(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
png(glue("Results/FRP_max_{lake_name}.png"), width = 900, height = 600)
plot(max_FRP_plot)  
dev.off()
}

Khamra_max_FRP_plot <- max_FRP_over_years(lake_data = Khamra_FRP_max, 
                                          FRP_max_years_data = fire_years_max_Khamra,
                                          lake_name = "Khamra")
Satagy_max_FRP_plot <- max_FRP_over_years(lake_data = Satagay_FRP_max, 
                                          FRP_max_years_data = fire_years_max_Satagay,
                                          lake_name = "Satagay")


# Sort the entire data to FRP # 
MODIS_Khamra_data       <- Khamra_MODIS_point
MODIS_Khamra_data$Year  <- as.factor(MODIS_Khamra_data$Year)
Khamra_sort_FRP         <- MODIS_Khamra_data %>% 
                           arrange(desc(FRP))
head(MODIS_Khamra_data)


# Filtering the single years and sort the FRP values from max to min # 
filter_years_sort_FRP   <- function(data_MODIs_point, year){
                           data_MODIs_point     %>% 
                           group_by(Year)       %>% 
                           filter(Year == year) %>% 
                           arrange(desc(FRP))}

Khamra_2001 <- filter_years_sort_FRP(data_MODIs_point = MODIS_Khamra_data, year = "2001")
Khamra_2007 <- filter_years_sort_FRP(data_MODIs_point = MODIS_Khamra_data, year = "2007")
Khamra_2014 <- filter_years_sort_FRP(data_MODIs_point = MODIS_Khamra_data, year = "2014")

head(Khamra_2001)
head(Khamra_2007)
head(Khamra_2014)




# Calculation the sum of FRP #
sum_Khamra  <- sum(Khamra_MODIS_point$FRP, na.rm = T)
sum_Satagay <- sum(Satagay_MODIS_point$FRP, na.rm = T)

data_sum <- data.frame(name   = c("Khamra", "Satagay"), 
                       values = c(sum_Khamra, sum_Satagay))

png(glue("Results/boxplot_sum_FRP_lakes.png"), width = 800, height = 800)
ggplot() +
  geom_bar(data = data_sum, mapping = aes(x = name, y = values), 
           stat = 'identity',  width = 0.6, colour = "black", fill = "firebrick") +
  theme_minimal()+
  labs(title = "The amount fire radiative power")+
  xlab("Lake") +
  ylab("sum of FRP") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))
dev.off()


#####
# Hansen data #
##### 

hansen_60_110_loss  <- raster("Data/HANSEN/Hansen_GFC-2020-v1.8_lossyear_60N_110E.tif")
hansen_70_110_loss  <- raster("Data/HANSEN/Hansen_GFC-2020-v1.8_lossyear_70N_110E.tif")

# Merging of dataset #
hansen_merge_raster <- function(map_60, map_70, extent){
        crop_map_60 <- crop(map_60, extent)
        crop_map_70 <- crop(map_70, extent)
        merge_map   <- merge(crop_map_60, crop_map_70, tolerance = 0.5, overlap=TRUE, ext = NULL)}

# Khamra #
hansen_loss_Khamra   <- hansen_merge_raster(map_60 = hansen_60_110_loss, map_70 = hansen_70_110_loss, extent = c(111.16, 114.80, 59.08, 60.9))
values(hansen_loss_Khamra)[values(hansen_loss_Khamra) <= 0] = NA
# save("hansen_loss_Khamrae.tif") # convert raster to polygon in QGIS because of the amount of ncell
writeRaster(hansen_loss_Khamra, "Data/HANSEN/hansen_loss_Khamra.tif", overwrite=TRUE) 

# Satagay #
extent_Satagay       <- extent(116.01, 120, 62.17, 63.98)
hansen_loss_Satagay  <- crop(hansen_70_110_loss, extent_Satagay)
values(hansen_loss_Satagay)[values(hansen_loss_Satagay) <= 0] = NA
# save("hansen_loss_Khamra.tif") # convert raster to polygon in QGIS because of the amount of ncell
writeRaster(hansen_loss_Satagay, "Data/HANSEN/hansen_loss_Satagay.tif", overwrite=TRUE)            

# Load shape file data #
hansen_loss_Khamra_polygon  <- read_sf("Data/HANSEN/hansen_loss_Khamra_polygon.shp") %>% 
                               st_transform(4326)
hansen_loss_Satagay_polygon <- read_sf("Data/HANSEN/hansen_loss_Satagay_polygon.shp") %>% 
                               st_transform(4326)

# change DN to years # 
names(hansen_loss_Khamra_polygon)[names(hansen_loss_Khamra_polygon) == "DN"]   <- "year"
hansen_loss_Khamra_polygon$year <- as.numeric(hansen_loss_Khamra_polygon$year)   + 2000
names(hansen_loss_Satagay_polygon)[names(hansen_loss_Satagay_polygon) == "DN"] <- "year"
hansen_loss_Satagay_polygon$year <- as.numeric(hansen_loss_Satagay_polygon$year) + 2000

# Joining Hansen and Modis data # 
sf::sf_use_s2(FALSE)
join_point_Khamra_polygon  <- st_join(x    = Khamra_MODIS_point,  y = hansen_loss_Khamra_polygon, 
                                      join = st_intersects)
join_point_Satagay_polygon <- st_join(x    = Satagay_MODIS_point, y = hansen_loss_Satagay_polygon, 
                                      join = st_intersects)

join_point_Khamra_polygon$coordinates  <- st_coordinates(join_point_Khamra_polygon)
join_point_Satagay_polygon$coordinates <- st_coordinates(join_point_Satagay_polygon)
save(join_point_Khamra_polygon,  file   = "Data/HANSEN/join_point_Khamra_polygon.Rdata")
save(join_point_Satagay_polygon, file   = "Data/HANSEN/join_point_Satagay_polygon.Rdata")

points_lakes <- list(join_point_Khamra_polygon, join_point_Satagay_polygon)

#####
# Calcultion of fire years #
#####

#for(i in 1:length(points_lakes)) {
    data       <- join_point_Khamra_polygon
    data$Year  <- as.factor(data$Year)
    median_FRP <- data %>%
                  group_by(Year) %>%
                  summarise("med_FRP" = median(FRP, na.rm = TRUE))
    
    # Creating boxplot #
    boxplot <- ggplot(data = data, mapping = aes(x = Year, y = FRP))+
               stat_boxplot(geom = 'errorbar')+
               geom_boxplot(outlier.size = 1) +
               geom_text(data = median_FRP, aes(y = med_FRP, label = round(med_FRP, 2)), size = 2) +
               ylim(c(0,100))  +
               theme_minimal() +
               scale_fill_manual(values= c("#69b3a2", "grey")) +
               scale_alpha_manual(values= c(1,0.1)) +
               theme(legend.position = "none") +
               xlab("Year") +
               ylab("The median of fire radiative power (FRP)") +
               theme(plot.title = element_text(size = 16, hjust = 0.5))
    
    # Creating barplot #
    MODIS_year_median <- data %>%
                         mutate(date   = as.Date(TIME)) %>%
                         mutate(year_2 = format(date, '%Y')) %>%
                         group_by(year_2) %>%
                         count()
    
    barplot <- ggplot() +
               geom_bar(data = MODIS_year_median, mapping = aes(x = year, y = n),
                        stat = 'identity', width = 0.6, colour = "black", fill = "white") +
               theme_minimal()+
               labs(title = glue("Lakes_inter{}")) +
               xlab("Year") +
               ylab("The amount pixels of FRP") +
               theme(plot.title = element_text(size = 20, hjust = 0.5))
            
               png(glue("Results/MODIS_Plot{i}.png"), width = 1200, height = 1200)
               fig_plots <- ggarrange(boxplot, barplot,nrow = 2, ncol = 1)
               print(annotate_figure(fig_plots, 
                                     top   = text_grob(glue("The averrage fire radiative power"), 
                                     vjust = 0.8, hjust = 0.5, face = "bold", size = 30)))
               dev.off()
               }









#####
# Plotting #
#####
plot(join_point_Khamra_polygon$geometry, col = "firebrick")
plot(Khamra$geometry, col = "blue", add = T)

plot(join_point_Satagay_polygon$geometry, col = "firebrick")
plot(Satagay$geometry, col = "blue", add = T)
