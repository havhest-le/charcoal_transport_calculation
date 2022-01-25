rm(list = ls(all= TRUE))


library(sp)
library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(ncdf4)
library(rnaturalearth)
library(tidyverse)
library(glue)
library(mapview)
library(units)

# Map
map     <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
ext     <- extent(c(103.82, 180, 50.07, 80.56))


#####
# Input data 
#####

load(file = "Results/lakes_buf.RData")
Elgy <- lakes_buf[["Elgy"]]
Kham <- lakes_buf[["Khamra"]]
Ill <- lakes_buf[["Ill"]]

# MODIS
MODIS <- read.csv2("N:/bioing/data/Projekte/Fires/East Siberia/R_Data/MODIS_east_siberia.txt")
MODIS <- select(.data = MODIS,"FID", "LATITUDE","LONGITUDE", "ACQ_DATE", "FRP")
names(MODIS)[names(MODIS) == "ACQ_DATE"] <- "TIME"
MODIS$TIME <- as.POSIXct(x = MODIS$TIME, format = '%d.%m.%Y %H:%M:%S')
MODIS$FRP <- as.numeric(MODIS$FRP)
MODIS$Year <- as.numeric(format(MODIS$TIME, "%Y"))
head(MODIS)

xy <- MODIS[,c(2,3)]

MODI_point <- st_as_sf(x = MODIS, 
                        coords = c("LONGITUDE", "LATITUDE"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


# Elgy
Inter_Elgy <- st_intersection(x = MODI_point, y = Elgy)
plot(Elgy)
plot(Inter_Elgy$geometry, add = T)

# Khamra
Inter_Kham <- st_intersection(x = MODI_point, y = Kham)
plot(Kham)
plot(Inter_Kham$geometry, add = T)

# Ill
Inter_Ill <- st_intersection(x = MODI_point, y = Ill)
plot(Ill)
plot(Inter_Ill$geometry, add = T)


Lakes_inter <- list(Elgy = Inter_Elgy, Kham = Inter_Kham, Ill = Inter_Ill)

#####
# Plotting
#####

for(i in 1:length(Lakes_inter)) {
  w <- Lakes_inter[[i]]
  
  w$Year <- as.factor(w$Year)
  meds <- w %>%
    group_by(Year) %>%
    summarise("med_FRP" = median(FRP, na.rm = TRUE))
  
  # boxplot
  boxplot <- ggplot(data = w, mapping = aes(x = Year, y = FRP))+
    stat_boxplot(geom = 'errorbar')+
    geom_boxplot(outlier.size = 1)+
    geom_text(data = meds, aes(y = med_FRP, label = round(med_FRP,2)),size = 2) +
    ylim(c(0,100))+
    theme_minimal()+
    scale_fill_manual(values=c("#69b3a2", "grey")) +
    scale_alpha_manual(values=c(1,0.1)) +
    theme(legend.position = "none") +
    xlab("Year") +
    ylab("The median of fire radiative power (FRP)") +
    theme(plot.title = element_text(size = 16, hjust = 0.5))
  
  # Barplot
  MODIS_year_median <- w %>%
    mutate(date = as.Date(TIME)) %>%
    mutate(year = format(date, '%Y')) %>%
    group_by(year) %>%
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
  fig_plots <- ggarrange(boxplot,barplot,nrow = 2, ncol = 1)
  print(annotate_figure(fig_plots, top = text_grob(glue("The averrage fire radiative power"), vjust = 0.8, hjust = 0.5, face = "bold", size = 30)))
  dev.off()
}

# sum of FRP 

# for(i in 1:length(Lakes_inter)) {
#   w <- Lakes_inter[[i]]
#   sum <- sum(w$FRP, na.rm = FALSE)
# }

sum_Elgy <- sum(Inter_Elgy$FRP, na.rm = T)
sum_Kham <- sum(Inter_Kham$FRP, na.rm = T)
sum_Ill  <- sum(Inter_Ill$FRP, na.rm = T)

data_sum <- data.frame(name = c("Kham", "Elgy", "Ill"), values = c(sum_Kham, sum_Elgy, sum_Ill))


  ggplot() +
    geom_bar(data = data_sum, mapping = aes(x = name, y = values), 
             stat = 'identity',  width = 0.6, colour = "black", fill = "white") +
    theme_minimal()+
    labs(title = "The amount fire radiative power")+
    xlab("Lake") +
    ylab("sum of FRP") +
    theme(plot.title = element_text(size = 20, hjust = 0.5))
  


#####
# For Elgy
#####
Inter_Elgy$Year <- as.factor(Inter_Elgy$Year)
meds <- Inter_Elgy %>%
  group_by(Year) %>%
  summarise("med_FRP" = median(FRP, na.rm = TRUE))

# boxplot
ggplot(data = Inter_Elgy, mapping = aes(x = Year, y = FRP))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(outlier.size = -2)+
  #geom_text(data = meds, aes(y = med_FRP, label = round(med_FRP,2)),size = 2) +
  ylim(c(0,70))+
  theme_minimal()+
  scale_fill_manual(values=c("#69b3a2", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme(legend.position = "none") +
  labs(title = "The averrage fire radiative power")+
  xlab("Year") +
  ylab("The median of fire radiative power (FRP)") +
  theme(plot.title = element_text(size = 16, hjust = 0.5))

# Barplot
MODIS_year_median <- Inter_Elgy %>%
  mutate(date = as.Date(TIME)) %>%
  mutate(year = format(date, '%Y')) %>%
  group_by(year) %>%
  count()

ggplot() +
  geom_bar(data = MODIS_year_median, mapping = aes(x = year, y = n),
           stat = 'identity', width = 0.6, colour = "black", fill = "white") +
  theme_minimal()+
  labs(title = "The averrage fire radiative power")+
  xlab("Year") +
  ylab("The amount pixels of FRP") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))


#####
# For Khamra
#####

Inter_Kham$Year <- as.factor(Inter_Kham$Year)
meds <- Inter_Kham %>%
  group_by(Year) %>%
  summarise("med_FRP" = median(FRP, na.rm = TRUE))

# boxplot
ggplot(data = Inter_Kham, mapping = aes(x = Year, y = FRP))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(outlier.size = -2)+
  #geom_text(data = meds, aes(y = med_FRP, label = round(med_FRP,2)),size = 2) +
  ylim(c(0,70))+
  theme_minimal()+
  scale_fill_manual(values=c("#69b3a2", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme(legend.position = "none") +
  labs(title = "The averrage fire radiative power")+
  xlab("Year") +
  ylab("The median of fire radiative power (FRP)") +
  theme(plot.title = element_text(size = 16, hjust = 0.5))

# Barplot
MODIS_year_median <- Inter_Kham %>%
  mutate(date = as.Date(TIME)) %>%
  mutate(year = format(date, '%Y')) %>%
  group_by(year) %>%
  count()

ggplot() +
  geom_bar(data = MODIS_year_median, mapping = aes(x = year, y = n),
           stat = 'identity', width = 0.6, colour = "black", fill = "white") +
  theme_minimal()+
  labs(title = "The averrage fire radiative power")+
  xlab("Year") +
  ylab("The amount pixels of FRP") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))


#####
# For Ill
#####

Inter_Ill$Year <- as.factor(Inter_Ill$Year)
meds <- Inter_Ill %>%
  group_by(Year) %>%
  summarise("med_FRP" = median(FRP, na.rm = TRUE))

# boxplot
ggplot(data = Inter_Ill, mapping = aes(x = Year, y = FRP))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(outlier.size = -2)+
  #geom_text(data = meds, aes(y = med_FRP, label = round(med_FRP,2)),size = 2) +
  ylim(c(0,70))+
  theme_minimal()+
  scale_fill_manual(values=c("#69b3a2", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme(legend.position = "none") +
  labs(title = "The averrage fire radiative power")+
  xlab("Year") +
  ylab("The median of fire radiative power (FRP)") +
  theme(plot.title = element_text(size = 16, hjust = 0.5))

# Barplot
MODIS_year_median <- Inter_Ill %>%
  mutate(date = as.Date(TIME)) %>%
  mutate(year = format(date, '%Y')) %>%
  group_by(year) %>%
  count()

ggplot() +
  geom_bar(data = MODIS_year_median, mapping = aes(x = year, y = n),
           stat = 'identity', width = 0.6, colour = "black", fill = "white") +
  theme_minimal()+
  labs(title = "The averrage fire radiative power")+
  xlab("Year") +
  ylab("The amount pixels of FRP") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))
