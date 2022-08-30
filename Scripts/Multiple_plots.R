#####
## Script to visualize the results of the both study region in multiple plots ##
#####

## rm(list = ls(all= TRUE))

#####
## Load packages ###############################################################
#####

library(tidyverse)
library(sf)
library(terra)
library(stars)
library(dplyr)
library(glue)
library(base)
library(raster)

## Plotting packages
library(ggpubr)
library(gridExtra)
library(ggsn)
library(ggnewscale)
library(colorspace)
library(maps)
library(viridis)
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(MASS)
library(survival)
library(fitdistrplus)
library(geosphere)

#####
## Lake Khamra #################################################################
######
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

#####
## Lake Satagaay 2.0 ###########################################################
######
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


#####
## Calculated plume injection heights ##########################################
#####

H_function_Khamra <- get(load("Results/Simulation/data/H_function_Khamra.rda"))
H_function_Khamra$FRP[H_function_Khamra$FRP == 0] <- NA
H_function_Khamra$H[H_function_Khamra$H == 0]     <- NA
H_function_Khamra                                 <- na.omit(H_function_Khamra)
H_function_Khamra$FRP_MW                          <- H_function_Khamra$FRP/238845.8966275
#
H_function_Satagay <- get(load("Results/Simulation/data/H_function_Satagay.rda"))
H_function_Satagay$FRP[H_function_Satagay$FRP == 0] <- NA
H_function_Satagay$H[H_function_Satagay$H == 0]     <- NA
H_function_Satagay                                  <- na.omit(H_function_Satagay)
H_function_Satagay$FRP_MW                           <- H_function_Satagay$FRP/238845.8966275

#####
## Observed plume injection heights ############################################
#####

PIH_Khamra   <- get(load("Results/Simulation/PIH/data/PIH_points_Khamra.rda"))
PIH_Satagay  <- get(load("Results/Simulation/PIH/data/PIH_points_Satagay.rda"))


#####
## Visualization ###############################################################
#####

H_function_Khamra_2500  <- H_function_Khamra %>% filter(H <= 2000)
H_function_Satagay_2500 <- H_function_Satagay %>% filter(H <= 2000)

H_points_Khamra  <- st_as_sf(H_function_Khamra_2500)
H_points_Satagay <- st_as_sf(H_function_Satagay_2500)


H_plot_1 <- ggplot(Khamra_buf_extent) +
            geom_sf(data = H_points_Khamra, aes(color = H), size = 2, na.rm = T, 
                    show.legend = T) +
            scale_color_distiller(palette = "Greys", direction = 1, limit = c(1,2000), 
                                  breaks = c(1,1000,2000))+
            geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.4) + 
            scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                              guide = guide_legend(override.aes = list(linetype = "blank",  
                                                                       shape = NA))) +  theme_minimal() +
            theme_minimal() +
            new_scale("fill")+
            geom_sf(data = Khamra, aes(fill = "Lake"), alpha = 0.9) +
            scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                              guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
            xlab("") +
            ylab("") +
            theme_minimal()+
            labs(title = "Study area Lake Khamra",
                 color = "Calculated PIH [m]") +
            theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
                  legend.title     = element_text(size = 20, vjust = 1),
                  legend.text      = element_text(size = 20),
                  axis.text        = element_text(size = 16),
                  legend.key.width = unit(1, 'cm'))

H_plot_2 <- ggplot(Satagay_buf_extent) +
            geom_sf(data = H_points_Satagay, aes(color = H), size = 2, na.rm = T, 
                    show.legend = F) +
            scale_color_distiller(palette = "Greys", direction = 1, limit = c(1,2000)) +
            geom_sf(data = Satagay_catchment, fill = "mediumorchid1", alpha = 0.8) +
            theme_minimal() +
            new_scale("fill")+
            geom_sf(data = Satagay, fill = "turquoise1", alpha = 0.9) +
            xlab("") +
            ylab("") +
            labs(title = "Study area Lake Satagay",
                 color = "Calculated PIH [m]", shape = "Lake") +
            theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
                  legend.title     = element_text(size = 20, vjust = 1),
                  legend.text      = element_text(size = 20),
                  axis.text        = element_text(size = 16),
                  legend.key.width = unit(1, 'cm'))


PIH_plot_3 <- ggplot(Khamra_buf_extent) +
              geom_sf(data = PIH_points_Khamra, aes(color = `50%`), size = 2, na.rm = T, 
                      show.legend = T) +
              scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000), 
                                    breaks = c(1,1000,2000))+
            
              geom_sf(data = Khamra, aes(fill = "Lake"), alpha = 0.9) +
              scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                                guide = guide_legend(override.aes = list(linetype = "blank", 
                                                                         shape = NA))) +
              new_scale("fill")+
               geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.4) + 
              scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                                guide = guide_legend(override.aes = list(linetype = "blank",  
                                                                         shape = NA))) +  theme_minimal() +
              theme_minimal() +
              xlab("") +
              ylab("") +
              labs(color = "Measured PIH [m] | Q[0.5]") +
              theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
                    legend.title     = element_text(size = 20, vjust = 1),
                    legend.text      = element_text(size = 20),
                    axis.text        = element_text(size = 16),
                    legend.key.width = unit(1, 'cm'))


PIH_plot_4 <- ggplot(Satagay_buf_extent) +
              geom_sf(data = PIH_points_Satagay, aes(color = `50%`), size = 2, na.rm = T, 
                      show.legend = F) +
              scale_color_distiller(palette = "Greys", direction = 1, limits = c(1,2000), 
                                    breaks = c(1,1000,2000))+
              geom_sf(data = Satagay_catchment, fill = "mediumorchid1", alpha = 0.8) + 
              theme_minimal() +
              new_scale("fill")+
              geom_sf(data = Satagay, fill = "turquoise1", alpha = 0.9) +
              xlab("") +
              ylab("") +
              labs(color = "Measured PIH [m] | Q[0.5]", shape = "Lake") +
              theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
                    legend.title     = element_text(size = 20, vjust = 1),
                    legend.text      = element_text(size = 20),
                    axis.text        = element_text(size = 16),
                    legend.key.width = unit(1, 'cm'))


## Multiple plot
library(ggpubr)

figure_1   <- ggarrange(H_plot_1, H_plot_2, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")

figure_2   <- ggarrange(PIH_plot_3, PIH_plot_4, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")



## Saving as png 
png(glue("Results/Simulation/H/Calculated_PIH_Map.png"), width = 1400, height = 800)
plot(figure_1)
dev.off()

png(glue("Results/Simulation/PIH/Measured_PIH_Map.png"), width = 1400, height = 800)
plot(figure_2)
dev.off()



#####
## Multiple plots of the link between FRP, PIH and wind speed ##################
#####

H_function_Khamra$FRP_MW <- H_function_Khamra$FRP/238845.8966275
plot_link_1 <- ggplot(data = H_function_Khamra,
                      aes(x = FRP_MW, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10), limit = c(1,20), guide = "none")+
  labs(title = "Study area of Lake Khamra")+
  xlab("FRP [MW]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 300))+
  ylim(c(0, 2000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -2),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))

H_function_Satagay$FRP_MW    <- H_function_Satagay$FRP/238845.8966275
plot_link_2 <- ggplot(data = H_function_Satagay,
                      aes(x = FRP_MW, y = H, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10),limit = c(1,20), guide = "none")+ 
  #breaks = c(2,4,6,8,10,12,14,16)) +
  labs(title = "Study area of Lake Satagay")+
  xlab("FRP [MW]") +
  ylab("Calculated PIH [m]")+
  xlim(c(0, 300))+
  ylim(c(0, 2000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -2),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))

## save
tibble_frp_pih_wind_Khamra   <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Khamra.rda"))
tibble_frp_pih_wind_Satagay  <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))



plot_link_3  <- ggplot(data = tibble_frp_pih_wind_Khamra,
                       aes(x = FRP, y = `50%`, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10), limit = c(1,20), guide = "none")+ 
  xlab("FRP [MW]") +
  ylab("Measured PIH [m] |  Q[0.5]")+
  xlim(c(0, 300))+
  ylim(c(0, 2000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))

plot_link_4 <- ggplot(data = tibble_frp_pih_wind_Satagay,
                      aes(x = FRP, y = `50%`, color = w_spd))+
  geom_point() +
  scale_color_gradientn(colours = rainbow(10), limit = c(1,20))+ 
  labs(color = "spd [m/s]")+
  xlab("FRP [MW]") +
  ylab("Measured PIH [m] |  Q[0.5]")+
  xlim(c(0, 300))+
  ylim(c(0, 2000))+
  theme_minimal() +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -4),
        axis.title    = element_text(size = 20, hjust = 0.5),
        legend.title  = element_text(size = 20, vjust = 0.5),
        legend.text   = element_text(size = 18, vjust = 0.75),
        axis.text     = element_text(size = 16, vjust = 0.75))


## Save
library(ggpubr)
library(gridExtra)

png(glue("Results/Simulation/Linkage_PIH_FRP_spd_calculated.png"), width = 1000, height = 900)
g <- ggarrange(plot_link_1, plot_link_2, plot_link_3, plot_link_4, 
               ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
plot(g)
dev.off()



#####
## Linear regression ###########################################################
#####

## Khamra | Calculated PIH #####################################################
lm_FRP_H_Khamra_MW  <- lm(H ~ FRP_MW, data = H_function_Khamra) # H ~ FRP_MW
coeff_Khamra_MW     <- coefficients(lm_FRP_H_Khamra_MW)
## (Intercept)             FRP 
## 742.942118         8.046233 
eq_K_MW <- paste0("y = ", round(coeff_Khamra_MW[1],2), " + ", round(coeff_Khamra_MW[2], 2), "*x")

## Satagay | Calculated PIH ####################################################
lm_FRP_H_Satagay_MW  <- lm(H ~ FRP_MW, data = H_function_Satagay)
coeff_satagay_MW     <- coefficients(lm_FRP_H_Satagay_MW)
## (Intercept)        FRP 
## 551.916070    3.627762 
eq_S_MW <- paste0("y = ", round(coeff_satagay_MW[1],2), " + ", round(coeff_satagay_MW[2], 2), "*x")


## Khamra | Meausred PIH #######################################################
tibble_frp_pih_wind_Khamra         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Khamra.rda"))
## Convert MW to cal/s
tibble_frp_pih_wind_Khamra$FRP_cal <- tibble_frp_pih_wind_Khamra$FRP*238845.8966275
lm_FRP_PIH_Khamra_MW  <- lm(`50%` ~ FRP, data = tibble_frp_pih_wind_Khamra)
coeff_Khamra_PIH_MW   <- coefficients(lm_FRP_PIH_Khamra_MW)
## (Intercept)             FRP 
## 318.45117685   0.09479513
eq_K_PIH_MW <- paste0("y = ", round(coeff_Khamra_PIH_MW[1],2), " + ", round(coeff_satagay_MW[2], 2), "*x")


## Satagay | Meausred PIH ######################################################
tibble_frp_pih_wind_Satagay         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))
## Convert MW to cal/s
tibble_frp_pih_wind_Satagay$FRP_cal <- tibble_frp_pih_wind_Satagay$FRP*238845.8966275
lm_FRP_PIH_Satagay_MW  <- lm(`50%` ~ FRP, data = tibble_frp_pih_wind_Satagay)
coeff_satagay_PIH_MW   <- coefficients(lm_FRP_PIH_Satagay_MW)
## (Intercept)           FRP 
## 370.29409101  -0.03533874
eq_S_PIH_MW <- paste0("y = ", round(coeff_satagay_PIH_MW[1],2), " ", round(coeff_satagay_PIH_MW[2], 2), "*x")



png(glue("Results/Simulation/lm_PIH_H.png"), width = 1000, height = 800)
layout(matrix(c(1,2,3,4),nrow = 2, ncol = 2, byrow = TRUE))

plot(x = H_function_Khamra$FRP_MW, y = H_function_Khamra$H,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.75, cex.axis = 1.75)
abline(lm_FRP_H_Khamra_MW, col = "red", lwd = 2)
mtext(mtext(eq_K_MW, 3, line=-23, col = "red", cex = 1.5))

plot(x = H_function_Satagay$FRP_MW, y = H_function_Satagay$H,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "", col = "grey40", cex.lab = 1.75, cex.axis = 1.75)
abline(lm_FRP_H_Satagay_MW, col = "red", lwd = 2)
mtext(mtext(eq_S_MW, 3, line=-23, col = "red", cex = 1.5))

plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`50%`,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.5]", col = "lightblue3", cex.lab = 1.75, cex.axis = 1.75)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
mtext(mtext(eq_K_PIH_MW, 3, line=-3, col = "red", cex = 1.5))

plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`50%`,
     xlim = c(0, 300), ylim = c(0, 2000),
     xlab = "FRP [MW]",
     ylab = "", col = "lightblue3", cex.lab = 1.75, cex.axis = 1.75)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2, )
mtext(mtext(eq_S_PIH_MW, 3, line=-3, col = "red", cex = 1.5))

dev.off()


#####
## Multiple regression #########################################################
#####

multi_Satagay_H_MW    <- lm(H ~ FRP_MW + w_spd, data = H_function_Satagay)
summary(multi_Satagay_H_MW) # H ~ FRP_MW + w_spd
mul_coeff_Satagay_MW  <- coefficients(multi_Satagay_H_MW)

multi_Khamra_H_MW    <- lm(H ~ FRP_MW + w_spd, data = H_function_Khamra)
mul_coeff_Khamra_MW  <- coefficients(multi_Khamra_H_MW)

multi_Satagay_PIH_MW     <- lm(`50%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Satagay)
mul_coeff_Satagay_PIH_MW <- coefficients(multi_Satagay_PIH_MW)

multi_Khamra_PIH_MW      <- lm(`50%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Khamra)
mul_coeff_Khamra_PIH_MW  <- coefficients(multi_Khamra_PIH_MW)

#####
## 3-D plot ####################################################################
#####


library(scatterplot3d)
png(glue("Results/Simulation/mul_PIH_H.png"), width = 1100, height = 800)
layout(matrix(c(1,2,3,4),nrow = 2, ncol = 2, byrow = TRUE))

scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP_MW, z = H_function_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), color = "grey40", 
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1.2, cex.lab = 1.6)

scatterplot3d(y = H_function_Satagay$H, x = H_function_Satagay$FRP_MW, z = H_function_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), color = "grey40", 
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "", 
              type = "h", font.lab = NULL, cex.axis = 1.2, cex.lab = 1.6)

scatterplot3d(y = tibble_frp_pih_wind_Khamra$`50%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.5]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1.2, cex.lab = 1.6)

scatterplot3d(y = tibble_frp_pih_wind_Satagay$`50%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.5]", zlab = "", 
              type = "h", font.lab = NULL, cex.axis = 1.2, cex.lab = 1.6)
dev.off()



#####
## Backward simulation results #################################################
#####
source_H_K <- ggplot(Khamra_buf_extent) +
  geom_sf(data = dist_data_crop_Khamra, aes(colour  = "Locations of\np-source areas\n(calculated PIH)"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(calculated PIH)" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Khamra_VR, aes(color = FRP_MW), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(color = "FRP [MW]") +
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, aes(fill = "Catchment"), alpha = 0.4) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA, alpha = 0.3))) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[1]], y = lat[[1]], col = "Lake"), data = data_coord_lakes, size = 3, stroke = 1.5) +
  scale_color_manual(values = c("Lake" = "turquoise1"),  name = NULL) +
  xlab("") +
  ylab("") +
  labs(title = "Study area of Lake Khamra")+
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
        legend.title     = element_text(size = 20, vjust = 1),
        legend.text      = element_text(size = 20),
        axis.text        = element_text(size = 16),
        legend.key.width = unit(1, 'cm'))


source_H_S <- ggplot(Satagay_buf_extent) +
  geom_sf(data = dist_data_crop_Satagay, aes(colour  = "Locations of\np-source areas\n(calculated PIH)"), show.legend = F) +
  scale_colour_manual(values = c("Locations of\np-source areas\n(calculated PIH)" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Satagay_VR, aes(color = FRP_MW), size = 2,  alpha = 0.5, na.rm = T, show.legend = F) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(color = "FRP [MW]") +
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, fill = "mediumorchid1", alpha = 0.8) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[2]], y = lat[[2]],
                           shape = location[[2]]), data = data_coord_lakes, colour = "turquoise1",
             size = 2, stroke = 1.5, show.legend = F) +
  xlab("") +
  ylab("") +
  labs(title = "Study area of Lake Satagay") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
        legend.title     = element_text(size = 20, vjust = 1),
        legend.text      = element_text(size = 20),
        axis.text        = element_text(size = 16),
        legend.key.width = unit(1, 'cm'))


source_PIH_K <- ggplot(Khamra_buf_extent) +
  geom_sf(data = dist_data_crop_Khamra_PIH_50, aes(colour  = "Locations of\np-source areas\n(measured PIH Q[0.5])"), show.legend = F) +
  scale_colour_manual(values = c("Locations of\np-source areas\n(measured PIH Q[0.5])" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Khamra_PIH_50, aes(color = FRP), size = 2,  alpha = 0.5, na.rm = T, show.legend = F) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(color = "FRP [MW]") +
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, fill = "mediumorchid1", alpha = 0.5) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[1]], y = lat[[1]],
                           shape = location[[1]]), data = data_coord_lakes, colour = "turquoise1",
             size = 3, stroke = 1.5, show.legend = F) +
  xlab("") +
  ylab("") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
        legend.title     = element_text(size = 20, vjust = 1),
        legend.text      = element_text(size = 20),
        axis.text        = element_text(size = 16),
        legend.key.width = unit(1, 'cm'))


source_PIH_S <- ggplot(Satagay_buf_extent) +
  geom_sf(data = dist_data_crop_Satagay_PIH_50, aes(colour  = "Locations of\np-source areas\n(measured PIH Q[0.5])"), show.legend = "point") +
  scale_colour_manual(values = c("Locations of\np-source areas\n(measured PIH Q[0.5])" = "darkslategrey"), name = NULL,
                      guide = guide_legend(override.aes = list(linetype = c("blank"),
                                                               shape = 16))) + # do not plot the color in legend
  new_scale("colour") +
  geom_sf(data = MODIS_Satagay_PIH_50, aes(color = FRP), size = 2,  alpha = 0.5, na.rm = T) +
  scale_color_continuous_sequential(palette = "Heat") +
  labs(color = "FRP [MW]") +
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.8) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA, alpha = 0.3))) +
  theme_minimal() +
  new_scale("colour") +
  geom_point(mapping = aes(x = lon[[2]], y = lat[[2]], col = "Lake"), data = data_coord_lakes, size = 3, stroke = 1.5) +
  scale_color_manual(values = c("Lake" = "turquoise1"),  name = NULL) +
  xlab("") +
  ylab("") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -3),
        legend.title     = element_text(size = 20, vjust = 1),
        legend.text      = element_text(size = 20),
        axis.text        = element_text(size = 16),
        legend.key.width = unit(1, 'cm'))




## Multiple plot
figure_1   <- ggarrange(source_H_K, source_H_S, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")

figure_2   <- ggarrange(source_PIH_K, source_PIH_S, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")



## Saving as png 
png(glue("Results/Simulation/p_source_areas_H.png"), width = 1400, height = 800)
plot(figure_1)
dev.off()

png(glue("Results/Simulation/p_source_areas_PIH.png"), width = 1400, height = 800)
plot(figure_2)
dev.off()



#####
## Multiple plot of potential source areas #####################################
#####
raster_H_K <- ggplot(Khamra_buf_extent) +
  geom_tile(data = df_r_K_na_VR , aes(x = x, y = y, fill = normalized), show.legend = F) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\np-source areas\n(calculated PIH)")+
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, fill = "mediumorchid1", alpha = 0.2) +
  new_scale("fill") +
  geom_sf(data = FRP_Khamra_VR, mapping = aes(fill = "firebrick2"), alpha = 0.2, colour = NA, show.legend = F) +
  new_scale("fill")+
  geom_sf(data = Khamra, fill = "turquoise1", alpha = 0.9) + 
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(title = "Study area Lake Khamra",shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = 5),
        legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 18))

raster_H_S <- ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_na_VR , aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99)), breaks = c(0,0.5,1.0)) +
  labs(fill = "Probability of\np-source areas\n(calculated PIH)")+
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.4) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill") +
  geom_sf(data = FRP_Satagay_VR, mapping = aes(fill = "Burned area"), alpha = 0.2, colour = NA) +
  scale_fill_manual(values = c("Burned area" = "firebrick2"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(title = "Study area Lake Satagay",shape = "Lake") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -3),
        legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 18))

raster_PIH_K <- ggplot(Khamra_buf_extent) +
  geom_tile(data = df_r_K_PIH_50_na , aes(x = x, y = y, fill = normalized), show.legend = F) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\np-source areas\n(measured PIH Q[0.5])")+
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, fill = "mediumorchid1", alpha = 0.2) +
  new_scale("fill") +
  geom_sf(data = FRP_Khamra_PIH_50, fill = "firebrick2", alpha = 0.2, colour = NA) +
  new_scale("fill")+
  geom_sf(data = Khamra, fill = "turquoise1", alpha = 0.9) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = 2),
        legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 18))

raster_PIH_S <- ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_PIH_50_na , aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99)), breaks = c(0,0.5,1.0)) +
  labs(fill = "Probability of\np-source areas\n(measured PIH Q[0.5])")+
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.4) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill") +
  geom_sf(data = FRP_Satagay_PIH_50, mapping = aes(fill = "Burned area"), alpha = 0.2, colour = NA) +
  scale_fill_manual(values = c("Burned area" = "firebrick2"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake"), alpha = 0.9) +
  scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(plot.title    = element_text(size = 26, hjust = 0.5, vjust = -2),
        legend.title  = element_text(size = 22, vjust = 0.5),
        legend.text   = element_text(size = 22, vjust = 0.75),
        axis.text     = element_text(size = 18))


## Multiple plot
figure_1   <- ggarrange(raster_H_K, raster_H_S, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")

figure_2   <- ggarrange(raster_PIH_K, raster_PIH_S, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend =T, legend = "bottom")



## Saving as png
png(glue("Results/Simulation/p_source_areas_H_raster.png"), width = 1100, height = 1200)
plot(figure_1)
dev.off()

png(glue("Results/Simulation/p_source_areas_PIH_raster.png"), width = 1100, height = 1200)
plot(figure_2)
dev.off()


######
## Actual charcoal source areas ################################################
######

raster_H_K <- ggplot(Khamra_buf_extent) +
  geom_tile(data = df_r_K_VR_mask, aes(x = x, y = y, fill = normalized), show.legend = F) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\nactual source areas\n(calculated PIH)")+
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, fill = "mediumorchid1", alpha = 0.2) +
  new_scale("fill")+
  geom_sf(data = Khamra, fill = "turquoise1") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(title = "Study area of Lake Khamra", shape = "Lake") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -4),
        legend.title     = element_text(size = 24, vjust = 1),
        legend.text      = element_text(size = 24),
        axis.text        = element_text(size = 18),
        legend.key.width = unit(1, 'cm'))


raster_H_S <- ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_VR_mask, aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99)), breaks = c(0,0.5,1.0)) +
  labs(fill = "Probability of\nactual source\nareas\n(calculated PIH)")+
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.2) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake")) +
  scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(title = "Study area of Lake Satagay",shape = "Lake") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = 1.8),
        legend.title     = element_text(size = 24, vjust = 1),
        legend.text      = element_text(size = 24),
        axis.text        = element_text(size = 18),
        legend.key.width = unit(1, 'cm'))


raster_PIH_K <- ggplot(Khamra_buf_extent) +
  geom_tile(data = df_r_K_mask_PIH_50, aes(x = x, y = y, fill = normalized), show.legend = F) +
  scale_fill_gradientn(colours = rev(viridis::mako(99))) +
  labs(fill = "Probability of\nactual source area\n(measured PIH Q[0.5])")+
  new_scale("fill") +
  geom_sf(data = Khamra_catchment, fill = "mediumorchid1", alpha = 0.2) +
  new_scale("fill")+
  geom_sf(data = Khamra, fill = "turquoise1") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -4),
        legend.title     = element_text(size = 24, vjust = 1),
        legend.text      = element_text(size = 24),
        axis.text        = element_text(size = 18),
        legend.key.width = unit(1, 'cm'))

raster_PIH_S <- ggplot(Satagay_buf_extent) +
  geom_tile(data = df_r_S_mask_PIH_50, aes(x = x, y = y, fill = normalized)) +
  scale_fill_gradientn(colours = rev(viridis::mako(99)), breaks = c(0,0.5,1.0)) +
  labs(fill = "Probability of\nactual source areas\n(measured PIH Q[0.5])")+
  new_scale("fill") +
  geom_sf(data = Satagay_catchment, aes(fill = "Catchment"), alpha = 0.2) +
  scale_fill_manual(values = c("Catchment" = "mediumorchid1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  
  new_scale("fill")+
  geom_sf(data = Satagay, aes(fill = "Lake")) +
  scale_fill_manual(values = c("Lake" = "turquoise1"), name = NULL,
                    guide = guide_legend(override.aes = list(linetype = "blank",  shape = NA))) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(shape = "Lake") +
  theme(plot.title       = element_text(size = 28, hjust = 0.5, vjust = -4),
        legend.title     = element_text(size = 24, vjust = 1),
        legend.text      = element_text(size = 24),
        axis.text        = element_text(size = 18),
        legend.key.width = unit(1, 'cm'))


## Multiple plot
figure_1   <- ggarrange(raster_H_K, raster_H_S, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend = T, legend = "bottom")

figure_2   <- ggarrange(raster_PIH_K, raster_PIH_S, nrow = 1, ncol = 2, 
                        widths = c(1,1),heights = c(1,1),
                        common.legend = T, legend = "bottom")


## Saving as png
png(glue("Results/Simulation/p_source_areas_H_raster_absolut.png"), width = 1200, height = 1200)
plot(figure_1)
dev.off()

png(glue("Results/Simulation/p_source_areas_PIH_raster_absolut.png"), width = 1200, height = 1200)
plot(figure_2)
dev.off()



#####
## Dispersal distance of charcoal ##############################################
#####

dist_data_Khamra_PIH_20     <- get(load("Results/Simulation/data/dist_data_Khamra_PIH_20.rda"))
dist_data_Satagay_PIH_20    <- get(load("Results/Simulation/data/dist_data_Satagay_PIH_20.rda"))
dist_data_Khamra_PIH_50     <- get(load("Results/Simulation/data/dist_data_Khamra_PIH_50.rda"))
dist_data_Satagay_PIH_50    <- get(load("Results/Simulation/data/dist_data_Satagay_PIH_50.rda"))
dist_data_Khamra_PIH_80     <- get(load("Results/Simulation/data/dist_data_Khamra_PIH_80.rda"))
dist_data_Satagay_PIH_80    <- get(load("Results/Simulation/data/dist_data_Satagay_PIH_80.rda"))
dist_data_Khamra_PIH_range  <- get(load("Results/Simulation/data/dist_data_Khamra_PIH_range.rda"))
dist_data_Satagay_PIH_range <- get(load("Results/Simulation/data/dist_data_Satagay_PIH_range.rda"))


#####
##  Lake Khamra ################################################################
#####

hist_dist_K_20 <- hist(dist_data_Khamra_PIH_20$dist[dist_data_Khamra_PIH_20$lake==1 & dist_data_Khamra_PIH_20$dist<400000]/1000, breaks = 100,
                       xlim = c(0,20), ylim = c(0, 6000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
hist_dist_K_50 <- hist(dist_data_Khamra_PIH_50$dist[dist_data_Khamra_PIH_50$lake==1 & dist_data_Khamra_PIH_50$dist<400000]/1000, breaks = 100,
                       xlim = c(0,20), ylim = c(0, 6000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
hist_dist_K_80 <- hist(dist_data_Khamra_PIH_80$dist[dist_data_Khamra_PIH_80$lake==1 & dist_data_Khamra_PIH_80$dist<400000]/1000, breaks = 100,
                       xlim = c(0,20), ylim = c(0, 6000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
hist_dist_K_range <- hist(dist_data_Khamra_PIH_range$dist[dist_data_Khamra_PIH_range$lake==1 & dist_data_Khamra_PIH_range$dist<400000]/1000, breaks = 100,
                          xlim = c(0,20), ylim = c(0, 6000), main = NULL,
                          xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")

legend("topright", legend=c("Measured PIH | Q[0.5]"), bty = "n", cex = 1.5)


png(glue("Results/Simulation/PIH/Khamra/boxplot_hist_travel_distances_Khamra_entire.png"), width = 1000, height = 800)
layout(matrix(c(1,2,3,4,5,5,5,5),nrow = 2, ncol = 4, byrow = TRUE))
boxplot(x = dist_data_Khamra_PIH_20$dist/1000, outline = F, 
        ylab = "Δx [km]", cex.main=1.5, cex.lab=1.8, cex.axis=1.8, 
        col = rgb(0,0,1,1/4), alpha = 0.1, ylim = c(0,10))
boxplot(x = dist_data_Khamra_PIH_50$dist/1000, outline = F, 
        cex.main=1.5, cex.lab=1.8, cex.axis=1.8, 
        col = "cadetblue1", ylim = c(0,10))
boxplot(x = dist_data_Khamra_PIH_80$dist/1000, outline = F, 
        cex.main=1.5, cex.lab=1.8, cex.axis=1.8, 
        col = "plum2", ylim = c(0,10))
boxplot(x = dist_data_Khamra_PIH_range$dist/1000, outline = F, 
        ylim = c(0,10),cex.main=1.2, cex.lab=1.8, cex.axis=1.8,
        col = "khaki")


## hist
plot(hist_dist_K_20, col = rgb(0,0,1,1/4), main = NULL, xlab = "Δx [km]", cex.lab = 1.7, xlim = c(0, 20), ylim = c(0,6000))
plot(hist_dist_K_50, col = "cadetblue1", add = T, main = NULL, xlab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5,  xlim = c(0, 20), ylim = c(0,6000))
plot(hist_dist_K_80, col = "plum2", add = T, main = NULL, xlab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5,  xlim = c(0, 20), ylim = c(0,6000))
plot(hist_dist_K_range, col = "khaki", add= T, main = NULL, xlab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5,  xlim = c(0, 20), ylim = c(0,6000))

legend("topright", legend = c("Q[0.2]","Q[0.5]", "Q[0.8]","Range between\nQ[0.2] and Q[0.8]"), 
       col=c(rgb(0,0,1,1/4), "cadetblue1", "plum2", "khaki"), 
       pt.cex=3, pch=15, cex = 1.8, title = "Measured PIH [m]", bty = "n")

mtext("Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density) into the atmosphere of study area Lake Khamra",
      line = -3.5 ,outer = TRUE, cex = 1.4, font = 2, adj = 0.6)
dev.off()


#####
##  Lake Satagay ###############################################################
#####

hist_dist_S_20 <- hist(dist_data_Satagay_PIH_20$dist[dist_data_Satagay_PIH_20$lake==2 & dist_data_Satagay_PIH_20$dist<400000]/1000, breaks = 50,
                       xlim = c(0,20), ylim = c(0, 3000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
hist_dist_S_50 <- hist(dist_data_Satagay_PIH_50$dist[dist_data_Satagay_PIH_50$lake==2 & dist_data_Satagay_PIH_50$dist<400000]/1000, breaks = 100,
                       xlim = c(0,20), ylim = c(0,3000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
hist_dist_S_80 <- hist(dist_data_Satagay_PIH_80$dist[dist_data_Satagay_PIH_80$lake==2 & dist_data_Satagay_PIH_80$dist<400000]/1000, breaks = 200,
                       xlim = c(0,20), ylim = c(0, 3000), main = NULL,
                       xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")
hist_dist_S_range <- hist(dist_data_Satagay_PIH_range$dist[dist_data_Satagay_PIH_range$lake==2 & dist_data_Satagay_PIH_range$dist<400000]/1000, breaks = 200,
                          xlim = c(0,20), ylim = c(0, 3000), main = NULL,
                          xlab = "Δx [km]", cex.main=1.7, cex.lab=1.5, cex.axis=1, col = "lightblue3")

legend("topright", legend=c("Measured PIH | Q[0.5]"), bty = "n", cex = 1.5)


png(glue("Results/Simulation/PIH/Satagay/boxplot_hist_travel_distances_K_S_entire.png"), width = 1000, height = 800)
layout(matrix(c(1,2,3,4,5,5,5,5),nrow = 2, ncol = 4, byrow = TRUE))
boxplot(x = dist_data_Satagay_PIH_20$dist/1000, outline = F, 
        ylab = "Δx [km]", cex.main=1.5, cex.lab=1.8, cex.axis=1.8, 
        col = rgb(0,0,1,1/4), alpha = 0.1, ylim = c(0,12))
boxplot(x = dist_data_Satagay_PIH_50$dist/1000, outline = F, 
        cex.main=1.5, cex.lab=1.8, cex.axis=1.8, 
        col = "cadetblue1", ylim = c(0,12))
boxplot(x = dist_data_Satagay_PIH_80$dist/1000, outline = F, 
        cex.main=1.5, cex.lab=1.8, cex.axis=1.8, 
        col = "plum2", ylim = c(0,12))
boxplot(x = dist_data_Satagay_PIH_range$dist/1000, outline = F, 
        ylim = c(0,12),cex.main=1.2, cex.lab=1.8, cex.axis=1.8,
        col = "khaki")

plot(hist_dist_S_20, col = rgb(0,0,1,1/4), main = NULL, xlab = "Δx [km]", cex.lab = 1.7, xlim = c(0, 20), ylim = c(0,6000))
plot(hist_dist_S_50, col = "cadetblue1", add = T, main = NULL, xlab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5,  xlim = c(0, 20), ylim = c(0,6000))
plot(hist_dist_S_80, col = "plum2", add = T, main = NULL, xlab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5,  xlim = c(0, 20), ylim = c(0,6000))
plot(hist_dist_S_range, col = "khaki", add= T, main = NULL, xlab = "Δx [km]", cex.lab = 1.7, cex.axis = 1.5,  xlim = c(0, 20), ylim = c(0,6000))


legend("topright", legend = c("Q[0.2]","Q[0.5]", "Q[0.8]","Range between\nQ[0.2] and Q[0.8]"), 
       col=c(rgb(0,0,1,1/4), "cadetblue1", "plum2", "khaki"), 
       pt.cex=3, pch=15, cex = 1.8, title = "Measured PIH [m]", bty = "n")

mtext("Travel distance of spherical particle (150-500 μm | 0.5 g/cm³ density) into the atmosphere of study area Lake Satagay",
      line = -3.5 ,outer = TRUE, cex = 1.4, font = 2, adj = 0.6)
dev.off()



