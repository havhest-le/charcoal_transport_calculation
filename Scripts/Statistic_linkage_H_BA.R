#####
## Statistical analysis of the linkage between FRP and wind speed, as a ########
## function of calculated PIH based on Vachula & Richter (2018) ################
#####

rm(list = ls(all= TRUE))

#####
## Load packages ###############################################################
#####

library(tidyverse)
library(sf)
library(terra)
library(stars)
library(dplyr)
library(glue)

## Plotting
library(maps)
library(viridis)
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(ggnewscale)
library(MASS)
library(survival)
library(fitdistrplus)

#####
## Statistics ################################################################## 
#####

#####
## Lake Khamra #################################################################
######
## Load and clean data #########################################################
H_function_Khamra <- get(load("Results/Simulation/data/H_function_Khamra.rda"))
H_function_Khamra$FRP[H_function_Khamra$FRP == 0] <- NA
H_function_Khamra$H[H_function_Khamra$H == 0]     <- NA
H_function_Khamra                                 <- na.omit(H_function_Khamra)
H_function_Khamra$FRP_MW                          <- H_function_Khamra$FRP/238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

## Pearson correlation
cor(H_function_Khamra$FRP_MW, H_function_Khamra$H,  method = "pearson", use = "complete.obs")
## The correlation between calculated PIH (H) and FRP is 64,45 %

#####
## Linear regression ###########################################################
#####

lm_FRP_H_Khamra_MW  <- lm(H ~ FRP_MW, data = H_function_Khamra)
coeff_Khamra_MW     <- coefficients(lm_FRP_H_Khamra_MW)
## (Intercept)             FRP 
## 742.942118         8.046233 
eq_K_MW <- paste0("y = ", round(coeff_Khamra_MW[1],2), " + ", round(coeff_Khamra_MW[2], 2), "*x")
## [1] "y = 8.05*x742.94"

## Plotting ####################################################################

png(glue("Results/Simulation/H/Khamra/lm/lm_FRP_H_Khamra_MW.png"), width = 1000, height = 800)
plot(x = H_function_Khamra$FRP_MW, y = H_function_Khamra$H,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Khamra | {eq_K_MW}"), cex.main = 1.5)
title(sub = "R² = 0.4154 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
abline(lm_FRP_H_Khamra_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/H/Khamra/lm/lm_FRP_H_Khamra_1000_MW.png"), width = 1000, height = 800)
plot(x = H_function_Khamra$FRP_MW, y = H_function_Khamra$H, 
     xlim = c(0, 300), ylim = c(0,1000), 
     xlab = "FRP [MW]", 
     ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Khamra | {eq_K_MW}"), cex.main = 1.5)
title(sub = "R² = 0.4154 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
abline(lm_FRP_H_Khamra_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/H/Khamra/lm/lm_FRP_H_Khamra_4000_MW.png"), width = 1000, height = 800)
plot(x = H_function_Khamra$FRP_MW, y = H_function_Khamra$H, 
     xlim = c(0, 300), ylim = c(0,4000), 
     xlab = "FRP [MW]", 
     ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Khamra | {eq_K_MW}"), cex.main = 1.5)
title(sub = "R² = 0.4154 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
abline(lm_FRP_H_Khamra_MW, col = "red", lwd = 2)
dev.off()

## Summary of the data #########################################################

summary(lm_FRP_H_Khamra_MW)
## There is a significant positive relationship between calculated PIH and FRP
## F-statistic:  3418 on 1 and 4810 DF,  p-value: < 2.2e-16
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.4154 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the variable $FRP$ increases by 1, $\widehat{PIH}$ increases by 8.05 units.
## (high significant, <2e-16 ***)


#####
## Multiple regression #########################################################
#####

multi_Khamra_H_MW    <- lm(H ~ FRP_MW + w_spd, data = H_function_Khamra)
mul_coeff_Khamra_MW  <- coefficients(multi_Khamra_H_MW)
## (Intercept)      FRP_MW       w_spd 
## 1719.308564    8.078606 -215.516521 

## Summary of data #############################################################

summary(multi_Khamra_H_MW)
## F-statistic:  2982 on 2 and 4809 DF,  p-value: < 2.2e-16
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.5536 (Variance clarification)
## If the FRP increases by 1, calculated PIH increases by 8.0786 units 
## (high significant, <2e-16 ***)
## If the wind speed increases by 1, calculated PIH decreases by 215.5165 units
## (high significant, <2e-16 ***)
eq_multi_H_K <- paste0("y = ", round(mul_coeff_Khamra_MW[1],2), " + ", round(mul_coeff_Khamra_MW[2], 2), "*x1 ", round(mul_coeff_Khamra_MW[3], 2), "*x2")
## "y = 1719.31 + 8.08*x1 -215.52*x2"


## Plotting ####################################################################

library(scatterplot3d)
png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP_MW, z = H_function_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), color = "grey40", 
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.6)
dev.off()


png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra_points_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP_MW, z = H_function_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), color = "grey40",  
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.6)
dev.off()

png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra_4000_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP_MW, z = H_function_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 4000), color = "grey40",  
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.6)
dev.off()

png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra_points_4000_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP_MW, z = H_function_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 4000), color = "grey40", 
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.6)
dev.off()


# #####
# ## Calculation with FRP [cal/s] ################################################
# #####
# 
# ## Correlation
# cor(H_function_Khamra$FRP, H_function_Khamra$H,  method = "pearson", use = "complete.obs")
# ## The correlation between calculated PIH (H) and FRP is 64,45 % 
# 
# #####
# ## Linear regression ###########################################################
# #####
# 
# lm_FRP_H_Khamra  <- lm(H ~ FRP, data = H_function_Khamra)
# coeff_Khamra     <- coefficients(lm_FRP_H_Khamra)
# ## (Intercept)             FRP 
# ## 7.429421e+02   3.368797e-05 
# eq_K <- paste0("y = ", round(coeff_Khamra[2], 5), "*x", round(coeff_Khamra[1],1))
# ## [1] "y = 3e-05*x742.9"
# 
# ## Plotting ###################################################################
# 
# png(glue("Results/Simulation/H/Khamra/lm/lm_FRP_H_Khamra.png"), width = 1000, height = 800)
# plot(x = H_function_Khamra$FRP, y = H_function_Khamra$H,
#      xlim = c(0, 71653768.9882488), ylim = c(0,10000),
#      xlab = "FRP [cal/s]",
#      ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
# title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Khamra | {eq_K}"), cex.main = 1.5)
# title(sub = "R² = 0.4154 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
# abline(lm_FRP_H_Khamra, col = "red", lwd = 2)
# dev.off()
# 
# png(glue("Results/Simulation/H/Khamra/lm/lm_FRP_H_Khamra_4000.png"), width = 1000, height = 800)
# plot(x = H_function_Khamra$FRP, y = H_function_Khamra$H, 
#      xlim = c(0, 71653768.9882488), ylim = c(0,4000), 
#      xlab = "FRP [cal/s]", 
#      ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
# title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Khamra | {eq_K}"), cex.main = 1.5)
# title(sub = "R² = 0.4154 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
# abline(lm_FRP_H_Khamra, col = "red", lwd = 2)
# dev.off()
# 
# png(glue("Results/Simulation/H/Khamra/lm/lm_FRP_H_Khamra_2000.png"), width = 1000, height = 800)
# plot(x = H_function_Khamra$FRP, y = H_function_Khamra$H, 
#      xlim = c(0, 71653768.9882488), ylim = c(0,2000), 
#      xlab = "FRP [cal/s]", 
#      ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
# title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Khamra | {eq_K}"), cex.main = 1.5)
# title(sub = "R² = 0.4154 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
# abline(lm_FRP_H_Khamra, col = "red", lwd = 2)
# dev.off()
# 
# ## Summary of the data #########################################################
# 
# summary(lm_FRP_H_Khamra)
# ## There is a significant positive relationship between calculated PIH and FRP
# ## F-statistic:  3418 on 1 and 4810 DF,  p-value: < 2.2e-16
# ## Null hypothesis rejected, because p-value is < 0.05
# ## Multiple R-squared:  0.4154 (Variance clarification) is not a good or stable 
# ## value. 1 would be perfect to explain the model with the variance of the 
# ## dependent variable. 
# ## If the FRP is increasing by 1 = calculated PIH is increasing by 3.369e-05 units 
# ## (high significant, <2e-16 ***)
# 
# #####
# ## Multiple regression #########################################################
# #####
# 
# multi_Khamra_H    <- lm(H ~ FRP + w_spd, data = H_function_Khamra)
# mul_coeff_Khamra  <- coefficients(multi_Khamra_H)
# ## (Intercept)            FRP         w_spd 
# ## 1.719309e+03  3.382351e-05 -2.155165e+02 
# 
# ## Summary of data #############################################################
# 
# summary(multi_Khamra_H)
# ## F-statistic:  2982 on 2 and 4809 DF,  p-value: < 2.2e-16 
# ## Null hypothesis rejected, because p-value is < 0.05
# ## Multiple R-squared:  0.5536 (Variance clarification)
# ## If the FRP is increasing by 1 = calculated PIH is increasing by 3.382e-05 units 
# ## (high significant, <2e-16 ***)
# ## If the wind speed is increasing by 1 = calculated PIH is decreasing by 2.155e+02 units
# ## (high significant, <2e-16 ***)
# 
# ## Plotting ####################################################################
# 
# library(scatterplot3d)
# png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 4000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
#       cex.main = 1.4)
# dev.off()
# 
# png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra_points.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 4000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
#       cex.main = 1.4)
# dev.off()
# 
# png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra_2000.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 2000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
#       cex.main = 1.4)
# dev.off()
# 
# png(glue("Results/Simulation/H/Khamra/mul/Multiple_reg_Khamra_points_2000.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 2000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5536 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
#       cex.main = 1.4)
# dev.off()


## Which variable has the greatest impact on the plume injection height? 
H_function_Khamra$H_Z     <- scale(H_function_Khamra$H)
#H_function_Khamra$FRP_Z   <- scale(H_function_Khamra$FRP)
H_function_Khamra$w_spd_Z <- scale(H_function_Khamra$w_spd)
H_function_Khamra$FRP_Z_MW   <- scale(H_function_Khamra$FRP_MW)

# # For FRP [cal/s]
# Z_model_H_Khamra <- lm(H_Z ~ FRP_Z + w_spd_Z, data = H_function_Khamra)
# summary(Z_model_H_Khamra) ## Comparable values
# ## F-statistic:  2982 on 2 and 4809 DF,  p-value: < 2.2e-16 (identical values)
# ## Multiple R-squared:  0.5536 (identical values)
# ## The fire intensity has a greater impact on the height of the smoke plume during
# ## a fire event (value is greater)

## For FRP [MW]
Z_model_H_Khamra_MW <- lm(H_Z ~ FRP_Z_MW + w_spd_Z, data = H_function_Khamra)
summary(Z_model_H_Khamra_MW) ## Comparable values
## F-statistic:  2982 on 2 and 4809 DF,  p-value: < 2.2e-16 (identical values)
## Multiple R-squared:  0.5536 (identical values)
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)


#####
## Lake Satagay ################################################################
#####

## Load and clean data #########################################################
H_function_Satagay <- get(load("Results/Simulation/data/H_function_Satagay.rda"))
H_function_Satagay$FRP[H_function_Satagay$FRP == 0] <- NA
H_function_Satagay$H[H_function_Satagay$H == 0]     <- NA
H_function_Satagay                                  <- na.omit(H_function_Satagay)
H_function_Satagay$FRP_MW                           <- H_function_Satagay$FRP/238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

## Correlation
cor(H_function_Satagay$H, H_function_Satagay$FRP_MW,  method = "pearson", use = "complete.obs")
## The correlation between calculated PIH (H) and FRP is 53,68 % 

#####
## Linear regression ###########################################################
#####

lm_FRP_H_Satagay_MW  <- lm(H ~ FRP_MW, data = H_function_Satagay)
coeff_satagay_MW     <- coefficients(lm_FRP_H_Satagay_MW)
## (Intercept)        FRP 
## 551.916070    3.627762 
eq_S_MW <- paste0("y = ", round(coeff_satagay_MW[1],2), " + ", round(coeff_satagay_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/H/Satagay/lm/lm_FRP_H_Satagay_MW.png"), width = 1000, height = 800)
plot(x = H_function_Satagay$FRP_MW, y = H_function_Satagay$H,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_MW}"), cex.main = 1.5)
title(sub = "R² = 0.2882 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
abline(lm_FRP_H_Satagay_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/H/Satagay/lm/lm_FRP_H_Satagay_4000_MW.png"), width = 1000, height = 800)
plot(x = H_function_Satagay$FRP_MW, y = H_function_Satagay$H,
     xlim = c(0, 300), ylim = c(0,4000),
     xlab = "FRP [MW]",
     ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_MW}"), cex.main = 1.5)
title(sub = "R² = 0.2882 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
abline(lm_FRP_H_Satagay_MW, col = "red", lwd = 2)
dev.off()


## Summary of data #############################################################

summary(lm_FRP_H_Satagay_MW)
## There is a significant positive relationship between calculated PIH and FRP
## F-statistic:  8392 on 1 and 20730 DF,  p-value: < 2.2e-16
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.2882 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, calculated PIH increases by 3.6278 units 
## (high significant, <2e-16 ***)

#####
## Multiple regression #########################################################
#####

multi_Satagay_H_MW    <- lm(H ~ FRP_MW + w_spd, data = H_function_Satagay)
mul_coeff_Satagay_MW  <- coefficients(multi_Satagay_H_MW)
## (Intercept)      FRP_MW       w_spd 
## 1736.358290    3.783258 -171.674485 

## Summary of data #############################################################

summary(multi_Satagay_H_MW)
## F-statistic: 1.061e+04 on 2 and 20729 DF,  p-value: < 2.2e-16
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.5059 (Variance clarification)
## If the FRP increases by 1, calculated PIH increases by 3.78326 units 
## (high significant, <2e-16 ***)
## If the wind speed increases by 1, calculated PIH decreases by 171.67449 units
## (high significant, <2e-16 ***)
eq_multi_H_S <- paste0("y = ", round(mul_coeff_Satagay_MW[1],2), " + ", round(mul_coeff_Satagay_MW[2], 2), "*x1 ", round(mul_coeff_Satagay_MW[3], 2), "*x2")
## "y = 1736.36 + 3.78*x1 -171.67*x2"

## Plotting ####################################################################

png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Satagay$H, x = H_function_Satagay$FRP_MW, z = H_function_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), color = "grey40", 
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.6)
dev.off()

png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay_points_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Satagay$H, x = H_function_Satagay$FRP_MW, z = H_function_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), color = "grey40", 
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.6)
dev.off()

png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay_4000_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Satagay$H, x = H_function_Satagay$FRP_MW, z = H_function_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 4000), color = "grey40",  
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.6)
dev.off()

png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay_points_4000_MW.png"), width = 1000, height = 800)
scatterplot3d(y = H_function_Satagay$H, x = H_function_Satagay$FRP_MW, z = H_function_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 4000), color = "grey40", 
              xlab = "FRP [MW]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.6)
dev.off()


# #####
# ## Calculation with FRP [cal/s] #################################################
# #####
# 
# ## Correlation
# cor(H_function_Satagay$H, H_function_Satagay$FRP,  method = "pearson", use = "complete.obs")
# ## The correlation between calculated PIH (H) and FRP is 53,68 % 
# 
# #####
# ## Linear regression ###########################################################
# #####
# 
# lm_FRP_H_Satagay  <- lm(H ~ FRP, data = H_function_Satagay)
# coeff_satagay     <- coefficients(lm_FRP_H_Satagay)
# ## (Intercept)           FRP 
# ## 5.519161e+02 1.518871e-05 
# eq_S <- paste0("y = ", round(coeff_satagay[2], 5), "*x", round(coeff_satagay[1],2))  
# ## [1] "y = 0*x551.92"
# 
# ## Plotting ####################################################################
# 
# png(glue("Results/Simulation/H/Satagay/lm/lm_FRP_H_Satagay.png"), width = 1000, height = 800)
# plot(x = H_function_Satagay$FRP, y = H_function_Satagay$H,
#      xlim = c(0, 71653768.9882488), ylim = c(0,10000),
#      xlab = "FRP [cal/s]",
#      ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
# title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S}"), cex.main = 1.5)
# title(sub = "R² = 0.2882 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
# abline(lm_FRP_H_Satagay, col = "red", lwd = 2)
# dev.off()
# 
# png(glue("Results/Simulation/H/Satagay/lm/lm_FRP_H_Satagay_4000.png"), width = 1000, height = 800)
# plot(x = H_function_Satagay$FRP, y = H_function_Satagay$H,
#      xlim = c(0, 71653768.9882488), ylim = c(0,4000),
#      xlab = "FRP [cal/s]",
#      ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
# title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S}"), cex.main = 1.5)
# title(sub = "R² = 0.2882 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
# abline(lm_FRP_H_Satagay, col = "red", lwd = 2)
# dev.off()
# 
# png(glue("Results/Simulation/H/Satagay/lm/lm_FRP_H_Satagay_2000.png"), width = 1000, height = 800)
# plot(x = H_function_Satagay$FRP, y = H_function_Satagay$H,
#      xlim = c(0, 71653768.9882488), ylim = c(0,2000),
#      xlab = "FRP [cal/s]",
#      ylab = "Calculated PIH [m]", col = "grey40", cex.lab = 1.5)
# title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S}"), cex.main = 1.5)
# title(sub = "R² = 0.2882 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5)
# abline(lm_FRP_H_Satagay, col = "red", lwd = 2)
# dev.off()
# 
# ## Summary of data #############################################################
# 
# summary(lm_FRP_H_Satagay)
# ## There is a significant positive relationship between calculated PIH and FRP
# ## F-statistic:  8392 on 1 and 20730 DF,  p-value: < 2.2e-16
# ## Null hypothesis rejected, because p-value is < 0.05
# ## Multiple R-squared:  0.2882 (Variance clarification) is not a good or stable 
# ## value. 1 would be perfect to explain the model with the variance of the 
# ## dependent variable. 
# ## If the FRP is increasing by 1 = calculated PIH is increasing by 1.519e-05 units 
# ## (high significant, <2e-16 ***)
# 
# 
# #####
# ## Multiple regression #########################################################
# #####
# 
# multi_Satagay_H    <- lm(H ~ FRP + w_spd, data = H_function_Satagay)
# mul_coeff_Satagay  <- coefficients(multi_Satagay_H)
# ## (Intercept)         FRP_MW         w_spd 
# ## 1.736358e+03  1.583975e-05 -1.716745e+02  
# 
# ## Summary of data #############################################################
# 
# summary(multi_Satagay_H)
# ## F-statistic: 1.061e+04 on 2 and 20729 DF,  p-value: < 2.2e-16
# ## Null hypothesis rejected, because p-value is < 0.05
# ## Multiple R-squared:  0.5059 (Variance clarification)
# ## If the FRP is increasing by 1 = calculated PIH is increasing by 1.584e-05 units 
# ## (high significant, <2e-16 ***)
# ## If the wind speed is increasing by 1 = calculated PIH is decreasing by 1.717e+02 units
# ## (high significant, <2e-16 ***)
# 
# 
# ## Plotting ####################################################################
# 
# library(scatterplot3d)
# 
# png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 4000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
#       cex.main = 1.4)
# dev.off()
# 
# 
# png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay_points.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 4000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
#       cex.main = 1.4)
# dev.off()
# 
# png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay_2000.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 2000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
#       cex.main = 1.4)
# dev.off()
# 
# png(glue("Results/Simulation/H/Satagay/mul/Multiple_reg_Satagay_points_2000.png"), width = 1000, height = 800)
# scatterplot3d(y = H_function_Khamra$H, x = H_function_Khamra$FRP, z = H_function_Khamra$w_spd, 
#               xlim = c(0, 71653768.9882488), ylim = c(0, 2000), color = "grey40", 
#               xlab = "FRP [cal/s]", ylab = "Calculated PIH [m]", zlab = "spd [m/s]", 
#               type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
# title(sub = "R² = 0.5059 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
#       main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
#       cex.main = 1.4)
# dev.off()
# 

## Which variable has the greatest impact on the plume injection height? 
H_function_Satagay$H_Z      <- scale(H_function_Satagay$H)
H_function_Satagay$FRP_Z    <- scale(H_function_Satagay$FRP)
H_function_Satagay$w_spd_Z  <- scale(H_function_Satagay$w_spd)
H_function_Satagay$FRP_Z_MW <- scale(H_function_Satagay$FRP_MW)

# ## For FRP [cal/s]
# Z_model_H_Satagay <- lm(H_Z ~ FRP_Z + w_spd_Z, data = H_function_Satagay)
# summary(Z_model_H_Satagay) ## Comparable values
# ## F-statistic: 1.061e+04 on 2 and 20729 DF,  p-value: < 2.2e-16 (identical values)
# ## Multiple R-squared:  0.5059 (identical values)
# ## The fire intensity has a greater impact on the height of the smoke plume during 
# ## a fire event (value is greater)

## For FRP [MW]
Z_model_H_Satagay_MW <- lm(H_Z ~ FRP_Z_MW + w_spd_Z, data = H_function_Satagay)
summary(Z_model_H_Satagay_MW) ## Comparable values
## F-statistic: 1.061e+04 on 2 and 20729 DF,  p-value: < 2.2e-16 (identical values)
## Multiple R-squared:  0.5059 (identical values)
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)



