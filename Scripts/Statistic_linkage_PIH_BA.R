#####
## Statistical anaylsis of the linkage between FRP and wind speed, as a ########
## function of measured PIH by MODIS ###########################################
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
library(scatterplot3d)

#####
## Statistics ################################################################## 
#####

#####
## Lake Khamra #################################################################
######
## Load and clean data #########################################################
tibble_frp_pih_wind_Khamra         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Khamra.rda"))
## Convert MW to cal/s
tibble_frp_pih_wind_Khamra$FRP_cal <- tibble_frp_pih_wind_Khamra$FRP*238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

## Correlation
cor_025 <- cor(tibble_frp_pih_wind_Khamra$FRP, tibble_frp_pih_wind_Khamra$`2.5%`,  method = "pearson", use = "complete.obs")
## The correlation between measured PIH (H) and FRP is -1.43 % 
cor_20 <- cor(tibble_frp_pih_wind_Khamra$FRP, tibble_frp_pih_wind_Khamra$`20%`,  method = "pearson", use = "complete.obs")
## The correlation between measured PIH (H) and FRP is 2,22 % 
cor_50 <- cor(tibble_frp_pih_wind_Khamra$FRP, tibble_frp_pih_wind_Khamra$`50%`,  method = "pearson", use = "complete.obs")
## The correlation between measured PIH (H) and FRP is 6,51 % 
cor_80 <- cor(tibble_frp_pih_wind_Khamra$FRP, tibble_frp_pih_wind_Khamra$`80%`,  method = "pearson", use = "complete.obs")
## The correlation between measured PIH (H) and FRP is 7,19 % 
cor_97 <- cor(tibble_frp_pih_wind_Khamra$FRP, tibble_frp_pih_wind_Khamra$`97.5%`,  method = "pearson", use = "complete.obs")
## The correlation between measured PIH (H) and FRP is 6,95 % 


## Q[0.025] ######################################################################

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Khamra_MW  <- lm(`2.5%` ~ FRP, data = tibble_frp_pih_wind_Khamra)
coeff_Khamra_PIH_MW   <- coefficients(lm_FRP_PIH_Khamra_MW)
## (Intercept)             FRP 
## 51.70744473      -0.01174436
eq_K_PIH_MW <- paste0("y = ",round(coeff_Khamra_PIH_MW[1],2), " + ", round(coeff_Khamra_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Khamra/2.5/lm/lm_FRP_PIH_Khamra_MW_2_5.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`2.5%`,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.025]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.0002248 | p-value = 0.3841", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/2.5/lm/lm_FRP_PIH_Khamra_1000_MW_2_5.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`80%`,
     xlim = c(0, 300), ylim = c(0,1000), 
     xlab = "FRP [MW]", 
     ylab = "Measured PIH [m] | Q[0.025]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.0002248 | p-value = 0.3841", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()


## Summary of the data #########################################################

summary(lm_FRP_PIH_Khamra_MW)
## There is a significant positive relationship between measured PIH and FRP
## F-statistic: 0.7577 on 1 and 3370 DF,  p-value: 0.3841
## Multiple R-squared:  0.0002248 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increase by 1, measured PIH decrease by 0.01174 units 
## (0.384)


#####
## Multiple regression #########################################################
#####

multi_Khamra_PIH_MW      <- lm(`2.5%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Khamra)
mul_coeff_Khamra_PIH_MW  <- coefficients(multi_Khamra_PIH_MW)
## (Intercept)        FRP_MW        w_spd 
## 65.53608625   -0.01052091   -3.31264305 

## Summary of data #############################################################

summary(multi_Khamra_PIH_MW)
## F-statistic: 12.28 on 2 and 3369 DF,  p-value: 4.864e-06
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.007236 (Variance clarification)
## If the FRP increases by 1 = measured PIH decreases by 0.01052 units (0.434)
## If the wind speed increases by 1, measured PIH decreases by 3.31264 units
## (high significant, 1.12e-06 ***)

## Plotting ####################################################################
png(glue("Results/Simulation/PIH/Khamra/2.5/mul/Multiple_reg_Khamra_PIH_MW_2_5.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`2.5%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20),color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.025]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.007236 | p-value = 4.864e-06", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/2.5/mul/Multiple_reg_Khamra_PIH_MW_points_2_5.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`2.5%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.025]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.007236 | p-value = 4.864e-06", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Khamra$PIH_Z   <- scale(tibble_frp_pih_wind_Khamra$`2.5%`)
tibble_frp_pih_wind_Khamra$FRP_Z   <- scale(tibble_frp_pih_wind_Khamra$FRP)
tibble_frp_pih_wind_Khamra$w_spd_Z <- scale(tibble_frp_pih_wind_Khamra$w_spd)

## For FRP [MW]
Z_model_PIH_Khamra_MW <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Khamra)
summary(Z_model_PIH_Khamra_MW) ## Comparable values
## F-statistic: 12.28 on 2 and 3369 DF,  p-value: 4.864e-06
## Multiple R-squared:  0.007236 (identical values)
## FRP_Z       -1.343e-02  1.717e-02  -0.782    0.434    
## w_spd_Z     -8.375e-02  1.717e-02  -4.878 1.12e-06 ***

## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)???


#####
## Lake Satagay ################################################################
#####

## Load and clean data #########################################################
tibble_frp_pih_wind_Satagay         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))

## Pearson correlation
cor_50 <- cor(tibble_frp_pih_wind_Satagay$FRP, tibble_frp_pih_wind_Satagay$`50%`,  method = "pearson", use = "complete.obs")
## The correlation between measured PIH (H) and FRP is -2.8 % 


## Convert MW to cal/s
tibble_frp_pih_wind_Satagay$FRP_cal <- tibble_frp_pih_wind_Satagay$FRP*238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Satagay_MW  <- lm(`2.5%` ~ FRP, data = tibble_frp_pih_wind_Satagay)
coeff_satagay_PIH_MW   <- coefficients(lm_FRP_PIH_Satagay_MW)
## (Intercept)             FRP 
## 51.77628967      -0.02762737
eq_S_PIH_MW <- paste0("y = ",round(coeff_satagay_PIH_MW[1],2), " + ", round(coeff_satagay_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/2.5/lm/lm_FRP_PIH_Satagay_MW_2_5.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`2.5%`,
     xlim = c(0, 300), ylim = c(0, 2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.025]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.003172 | p-value = 2.224e-13", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()


png(glue("Results/Simulation/PIH/Satagay/2.5/lm/lm_FRP_PIH_Satagay_MW_1000_2_5.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`2.5%`,
     xlim = c(0, 300), ylim = c(0, 1000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.025]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.003172 | p-value = 2.224e-13", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()

## Summary of data #############################################################

summary(lm_FRP_PIH_Satagay_MW)
## There is a significant positive relationship between calculated PIH and FRP
## F-statistic: 53.88 on 1 and 16936 DF,  p-value: 2.224e-13
## Multiple R-squared:  0.003172 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, calculated PIH decreases by 0.027627 units 
## (high significant, 2.22e-13 ***)

#####
## Multiple regression #########################################################
#####

multi_Satagay_PIH_MW    <- lm(`2.5%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Satagay)
mul_coeff_Satagay_PIH_MW  <- coefficients(multi_Satagay_PIH_MW)
## (Intercept)          FRP_MW         w_spd 
## 60.1926751       -0.0265382    -1.2450753   

## Summary of data #############################################################

summary(multi_Satagay_PIH_MW)
## F-statistic: 43.23 on 2 and 16935 DF,  p-value: < 2.2e-16
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.00508 (Variance clarification)
## If the FRP increases by 1, calculated PIH decreases by 0.026538   units 
## (high significant,  1.88e-12 ***)
## If the wind speed increases by 1, calculated PIH decreases by 1.245075 units
## (high significant,  1.22e-08 ***)

## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/2.5/mul/Multiple_reg_Satagay_MW_2_5.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`2.5%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.025]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.00508 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Satagay/2.5/mul/Multiple_reg_Satagay_points_MW_2_5.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`2.5%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.025]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.00508 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Satagay$PIH_Z    <- scale(tibble_frp_pih_wind_Satagay$`2.5%`)
tibble_frp_pih_wind_Satagay$FRP_Z    <- scale(tibble_frp_pih_wind_Satagay$FRP)
tibble_frp_pih_wind_Satagay$w_spd_Z  <- scale(tibble_frp_pih_wind_Satagay$w_spd)


Z_model_PIH_Satagay <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Satagay)
summary(Z_model_PIH_Satagay) ## Comparable values
## F-statistic: 43.23 on 2 and 16935 DF,  p-value: < 2.2e-16 (identical values)
## Multiple R-squared:  0.00508 (identical values)

## ??
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)



## Q[0.2] ######################################################################

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Khamra_MW  <- lm(`20%` ~ FRP, data = tibble_frp_pih_wind_Khamra)
coeff_Khamra_PIH_MW   <- coefficients(lm_FRP_PIH_Khamra_MW)
## (Intercept)             FRP 
## 163.35851539     0.02467838 
eq_K_PIH_MW <- paste0("y = ",round(coeff_Khamra_PIH_MW[1],2), " + ", round(coeff_Khamra_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Khamra/20/lm/lm_FRP_PIH_Khamra_MW_20.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`20%`,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.2]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.0004949 | p-value = 0.1965", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/20/lm/lm_FRP_PIH_Khamra_1000_MW_20.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`20%`,
     xlim = c(0, 300), ylim = c(0,1000), 
     xlab = "FRP [MW]", 
     ylab = "Measured PIH [m] | Q[0.2]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.0004949 | p-value = 0.1965", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()


## Summary of the data #########################################################

summary(lm_FRP_PIH_Khamra_MW)
## There is a significant positive relationship between measured PIH and FRP
## F-statistic: 1.669 on 1 and 3370 DF,  p-value: 0.1965
## Multiple R-squared:  0.0004949 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, measured PIH increases by 0.02468 units 
## (0.197)


#####
## Multiple regression #########################################################
#####

multi_Khamra_PIH_MW      <- lm(`20%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Khamra)
mul_coeff_Khamra_PIH_MW  <- coefficients(multi_Khamra_PIH_MW)
## (Intercept)        FRP_MW        w_spd 
## 201.37611823   0.02804188  -9.10709470  

## Summary of data #############################################################

summary(multi_Khamra_PIH_MW)
## F-statistic:  46.6 on 2 and 3369 DF,  p-value: < 2.2e-16
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.02692 (Variance clarification)
## If the FRP increases by 1, measured PIH increases by 0.02804 units (0.137)
## If the wind speed increases by 1, measured PIH decreases by 9.10709 units
## (high significant, <2e-16 ***)

## Plotting ####################################################################

library(scatterplot3d)
png(glue("Results/Simulation/PIH/Khamra/20/mul/Multiple_reg_Khamra_PIH_MW_20.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`20%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 1500), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.2]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.02692 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/20/mul/Multiple_reg_Khamra_PIH_MW_points_20.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`20%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.2]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.02692 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Khamra$PIH_Z   <- scale(tibble_frp_pih_wind_Khamra$`20%`)
tibble_frp_pih_wind_Khamra$FRP_Z   <- scale(tibble_frp_pih_wind_Khamra$FRP)
tibble_frp_pih_wind_Khamra$w_spd_Z <- scale(tibble_frp_pih_wind_Khamra$w_spd)

## For FRP [MW]
Z_model_PIH_Khamra_MW <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Khamra)
summary(Z_model_PIH_Khamra_MW) ## Comparable values
## F-statistic:  46.6 on 2 and 3369 DF,  p-value: < 2.2e-16
## Multiple R-squared:  0.02692 (identical values)
## FRP_Z        2.528e-02  1.700e-02   1.487    0.137    
## w_spd_Z     -1.626e-01  1.700e-02  -9.565   <2e-16 ***
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)


#####
## Lake Satagay ################################################################
#####

## Load and clean data #########################################################
tibble_frp_pih_wind_Satagay         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))

## Convert MW to cal/s
tibble_frp_pih_wind_Satagay$FRP_cal <- tibble_frp_pih_wind_Satagay$FRP*238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Satagay_MW  <- lm(`20%` ~ FRP, data = tibble_frp_pih_wind_Satagay)
coeff_satagay_PIH_MW   <- coefficients(lm_FRP_PIH_Satagay_MW)
## (Intercept)           FRP 
## 192.40534632  -0.03689595  
eq_S_PIH_MW <- paste0("y = ",round(coeff_satagay_PIH_MW[1],2), " + ", round(coeff_satagay_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/20/lm/lm_FRP_PIH_Satagay_MW_20.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`20%`,
     xlim = c(0, 300), ylim = c(0, 2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.2]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.001914 | p-value = 1.224e-08", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()


png(glue("Results/Simulation/PIH/Satagay/20/lm/lm_FRP_PIH_Satagay_MW_1000_20.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`20%`,
     xlim = c(0, 300), ylim = c(0, 1000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.2]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.001914 | p-value = 1.224e-08", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()

## Summary of data #############################################################

summary(lm_FRP_PIH_Satagay_MW)
## There is a significant positive relationship between calculated PIH and FRP
## F-statistic: 32.48 on 1 and 16936 DF,  p-value: 1.224e-08
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.001914 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, calculated PIH decreases by 0.036896  units 
## (high significant, 1.22e-08 ***)

#####
## Multiple regression #########################################################
#####

multi_Satagay_PIH_MW    <- lm(`20%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Satagay)
mul_coeff_Satagay_PIH_MW  <- coefficients(multi_Satagay_PIH_MW)
## (Intercept)        FRP_MW        w_spd 
## 204.68429891  -0.03530692  -1.81648290 

## Summary of data #############################################################

summary(multi_Satagay_PIH_MW)
## F-statistic: 27.94 on 2 and 16935 DF,  p-value: 7.698e-13
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.003289 (Variance clarification)
## If the FRP increases by 1, calculated PIH decreases by 0.035307 units 
## (high significant, 5.10e-08 ***)
## If the wind speed increases by 1, calculated PIH decreases by 1.816483 units
## (high significant, 1.36e-06 ***)

## Plotting ####################################################################


png(glue("Results/Simulation/PIH/Satagay/20/mul/Multiple_reg_Satagay_MW_20.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`20%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.2]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.003289 | p-value = 7.698e-13", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Satagay/20/mul/Multiple_reg_Satagay_points_MW_20.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`20%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.2]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.003289 | p-value = 7.698e-13", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Satagay$PIH_Z    <- scale(tibble_frp_pih_wind_Satagay$`20%`)
tibble_frp_pih_wind_Satagay$FRP_Z    <- scale(tibble_frp_pih_wind_Satagay$FRP)
tibble_frp_pih_wind_Satagay$w_spd_Z  <- scale(tibble_frp_pih_wind_Satagay$w_spd)


Z_model_PIH_Satagay <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Satagay)
summary(Z_model_PIH_Satagay) ## Comparable values
## F-statistic: 27.94 on 2 and 16935 DF,  p-value: 7.698e-13 (identical values)
## Multiple R-squared:  0.003289 (identical values)
## FRP_Z       -4.187e-02  7.682e-03  -5.450 5.10e-08 ***
## w_spd_Z     -3.712e-02  7.682e-03  -4.833 1.36e-06 ***

## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)???


## Q[0.5] ######################################################################

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Khamra_MW  <- lm(`50%` ~ FRP, data = tibble_frp_pih_wind_Khamra)
coeff_Khamra_PIH_MW   <- coefficients(lm_FRP_PIH_Khamra_MW)
## (Intercept)             FRP 
## 318.45117685     0.09479513 
eq_K_PIH_MW <- paste0("y = ",round(coeff_Khamra_PIH_MW[1],2), " + ", round(coeff_Khamra_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Khamra/50/lm/lm_FRP_PIH_Khamra_MW_50.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`50%`,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.5]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.00424 | p-value = 0.0001545", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/50/lm/lm_FRP_PIH_Khamra_1000_MW_50.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`50%`,
     xlim = c(0, 300), ylim = c(0,1000), 
     xlab = "FRP [MW]", 
     ylab = "Measured PIH [m] | Q[0.5]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.00424 | p-value = 0.0001545", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()


## Summary of the data #########################################################

summary(lm_FRP_PIH_Khamra_MW)
## There is a significant positive relationship between measured PIH and FRP
## F-statistic: 14.35 on 1 and 3370 DF,  p-value: 0.0001545
## Null hypothesis accepted, because p-value is < 0.05
## Multiple R-squared:  0.00424 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, measured PIH increases by 0.09480 units 
## (high significant, 0.000155 ***)


#####
## Multiple regression #########################################################
#####

multi_Khamra_PIH_MW      <- lm(`50%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Khamra)
mul_coeff_Khamra_PIH_MW  <- coefficients(multi_Khamra_PIH_MW)
## (Intercept)        FRP_MW        w_spd 
## 343.01853576   0.09696866  -5.88509657 

## Summary of data #############################################################

summary(multi_Khamra_PIH_MW)
## F-statistic: 18.13 on 2 and 3369 DF,  p-value: 1.477e-08
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.01065 (Variance clarification)
## If the FRP increases by 1, measured PIH increases by 0.09697 units 
## (high significant, 0.000104 ***)
## If the wind speed increases by 1, measured PIH decreases by 5.88510 units
## (high significant, 3.12e-06 ***)
eq_multi_PIH_K <- paste0("y = ", round(mul_coeff_Khamra_PIH_MW[1],2), " + ", round(mul_coeff_Khamra_PIH_MW[2], 2), "*x1 ", round(mul_coeff_Khamra_PIH_MW[3], 2), "*x2")
## "y = 343.02 + 0.1*x1 -5.89*x2" 


## Plotting ####################################################################

library(scatterplot3d)
png(glue("Results/Simulation/PIH/Khamra/50/mul/Multiple_reg_Khamra_PIH_MW_50.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`50%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.5]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.01065 | p-value = 1.477e-08", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/50/mul/Multiple_reg_Khamra_PIH_MW_points_50.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`50%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.5]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.01065 | p-value = 1.477e-08", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Khamra$PIH_Z   <- scale(tibble_frp_pih_wind_Khamra$`50%`)
tibble_frp_pih_wind_Khamra$FRP_Z   <- scale(tibble_frp_pih_wind_Khamra$FRP)
tibble_frp_pih_wind_Khamra$w_spd_Z <- scale(tibble_frp_pih_wind_Khamra$w_spd)

## For FRP [MW]
Z_model_PIH_Khamra_MW <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Khamra)
summary(Z_model_PIH_Khamra_MW) ## Comparable values
## F-statistic: 18.13 on 2 and 3369 DF,  p-value: 1.477e-08
## Multiple R-squared:  0.01065 (identical values)
## FRP_Z        6.661e-02  1.714e-02   3.886 0.000104 ***
## w_spd_Z     -8.006e-02  1.714e-02  -4.671 3.12e-06 ***
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)


#####
## Lake Satagay ################################################################
#####

## Load and clean data #########################################################
tibble_frp_pih_wind_Satagay         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))

## Convert MW to cal/s
tibble_frp_pih_wind_Satagay$FRP_cal <- tibble_frp_pih_wind_Satagay$FRP*238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Satagay_MW  <- lm(`50%` ~ FRP, data = tibble_frp_pih_wind_Satagay)
coeff_satagay_PIH_MW   <- coefficients(lm_FRP_PIH_Satagay_MW)
## (Intercept)           FRP 
## 370.29409101  -0.03533874 
eq_S_PIH_MW <- paste0("y = ",round(coeff_satagay_PIH_MW[1],2), " + ", round(coeff_satagay_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/50/lm/lm_FRP_PIH_Satagay_MW_50.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`50%`,
     xlim = c(0, 300), ylim = c(0, 2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.5]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.0007987 | p-value = 0.0002346", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()


png(glue("Results/Simulation/PIH/Satagay/50/lm/lm_FRP_PIH_Satagay_MW_1000_50.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`50%`,
     xlim = c(0, 300), ylim = c(0, 1000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.5]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.0007987 | p-value = 0.0002346", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()

## Summary of data #############################################################

summary(lm_FRP_PIH_Satagay_MW)
## There is a significant positive relationship between calculated PIH and FRP
## F-statistic: 13.54 on 1 and 16936 DF,  p-value: 0.0002346
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.0007987 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, calculated PIH decreases by 0.035339 units 
## (high significant, 0.000235 ***)

#####
## Multiple regression #########################################################
#####

multi_Satagay_PIH_MW    <- lm(`50%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Satagay)
mul_coeff_Satagay_PIH_MW  <- coefficients(multi_Satagay_PIH_MW)
## (Intercept)        FRP_MW        w_spd 
## 380.99637223  -0.03395374  -1.58323853 

## Summary of data #############################################################

summary(multi_Satagay_PIH_MW)
## F-statistic:  10.8 on 2 and 16935 DF,  p-value: 2.058e-05
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.001274 (Variance clarification)
## If the FRP increases by 1, calculated PIH decreases by 0.033954 units 
## (high significant, 0.000415 ***)
## If the wind speed increases by 1, calculated PIH decreases by 1.583239 units
## (significant, 0.004548**)
eq_multi_PIH_S <- paste0("y = ", round(mul_coeff_Satagay_PIH_MW[1],2), " + ", round(mul_coeff_Satagay_PIH_MW[2], 2), "*x1 ", round(mul_coeff_Satagay_PIH_MW[3], 2), "*x2")
## "y = 381 + -0.03*x1 -1.58*x2"

## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/50/mul/Multiple_reg_Satagay_MW_50.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`50%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.5]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.001274 | p-value = 2.058e-05", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Satagay/50/mul/Multiple_reg_Satagay_points_MW_50.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`50%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.5]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.001274 | p-value = 2.058e-05", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Satagay$PIH_Z    <- scale(tibble_frp_pih_wind_Satagay$`50%`)
tibble_frp_pih_wind_Satagay$FRP_Z    <- scale(tibble_frp_pih_wind_Satagay$FRP)
tibble_frp_pih_wind_Satagay$w_spd_Z  <- scale(tibble_frp_pih_wind_Satagay$w_spd)


Z_model_PIH_Satagay <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Satagay)
summary(Z_model_PIH_Satagay) ## Comparable values
## F-statistic:  10.8 on 2 and 16935 DF,  p-value: 2.058e-05 (identical values)
## Multiple R-squared:  0.001274 (identical values)
## FRP_Z       -2.715e-02  7.689e-03  -3.531 0.000415 ***
## w_spd_Z     -2.182e-02  7.689e-03  -2.838 0.004548 ** 
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)???


## Q[0.8] ######################################################################

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Khamra_MW  <- lm(`80%` ~ FRP, data = tibble_frp_pih_wind_Khamra)
coeff_Khamra_PIH_MW   <- coefficients(lm_FRP_PIH_Khamra_MW)
## (Intercept)             FRP 
## 518.2157655        0.1418887
eq_K_PIH_MW <- paste0("y = ",round(coeff_Khamra_PIH_MW[1],2), " + ", round(coeff_Khamra_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Khamra/80/lm/lm_FRP_PIH_Khamra_MW_80.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`80%`,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.8]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.005184 | p-value = 2.854e-05", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/80/lm/lm_FRP_PIH_Khamra_1000_MW_80.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`80%`,
     xlim = c(0, 300), ylim = c(0,1000), 
     xlab = "FRP [MW]", 
     ylab = "Measured PIH [m] | Q[0.8]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.005184 | p-value = 2.854e-05", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()


## Summary of the data #########################################################

summary(lm_FRP_PIH_Khamra_MW)
## There is a significant positive relationship between measured PIH and FRP
## F-statistic: 17.56 on 1 and 3370 DF,  p-value: 2.854e-05
## Multiple R-squared:  0.005184 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, measured PIH increases by 0.14189 units 
## (high significant, 2.85e-05 ***)


#####
## Multiple regression #########################################################
#####

multi_Khamra_PIH_MW      <- lm(`80%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Khamra)
mul_coeff_Khamra_PIH_MW  <- coefficients(multi_Khamra_PIH_MW)
## (Intercept)        FRP_MW        w_spd 
## 558.2814046     0.1454334   -9.5977006  

## Summary of data #############################################################

summary(multi_Khamra_PIH_MW)
## F-statistic: 24.76 on 2 and 3369 DF,  p-value: 2.123e-11
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.01448 (Variance clarification)
## If the FRP increases by 1, measured PIH increases by 0.14543 units 
## (high significant, 1.65e-05 ***)
## If the wind speed increases by 1, measured PIH decreases by 9.59770 units
## (high significant, 1.86e-08 ***)

## Plotting ####################################################################

library(scatterplot3d)
png(glue("Results/Simulation/PIH/Khamra/80/mul/Multiple_reg_Khamra_PIH_MW_80.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`80%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.8]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.01448 | p-value = 2.123e", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/80/mul/Multiple_reg_Khamra_PIH_MW_points_80.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`80%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.8]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.01448 | p-value = 2.123e", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Khamra$PIH_Z   <- scale(tibble_frp_pih_wind_Khamra$`80%`)
tibble_frp_pih_wind_Khamra$FRP_Z   <- scale(tibble_frp_pih_wind_Khamra$FRP)
tibble_frp_pih_wind_Khamra$w_spd_Z <- scale(tibble_frp_pih_wind_Khamra$w_spd)

## For FRP [MW]
Z_model_PIH_Khamra_MW <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Khamra)
summary(Z_model_PIH_Khamra_MW) ## Comparable values
## F-statistic: 24.76 on 2 and 3369 DF,  p-value: 2.123e-11
## Multiple R-squared:  0.01448 (identical values)
## FRP_Z        7.380e-02  1.711e-02   4.314 1.65e-05 ***
## w_spd_Z     -9.645e-02  1.711e-02  -5.638 1.86e-08 ***
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)


#####
## Lake Satagay ################################################################
#####

## Load and clean data #########################################################
tibble_frp_pih_wind_Satagay         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))

## Convert MW to cal/s
tibble_frp_pih_wind_Satagay$FRP_cal <- tibble_frp_pih_wind_Satagay$FRP*238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Satagay_MW  <- lm(`80%` ~ FRP, data = tibble_frp_pih_wind_Satagay)
coeff_satagay_PIH_MW   <- coefficients(lm_FRP_PIH_Satagay_MW)
## (Intercept)             FRP 
## 597.793051704  -0.009205707
eq_S_PIH_MW <- paste0("y = ",round(coeff_satagay_PIH_MW[1],2), " + ", round(coeff_satagay_PIH_MW[2], 2), "*x")
 

## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/80/lm/lm_FRP_PIH_Satagay_MW_80.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`80%`,
     xlim = c(0, 300), ylim = c(0, 2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.8]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 3.219e-05 | p-value = 0.4603", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()


png(glue("Results/Simulation/PIH/Satagay/80/lm/lm_FRP_PIH_Satagay_MW_1000_80.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`80%`,
     xlim = c(0, 300), ylim = c(0, 1000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.8]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 3.219e-05 | p-value = 0.4603", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()

## Summary of data #############################################################

summary(lm_FRP_PIH_Satagay_MW)
## There is a significant positive relationship between calculated PIH and FRP
## F-statistic: 0.5452 on 1 and 16936 DF,  p-value: 0.4603
## Multiple R-squared:  3.219e-05 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, calculated PIH decreases by 0.009206 units 
## (0.46)

#####
## Multiple regression #########################################################
#####

multi_Satagay_PIH_MW    <- lm(`80%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Satagay)
mul_coeff_Satagay_PIH_MW  <- coefficients(multi_Satagay_PIH_MW)
## (Intercept)          FRP_MW         w_spd 
## 630.275795501  -0.005002076  -4.805324257  

## Summary of data #############################################################

summary(multi_Satagay_PIH_MW)
## F-statistic: 22.33 on 2 and 16935 DF,  p-value: 2.056e-10
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.002631 (Variance clarification)
## If the FRP increases by 1, calculated PIH decreases by 0.005002 units 
## (0.688)
## If the wind speed increases by 1, calculated PIH decreases by 4.805324 units
## (high significant, 3.18e-11 ***)

## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/80/mul/Multiple_reg_Satagay_MW_80.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`80%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.8]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.002631 | p-value = 2.056e-10", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Satagay/80/mul/Multiple_reg_Satagay_points_MW_80.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`80%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.8]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.002631 | p-value = 2.056e-10", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Satagay$PIH_Z    <- scale(tibble_frp_pih_wind_Satagay$`80%`)
tibble_frp_pih_wind_Satagay$FRP_Z    <- scale(tibble_frp_pih_wind_Satagay$FRP)
tibble_frp_pih_wind_Satagay$w_spd_Z  <- scale(tibble_frp_pih_wind_Satagay$w_spd)


Z_model_PIH_Satagay <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Satagay)
summary(Z_model_PIH_Satagay) ## Comparable values
## F-statistic: 22.33 on 2 and 16935 DF,  p-value: 2.056e-10(identical values)
## Multiple R-squared:  0.002631 (identical values)
## FRP_Z       -3.083e-03  7.684e-03  -0.401    0.688    
## w_spd_Z     -5.104e-02  7.684e-03  -6.642 3.18e-11 ***

## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater) ????


## Q[0.975] ######################################################################

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Khamra_MW  <- lm(`97.5%` ~ FRP, data = tibble_frp_pih_wind_Khamra)
coeff_Khamra_PIH_MW   <- coefficients(lm_FRP_PIH_Khamra_MW)
## (Intercept)             FRP 
## 843.763319         0.209933 
eq_K_PIH_MW <- paste0("y = ",round(coeff_Khamra_PIH_MW[1],2), " + ", round(coeff_Khamra_PIH_MW[2], 2), "*x")


## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Khamra/97.5/lm/lm_FRP_PIH_Khamra_MW_97.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`97.5%`,
     xlim = c(0, 300), ylim = c(0,2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.975]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.004701 | p-value = 6.759e-05", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/97.5/lm/lm_FRP_PIH_Khamra_1000_MW_97.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Khamra$FRP, y = tibble_frp_pih_wind_Khamra$`97.5%`,
     xlim = c(0, 300), ylim = c(0,1000), 
     xlab = "FRP [MW]", 
     ylab = "Measured PIH [m] | Q[0.975]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between measured PIH and FRP for study area Lake Khamra | {eq_K_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.004701 | p-value = 6.759e-05", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Khamra_MW, col = "red", lwd = 2)
dev.off()


## Summary of the data #########################################################

summary(lm_FRP_PIH_Khamra_MW)
## There is a significant positive relationship between measured PIH and FRP
## F-statistic: 15.92 on 1 and 3370 DF,  p-value: 6.759e-05
## Multiple R-squared:  0.004701 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, measured PIH increases by 0.20993 units 
## (high significant, 6.76e-05 ***)


#####
## Multiple regression #########################################################
#####

multi_Khamra_PIH_MW      <- lm(`97.5%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Khamra)
mul_coeff_Khamra_PIH_MW  <- coefficients(multi_Khamra_PIH_MW)
## (Intercept)        FRP_MW        w_spd 
## 904.2957343      0.2152884 -14.5005049  

## Summary of data #############################################################

summary(multi_Khamra_PIH_MW)
## F-statistic: 23.04 on 2 and 3369 DF,  p-value: 1.151e-10
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.01349 (Variance clarification)
## If the FRP increases by 1, measured PIH increases by 0.2153 units 
## (high significant, 4.08e-05 ***)
## If the wind speed increases by 1, measured PIH decreases by 14.5005 units
## (high significant, 4.57e-08 ***)

## Plotting ####################################################################

library(scatterplot3d)
png(glue("Results/Simulation/PIH/Khamra/97.5/mul/Multiple_reg_Khamra_PIH_MW_97.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`97.5%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.975]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.01349 | p-value = 1.151e-10", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Khamra/97.5/mul/Multiple_reg_Khamra_PIH_MW_points_97.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Khamra$`97.5%`, x = tibble_frp_pih_wind_Khamra$FRP, z = tibble_frp_pih_wind_Khamra$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.975]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.01349 | p-value = 1.151e-10", adj = 1, cex.sub = 1.5, 
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Khamra"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Khamra$PIH_Z   <- scale(tibble_frp_pih_wind_Khamra$`97.5%`)
tibble_frp_pih_wind_Khamra$FRP_Z   <- scale(tibble_frp_pih_wind_Khamra$FRP)
tibble_frp_pih_wind_Khamra$w_spd_Z <- scale(tibble_frp_pih_wind_Khamra$w_spd)

## For FRP [MW]
Z_model_PIH_Khamra_MW <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Khamra)
summary(Z_model_PIH_Khamra_MW) ## Comparable values
## F-statistic: 23.04 on 2 and 3369 DF,  p-value: 1.151e-10
## Multiple R-squared:  0.01349 (identical values)
## FRP_Z        7.031e-02  1.711e-02   4.108 4.08e-05 ***
## w_spd_Z     -9.379e-02  1.711e-02  -5.480 4.57e-08 ***
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)


#####
## Lake Satagay ################################################################
#####

## Load and clean data #########################################################
tibble_frp_pih_wind_Satagay         <- get(load("Results/Simulation/data/tibble_frp_pih_wind_Satagay.rda"))

## Convert MW to cal/s
tibble_frp_pih_wind_Satagay$FRP_cal <- tibble_frp_pih_wind_Satagay$FRP*238845.8966275

#####
## Calculation with FRP [MW] ###################################################
#####

#####
## Linear regression ###########################################################
#####

lm_FRP_PIH_Satagay_MW  <- lm(`97.5%` ~ FRP, data = tibble_frp_pih_wind_Satagay)
coeff_satagay_PIH_MW   <- coefficients(lm_FRP_PIH_Satagay_MW)
## (Intercept)             FRP 
## 960.3124155       0.1077719 
eq_S_PIH_MW <- paste0("y = ",round(coeff_satagay_PIH_MW[1],2), " + ", round(coeff_satagay_PIH_MW[2], 2), "*x")

## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/97.5/lm/lm_FRP_PIH_Satagay_MW_97.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`97.5%`,
     xlim = c(0, 300), ylim = c(0, 2000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.975]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.002487 | p-value = 8.34e-11", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()


png(glue("Results/Simulation/PIH/Satagay/97.5/lm/lm_FRP_PIH_Satagay_MW_1000_97.png"), width = 1000, height = 800)
plot(x = tibble_frp_pih_wind_Satagay$FRP, y = tibble_frp_pih_wind_Satagay$`97.5%`,
     xlim = c(0, 300), ylim = c(0, 1000),
     xlab = "FRP [MW]",
     ylab = "Measured PIH [m] | Q[0.975]", col = "lightblue3", cex.lab = 1.5)
title(main = glue("The linear regression between calculated PIH and FRP for study area Lake Satagay | {eq_S_PIH_MW}"), cex.main = 1.5)
title(sub = "R² = 0.002487 | p-value = 8.34e-11", adj = 1, cex.sub = 1.5)
abline(lm_FRP_PIH_Satagay_MW, col = "red", lwd = 2)
dev.off()

## Summary of data #############################################################

summary(lm_FRP_PIH_Satagay_MW)
## There is a significant positive relationship between calculated PIH and FRP
## F-statistic: 42.23 on 1 and 16936 DF,  p-value: 8.34e-11
## Multiple R-squared:  0.002487 (Variance clarification) is not a good or stable 
## value. 1 would be perfect to explain the model with the variance of the 
## dependent variable. 
## If the FRP increases by 1, calculated PIH decreases by 0.10777 units 
## (high significant, 8.34e-11 ***)

#####
## Multiple regression #########################################################
#####

multi_Satagay_PIH_MW    <- lm(`97.5%` ~ FRP + w_spd, data = tibble_frp_pih_wind_Satagay)
mul_coeff_Satagay_PIH_MW  <- coefficients(multi_Satagay_PIH_MW)
## (Intercept)          FRP_MW          w_spd 
## 1016.6883455       0.1150676    -8.3399551  

## Summary of data #############################################################

summary(multi_Satagay_PIH_MW)
## F-statistic: 58.83 on 2 and 16935 DF,  p-value: < 2.2e-16
## Null hypothesis rejected, because p-value is < 0.05
## Multiple R-squared:  0.0069 (Variance clarification)
## If the FRP increases by 1, calculated PIH increases by 0.11507 units 
## (high significant, 3.93e-12 ***)
## If the wind speed increases by 1, calculated PIH decreases by 8.33996 units
## (high significant, < 2e-16 ***)

## Plotting ####################################################################

png(glue("Results/Simulation/PIH/Satagay/97.5/mul/Multiple_reg_Satagay_MW_97.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`97.5%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.975]", zlab = "spd [m/s]", 
              type = "h", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.0069 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()

png(glue("Results/Simulation/PIH/Satagay/97.5/mul/Multiple_reg_Satagay_points_MW_97.png"), width = 1000, height = 800)
scatterplot3d(y = tibble_frp_pih_wind_Satagay$`97.5%`, x = tibble_frp_pih_wind_Satagay$FRP, z = tibble_frp_pih_wind_Satagay$w_spd, 
              xlim = c(0, 300), ylim = c(0, 2000), zlim = c(0,20), color = "lightblue3", 
              xlab = "FRP [MW]", ylab = "Measured PIH [m] | Q[0.975]", zlab = "spd [m/s]", 
              type = "p", font.lab = NULL, cex.axis = 1, cex.lab = 1.5)
title(sub = "R² = 0.0069 | p-value < 2.2e-16", adj = 1, cex.sub = 1.5,
      main = glue("The multiple regression between wind speed, FRP and calculated PIH for study area Lake Satagay"),
      cex.main = 1.4)
dev.off()


## Which variable has the greatest impact on the plume injection height? 
tibble_frp_pih_wind_Satagay$PIH_Z    <- scale(tibble_frp_pih_wind_Satagay$`97.5%`)
tibble_frp_pih_wind_Satagay$FRP_Z    <- scale(tibble_frp_pih_wind_Satagay$FRP)
tibble_frp_pih_wind_Satagay$w_spd_Z  <- scale(tibble_frp_pih_wind_Satagay$w_spd)


Z_model_PIH_Satagay <- lm(PIH_Z ~ FRP_Z + w_spd_Z, data = tibble_frp_pih_wind_Satagay)
summary(Z_model_PIH_Satagay) ## Comparable values
## F-statistic: 58.83 on 2 and 16935 DF,  p-value: < 2.2e-16 (identical values)
## Multiple R-squared:  0.0069 (identical values)
##FRP_Z        5.325e-02  7.668e-03   6.945 3.93e-12 ***
##  w_spd_Z     -6.651e-02  7.668e-03  -8.675  < 2e-16 ***
## ??
## The fire intensity has a greater impact on the height of the smoke plume during 
## a fire event (value is greater)


