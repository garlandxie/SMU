##### Plots for SMU_FDPD ########

### Libraries #####
library(ggplot2)
library(data.table)
library(dplyr)


### Metadata #####

## lm_rev_mods: list of "re-evaluated" linear regression models from MuMIn 
# $org - lm(org ~ scale(RaoQ))
# $P - lm(P ~ scale(CWM_H))
# $K - lm(K ~ scale(CWM_H) + scale(CWM_LA) + scale(RaoQ) + scale(CWM_SLA))
# $NO3 - lm(NO3 ~ scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA))
# $totbio - lm(totbio ~ scale(CWM_LA) + scale(CWM_SLA))
# $strat10 - lm(strat10 ~ scale(CWM_H) + scale(mPD) + scale(PD))
# nothing significant for crat10 from MuMIn
# hits10 - lm(hits10 ~ scale(CWM_LA) + scale(CWM_SLA) + scale(RaoQ))
# nothing significant for ind1 from MuMIn
# ind2 ~ lm(ind2 ~ scale(CWM_H) + scale(CWM_LA))

### Mean values of diversity treatments ####
fnn <- data.table(fnn) 
fnn_agg <- as.data.frame(fnn[, j= list(mean(mPD, na.rm = TRUE),
                                       mean(PD, na.rm = TRUE), 
                                       mean(FD, na.rm = TRUE),
                                       mean(RaoQ, na.rm = TRUE),
                                       sd(mPD, na.rm = TRUE), 
                                       sd(PD, na.rm = TRUE), 
                                       sd(FD, na.rm = TRUE),
                                       sd(RaoQ, na.rm = TRUE)),
                                       by = treat])
names(fnn_agg) <- c("treat", "mean.MPD", "mean.PD", "mean.FD", "mean.RaoQ", "sd.MPD", "sd.PD", "sd.FD", "sd.RaoQ")
           

### Multi-panel figures - CWM.H ####

## CWM_H 
## try to write a custom function next time
par(mfrow = c(2,3))
par(mar=c(5.1,5,4.1,2.1))
    
termplot(lm_rev_mods$NO3, 
         term = "scale(CWM_H)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Nitrate-N")
termplot(lm_rev_mods$K, 
         term = "scale(CWM_H)", 
         cex.lab = 1.5,
         pch = 16,
         partial.resid = T, 
         xlab = "",
         ylab = "Potassium")
termplot(lm_rev_mods$ind2,  
         term = "scale(CWM_H)", 
         cex.lab = 1.5,
         pch = 16,
         partial.resid = T, 
         xlab = "",
         ylab = "Multifunctionality index 2")
termplot(lm_rev_mods$P, 
         term = "scale(CWM_H)", 
         cex.lab = 1.5,
         pch = 16,
         partial.resid = T, 
         xlab = "",
         ylab = "Phosphorus")
termplot(lm_rev_mods$strat10, 
         term = "scale(CWM_H)", 
         cex.lab = 1.5,
         pch = 16,
         partial.resid = T, 
         xlab = "",
         ylab = "Soil Temperature Index")

# CWM_LA
par(mfrow = c(2,2))
par(mar=c(5.1,5,4.1,2.1))

termplot(lm_rev_mods$NO3, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Nitrate-N")
termplot(lm_rev_mods$hits10, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Canopy Density")
termplot(lm_rev_mods$K, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Potassium")
termplot(lm_rev_mods$ind2, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Multifunctionality Index 2")

# CWM_SLA 
par(mfrow = c(1,3))
par(mar=c(5.1,5,4.1,2.1))

termplot(lm_rev_mods$NO3, 
         term = "scale(CWM_SLA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Nitrate-N")
termplot(lm_rev_mods$hits10, 
         term = "scale(CWM_SLA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Canopy Density")
termplot(lm_rev_mods$totbio, 
         term = "scale(CWM_SLA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Aboveground Biomass")

# RaoQ
par(mfrow = c(1,3))
par(mar=c(5.1,5,4.1,2.1))

termplot(lm_rev_mods$hits10, 
         term = "scale(RaoQ)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Canopy Density")
termplot(lm_rev_mods$org, 
         term = "scale(RaoQ)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Soil Organic Matter")
termplot(lm_rev_mods$K, 
         term = "scale(RaoQ)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Substrate Potassium Content")

# PD
par(mfrow = c(1,1))
par(mar=  c(5.1,4.1,4.1,2.1))
termplot(lm_rev_mods$strat10, 
         term = "scale(PD)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Faith's PD")
