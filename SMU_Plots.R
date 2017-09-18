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
# create two horizontal panels that contains three plots 
par(mfrow = c(2,3))
# adjust plot margins
par(mar=c(5.1,5.3,4.1,2.1))
# adjust outer bottom margin to include CWM_H text
par(oma = c(4,0,0,0))
#plots
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

# add "CWM_H" to the bottom outer figure margin of the plot
# outer = T - use the figure margins for text
# side = 1 - bottom margin
mtext("CWM_H", outer = T, side = 1)

# CWM_LA
# create two horizontal panels that contains three plots 
par(mfrow = c(2,3))
# adjust plot margins
par(mar=c(5.1,5,4.1,2.1))
# adjust outer bottom margin to include CWM_H text
par(oma = c(4,0,0,0))

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
termplot(lm_rev_mods$totbio, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Aboveground biomass")

# add "CWM_LA" to the bottom outer figure margin of the plot
# outer = T - use the figure margins for text
# side = 1 - bottom margin
mtext("CWM_LA", outer = T, side = 1)


# CWM_SLA 
# create two horizontal panels that contains two plots 
par(mfrow = c(2,2))
# adjust plot margins
par(mar=c(5.1,5,4.1,2.1))
# adjust outer bottom margin to include CWM_H text
par(oma = c(4,0,0,0))

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
termplot(lm_rev_mods$K, 
         term = "scale(CWM_SLA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Potassium")

# add "CWM_SLA" to the bottom outer figure margin of the plot
# outer = T - use the figure margins for text
# side = 1 - bottom margin
mtext("CWM_SLA", outer = T, side = 1)

# diversity metrics
par(mfrow = c(2,3))
par(mar=c(5.1,5,4.1,2.1))
# RaoQ
termplot(lm_rev_mods$hits10, 
         term = "scale(RaoQ)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "Rao's Q",
         ylab = "Canopy Density")
termplot(lm_rev_mods$org, 
         term = "scale(RaoQ)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "Rao's Q",
         ylab = "Soil Organic Matter")
termplot(lm_rev_mods$K, 
         term = "scale(RaoQ)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "Rao's Q",
         ylab = "Substrate Potassium Content")
# PD 
termplot(lm_rev_mods$strat10, 
         term = "scale(PD)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "PD",
         ylab = "Faith's PD")
# mPD
termplot(lm_rev_mods$strat10, 
         term = "scale(mPD)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "mPD",
         ylab = "Soil Tmperature Index")

