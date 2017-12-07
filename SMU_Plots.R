##### Plots for SMU_FDPD ########

#### Libraries ####
library(ggplot2)
library(data.table)

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

### Calculate aggregate values of diversity treatments ####

# data.table - syntax can be used to aggregate rows 
fnn <- data.table(fnn) 

# calculate mean and std values for each treatment; relabel columns 
fnn_agg <- as.data.frame(fnn[, j= list(mean.MPD = mean(mPD, na.rm = TRUE),
                                       mean.PD = mean(PD, na.rm = TRUE), 
                                       mean.FR = mean(FD, na.rm = TRUE),
                                       mean.RaoQ = mean(RaoQ, na.rm = TRUE),
                                       sd.MPD = sd(mPD, na.rm = TRUE), 
                                       sd.PD = sd(PD, na.rm = TRUE), 
                                       sd.FR = sd(FD, na.rm = TRUE),
                                       sd.RaoQ = sd(RaoQ, na.rm = TRUE)),
                                       by = treat])

# reorder factor levels for treatment; helps clean up the supplementary graphs
fnn_agg$treat <- factor(fnn_agg$treat, levels = c("d", "g", "s", "t", 
                                                  "dgs", "dgt", "dgc", "dts", "dct", "dcs",
                                                  "gcs", "gct", "gts", 
                                                  "cts",
                                                  "cdgst"))

### Mean RaoQ, mixtures only #####

# set up data
ggplot(data = fnn_agg, aes(x = factor(treat), y = mean.RaoQ), ylim = 0.15) + 
  # bar graph layer
  geom_bar(position = position_dodge(), stat = "identity") +
  # error bar layer
  geom_errorbar(aes(ymin = mean.RaoQ, ymax = mean.RaoQ + sd.RaoQ), position = position_dodge(.9)) +
  # axis titles
  labs(x = "Treatment", y = "Mean RaoQ") +
  # font size and margins for axis titles
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, b = 15)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 15)),
        axis.text = element_text(size = 12))

### Mean PD, mixtures only ####

# set up data
ggplot(data = fnn_agg, aes(x = factor(treat), y = mean.PD), ylim = 0.15) + 
  # bar graph layer
  geom_bar(position = position_dodge(), stat = "identity") +
  # error bar layer
  geom_errorbar(aes(ymin = mean.PD, ymax = mean.PD + sd.PD), position = position_dodge(.9)) +
  # axis titles
  labs(x = "Treatment", y = "Mean Faith's PD") +
  # font size and margins for axis titles
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, b = 15)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 15)),
        axis.text = element_text(size = 12))

### Mean MPD, mixtures only ####

# set up data
ggplot(data = fnn_agg, aes(x = factor(treat), y = mean.MPD), ylim = 0.15) + 
  # bar graph layer
  geom_bar(position = position_dodge(), stat = "identity") +
  # error bar layer
  geom_errorbar(aes(ymin = mean.MPD, ymax = mean.MPD + sd.MPD), position = position_dodge(.9)) +
  # axis titles
  labs(x = "Treatment", y = "Mean MPD") +
  # font size and margins for axis titles
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, b = 15)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 15)),
        axis.text = element_text(size = 12))

### Mean FR, mixtures only #####

# set up data
ggplot(data = fnn_agg, aes(x = factor(treat), y = mean.FR), ylim = 0.15) + 
  # bar graph layer
  geom_bar(position = position_dodge(), stat = "identity") +
  # error bar layer
  geom_errorbar(aes(ymin = mean.FR, ymax = mean.FR + sd.FR), position = position_dodge(.9)) +
  # axis titles
  labs(x = "Treatment", y = "Mean dendrogram-based FR") +
  # font size and margins for axis titles
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, b = 15)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 15)),
        axis.text = element_text(size = 12))

### Multi-panel figures - CWM_H ####

## CWM_H 
## try to write a custom function next time
# create two horizontal panels that contains three plots 
par(mfrow = c(2,2))
# adjust plot margins
par(mar=c(5.1,5.3,4.1,2.1))
# adjust outer bottom margin to include CWM_H text
par(oma = c(4,0,0,0))
#plots

termplot(lm_rev_mods$NO3, 
         term = "scale(CWM_H)",
         cex.lab = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Substrate Nitrate-N")

# plot sublabels
text(x = 3.6, y = 3.6, labels = "a)", pos = 1, cex = 1.5)

#termplot(lm_rev_mods$K, 
#         term = "scale(CWM_H)", 
#         cex.lab = 1.5,
#         pch = 16,
#         partial.resid = T, 
#         xlab = "",
#         ylab = "Substrate Potassium")

termplot(lm_rev_mods$ind2,  
         term = "scale(CWM_H)", 
         cex.lab = 1.5,
         pch = 16,
         partial.resid = T, 
         xlab = "",
         ylab = "Multifunctionality index 2")

# plot sublabels
text(x = 3.6, y = 7.8, labels = "b)", pos = 1, cex = 1.5)

termplot(lm_rev_mods$P, 
         term = "scale(CWM_H)", 
         cex.lab = 1.5,
         pch = 16,
         partial.resid = T, 
         xlab = "",
         ylab = "Substrate Phosphorus")

# plot sublabels
text(x = 3.6, y = 75, labels = "c)", pos = 1, cex = 1.5)

termplot(lm_rev_mods$strat10, 
         term = "scale(CWM_H)", 
         cex.lab = 1.5,
         pch = 16,
         partial.resid = T, 
         xlab = "",
         ylab = "Substrate Temperature Index")

# plot sublabels
text(x = 3.6, y = 35, labels = "d)", pos = 1, cex = 1.5)


# add "CWM_H" to the bottom outer figure margin of the plot
# outer = T - use the figure margins for text
# side = 1 - bottom margin
mtext("CWM H", outer = T, side = 1)


### Multi-panel figures - CWM_LA #####

# CWM_LA
# create two horizontal panels that contains three plots 
par(mfrow = c(2,2))
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
         ylab = "Substrate Nitrate-N Content")

# plot sublabels
text(x = 2, y = 7.5, labels = "a)", pos = 1, cex = 1.5)

termplot(lm_rev_mods$hits10, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Canopy Density")

# plot sublabels
text(x = 2, y = 225, labels = "b)", pos = 1, cex = 1.5)

## removing this plot from final manuscript
#termplot(lm_rev_mods$K, 
#         term = "scale(CWM_LA)",
#         cex.lab = 1.5,
#         cex.axis = 1.5,
#         pch = 16, 
#         partial.resid = T, 
#         xlab = "",
#         ylab = "Substrate Potassium")

termplot(lm_rev_mods$ind2, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Multifunctionality Index 2")

# plot sublabels
text(x = 2, y = 10, labels = "c)", pos = 1, cex = 1.5)

termplot(lm_rev_mods$totbio, 
         term = "scale(CWM_LA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Aboveground biomass")

# plot sublabels
text(x = 2, y = 40, labels = "d)", pos = 1, cex = 1.5)

# add "CWM_LA" to the bottom outer figure margin of the plot
# outer = T - use the figure margins for text
# side = 1 - bottom margin
mtext("CWM LA", outer = T, side = 1)

### Multi-panel figures - CWM_SLA ####

# CWM_SLA 
# create two horizontal panels that contains two plots 
par(mfrow = c(1,3))
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
         ylab = "Substrate Nitrate-N")

# plot sublabels
text(x = 2.3, y = 3.7, labels = "a)", pos = 1, cex = 1.5)

termplot(lm_rev_mods$hits10, 
         term = "scale(CWM_SLA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Canopy Density")

# plot sublabels
text(x = 2.3, y = 190, labels = "b)", pos = 1, cex = 1.5)


termplot(lm_rev_mods$totbio, 
         term = "scale(CWM_SLA)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "",
         ylab = "Aboveground Biomass")

# plot sublabels
text(x = 2.3, y = 45, labels = "c)", pos = 1, cex = 1.5)


# removing this figure from manuscript
#termplot(lm_rev_mods$K, 
#         term = "scale(CWM_SLA)",
#         cex.lab = 1.5,
#         cex.axis = 1.5,
#         pch = 16, 
#         partial.resid = T, 
#         xlab = "",
#         ylab = "Substrate Potassium")

# add "CWM_SLA" to the bottom outer figure margin of the plot
# outer = T - use the figure margins for text
# side = 1 - bottom margin
mtext("CWM SLA", outer = T, side = 1)

### Multi-panel figures - Diversity Metrics #####

# diversity metrics
par(mfrow = c(2,2))
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

# plot sublabels
text(x = 0.11, y = 175, labels = "a)", pos = 1, cex = 1.5)

termplot(lm_rev_mods$org, 
         term = "scale(RaoQ)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "Rao's Q",
         ylab = "Substrate Organic Matter")

# plot sublabels
text(x = 0.11, y = 4.5, labels = "b)", pos = 1, cex = 1.5)


#termplot(lm_rev_mods$K, 
#         term = "scale(RaoQ)",
#         cex.lab = 1.5,
#         cex.axis = 1.5,
#         pch = 16, 
#         partial.resid = T, 
#         xlab = "Rao's Q",
#         ylab = "Substrate Potassium")

# PD 
termplot(lm_rev_mods$strat10, 
         term = "scale(PD)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "PD",
         ylab = "Substrate Temperature Index")

# plot sublabels
text(x = 800, y = 50, labels = "c)", pos = 1, cex = 1.5)

# mPD
termplot(lm_rev_mods$strat10, 
         term = "scale(mPD)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 16, 
         partial.resid = T, 
         xlab = "MPD",
         ylab = "Substrate Temperature Index")

# plot sublabels
text(x = 210, y = 45, labels = "d)", pos = 1, cex = 1.5)
