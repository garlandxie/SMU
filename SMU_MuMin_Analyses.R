########## SMU_MuMin_Analyses.R #########

### Import libraries ####
library(MuMIn) # multimodel inference
library(nlme) # mixed effect models, gives p-values 

### Metadata #####

## fnn 
# mod: ID for each green roof module (excludes monocultures since we couldn't calculate MPD nor Rao's Q)
# treat: experimental treatments 
# d, g, s, t (life forms, 3 spp)
# cts, dcs, dct, dgc, dgs, dgt, dts, gcs, gts (combos of lifeforms, 9 spp)
# cdgst (all spp, 15 spp)
# nosp - planted species richness (nosp)
# mPD - mean pairwise phylogenetic distance (phylo diversity measure)
# PD - faith's PD as phylogenetic richness (phylo diversity measure)
# SR - realized species richness (species richness during 4th year of experiment)
# RaoQ - Rao's quadratic entropy, meant to represent a form of variance in functional trait space (functional diversity)
# FD - functional richness calculated from a functional dendogram (incoporates PCA measures to reduce redudancy) 
# CWM_H - community weighted mean trait value for plant height
# CWM_LA - community weighted mean trait value for leaf area
# CWM_SLA - community weighted mean trait vaue for specific leaf area
# CWM_LDMC - community weighted mean trait value for specific dry leaf matter content
# org - service indicator for soil organic matter
# P - service indicator for substrate phosphorus content 
# K - service indicator for potassium phosphorus content
# NO3 - service indicator for substrate nitrate content
# totbio - service indicator for total community biomass
# strat10 - service indicator for soil temperature index
# crat10 - service indicator for stormwater capture index 
# ind 1 - multifunctionality index 1 (see below for more details)
# ind 2 - multifunctionality index 2 (see below for more details)

## avg_models: list of averaging R objects from MuMIn package
# $org - summary MuMIn output for substrate phosphorus content
# $K - summary MuMIn output for substrate potassium content
# $NO3 - summary MuMIn output for substrate potassium content
# $totbio - summary MuMIn output for aboveground biomass
# $strat10 - summary MuMIn output for soil temperature index
# $crat10 - summary MuMIn output for stormwater capture index
# $ind1 - summary MuMIn output for multifunctionality index 1 
# $ind2 - summary MuMIn output for multifunctionality index 2 

## dredge_mods: list of "model.selection" "data.frame" objects using MuMIn package
# $org - summary of all possible candidate models for substrate phosphorus content 
# $K - summary of all possible candidate models for substrate potassium content
# $NO3 - summary of all possible candidate models for substrate nitrate content 
# $totbio - summary of all possible candidate models for aboveground biomass
# $strat10 - summary of all possible candidate models for soil temperature index
# $crat10 - summary of all possible candidate models for stormwater capture index
# $ind1 - summary of all possible candidate models for multifunctionality index 1
# $ind2 - summary of all possible candidate models for multifunctionality index 2 

### Import files ####

# imported from Google Drive
load("C:/Users/garla/Google Drive/Research Projects/SMU_12sp_Phylogeny/SMU_PD-FD_Data/SMU_MuMin_Analyses.RData")

### Calculating multifunctionality indices ####

# multifunctionality index 1: includes both stormwater capture and soil temperature index
# note: you should center both service indicators first using scale()
# note: adding ind1 as a separate vector to fnn (reduces clutter in global R environment)
fnn$ind1 <- (scale(fnn$crat10) + scale(-1*(fnn$strat10)))

# multifunctionality index 2: includes all of the service indicators (see metadat for more info)
# note: you should center both service indicators first using scale()
# note: adding ind1 as a separate vector to fnn (reduces clutter in global R environment)
fnn$ind2 <- scale(fnn$crat10) + scale(-1*(fnn$strat10)) + scale(-1*(fnn$P)) + scale(-1*(fnn$NO3)) + scale(-1*(fnn$K)) + scale(fnn$totbio) + scale(fnn$hits10)+ scale(fnn$org)

# Reorganize dataset to make it more neat (e.g. group service indicators, then group diversity metrics together
# note: better to organize by column names as opposed to vector indices (more clear to other readers)
fnn <- fnn[, c("mod", "treat", "nosp", # treatments
               "mPD", "PD", "SR", "RaoQ", "FD" , # phylo and functional diversity metrics
               "CWM_H", "CWM_LA", "CWM_SLA", "CWM_LDMC", # CWM values
               "org", "P", "K", "NO3", "totbio","strat10", "crat10", "hits10", "ind1", "ind2") # service indicators
          ]
### Multimodel inference ####

## Mixed-effect models ##
# Create a large list of different lm models consisting of service indicators and diversity measures
# reason: easy to access and better management(?) of R workspace
# can't seem to use lapply() to shorten code: some bug associated with lme(); not elegant but also not a huge deal 
lm_models <- list(lme(org ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~1|scale(SR), method = "ML", data = fnn),
                  lme(P ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML", data = fnn),
                  lme(K ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML", data = fnn),
                  lme(NO3 ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML", data = fnn),
                  lme(totbio ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML",data = fnn),
                  lme(strat10 ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML", data = fnn),
                  lme(crat10 ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML", data = fnn),
                  lme(hits10 ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML", data = fnn),
                  lme(ind1 ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML", data = fnn),
                  lme(ind2 ~ scale(mPD) + scale(PD) + scale(RaoQ) + scale(FD) + scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), random = ~ 1|scale(SR), method = "ML",data = fnn)
                  )

# create names for each list object (lme model) for faster indexing
names(lm_models) <- c("org", "P", "K", "NO3", "totbio", "strat10", "crat10", "hits10", "ind1", "ind2")

## Model dredging ##
# Create a large list of dredge MuMin models using lapply and indexing from lm_models
# this assumes that order of service indicators is consistent 
dredge_mods <- lapply(1:length(lm_models), function(x) dredge(lm_models[[x]]))
                      
# create names for each list object
names(dredge_mods) <- c("org", "P", "K", "NO3", "totbio", "strat10", "crat10", "hits10", "ind1", "ind2")

## Model averaging ##
# Create a large list of model averaging outputs for each service indicator using lapply and indexing from dredge_mods
# combine all models based on Xingfeng's recommendation
avg_models <- lapply(1:length(dredge_mods), function(x) model.avg(dredge_mods[[x]], beta = TRUE))

# create names for each list object
names(avg_models) <- c("org", "P", "K", "NO3", "totbio", "strat10", "crat10", "hits10", "ind1", "ind2")

### Summary outputs for MuMin models ####
summary(avg_models$org)

#Model-averaged coefficients:  
#  (full average) 
#Estimate Std. Error Adjusted SE z value Pr(>|z|)  
#(Intercept)     0.00000    0.00000     0.00000      NA       NA  
#scale(RaoQ)     0.34752    0.17378     0.17537   1.982   0.0475 *
#scale(mPD)     -0.03644    0.10612     0.10725   0.340   0.7341  
#scale(FD)      -0.01778    0.11201     0.11354   0.157   0.8756  
#scale(PD)      -0.01541    0.10758     0.10912   0.141   0.8877  
#scale(CWM_H)    0.02586    0.10819     0.10936   0.236   0.8130  
#scale(CWM_SLA) -0.01527    0.06713     0.06796   0.225   0.8222  
#scale(CWM_LA)   0.01201    0.08692     0.08810   0.136   0.8915  

#(conditional average) 
#Estimate Std. Error Adjusted SE z value Pr(>|z|)   
#(Intercept)     0.00000    0.00000     0.00000      NA       NA   
#scale(RaoQ)     0.38017    0.14362     0.14571   2.609  0.00908 **
#scale(mPD)     -0.11574    0.16307     0.16540   0.700  0.48410   
#scale(FD)      -0.06474    0.20651     0.20952   0.309  0.75733   
#scale(PD)      -0.05752    0.20196     0.20500   0.281  0.77903   
#scale(CWM_H)    0.09048    0.18736     0.18972   0.477  0.63341   
#scale(CWM_SLA) -0.05578    0.11916     0.12087   0.461  0.64449   
#scale(CWM_LA)   0.04650    0.16623     0.16862   0.276  0.78275   
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#  scale(RaoQ) scale(mPD) scale(CWM_H) scale(FD) scale(CWM_SLA) scale(PD) scale(CWM_LA)
#Importance:          0.91        0.31       0.29         0.27      0.27           0.27      0.26         
#N containing models:   64          64         64           64        64             64        64 

summary(avg_models$P)
#Model-averaged coefficients:  
#  (full average) 
#              Estimate Std. Error Adjusted SE z value Pr(>|z|)
#(Intercept)     0.00000    0.00000     0.00000      NA       NA
#scale(CWM_H)    0.49294    0.34855     0.35029   1.407    0.159
#scale(RaoQ)    -0.22744    0.24049     0.24180   0.941    0.347
#scale(CWM_LA)  -0.17984    0.25638     0.25795   0.697    0.486
#scale(CWM_SLA)  0.07644    0.13867     0.13950   0.548    0.584
#scale(mPD)      0.03679    0.13164     0.13271   0.277    0.782
#scale(PD)       0.07343    0.21015     0.21158   0.347    0.729
#scale(FD)      -0.06259    0.20751     0.20895   0.300    0.765

#(conditional average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)  
#(Intercept)      0.0000     0.0000      0.0000      NA       NA  
#scale(CWM_H)     0.5588     0.3177      0.3198   1.747   0.0806 .
#scale(RaoQ)     -0.3485     0.2155      0.2177   1.601   0.1094  
#scale(CWM_LA)   -0.3385     0.2646      0.2674   1.266   0.2056  
#scale(CWM_SLA)   0.1785     0.1634      0.1650   1.082   0.2793  
#scale(mPD)       0.1147     0.2123      0.2144   0.535   0.5928  
#scale(PD)        0.2286     0.3194      0.3223   0.709   0.4781  
#scale(FD)       -0.2018     0.3328      0.3357   0.601   0.5477  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#                  scale(CWM_H) scale(RaoQ) scale(CWM_LA) scale(CWM_SLA) scale(PD) scale(mPD) scale(FD)
#Importance:          0.88         0.65        0.53          0.43           0.32      0.32       0.31     
#N containing models:   64           64          64            64             64        64         64  

summary(avg_models$K)

# Model-averaged coefficients:  
#  (full average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)
#(Intercept)     0.00000    0.00000     0.00000      NA       NA
#scale(CWM_H)   -0.49355    0.38725     0.38876   1.270    0.204
#scale(CWM_LA)   0.31553    0.29526     0.29684   1.063    0.288
#scale(mPD)     -0.06704    0.16769     0.16874   0.397    0.691
#scale(PD)       0.16652    0.24100     0.24265   0.686    0.493
#scale(RaoQ)     0.23180    0.23665     0.23784   0.975    0.330
#scale(CWM_SLA)  0.14102    0.16061     0.16140   0.874    0.382
#scale(FD)      -0.01313    0.19689     0.19861   0.066    0.947

#(conditional average) 
#Estimate Std. Error Adjusted SE z value Pr(>|z|)  
#(Intercept)     0.00000    0.00000     0.00000      NA       NA  
#scale(CWM_H)   -0.63311    0.32250     0.32483   1.949   0.0513 .
#scale(CWM_LA)   0.45802    0.24755     0.25028   1.830   0.0672 .
#scale(mPD)     -0.17139    0.23240     0.23433   0.731   0.4645  
#scale(PD)       0.30880    0.25253     0.25545   1.209   0.2267  
#scale(RaoQ)     0.35628    0.20427     0.20638   1.726   0.0843 .
#scale(CWM_SLA)  0.23907    0.14245     0.14395   1.661   0.0967 .
#scale(FD)      -0.03757    0.33161     0.33453   0.112   0.9106  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#                  scale(CWM_H) scale(CWM_LA) scale(RaoQ) scale(CWM_SLA) scale(PD) scale(mPD) scale(FD)
#Importance:          0.78         0.69          0.65        0.59           0.54      0.39       0.35     
#N containing models:   64           64            64          64             64        64         64

summary(avg_models$NO3)

#Model-averaged coefficients:  
#  (full average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)     0.00000    0.00000     0.00000      NA       NA    
#scale(CWM_H)    1.03370    0.25116     0.25399   4.070  4.7e-05 ***
#scale(CWM_LA)  -1.28153    0.19223     0.19510   6.569  < 2e-16 ***
#scale(CWM_SLA)  0.38187    0.13241     0.13376   2.855  0.00431 ** 
#scale(mPD)      0.12279    0.17542     0.17643   0.696  0.48643    
#scale(RaoQ)    -0.18983    0.20796     0.20920   0.907  0.36419    
#scale(PD)      -0.02385    0.11029     0.11160   0.214  0.83076    
#scale(FD)      -0.03185    0.11104     0.11234   0.284  0.77679    

# (conditional average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)     0.00000    0.00000     0.00000      NA       NA    
#scale(CWM_H)    1.03379    0.25099     0.25381   4.073 4.64e-05 ***
#scale(CWM_LA)  -1.28153    0.19223     0.19510   6.569  < 2e-16 ***
#scale(CWM_SLA)  0.38460    0.12887     0.13027   2.952  0.00315 ** 
#scale(mPD)      0.24212    0.17828     0.18023   1.343  0.17914    
#scale(RaoQ)    -0.30624    0.18471     0.18696   1.638  0.10143    
#scale(PD)      -0.08206    0.19255     0.19512   0.421  0.67406    
#scale(FD)      -0.10660    0.18248     0.18513   0.576  0.56474    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#                scale(CWM_LA) scale(CWM_H) scale(CWM_SLA) scale(RaoQ) scale(mPD) scale(FD) scale(PD)
#Importance:          1.00          1.00         0.99           0.62        0.51       0.30      0.29     
#N containing models:   64            64           64             64          64         64        64 

summary(avg_models$totbio)
# Model-averaged coefficients:  
#  (full average) 
#Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)     0.000000   0.000000    0.000000      NA       NA    
#scale(CWM_H)   -0.160510   0.201712    0.203160   0.790  0.42949    
#scale(CWM_LA)   0.596301   0.184642    0.186385   3.199  0.00138 ** 
#scale(CWM_SLA)  0.458606   0.098116    0.099428   4.612    4e-06 ***
#scale(mPD)     -0.016416   0.071075    0.071947   0.228  0.81952    
#scale(PD)      -0.014127   0.085048    0.086314   0.164  0.86999    
#scale(RaoQ)     0.014436   0.083736    0.084652   0.171  0.86459    
#scale(FD)      -0.007221   0.082996    0.084257   0.086  0.93170    

#(conditional average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)     0.00000    0.00000     0.00000      NA       NA    
#scale(CWM_H)   -0.29337    0.18812     0.19095   1.536  0.12445    
#scale(CWM_LA)   0.59666    0.18412     0.18587   3.210  0.00133 ** 
#scale(CWM_SLA)  0.45952    0.09605     0.09739   4.718  2.4e-06 ***
#scale(mPD)     -0.06217    0.12762     0.12946   0.480  0.63105    
#scale(PD)      -0.05568    0.16185     0.16447   0.339  0.73495    
#scale(RaoQ)     0.05359    0.15470     0.15654   0.342  0.73209    
#scale(FD)      -0.02895    0.16428     0.16684   0.174  0.86223    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#                 scale(CWM_LA) scale(CWM_SLA) scale(CWM_H) scale(RaoQ) scale(mPD) scale(PD) scale(FD)
#Importance:          1.00          1.00           0.55         0.27        0.26       0.25      0.25     
#N containing models:   64            64             64           64          64         64        64 

summary(avg_models$strat10)
#Model-averaged coefficients:  
#  (full average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)
#(Intercept)     0.00000    0.00000     0.00000      NA       NA
#scale(CWM_H)    0.46554    0.28589     0.28767   1.618    0.106
#scale(mPD)      0.24328    0.23360     0.23500   1.035    0.301
#scale(PD)      -0.47504    0.28869     0.29072   1.634    0.102
#scale(RaoQ)    -0.14405    0.21053     0.21174   0.680    0.496
#scale(CWM_LA)  -0.13689    0.21963     0.22113   0.619    0.536
#scale(CWM_SLA)  0.02198    0.09766     0.09846   0.223    0.823
#scale(FD)      -0.03529    0.22478     0.22653   0.156    0.876

# (conditional average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)  
#(Intercept)     0.00000    0.00000     0.00000      NA       NA  
#scale(CWM_H)    0.50156    0.26457     0.26663   1.881   0.0600 .
#scale(mPD)      0.35617    0.19920     0.20160   1.767   0.0773 .
#scale(PD)      -0.55012    0.23497     0.23785   2.313   0.0207 *
#scale(RaoQ)    -0.28317    0.21848     0.22077   1.283   0.1996  
#scale(CWM_LA)  -0.28187    0.24178     0.24459   1.152   0.2491  
#scale(CWM_SLA)  0.07034    0.16469     0.16619   0.423   0.6721  
#scale(FD)      -0.10176    0.37272     0.37576   0.271   0.7865  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#                scale(CWM_H) scale(PD) scale(mPD) scale(RaoQ) scale(CWM_LA) scale(FD) scale(CWM_SLA)
#Importance:          0.93         0.86      0.68       0.51        0.49          0.35      0.31          
#N containing models:   64           64        64         64          64            64        64   

summary(avg_models$crat10)
#Model-averaged coefficients:  
#  (full average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)
#(Intercept)     0.0000000  0.0000000   0.0000000      NA       NA
#scale(CWM_H)    0.1542155  0.1843150   0.1855544   0.831    0.406
#scale(CWM_LA)   0.1317712  0.1736562   0.1748414   0.754    0.451
#scale(RaoQ)    -0.0172786  0.0932008   0.0943372   0.183    0.855
#scale(FD)      -0.0275125  0.1228614   0.1242930   0.221    0.825
#scale(mPD)     -0.0020998  0.0736651   0.0747994   0.028    0.978
#scale(PD)       0.0124693  0.1207992   0.1222762   0.102    0.919
#scale(CWM_SLA)  0.0006861  0.0586786   0.0595416   0.012    0.991

#  (conditional average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)
#(Intercept)     0.000000   0.000000    0.000000      NA       NA
#scale(CWM_H)    0.268024   0.168936    0.171278   1.565    0.118
#scale(CWM_LA)   0.245950   0.167942    0.170222   1.445    0.148
#scale(RaoQ)    -0.062688   0.169317    0.171585   0.365    0.715
#scale(FD)      -0.102822   0.220614    0.223591   0.460    0.646
#scale(mPD)     -0.008567   0.148605    0.150899   0.057    0.955
#scale(PD)       0.049032   0.235772    0.238747   0.205    0.837
#scale(CWM_SLA)  0.002766   0.117806    0.119539   0.023    0.982

#Relative variable importance: 
#                  scale(CWM_H) scale(CWM_LA) scale(RaoQ) scale(FD) scale(PD) scale(CWM_SLA) scale(mPD)
# Importance:          0.58         0.54          0.28        0.27      0.25      0.25           0.25      
#N containing models:   64           64            64          64        64        64             64  

summary(avg_models$hits10)
#Model-averaged coefficients:  
#  (full average) 
#              Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)     0.00000    0.00000     0.00000      NA       NA    
#scale(CWM_H)    0.16432    0.23016     0.23150   0.710  0.47781    
#scale(CWM_LA)   0.55544    0.18991     0.19149   2.901  0.00372 ** 
#scale(CWM_SLA)  0.48828    0.10849     0.10971   4.451  8.6e-06 ***
#scale(RaoQ)    -0.24704    0.17445     0.17563   1.407  0.15956    
#scale(PD)       0.05251    0.14589     0.14702   0.357  0.72098    
#scale(FD)      -0.03059    0.13752     0.13866   0.221  0.82539    
#scale(mPD)     -0.04333    0.10828     0.10910   0.397  0.69127    

#(conditional average) 
#              Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)      0.0000     0.0000      0.0000      NA       NA    
#scale(CWM_H)     0.3217     0.2304      0.2330   1.381 0.167297    
#scale(CWM_LA)    0.5703     0.1690      0.1708   3.338 0.000843 ***
#scale(CWM_SLA)   0.4883     0.1085      0.1097   4.452  8.5e-06 ***
#scale(RaoQ)     -0.3063     0.1400      0.1418   2.160 0.030768 *  
#scale(PD)        0.1640     0.2195      0.2219   0.739 0.459844    
#scale(FD)       -0.1081     0.2417      0.2440   0.443 0.657852    
#scale(mPD)      -0.1234     0.1533      0.1550   0.796 0.425798    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#                 scale(CWM_SLA) scale(CWM_LA) scale(RaoQ) scale(CWM_H) scale(mPD) scale(PD) scale(FD)
#Importance:          1.00           0.97          0.81        0.51         0.35       0.32      0.28     
#N containing models:   64             64            64          64           64         64        64     

summary(avg_models$ind1)
#Model-averaged coefficients:  
#  (full average) 
#                Estimate Std. Error Adjusted SE z value Pr(>|z|)
#(Intercept)     0.00000    0.00000     0.00000      NA       NA
#scale(PD)       0.25693    0.25103     0.25290   1.016    0.310
#scale(mPD)     -0.06565    0.14214     0.14333   0.458    0.647
#scale(FD)       0.01338    0.20055     0.20235   0.066    0.947
#scale(CWM_LA)   0.06736    0.14879     0.14989   0.449    0.653
#scale(CWM_SLA)  0.01131    0.06291     0.06376   0.177    0.859
#scale(CWM_H)   -0.03826    0.13651     0.13766   0.278    0.781
#scale(RaoQ)     0.01011    0.08458     0.08569   0.118    0.906

#(conditional average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)
#(Intercept)     0.00000    0.00000     0.00000      NA       NA
#scale(PD)       0.36299    0.22479     0.22773   1.594    0.111
#scale(mPD)     -0.17395    0.18625     0.18866   0.922    0.357
#scale(FD)       0.03492    0.32287     0.32579   0.107    0.915
#scale(CWM_LA)   0.18501    0.19759     0.19985   0.926    0.355
#scale(CWM_SLA)  0.04345    0.11748     0.11923   0.364    0.716
#scale(CWM_H)   -0.12738    0.22513     0.22745   0.560    0.575
#scale(RaoQ)     0.03899    0.16268     0.16491   0.236    0.813

#Relative variable importance: 
#                 scale(PD) scale(FD) scale(mPD) scale(CWM_LA) scale(CWM_H) scale(CWM_SLA) scale(RaoQ)
#Importance:          0.71      0.38      0.38       0.36          0.30         0.26           0.26       
#N containing models:   64        64        64         64            64           64             64  

summary(avg_models$ind2)

#Model-averaged coefficients:  
#(full average) 
#              Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)     0.000000   0.000000    0.000000      NA       NA    
#scale(CWM_H)   -0.373057   0.269449    0.271393   1.375 0.169255    
#scale(CWM_LA)   0.780349   0.229692    0.231911   3.365 0.000766 ***
#scale(mPD)     -0.078894   0.139829    0.140780   0.560 0.575204    
#scale(RaoQ)     0.089924   0.153958    0.154984   0.580 0.561771    
#scale(CWM_SLA)  0.033578   0.087522    0.088257   0.380 0.703607    
#scale(PD)       0.007019   0.109474    0.110864   0.063 0.949517    
#scale(FD)       0.028077   0.114091    0.115424   0.243 0.807813    

#(conditional average) 
#               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#(Intercept)     0.00000    0.00000     0.00000      NA       NA    
#scale(CWM_H)   -0.46733    0.21655     0.21957   2.128 0.033305 *  
#scale(CWM_LA)   0.78065    0.22922     0.23145   3.373 0.000744 ***
#scale(mPD)     -0.18749    0.16157     0.16352   1.147 0.251563    
#scale(RaoQ)     0.20522    0.17445     0.17651   1.163 0.244964    
#scale(CWM_SLA)  0.09488    0.12581     0.12725   0.746 0.455927    
#scale(PD)       0.02648    0.21140     0.21411   0.124 0.901591    
#scale(FD)       0.09944    0.19750     0.20022   0.497 0.619436    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Relative variable importance: 
#                 scale(CWM_LA) scale(CWM_H) scale(RaoQ) scale(mPD) scale(CWM_SLA) scale(FD) scale(PD)
#Importance:          1.00          0.80         0.44        0.42       0.35           0.28      0.27     
#N containing models:   64            64           64          64         64             64        64 



### Linear model (re-evaluation) ####

# re-evaluate models using significant predictors from MuMIn analyses based on conditional averaging
# hard to use lapply() since there is no consistent set of predictors
lm_rev_mods <- list(lm(org ~ scale(RaoQ), data = fnn), #org
                    lm(P ~ scale(CWM_H), data = fnn), # P
                    lm(K ~ scale(CWM_H) + scale(CWM_LA) + scale(RaoQ) + scale(CWM_SLA), data = fnn), # K
                    lm(NO3 ~ scale(CWM_H) + scale(CWM_LA) + scale(CWM_SLA), data = fnn), # NO3
                    lm(totbio ~ scale(CWM_LA) + scale(CWM_SLA), data = fnn), # totbio
                    lm(strat10 ~ scale(CWM_H) + scale(mPD) + scale(PD), data = fnn), # strat10
                    # nothing significant for crat10 from MuMIn
                    lm(hits10 ~ scale(CWM_LA) + scale(CWM_SLA) + scale(RaoQ), data = fnn),
                    # nothing significant for ind1 from MuMIn
                    lm(ind2 ~ scale(CWM_H) + scale(CWM_LA), data = fnn)
                    )

names(lm_rev_mods) <- c("org", "P", "K", "NO3", "totbio", "strat10", "hits10", "ind2")

### Summary outputs for re-evaluated linear models

summary(lm_rev_mods$org)

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)   7.7276     0.1034  74.763  < 2e-16 ***
#  scale(RaoQ)   0.3502     0.1040   3.369  0.00114 ** 
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 0.9641 on 85 degrees of freedom
#Multiple R-squared:  0.1178,	Adjusted R-squared:  0.1074 
#F-statistic: 11.35 on 1 and 85 DF,  p-value: 0.001136

summary(lm_rev_mods$P)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   219.368      3.280   66.87   <2e-16 ***
#  scale(CWM_H)    6.697      3.299    2.03   0.0455 *  
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 30.6 on 85 degrees of freedom
# Multiple R-squared:  0.04623,	Adjusted R-squared:  0.03501 
# F-statistic:  4.12 on 1 and 85 DF,  p-value: 0.0455

## K

summary(lm_rev_mods$K)

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     219.046      3.984  54.978  < 2e-16 ***
#  scale(CWM_H)    -27.165      9.836  -2.762  0.00709 ** 
#  scale(CWM_LA)    18.964      8.498   2.232  0.02837 *  
#  scale(RaoQ)      12.945      5.928   2.184  0.03183 *  
#  scale(CWM_SLA)    5.167      4.762   1.085  0.28107    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 37.16 on 82 degrees of freedom
#Multiple R-squared:  0.1904,	Adjusted R-squared:  0.1509 
#F-statistic: 4.822 on 4 and 82 DF,  p-value: 0.001526

# model diagnostics
(plot(lm_rev_mods$K))


summary(lm_rev_mods$NO3)
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)     1.00345    0.08256  12.154  < 2e-16 ***
#  scale(CWM_H)    0.86113    0.17915   4.807 6.74e-06 ***
#  scale(CWM_LA)  -1.21365    0.17585  -6.901 9.50e-10 ***
#  scale(CWM_SLA)  0.29385    0.08929   3.291  0.00147 ** 
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 0.7701 on 83 degrees of freedom
#Multiple R-squared:  0.4105,	Adjusted R-squared:  0.3892 
#F-statistic: 19.26 on 3 and 83 DF,  p-value: 1.428e-09

summary(lm_rev_mods$totbio)
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)      42.034      2.005  20.962  < 2e-16 ***
#  scale(CWM_LA)    11.186      2.017   5.546 3.31e-07 ***
#  scale(CWM_SLA)   13.072      2.017   6.481 5.90e-09 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 18.7 on 84 degrees of freedom
# Multiple R-squared:  0.462,	Adjusted R-squared:  0.4492 
# F-statistic: 36.06 on 2 and 84 DF,  p-value: 4.942e-12

summary(lm_rev_mods$strat10)
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)    39.736      2.210  17.982  < 2e-16 ***
#  scale(CWM_H)    6.495      2.412   2.693 0.008568 ** 
#  scale(mPD)      5.604      3.445   1.626 0.107651    
#  scale(PD)     -13.126      3.487  -3.764 0.000311 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 20.61 on 83 degrees of freedom
# Multiple R-squared:  0.1881,	Adjusted R-squared:  0.1587 
# F-statistic: 6.409 on 3 and 83 DF,  p-value: 0.0005862

summary(lm_rev_mods$hits10)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)     139.092      5.890  23.615  < 2e-16 ***
#  scale(CWM_LA)    53.547      7.561   7.082 4.22e-10 ***
#  scale(CWM_SLA)   35.207      6.115   5.757 1.40e-07 ***
#  scale(RaoQ)     -19.594      7.703  -2.544   0.0128 *  
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 54.94 on 83 degrees of freedom
# Multiple R-squared:  0.5003,	Adjusted R-squared:  0.4823 
# F-statistic:  27.7 on 3 and 83 DF,  p-value: 1.642e-12

summary(lm_rev_mods$ind2)
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   -2.803e-16  3.226e-01   0.000    1.000    
#scale(CWM_H)  -1.496e+00  6.510e-01  -2.299    0.024 *  
#scale(CWM_LA)  3.033e+00  6.510e-01   4.660 1.18e-05 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 3.009 on 84 degrees of freedom
#Multiple R-squared:  0.2876,	Adjusted R-squared:  0.2706 
#F-statistic: 16.95 on 2 and 84 DF,  p-value: 6.526e-07


