################ SMU PD-FD ANALYSIS #####################

### Import libraries ####

library(broom) # cleaning data by taking lm model summaries into dataframes

### Metadata ####

### SMU_FDPD_MuMin.RData ###

## cwm_mod
# community weighted mean values for four leaf traits (LA, SLA, SDMC, H) of each green roof modules

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
# NO3 - service indicator for substrate nitrate content
# totbio - service indicator for total community biomass
# strat10 - service indicator for soil temperature index
# crat10 - service indicator for stormwater capture index 
# ind 1 - multifunctionality index 1 (see below for more details)
# ind 2 - multifunctionality index 2 (see below for more details)

## models_SFP
# a list of lm objects 
# $org: lm(org ~ FD + SR + PD)
# $K: lm(K ~ FD + SR + PD, data = fnn)
# $P: lm(P ~ FD + SR + PD, data = fnn)
# $NO3: lm(NO3 ~ FD + SR + PD, data = fnn)
# $totbio: lm(totbio ~ FD + SR + PD, data = fnn)
# $strat10: lm(strat10 ~ FD + SR + PD, data = fnn)
# $crat10: lm(crat10 ~ FD + SR + PD, data = fnn)
# $hits10: lm(hits10 ~ FD + SR + PD, data = fnn)
# $ind1: lm(ind1 ~ FD + SR + PD, data = fnn)
# $ind2: lm(ind2 ~ FD + SR + PD, data = fnn) 

### Import RData files ####

# imported from Google Drive
load("C:/Users/garla/Google Drive/Research Projects/SMU_12sp_Phylogeny/SMU_PD-FD_Data/SMU_Multi_Regress_SFP.RData")

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

### Multiple Linear Regression: FD, SR, PD and service indicators ####

# first, cram all of the possible multiple regression models into a list object
# reason: less memory space, easier access and a clean environment workspace

models_SFP <- list(lm(log(org) ~ FD + SR + PD, data = fnn), # org
                   lm(K ~ FD + SR + PD, data = fnn), # K
                   lm(P ~ FD + SR + PD, data = fnn), # P
                   lm(NO3 ~ FD + SR + PD, data = fnn), # NO3
                   lm(totbio ~ FD + SR + PD, data = fnn), #totbio
                   lm(strat10 ~ FD + SR + PD, data = fnn), # strat10
                   lm(crat10 ~ FD + SR + PD, data = fnn), # crat10
                   lm(hits10 ~ FD + SR + PD, data = fnn), # hits10
                   lm(ind1 ~ FD + SR + PD, data = fnn), # ind1
                   lm(ind2 ~ FD + SR + PD, data = fnn) # ind2
                   )
                   
# next, create names for each lm object in the list so that it can be indexed quickly
# all linear models have consistent set of predictor vars; using service indicators is fine 
# make sure it's named in the correct order to service indicators 

names(models_SFP) <- c("org", "K", "P", "NO3", "totbio", 
                   "strat10", "crat10", "hits10", "ind1", "ind2")
                   

## org - linear model input
summary(models_SFP$org)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -0.28664 -0.07719 -0.00563  0.05687  0.50639 

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  1.993e+00  5.271e-02  37.813   <2e-16 ***
#  FD          -1.039e-02  3.842e-02  -0.271    0.787    
#  SR           1.798e-02  3.152e-02   0.571    0.570    
#  PD           1.637e-05  2.263e-04   0.072    0.943    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 0.1247 on 83 degrees of freedom
#Multiple R-squared:  0.0275,	Adjusted R-squared:  -0.007649 
#F-statistic: 0.7824 on 3 and 83 DF,  p-value: 0.5071

## model diagnostics for multiple regression - lm_org2
(plot(models_SFP$org))

# Residuals vs Fitted: looks fine
# QQ plot: skewness on the right side of plot, a couple of outliers
# scale-location: some evidence of heterodaskicity (decrease then goes a bit flat)
# residuals vs leverage: high residual, low leverage 

## model validity - histogram of residuals from lm_org2 model
(hist(resid(models_SFP$org)))

# right-skewness, departure from normality

## K - linear model output
summary(models_SFP$K)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -67.780 -27.131   2.444  22.572 148.376 

#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 196.148983  16.958672  11.566   <2e-16 ***
#  FD            0.471022  12.362472   0.038    0.970    
#  SR            4.437733  10.140422   0.438    0.663    
#  PD           -0.002527   0.072821  -0.035    0.972    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 40.11 on 83 degrees of freedom
# Multiple R-squared:  0.04541,	Adjusted R-squared:  0.0109 
# F-statistic: 1.316 on 3 and 83 DF,  p-value: 0.2747

## model diagnostics for multiple regression - lm_K2
(plot(models_SFP$K))

# residual vs fitted: looks fine
# QQ plot: skewed distribution
# scale-location: sloped lowess curve, evidence of heteroskedasticity
# resid vs lev: looks OK?

## P - linear model output

summary(models_SFP$P)
# Residuals:
#  Min      1Q  Median      3Q     Max 
# -76.223 -22.532  -2.106  23.375  85.907 

#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  219.368      3.310   66.27   <2e-16 ***
#  scale(FD)     -7.946     15.001   -0.53   0.5977    
#  scale(SR)    -11.945     14.562   -0.82   0.4144    
#  scale(PD)     20.415      9.720    2.10   0.0387 *  
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 30.88 on 83 degrees of freedom
# Multiple R-squared:  0.05147,	Adjusted R-squared:  0.01718 
# F-statistic: 1.501 on 3 and 83 DF,  p-value: 0.2203

## model diagnostics for multiple regression - lm_P2
(plot(models_SFP$P))

## NO3 - linear model output

summary(models_SFP$NO3)
# Residuals:
#  Min      1Q  Median      3Q     Max 
# -1.2401 -0.4996 -0.1779  0.2398  6.1556 

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)   1.0034     0.1008   9.960 7.91e-16 ***
#  scale(FD)     1.1943     0.4565   2.616   0.0106 *  
#  scale(SR)    -1.0718     0.4432  -2.418   0.0178 *  
#  scale(PD)    -0.3299     0.2958  -1.115   0.2679    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.9397 on 83 degrees of freedom
# Multiple R-squared:  0.1221,	Adjusted R-squared:  0.09032 
# F-statistic: 3.846 on 3 and 83 DF,  p-value: 0.01246

# model diagnostics for multiple regression - lm_NO3
(plot(models_SFP$NO3))

# Residuals vs Fitted: possible evidence of nonconstant variance but it looks fine
# QQ plot: departure from normality; skewness
# scale-location: some evidence of heterodaskicity (decrease then goes a bit flat)
# residuals vs leverage: high residual, low leverage 

# NOTES:
# FD is a stronger predictor than SR (based on effect sizes)
# multicollinearity is still an issue though

## totbio - linear model output

summary(models_SFP$totbio)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -44.757 -21.481  -3.625  20.371  47.261 

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  42.0345     2.7144  15.486   <2e-16 ***
#  scale(FD)    -0.7427    12.2996  -0.060    0.952    
#  scale(SR)     8.3832    11.9401   0.702    0.485    
#  scale(PD)    -4.2394     7.9693  -0.532    0.596    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 25.32 on 83 degrees of freedom
# Multiple R-squared:  0.0259,	Adjusted R-squared:  -0.009307 
# F-statistic: 0.7357 on 3 and 83 DF,  p-value: 0.5337

## model diagnostics for multiple regression - lm_NO3
(plot(models_SFP$totbio))


## strat10 - linear model output

summary(models_SFP$strat10)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -45.187 -17.711  -0.831  18.002  39.030 

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  39.7356     2.3449  16.946   <2e-16 ***
#  scale(FD)     0.3144    10.6253   0.030    0.976    
#  scale(SR)    -4.2642    10.3148  -0.413    0.680    
#  scale(PD)    -2.7399     6.8845  -0.398    0.692    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 21.87 on 83 degrees of freedom
# Multiple R-squared:  0.08578,	Adjusted R-squared:  0.05274 
# F-statistic: 2.596 on 3 and 83 DF,  p-value: 0.05786

## model diagnostics for multiple regression - lm_strat10
(plot(lm_strat10_SFP)) 

# looks fine

## crat10 - linear model output

summary(models_SFP$crat10)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -0.18242 -0.08022 -0.02893  0.05051  0.31634 

# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  1.0399332  0.0467839  22.228   <2e-16 ***
#  FD          -0.0222292  0.0341044  -0.652    0.516    
#  SR          -0.0041769  0.0279744  -0.149    0.882    
#  PD           0.0002839  0.0002009   1.413    0.161    
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.1107 on 83 degrees of freedom
# Multiple R-squared:  0.02507,	Adjusted R-squared:  -0.01017 
# F-statistic: 0.7115 on 3 and 83 DF,  p-value: 0.5478

## model diagnostics for multiple regression - lm_crat10
(plot(models_SFP$crat10)) 

# looks reasonable


## hits10 - linear model output
summary(models_SFP$hits10)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -107.293  -53.857   -8.456   33.147  243.657 

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 115.7724    32.3396   3.580 0.000578 ***
#  FD          -12.0867    23.5748  -0.513 0.609528    
#  SR           -0.7718    19.3375  -0.040 0.968258    
#  PD            0.1766     0.1389   1.271 0.207127    
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 76.49 on 83 degrees of freedom
# Multiple R-squared:  0.03142,	Adjusted R-squared:  -0.003586 
# F-statistic: 0.8976 on 3 and 83 DF,  p-value: 0.4461

## model diagnostics for multiple regression - lm_hits10
(plot(models_SFP$hits10)) 


# residuals vs fitted: looks fine
# QQ plot: skewed distribution (might not be a big deal?)
# scale-location: sloped lowess curve, heteroscedastic 
# cook's distance: looks fine

## ind1 - linear model output

summary(models_SFP$ind1)

# Residuals:
#  Min       1Q   Median       3Q      Max 
#  -2.49449 -1.08294  0.03093  0.85238  3.07845 

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.887007   0.573487  -1.547    0.126
# FD          -0.210788   0.418059  -0.504    0.615
# SR           0.063784   0.342916   0.186    0.853
# PD           0.003282   0.002463   1.333    0.186

# Residual standard error: 1.356 on 83 degrees of freedom
# Multiple R-squared:  0.07507,	Adjusted R-squared:  0.04164 
# F-statistic: 2.246 on 3 and 83 DF,  p-value: 0.08905

## ind2 - linear model output
summary(models_SFP$ind2)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -13.4926  -2.1399  -0.0446   2.2294  10.1548 

# Coefficients:
#           Estimate Std. Error t value Pr(>|t|)
# (Intercept) -1.087068   1.473597  -0.738    0.463
# FD          -1.046183   1.074217  -0.974    0.333
# SR           1.032701   0.881136   1.172    0.245
# PD           0.002602   0.006328   0.411    0.682

# Residual standard error: 3.485 on 83 degrees of freedom
# Multiple R-squared:  0.05549,	Adjusted R-squared:  0.02135 
# F-statistic: 1.625 on 3 and 83 DF,  p-value: 0.1897

## model diagnostics for multiple regression - lm_ind2
(plot(models_SFP$ind2))

# everything looks reasonable

### Coercing model summary to dataframes ####

# since all of the lm model objects are in a single list, 
# you can create a dataframe of the model outputs by looping through each lm object, 
# then applying the tidy function to create the dataframe
all_coefs_SFP <- plyr::ldply(models_SFP, tidy, .id = "model")

# change the working directory to save your output
setwd("C:/Users/garla/Google Drive/Research Projects/SMU_12sp_Phylogeny/SMU_PD-FD_Data/SMU_CSV")

# coerce the data-frame into a csv file 
write.csv(all_coefs_SFP, "appendix_lm_SFP_table.csv")


