
################ SMU PD-FD ANALYSIS #####################

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

### Import RData files ####

# imported from Google Drive
load("C:\\Users\\garla\\Google Drive\\Research Projects\\SMU_12sp_Phylogeny\\SMU_PD-FD_Data\\SMU_FDPD_MuMin.RData")

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

## org - linear model input
lm_org2 <- lm(org ~ FD + SR + PD, data = fnn)

#  Residuals:
#  Min      1Q  Median      3Q     Max 
# -1.9916 -0.6233 -0.0954  0.4113  4.9797 

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)  7.4402343  0.4353786  17.089   <2e-16 ***
#  FD          -0.0403712  0.3173807  -0.127    0.899    
#  SR           0.1245129  0.2603342   0.478    0.634    
#  PD          -0.0002405  0.0018695  -0.129    0.898    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 1.03 on 83 degrees of freedom
# Multiple R-squared:  0.01724,	Adjusted R-squared:  -0.01828 
# F-statistic: 0.4855 on 3 and 83 DF,  p-value: 0.6933

## model diagnostics for multiple regression - lm_org2
(plot(lm_org2))

# Residuals vs Fitted: looks fine
# QQ plot: skewness on the right side of plot
# scale-location: some evidence of heterodaskicity (decrease then goes a bit flat)
# residuals vs leverage: high residual, low leverage 

## model validity - histogram of residuals from lm_org2 model
(hist(lm_org2$residuals))

# right-skewness

# K
lm_K2 <- lm(K ~ FD + SR + PD, data = fnn)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -67.780 -27.131   2.444  22.572 148.376 

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 196.148983  16.958672  11.566   <2e-16 ***
# FD            0.471022  12.362472   0.038    0.970    
# SR            4.437733  10.140422   0.438    0.663    
# PD           -0.002527   0.072821  -0.035    0.972    
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 40.11 on 83 degrees of freedom
# Multiple R-squared:  0.04541,	Adjusted R-squared:  0.0109 
# F-statistic: 1.316 on 3 and 83 DF,  p-value: 0.2747

# P 
lm_P2 <- lm(P ~ FD + SR + PD, data = fnn)

# Residuals:
#  Min      1Q  Median      3Q     Max 
#-76.223 -22.532  -2.106  23.375  85.907 

#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#  (Intercept) 217.47849   13.05540   16.66   <2e-16 ***
#  FD           -5.04100    9.51708   -0.53   0.5977    
#  SR           -6.40348    7.80646   -0.82   0.4144    
#  PD            0.11775    0.05606    2.10   0.0387 *  
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 30.88 on 83 degrees of freedom
# Multiple R-squared:  0.05147,	Adjusted R-squared:  0.01718 
# F-statistic: 1.501 on 3 and 83 DF,  p-value: 0.2203

# NO3
lm_NO3 <- lm(NO3 ~ FD + SR + PD, data = fnn)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -1.2401 -0.4996 -0.1779  0.2398  6.1556 

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)  
#  (Intercept)  0.949022   0.397329   2.389   0.0192 *
#  FD           0.757690   0.289644   2.616   0.0106 *
#  SR          -0.574550   0.237583  -2.418   0.0178 *
#  PD          -0.001903   0.001706  -1.115   0.2679  
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.9397 on 83 degrees of freedom
# Multiple R-squared:  0.1221,	Adjusted R-squared:  0.09032 
# F-statistic: 3.846 on 3 and 83 DF,  p-value: 0.01246

# NOTES: only super interesting result here - FD is a stronger predictor than SR (based on effect sizes)

# totbio
lm_totbio2 <- lm(totbio ~ FD + SR + PD, data = fnn)

# Residuals:
#  Min      1Q  Median      3Q     Max 
#-44.757 -21.481  -3.625  20.371  47.261 

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)   
# ( Intercept) 34.55796   10.70452   3.228  0.00178 **
#  FD          -0.47120    7.80334  -0.060  0.95199   
#  SR           4.49398    6.40076   0.702  0.48458   
#  PD          -0.02445    0.04597  -0.532  0.59617   
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 25.32 on 83 degrees of freedom
# Multiple R-squared:  0.0259,	Adjusted R-squared:  -0.009307 
# F-statistic: 0.7357 on 3 and 83 DF,  p-value: 0.5337

# strat10
lm_strat10_2 <- lm(strat10 ~ FD + SR + PD, data = fnn)

# Residuals:
#  Min      1Q  Median      3Q     Max 
# -45.187 -17.711  -0.831  18.002  39.030 

# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#  (Intercept) 57.83451    9.24734   6.254 1.66e-08 ***
#  FD           0.19944    6.74109   0.030    0.976    
#  SR          -2.28592    5.52944  -0.413    0.680    
#  PD          -0.01580    0.03971  -0.398    0.692    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 21.87 on 83 degrees of freedom
# Multiple R-squared:  0.08578,	Adjusted R-squared:  0.05274 
# F-statistic: 2.596 on 3 and 83 DF,  p-value: 0.05786

# crat10
lm_crat10_2 <- lm(crat10 ~ FD + SR + PD, data = fnn)

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

# hits10
lm_hits10_2 <- lm(hits10 ~ FD + SR + PD, data = fnn)

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

# ind1
lm_ind1_2 <- lm(ind1 ~ FD + SR + PD, data = fnn)

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

# ind2
lm_ind2_2 <- lm(ind2 ~ FD + SR + PD, data = fnn)

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
