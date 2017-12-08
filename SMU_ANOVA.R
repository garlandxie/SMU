## One-way ANOVA analyses ##

# Main question:
# does 15sp and 9sp treatment had higher MPD, FR, RaoQ and PD compared to 3sp?

summary(lm(mPD ~ factor(nosp), data = fnn))

# Residuals:
#   Min       1Q   Median       3Q      Max 
#-113.617  -30.741   -1.602   40.653   94.249 

#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)       27.31      11.74   2.326   0.0224 *  
#  factor(nosp)9     96.37      13.73   7.018 5.36e-10 ***
#  factor(nosp)15   124.49      16.19   7.691 2.51e-11 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 49.82 on 84 degrees of freedom
# Multiple R-squared:  0.445,	Adjusted R-squared:  0.4318 
# F-statistic: 33.68 on 2 and 84 DF,  p-value: 1.817e-11

summary(lm(RaoQ ~ factor(nosp), data = fnn))

# Residuals:
#  Min        1Q    Median        3Q       Max 
# -0.060440 -0.022428 -0.002835  0.022705  0.075163 

#Coefficients:
#           Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)    0.024180   0.007471   3.236  0.00173 ** 
#  factor(nosp)9  0.038758   0.008737   4.436 2.76e-05 ***
#  factor(nosp)15 0.058175   0.010299   5.649 2.15e-07 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.0317 on 84 degrees of freedom
# Multiple R-squared:  0.284,	Adjusted R-squared:  0.267 
# F-statistic: 16.66 on 2 and 84 DF,  p-value: 8.056e-07

summary(lm(FD ~ factor(nosp), data = fnn))

#  Residuals:
#  Min      1Q  Median      3Q     Max 
# -2.6357 -0.5745  0.1989  0.7497  1.5252 

#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)      2.8858     0.2289  12.607  < 2e-16 ***
#  factor(nosp)9    2.3145     0.2677   8.647 3.05e-13 ***
#  factor(nosp)15   3.7298     0.3155  11.821  < 2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.9712 on 84 degrees of freedom
# Multiple R-squared:  0.6292,	Adjusted R-squared:  0.6204 
# F-statistic: 71.27 on 2 and 84 DF,  p-value: < 2.2e-16

summary(lm(PD ~ factor(nosp), data = fnn))

#  Residuals:
#  Min       1Q   Median       3Q      Max 
# -176.411  -72.840   -0.955   82.728  183.668 

# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)      251.74      23.73  10.610  < 2e-16 ***
#  factor(nosp)9    266.55      27.74   9.607 3.56e-15 ***
#  factor(nosp)15   422.65      32.70  12.923  < 2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 100.7 on 84 degrees of freedom
# Multiple R-squared:  0.6707,	Adjusted R-squared:  0.6629 
# F-statistic: 85.56 on 2 and 84 DF,  p-value: < 2.2e-16
