# Riviera-Guitierrez 2017

# Libraries====
library(here)
library(MuMIn)

# Set wd====
dir<-here()

# Data====
riviera<-read.csv("Data/#215 Riviera Gutierrez 2017.csv")
str(riviera)

scatter.smooth(x=riviera$exploration, y=riviera$clutchsize, main="exploration ~ clutch size")  # scatterplot

linearMod <- lm(clutchsize ~ exploration, data=riviera)  # build linear regression model on full data
print(linearMod)
summary(linearMod)
r.squaredGLMM(linearMod)

#Call:
#  lm(formula = clutchsize ~ exploration, data = riviera)
#
#Residuals:
#  Min     1Q Median     3Q    Max 
#-3.415 -1.201 -0.244  0.699  4.670 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  7.98726    0.76863  10.392 6.19e-11 ***
#  exploration  0.01426    0.03087   0.462    0.648    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 1.717 on 27 degrees of freedom
#Multiple R-squared:  0.007842,	Adjusted R-squared:  -0.0289 
#F-statistic: 0.2134 on 1 and 27 DF,  p-value: 0.6478