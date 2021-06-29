# Penteriani 2002
# data from figure 1

# Libraries====
library(here)

# Set wd====
dir<-here()

# Data loading ====
penteriani<-read.csv("Data/penteriani 2002 Fig 1.csv")
str(penteriani)

scatter.smooth(x=penteriani$duration_vocalisation, y=penteriani$nb_fledged_young, main="vocalisation ~ productivity")  # scatterplot

linearMod <- lm(duration_vocalisation ~ nb_fledged_young, data=penteriani)  # build linear regression model on full data
print(linearMod)
summary(linearMod)

cor(penteriani$duration_vocalisation, penteriani$nb_fledged_young, method="spearman")
#0.4521959

#Call:
#lm(formula = duration_vocalisation ~ nb_fledged_young, data = penteriani)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-1059.1  -668.2  -554.0   -29.0  3535.2 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)         796.1      707.0   1.126    0.297
#nb_fledged_young    192.7      452.2   0.426    0.683

#Residual standard error: 1492 on 7 degrees of freedom
#Multiple R-squared:  0.02528,	Adjusted R-squared:  -0.114 
#F-statistic: 0.1816 on 1 and 7 DF,  p-value: 0.6828