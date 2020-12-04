# Yoshida 2016
# Data from figs 5 and 6

# libraries====
library(lme4)
library(MuMIn)
library(here)

# Set wd====
dir<-here()

# Figure 5====
yoshida1<-read.csv("Data/yoshida 2016 fig 5.csv")
str(yoshida1)

scatter.smooth(x=yoshida1$Boldness, y=yoshida1$Longevity, main="boldness ~ longevity")  # scatterplot

linearMod <- lm(Longevity ~ Boldness, data=yoshida1)  # build linear regression model on full data
print(linearMod)
summary(linearMod)
r.squaredGLMM(linearMod)

cor(yoshida1$Boldness,yoshida1$Longevity, method="spearman")

# Figure 6====
yoshida2<-read.csv("Data/yoshida 2016 fig 6.csv")
str(yoshida2)

scatter.smooth(x=yoshida2$Group_joining, y=yoshida2$Longevity, main="sociability ~ longevity")  # scatterplot

linearMod <- lm( Longevity~Group_joining, data=yoshida2)  # build linear regression model on full data
print(linearMod)
summary(linearMod)
r.squaredGLMM(linearMod)

cor(yoshida2$Group_joining,yoshida2$Longevity, method="spearman")