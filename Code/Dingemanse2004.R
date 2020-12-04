# Dingemanse 2004
# data extracted from figure for phenotypic estimate

# Load libraries====
library(here)
library(MuMIn)

# Set wd====
dir<-here()

#Data loading ====
exploration<-read.csv("Data/Dingemanse2004.csv")
str(exploration)

scatter.smooth(x=exploration$exploration, y=exploration$recruits, main="exploration ~ number of recruits")  # scatterplot

linearMod <- lm(recruits ~ exploration, data=exploration)  # build linear regression model on full data
print(linearMod)
summary(linearMod)
linearMod
r.squaredGLMM(linearMod)

cor(exploration$exploration,exploration$recruits, method="spearman")