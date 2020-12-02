# Belgrad 2018 Rec 27
# Phenotypic extractions

# Load libraries====
library(here)

# Set wd====
dir<-here()

# Load data====
Belgrad<-read.csv("Data/Belgrad2018.csv")

BelF<-subset(Belgrad, gender=="f")
cor.test(BelF$behavior, BelF$censor) #r=-0.05060336 (switch sign to reflect boldness)

BelM<-subset(Belgrad, gender=="m")
cor.test(BelM$behavior, BelM$censor) #r=-0.1379683 (switch sign to reflect boldness)
