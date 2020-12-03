# Belgrad 2018 Rec 27
# Phenotypic extractions

# Load libraries====
library(here)

# Set wd====
dir<-here()

# Load data====
Belgrad<-read.csv("Data/Belgrad2018.csv")

# note from authors:
# behaviour 1= extremely bold, 0= extremely shy
BelF<-subset(Belgrad, gender=="f")
cor.test(BelF$behavior, BelF$censor) # high values= high boldness, sign does not need to be switched

BelM<-subset(Belgrad, gender=="m")
cor.test(BelM$behavior, BelM$censor) # high values= high boldness, sign does not need to be switched
