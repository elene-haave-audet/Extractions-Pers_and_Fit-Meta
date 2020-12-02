# Belgrad 2016 Rec 26
# Phenotypic extractions

# Load libraries====
library(here)

# Set wd====
dir<-here()

Belgrad2016<-read.csv("Data/Belgrad2016.csv")

Bel2016F<-subset(Belgrad2016, gender=="f")
Bel2016M<-subset(Belgrad2016, gender=="m")

cor.test(Bel2016M$time_in_refuge, Bel2016M$survival_days) #r=0.08370364 (switch sign, because time in refuge is inverse of boldness)
cor.test(Bel2016F$time_in_refuge, Bel2016F$survival_days) #0.1312222 (switched sign to represent inverse of time spent hiding)
