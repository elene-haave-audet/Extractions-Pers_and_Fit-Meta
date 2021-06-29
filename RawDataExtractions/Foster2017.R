# Foster 2017 rec 84
# Phenotypic extractions

# Load libraries====
library(here)

# Set wd====
dir<-here()

# Load data=====
Fost<-read.csv("Data/Foster2017.csv")

cor.test(Fost$Survived,Fost$Average_Time_Out)
#r=-0.1420765 high time=low boldness: need to switch sign
cor.test(Fost$Survived,Fost$Average_Height_cm)
#r=-0.1337131 (switch sign to reflect boldness)