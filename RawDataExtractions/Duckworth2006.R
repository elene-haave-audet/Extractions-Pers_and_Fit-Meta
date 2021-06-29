# Duckworth 2006

# Set-up----
# Load libraries====
library(here)

# Set wd====
dir<-here()

data<-read.csv("Data/Duckworth2006.csv")

# Extract estimate----
m1<-glm(Repro~Agg, family="poisson", data=data)

cor.test(data$Agg, data$Repro)
#r=-0.7347