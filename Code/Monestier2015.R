# Monestier 2015, Rec 167

# Load libraries====
library(here)
library(lme4)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/Monestier2015.csv")
data$Proactivity<-data$ï..Proactivity

# Get estimates====
cor.test(data$Proactivity, data$Survival) #-0.2437

m1<-glm(Survival~Proactivity, data=data, family = "binomial")
summary(m1) #-0.6258, SE=0.3504
