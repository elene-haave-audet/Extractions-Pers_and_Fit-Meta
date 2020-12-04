# Niemela 2019 rec 328
# Phenotypic extractions

# Load libraries====
library(here)
library(lme4)
library(rptR)

# Set wd====
dir<-here()

# Load data====
Nim2019AH<-read.csv("Data/Niemela2019_High_Adult.csv")
Nim2019AL<-read.csv("Data/Niemela2019_Low_Adult.csv")
Nim2019JH<-read.csv("Data/Niemela2019_High_juvenile.csv")
Nim2019JL<-read.csv("Data/Niemela2019_Low_juvenile.csv")

hist(Nim2019AH$longevity) #use poisson glmer
hist(Nim2019AL$longevity)
hist(Nim2019JH$longevity)
hist(Nim2019JL$longevity)
hist(Nim2019AH$Exploration) #all quite zero inflated
hist(Nim2019AL$Exploration)
hist(Nim2019JH$Exploration)
hist(Nim2019JL$Exploration)

# Models====
m1<-glmer(longevity~scale(Exploration)+ (1|CricketID), family="poisson", data=Nim2019AH)
summary(m1)
rpt(Exploration~(1|CricketID), grname = "CricketID", data=Nim2019AH)
#r=0.467

m2<-glmer(longevity~scale(Exploration)+ (1|CricketID), family="poisson", data=Nim2019AL)
summary(m2)
rpt(Exploration~(1|CricketID), grname = "CricketID", data=Nim2019AL)
#r= 0.454

m3<-glmer(longevity~scale(Exploration)+ (1|CricketID), family="poisson", data=Nim2019JH)
summary(m3)
rpt(Exploration~(1|CricketID), grname = "CricketID", data=Nim2019JH)
#r=0.558

m4<-glmer(longevity~scale(Exploration)+ (1|CricketID), family="poisson", data=Nim2019JL)
summary(m4)
rpt(Exploration~(1|CricketID), grname = "CricketID", data=Nim2019JL)
#r=0.384