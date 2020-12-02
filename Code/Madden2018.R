# Madden 2018 Rec 147
# Phenotypic extractions

# Load libraries====
library(lme4)
library(here)

# Set wd====
dir<-here()

# Load data====
# Sex: 1=M 2=F Fate:0=dead after 60 days, 1=alive after 60 days
Madden<-read.csv("Data/Madden2018.csv")

# Models====
m7<-glm(Fate~Exploratory_behaviour, family=binomial, data=Madden)
summary(m7)

MaddenF<-subset(Madden, Sex==2)
MaddenM<-subset(Madden, Sex==1)

m8<-glm(Fate~Exploratory_behaviour, family=binomial, data=MaddenF)
summary(m8)
table(MaddenF$Fate)
cor.test(MaddenF$Fate, MaddenF$Exploratory_behaviour)

m9<-glm(Fate~Exploratory_behaviour, family=binomial, data=MaddenM)
summary(m9)
table(MaddenM$Fate)
cor.test(MaddenM$Exploratory_behaviour, MaddenM$Fate) #r=0.2602139