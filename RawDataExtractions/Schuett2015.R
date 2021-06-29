# Schuett 2015 Rec 226
# Phenotypic extractions

# Load libraries====
library(here)
library(lme4)

# Set wd====
dir<-here()

# Load data====
schuett<-read.csv("Data/Schuett2015.csv")
str(schuett)
hist(schuett$Lifespan) #BEAUTIFUL!

schuett_T<-subset(schuett, Treatment=="T")
hist(schuett_T$Lifespan)

schuett_C<-subset(schuett, Treatment=="C")
hist(schuett_C$Lifespan)

m1<-lmer(Lifespan~ BehavType + (1|ID), data=schuett_T)
plot(m1)
summary(m1)

m2<-lmer(Lifespan~ BehavType + (1|ID), data=schuett_C)
plot(m2)
summary(m2)