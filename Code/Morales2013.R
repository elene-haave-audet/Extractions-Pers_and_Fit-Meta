# Morales 2013 rec 325
# Phenotypic extractions

# Load libraries====
library(here)
library(lme4)

# Set wd====
dir<-here()

# Load data====
Moral<-read.csv("Data/Morales2013.csv")
MorF<-subset(Moral, SEX=="F")
MorM<-subset(Moral, SEX=="M")

hist(Moral$Dead) #poisson

# Models====
m1<-glm(Dead~Distance_flown, family="poisson", data=MorF)
summary(m1)
cor.test(Moral$Dead,Moral$Distance_flown)

m2<-glm(Dead~distance_walked, family="poisson", data=MorF)
summary(m2)

m3<-glm(Dead~resting_time, family="poisson", data=MorF)
summary(m3) #switch sign so high value represents high exploration

m4<-glm(Dead~walking_velocity, family="poisson", data=MorF)
summary(m4)

m5<-glm(Dead~free_fall, family="poisson", data=MorF)
summary(m5) 

m6<-glm(Dead~death_feigning, family="poisson", data=MorF)
summary(m6) #keep sign to represent antipredator behaviour

m7<-glm(Dead~body_righting, family="poisson", data=MorF)
summary(m7) #high value reflects long latency, so switch sign

m8<-glm(Dead~intra_interaction, family="poisson", data=MorF)
summary(m8) #high value reflects long latency, therefore switch sign

m9<-glm(Dead~inter_interaction, family="poisson", data=MorF)
summary(m9) #switch sign as above

m10<-glm(Dead~Wingbeat, family="poisson", data=MorF)
summary(m10)

m11<-glm(Dead~Distance_flown, family="poisson", data=MorM)
summary(m11)

m12<-glm(Dead~distance_walked, family="poisson", data=MorM)
summary(m12)

m13<-glm(Dead~resting_time, family="poisson", data=MorM)
summary(m13) #switch sign so high value represents high exploration

m14<-glm(Dead~walking_velocity, family="poisson", data=MorM)
summary(m14)

m15<-glm(Dead~free_fall, family="poisson", data=MorM)
summary(m15) 

m16<-glm(Dead~death_feigning, family="poisson", data=MorM)
summary(m16) #high value reflects long latency, so switch sign

m17<-glm(Dead~body_righting, family="poisson", data=MorM)
summary(m17) #high value reflects long latency, so switch sign

m18<-glm(Dead~intra_interaction, family="poisson", data=MorM)
summary(m18) #high value reflects long latency, therefore switch sign

m19<-glm(Dead~inter_interaction, family="poisson", data=MorM)
summary(m19) #switch sign as above

m20<-glm(Dead~Wingbeat, family="poisson", data=MorM)
summary(m20)