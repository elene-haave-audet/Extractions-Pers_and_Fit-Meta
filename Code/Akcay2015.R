# Akcay 2015 Rec #10
# Phenotypic extractions

# Load libraries====
library(here)
library(lme4)

# Set wd====
dir<-here()

# Load data====
Ak15<-read.csv("Data/Akcay2015.csv")
str(Ak15)

# Visualize & transform data====
hist(Ak15$Lsrate_t) #Poison
Ak15$Lsrate_tt<-log((Ak15$Lsrate_t)+0.5)
hist(Ak15$Lsrate_tt) #not particularly normal, but best I can do

hist(Ak15$Wwrate_t) #zero inflated Poison
Ak15$Wwrate_tt<-log((Ak15$Wwrate_t)+0.5)
hist(Ak15$Wwrate_tt) #slightly better, still zero inflated

hist(Ak15$Ssrate_t)
Ak15$Ssrate_tt<-log((Ak15$Ssrate_t)+0.5)
hist(Ak15$Ssrate_tt)

hist(Ak15$Flight_t)
Ak15$Flight_t<-sqrt(Ak15$Flight)
hist(Ak15$Flight_t) #looking normal-ish

hist(Ak15$TimeW5)
Ak15$TimeW5_t<-sqrt(Ak15$TimeW5)
hist(Ak15$TimeW5_t) #still looks bimodal

hist(Ak15$Closest)
Ak15$Closest_t<-log((Ak15$Closest)+0.1)
hist(Ak15$Closest_t) #not great, but what's new?

# Models====
m1<-glmer(survivedTo2014~ scale(Lsrate_tt) + (1|Bird), family=binomial, data=Ak15)
plot(m1)
summary(m1)
rpt(Lsrate~(1|Bird), grname = "Bird", data=Ak15)
#r=0.262

m2<-glmer(survivedTo2014~ scale(Wwrate_tt) + (1|Bird), family=binomial, data=Ak15)
plot(m2)
summary(m2)
rpt(Wwrate~(1|Bird), grname = "Bird", data=Ak15)
#r=0.508

m3<-glmer(survivedTo2014~ scale(Ssrate_tt) + (1|Bird), family=binomial, data=Ak15)
plot(m3)
summary(m3)
rpt(Ssrate~(1|Bird), grname = "Bird", data=Ak15)
#r=0.3

m4<-glmer(survivedTo2014~ scale(Flight_t) + (1|Bird), family=binomial, data=Ak15)
plot(m4)
summary(m4)
rpt(Flight~(1|Bird), grname = "Bird", data=Ak15)
#r=0.426

m5<-glmer(survivedTo2014~ scale(TimeW5_t) + (1|Bird), family=binomial, data=Ak15)
plot(m5)
summary(m5)
rpt(TimeW5~(1|Bird), grname = "Bird", data=Ak15)
#r=0.405

m6<-glmer(survivedTo2014~ scale(Closest_t) + (1|Bird), family=binomial, data=Ak15)
plot(m6)
summary(m6) #reverse sign to such that near distance reflects higher aggression
rpt(Closest~(1|Bird), grname = "Bird", data=Ak15)
#r=0.278