# Archie 2014 Rec 317
# Phenotypic extractions

# Load libraries====
library(rptR)
library(tidyverse)
library(here)

# Set wd====
dir<-here()

# Load data====
ArchieRep<-read.csv("Data/Archie2014Repeat.csv")

hist(ArchieRep$SCI.F)
shapiro.test(ArchieRep$SCI.F)
hist(ArchieRep$SCI.M)
shapiro.test(ArchieRep$SCI.M)

m4<-rpt(SCI.F~(1|Animal.identity), grname=c("Animal.identity"), data=ArchieRep)
m4 #R=0.343

m5<-rpt(SCI.M~(1|Animal.identity), grname=c("Animal.identity"), data=ArchieRep)
m5 #R=0.246


m6<-rpt(dominance.rank~(1|Animal.identity), grname=c("Animal.identity"), data=ArchieRep)
m6 #R=0.783

ArchieClean<-subset(ArchieRep, SCI.F!="NULL")

Archie2<-subset(ArchieRep, dominance.rank!="NULL")

Dom = Archie2 %>% 
  group_by(Animal.identity) %>% 
  summarize(AvgDom = mean(dominance.rank))


SCI.F = ArchieClean %>% 
  group_by(Animal.identity) %>% 
  summarize(AvgSCI.F = mean(SCI.F))

SCI.M = ArchieClean %>% 
  group_by(Animal.identity) %>% 
  summarize(AvgSCI.M = mean(SCI.M))

ArchieSurv<-subset(ArchieRep, Failure==1)

ArchieFull<-merge(SCI.F,ArchieSurv, by="Animal.identity")
ArchieFull<-merge(SCI.M, ArchieFull, by="Animal.identity")
ArchieDom<-merge(Dom, ArchieSurv, by="Animal.identity")

cor.test(ArchieDom$AvgDom,ArchieDom$Survival.time..cumulative.number.of.days.as.an.adult.)
#r=0.002217117

cor.test(ArchieFull$AvgSCI.M,ArchieFull$Survival.time..cumulative.number.of.days.as.an.adult.)
#r=0.1776592 SCI.M~SurvivalTime

cor.test(ArchieFull$AvgSCI.F,ArchieFull$Survival.time..cumulative.number.of.days.as.an.adult.)
#r=0.05012702 SCI.F~SurvivalTime