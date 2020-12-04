# Van Overveld 2015 rec 261
# Phenotypic estimates

# Load libraries
library(here)
library(lme4)
library(rptR)
library(MuMIn)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/VanOverveld.csv")
data.f<-subset(data, sex=="F")
data.m<-subset(data, sex=="M")
hist(data.f$age)

# Models====
m1<-lmer(age~Exploration_Score+(1|Bird_ID), data=data.f)
summary(m1)
r.squaredGLMM(m1)
rpt(Exploration_Score~(1|Bird_ID), grname = "Bird_ID", data = data.f)

m2<-lmer(age~Exploration_Score+(1|Bird_ID), data=data.m)
summary(m2)
r.squaredGLMM(m2)
rpt(Exploration_Score~(1|Bird_ID), grname = "Bird_ID", data=data.m)
