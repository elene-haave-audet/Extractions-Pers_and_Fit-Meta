# Boulton 2018 rec 41
# Phenotypic estimates

# Load libraries====
library(here)
library(lme4)
library(tidyverse)

# Set wd====
dir<-here()

# Load data====
survival<-read.csv("Data/BoultonSurvival.csv")
dominance<-read.csv("Data/BoultonDominance.csv")
activity<-read.csv("Data/BoultonActivity.csv")

# Males====
male<-subset(survival, Sex==2) #pair multiple dom obs w/ survival
male<-left_join(survival, dominance, by="Fish")
male<-filter(male,!is.na(DOMscore))

m1<-glmer(surv~ DOMscore+ (1|Fish), data=male, family="binomial")
summary(m1)
r.squaredGLMM(m1)
plot(male$surv~male$DOMscore)
abline(lm(male$surv~male$DOMscore))

m.act<-left_join(male, activity, by="Fish")

m2<-glmer(surv~log(Act)+(1|Fish), data=m.act, family="binomial")
summary(m2)
r.squaredGLMM(m2)

# Females====
female<-subset(survival, sex=1)
female<-left_join(survival, activity, by="Fish")

m3<-glmer(surv~log(Act)+(1|Fish), data=female, family="binomial")
summary(m3)
r.squaredGLMM(m3)