# Quinn 2009 rec 206
# Variance partitioning

# Load libraries====
library(here)
library(MCMCglmm)
library(lme4)
library(rptR)
library(MuMIn)
library(tidyverse)
library(lmerTest)


# Set wd====
dir<-here()

# Load data====
eb.data<-read.csv("Data/QuinnEB.csv")
breeding.data<-read.csv("Data/QuinnBreeding.csv")

# Prep data====
female.data<-subset(breeding.data, select=-c(Father))
female.data<-plyr::rename(female.data, c("Mother"="BirdID"))
female.data<-filter(female.data, !is.na(BirdID))
female.eb<-left_join(eb.data, female.data, by="BirdID")
female.eb<-filter(female.eb,!is.na(Clutch.size))  

# Prior====
prior.miw<-list(R=list(V=diag(2), nu=2.002), G=list(G1=list(V=diag(2), nu=2.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000)))

# Females====
hist(female.eb$Clutch.size) #beautifully normal!
hist(female.eb$hopsandflights) #zero inflated poisson
hist(female.eb$longevity_max.age.recorded) #poisson

rpt(hopsandflights~(1|BirdID), grname = "BirdID", data=female.eb)
#r=0.603

m.pheno.repro.f<-lmer(Clutch.size~scale(hopsandflights)+(1|BirdID), data=female.eb)
summary(m.pheno.repro.f) #neg sign
r.squaredGLMM(m.pheno.repro.f) #R2=0.00126
anova(m.pheno.repro.f)

m2.pheno.repro.f<-MCMCglmm(cbind(hopsandflights, Clutch.size) ~ (trait-1), random = ~us(trait):BirdID ,rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=female.eb, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m2.pheno.repro.f) #beautiful!

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m2.pheno.repro.f$VCV[,1:4]) #-0.09
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m2.pheno.repro.f$VCV[,5:8]) #0
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

m.pheno.surv.f<-glmer(longevity_max.age.recorded~scale(hopsandflights)+(1|BirdID), family="poisson", data=female.eb)
summary(m.pheno.surv.f) #pos sign
r.squaredGLMM(m.pheno.surv.f) #R2=0.0000047


# Males====
male.data<-subset(breeding.data, select=-c(Mother))
male.data<-plyr::rename(male.data, c("Father"="BirdID"))
male.data<-filter(male.data, !is.na(BirdID))
male.eb<-left_join(eb.data, male.data, by="BirdID")
male.eb<-filter(male.eb,!is.na(Clutch.size))  

hist(male.eb$Clutch.size) #beautifully normal!
hist(male.eb$hopsandflights) #zero inflated poisson
hist(male.eb$longevity_max.age.recorded) #poisson

rpt(hopsandflights~(1|BirdID), grname = "BirdID", data=male.eb)
#r=0.409

m.pheno.repro.m<-lmer(Clutch.size~scale(hopsandflights)+(1|BirdID), data=male.eb)
summary(m.pheno.repro.m) #neg sign
r.squaredGLMM(m.pheno.repro.m) #R2=0.00019709
anova(m.pheno.repro.m)

m2.pheno.repro.m<-MCMCglmm(cbind(hopsandflights, Clutch.size) ~ (trait-1), random = ~us(trait):BirdID ,rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=male.eb, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m2.pheno.repro.m) #beautiful!

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m2.pheno.repro.m$VCV[,1:4]) #-0.05
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m2.pheno.repro.m$VCV[,5:8]) #0
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

m.pheno.surv.m<-glmer(longevity_max.age.recorded~scale(hopsandflights)+(1|BirdID), family="poisson", data=male.eb)
summary(m.pheno.surv.m) #pos sign
r.squaredGLMM(m.pheno.surv.m) #R2=0.001360041