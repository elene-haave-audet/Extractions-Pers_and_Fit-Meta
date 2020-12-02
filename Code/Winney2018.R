#=================================================
# Variance partitioning for Winney 2018, Rec 280
# March 25, 2020 (I'm so bored of being at home)
#=================================================

# 1) set working directory to data location
setwd("G:/Shared drives/Personality & Fitness Meta-Analysis/R_Pers&FitExtractions/Data")

# 2) Load libraries
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(GeneNet)

# 3) Load data
exp<-read.csv("Winney2018_Exploration.csv")
bold<-read.csv("Winney2018_Boldness.csv")
#ped<-read.table("Winney2018_Pedigree.txt", header=TRUE)
#won't need to use predigree; use n offspring in box at time of assay as fitness proxy

# 4) Load priors
prior.miw<-list(R=list(V=diag(2), nu=2.002), G=list(G1=list(V=diag(2), nu=2.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000)))

# 5) model effect of neophilia (ne boldness) on n offspring
#study used latency from touch box to entering box 
#split males and females
bold.f<-subset(bold, sex=="F")
bold.m<-subset(bold, sex=="M")

#visualize data
##females
hist(bold.f$offspring) #poisson
hist(bold.f$logentrylatency) #normal; data normalized and * by -1, so that >values = >bold

sample_size_f<-bold.f %>% 
  summarize(n_distinct(asrid)) #n ind=122, obs= 320

##males
hist(bold.m$offspring) #poisson
hist(bold.m$logentrylatency) #normalish

sample_size_m<-bold.m %>% 
  summarize(n_distinct(asrid)) #n ind=99, obs=220

##MCMCglmmm
m.bold.f<-MCMCglmm(cbind(logentrylatency, offspring) ~ (trait-1), random = ~us(trait):asrid ,rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=bold.f, prior =prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m.bold.f)
rpt(logentrylatency~(1|asrid), grname = "asrid", data=bold.f) #r=0.26

m.bold.m<-MCMCglmm(cbind(logentrylatency, offspring) ~ (trait-1), random = ~us(trait):asrid ,rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=bold.m, prior =prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m.bold.m)
rpt(logentrylatency~(1|asrid), grname = "asrid", data=bold.m) #r=0.081

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
#females
c1 <- posterior.cor(m.bold.f$VCV[,1:4]) #among=0.01
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m.bold.f$VCV[,5:8]) #within=0.12
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

#males
c3 <- posterior.cor(m.bold.m$VCV[,1:4]) #among=-0.06
round(apply(c3,2,mean),2)
round(apply(c3,2, quantile, c(0.025, 0.975)),2)
c4 <- posterior.cor(m.bold.m$VCV[,5:8]) #within=0.02
round(apply(c4,2,mean),2)
round(apply(c4,2, quantile, c(0.025, 0.975)),2)

# 6) model effect of exploration on fitness
#use logfr.notdoubled
#tests conducted after breeding season; could link EB with fitness of next year's brood

#pull fitness data from boldness assay
exp.fit<-select(exp, asrid, winter, sex, logfr.notdoubled)
exp.fit<-rename(exp.fit, year=winter)
exp.fit$year<-exp.fit$year+1 #add 1 to year so that matches brood size of that season
exp.fit$year<-as.factor(exp.fit$year) #make factor so data frames can be merged
exp.fit$asrid<-as.factor(exp.fit$asrid)

bold.fit<-select(bold, asrid, year, meanoffspring) #take avg brood size in single breeding season (can have >1 brood/yr)
bold.fit<-bold.fit %>% 
  group_by(asrid, year) %>% slice(1) %>% ungroup() #extract single obs per id per yr
bold.fit$year<-as.factor(bold.fit$year)
bold.fit$asrid<-as.factor(bold.fit$asrid)

exp.fit<-full_join(exp.fit, bold.fit, by=c("asrid", "year"))
exp.fit<-subset(exp.fit, logfr.notdoubled!="NA") #remove NAs
exp.fit<-subset(exp.fit, meanoffspring!="NA")

#split males and females
exp.f<-subset(exp.fit, sex=="F")

n_ind_f<-exp.f %>% 
  summarise(n_distinct(asrid)) #n ind=48, obs=77

exp.m<-subset(exp.fit, sex=="M")

n_ind_m<-exp.m %>% 
  summarise(n_distinct(asrid)) #n ind=41, obs=68

#MCMCglmm
m.exp.f<-MCMCglmm(cbind(logfr.notdoubled, meanoffspring) ~ (trait-1), random = ~us(trait):asrid ,rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=exp.f, prior =prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m.exp.f)
rpt(logfr.notdoubled~(1|asrid), grname = "asrid", data=exp.f) #r=0.408

m.exp.m<-MCMCglmm(cbind(logfr.notdoubled, meanoffspring) ~ (trait-1), random = ~us(trait):asrid ,rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=exp.m, prior =prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m.exp.m)
rpt(logfr.notdoubled~(1|asrid), grname = "asrid", data = exp.m) #r=0

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
#females
c5 <- posterior.cor(m.exp.f$VCV[,1:4]) #among=0.07
round(apply(c5,2,mean),2)
round(apply(c5,2, quantile, c(0.025, 0.975)),2)
c6 <- posterior.cor(m.exp.f$VCV[,5:8]) #within=0.02
round(apply(c6,2,mean),2)
round(apply(c6,2, quantile, c(0.025, 0.975)),2)

#males
c7 <- posterior.cor(m.exp.m$VCV[,1:4]) #among=0.36
round(apply(c7,2,mean),2)
round(apply(c7,2, quantile, c(0.025, 0.975)),2)
c8 <- posterior.cor(m.exp.m$VCV[,5:8]) #within=0.02
round(apply(c8,2,mean),2)
round(apply(c8,2, quantile, c(0.025, 0.975)),2)
