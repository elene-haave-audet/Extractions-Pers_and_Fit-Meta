# Variance Partitioning for Grist et al 2017, Rec 97
# Elene Haave Audet, Mar 23, 2020
# Day 12 pandemic

# Load libraries====
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(here)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/Grist2017.csv")

# Visualize data====
hist(data$Breed_Success) #count data that's not typical poisson
  #Sex: Male=1 Female=0
  #migration: Migratory=1 Resident =0

plot(data$Breed_Success, data$Migratory_Strategy)

data.f<-subset(data, Sex==0)
data.m<-subset(data, Sex==1)

table(data$Sex)

summary.f=data.f %>% 
  summarize(n_distinct(Bird_ID)) #211 Ind, 292 Obs

summary.m=data.m %>% 
  summarize(n_distinct(Bird_ID)) #224 Ind, 353 Obs

# Load priors====
prior.miw<-list(R=list(V=diag(2), nu=2.002), G=list(G1=list(V=diag(2), nu=2.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000)))

# Run MCMC glmm====
m.f<-MCMCglmm(cbind(Migratory_Strategy, Breed_Success) ~ (trait-1), random = ~us(trait):Bird_ID ,rcov = ~us(trait):units, family = c("categorical", "gaussian"), data=data.f, prior =prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m.f)

m.m<-MCMCglmm(cbind(Migratory_Strategy, Breed_Success) ~ (trait-1), random = ~us(trait):Bird_ID ,rcov = ~us(trait):units, family = c("categorical", "gaussian"), data=data.m, prior =prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m.m)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
#females
c1 <- posterior.cor(m.f$VCV[,1:4]) #among=-0.26
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m.f$VCV[,5:8]) #within=-0.4
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

#males
c3 <- posterior.cor(m.m$VCV[,1:4]) #among=-0.39
round(apply(c3,2,mean),2)
round(apply(c3,2, quantile, c(0.025, 0.975)),2)
c4 <- posterior.cor(m.m$VCV[,5:8]) #within=-0.48
round(apply(c4,2,mean),2)
round(apply(c4,2, quantile, c(0.025, 0.975)),2)

# Estimate repeatability====
#females
rpt(Migratory_Strategy~(1|Bird_ID), grname = "Bird_ID", datatype = c("Binary"), data=data.f)
#r=0.518

#males
rpt(Migratory_Strategy~(1|Bird_ID), grname = "Bird_ID", datatype = c("Binary"), data=data.m)
#r= 0.275

# Estimate Pheno correlation====
#females
cor.test(data.f$Migratory_Strategy,data.f$Breed_Success)
#r=-0.1282292

#males
cor.test(data.m$Migratory_Strategy,data.m$Breed_Success)
#r=-0.1428272