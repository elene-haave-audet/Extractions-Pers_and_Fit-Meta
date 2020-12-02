setwd("G:/Shared drives/Personality & Fitness Meta-Analysis/R_Pers&FitExtractions/Data")

# 2) Load libraries
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)

#load priors for MCMCglmm

prior.miw<-list(R=list(V=diag(2), nu=2.002), G=list(G1=list(V=diag(2), nu=2.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000)))

prior.invgamma<-list(R=list(V=diag(2), nu=1.02), G=list(G1=list(V=diag(2), nu=1.02)))

#use when behaviour is binomial and survival is gaussian
Prior3 <- list(R=list(V=diag(2),nu=3, fix = 2),G=list(G1 =list(V = diag(2), nu=3, alpha.mu = c(0,0), alpha.V = diag(c(1000, 25^2)))))

#use when both response and predictor are binomial
Prior2 <- list(R = list(V=diag(2),nu=3),
               G=list(G1 =list(V = diag(2), nu=3, alpha.V = diag(2)*1000)))

#use when both response and predictor are Gaussian
prior1 <-list(R = list(V = diag(c(1,1),2,2), nu = 3, fix = 2),
              G=list(G1 =list(V = diag(2), nu=3, alpha.mu = c(0,0), alpha.V = diag(25^2,2))))

#use when behaviour is gaussian and survival is binomial
prior4 <- list(R = list(V=diag(2), nu=3,fix=2),
               G=list(G1 =list(V = diag(2), nu=3, alpha.mu = c(0,0), alpha.V=diag(c(25^2,1000)))))

#load data
MarshMSurv<-read.csv("Marsh_Male_Surv.csv")
MarshFSurv<-read.csv("Marsh_Fem_Surv.csv")
MarshFBS<-read.csv("Marsh_Female_BS.csv")
MarshMBS<-read.csv("Marsh_Male_BS.csv")
MarshMMG<-read.csv("Marsh_Male_MG.csv")

#keep last observation for each individual to get survival
count(MarshMSurv, indiv) #75 individuals

male_survival=MarshMSurv %>% 
  group_by(indiv) %>% 
  slice(c(n())) %>% 
  ungroup() #end with 75 observations

count(MarshFSurv, indiv) #39 individuals

female_survival=MarshFSurv %>% 
  group_by(indiv) %>% 
  slice(c(n())) %>% 
  ungroup() #end with 39 observations

#combine survival & behavioural data into a single df
maleBS=left_join(male_survival, MarshMBS, by="indiv")
maleMG=left_join(male_survival, MarshMMG, by="indiv")
femaleBS=left_join(female_survival, MarshFBS, by="indiv")

hist(maleBS$bs.freq) #poisson
hist(maleBS$bs.sess) #right skewed
hist(maleMG$mg.freq) #poisson
hist(maleMG$mg.sess)#poisson
hist(femaleBS$bs.freq) #poisson
hist(femaleBS$bs.sess) #almost normal

# partition variance for male babysitting sessions
m1<-MCMCglmm(cbind(bs.sess, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=maleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m1)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m1$VCV[,1:4]) #-0.05
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m1$VCV[,5:8]) #0
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

#male babysitting frequency
m2<-MCMCglmm(cbind(bs.freq, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("poisson", "categorical"), data=maleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m2)

c3 <- posterior.cor(m2$VCV[,1:4])
round(apply(c3,2,mean),2)
round(apply(c3,2, quantile, c(0.025, 0.975)),2)
c4 <- posterior.cor(m2$VCV[,5:8])
round(apply(c4,2,mean),2)
round(apply(c4,2, quantile, c(0.025, 0.975)),2)

#female babysitting session
m3<-MCMCglmm(cbind(bs.sess, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=femaleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m3)

c5 <- posterior.cor(m3$VCV[,1:4])
round(apply(c5,2,mean),2)
round(apply(c5,2, quantile, c(0.025, 0.975)),2)
c6 <- posterior.cor(m3$VCV[,5:8])
round(apply(c6,2,mean),2)
round(apply(c6,2, quantile, c(0.025, 0.975)),2)

#female babysiting frequency
m4<-MCMCglmm(cbind(bs.freq, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=femaleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m4)

c7 <- posterior.cor(m4$VCV[,1:4])
round(apply(c7,2,mean),2)
round(apply(c7,2, quantile, c(0.025, 0.975)),2)
c8 <- posterior.cor(m4$VCV[,5:8])
round(apply(c8,2,mean),2)
round(apply(c8,2, quantile, c(0.025, 0.975)),2)