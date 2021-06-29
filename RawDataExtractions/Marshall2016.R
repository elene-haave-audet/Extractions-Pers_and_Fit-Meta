# Marshall et al 2016, Rec 323

# Load libraries====
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(here)

# Set wd====
dir<-here()

# Load priors====

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

# Load data====
MarshMSurv<-read.csv("Data/Marsh_Male_Surv.csv")
MarshFSurv<-read.csv("Data/Marsh_Fem_Surv.csv")
MarshFBS<-read.csv("Data/Marsh_Female_BS.csv")
MarshMBS<-read.csv("Data/Marsh_Male_BS.csv")
MarshMMG<-read.csv("Data/Marsh_Male_MG.csv")

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

# Extract estimates====
# Male babysitting freq####

hist(maleBS$bs.freq) #poisson

#phenotypic
#remove NAs
maleBS_clean<-filter(maleBS, !is.na(bs.freq))

m1<-glmer(survival~bs.freq + (1|indiv), family="binomial", data=maleBS_clean)
plot(resid(m1))
summary(m1)
rpt(bs.freq~(1|indiv), grname = "indiv", data=maleBS, datatype = c("Poisson"))
#r=0.021
r.squaredGLMM(m1)

#among
m2<-MCMCglmm(cbind(bs.freq, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("poisson", "categorical"), data=maleBS, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m2)

c3 <- posterior.cor(m2$VCV[,1:4])
round(apply(c3,2,mean),2)
round(apply(c3,2, quantile, c(0.025, 0.975)),2)
c4 <- posterior.cor(m2$VCV[,5:8])
round(apply(c4,2,mean),2)
round(apply(c4,2, quantile, c(0.025, 0.975)),2)

# Female babysiting freq####
hist(femaleBS$bs.freq) #poisson

# Phenotypic
# remove NAs
femalebs_clean<-filter(femaleBS, !is.na(bs.freq))
m3<-glmer(survival~bs.freq + (1|indiv), family="binomial", data=femalebs_clean)
summary(m3)
rpt(bs.freq~(1|indiv), grname = "indiv", data=femaleBS, datatype = c("Poisson"))
#r=0.049
r.squaredGLMM(m3)

m4<-MCMCglmm(cbind(bs.freq, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("poisson", "categorical"), data=femaleBS, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m4)

c7 <- posterior.cor(m4$VCV[,1:4]) 
round(apply(c7,2,mean),2)
round(apply(c7,2, quantile, c(0.025, 0.975)),2)
c8 <- posterior.cor(m4$VCV[,5:8])
round(apply(c8,2,mean),2)
round(apply(c8,2, quantile, c(0.025, 0.975)),2)

# Male mate-guarding freq####

hist(maleMG$mg.freq) #poisson

maleMG<-filter(maleMG, !is.na(mg.freq))

# Phenotypic
m5<-glmer(survival~mg.freq + (1|indiv), family="binomial", data=maleMG)
summary(m5)
rpt(mg.freq~(1|indiv), grname = "indiv", data=maleMG, datatype = c("Poisson"))
#r=0.297
r.squaredGLMM(m5)

# Among
m6<-MCMCglmm(cbind(mg.freq, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("poisson", "categorical"), data=maleMG, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m6)

c9 <- posterior.cor(m6$VCV[,1:4]) 
round(apply(c9,2,mean),2)
round(apply(c9,2, quantile, c(0.025, 0.975)),2)
c10 <- posterior.cor(m6$VCV[,5:8])
round(apply(c10,2,mean),2)
round(apply(c10,2, quantile, c(0.025, 0.975)),2)
