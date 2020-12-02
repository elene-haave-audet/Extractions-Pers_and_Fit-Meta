#=======================================================
# Phenotypic extractions for Kralj-Fiser 2017, Rec 326
# ======================================================

# 1) set working directory to data location
setwd("G:/Shared drives/Personality & Fitness Meta-Analysis/R_Pers&FitExtractions/Data")

# 2) Load libraries
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(GeneNet)

# 2) Load data, female=1, male=2
data<-read.csv("KF2017.csv")

# 3)Pull second observation for each individual into separate dataframe
data.r1<-select(data, ID, sex,Novel_env_time._to.stop1,Novel_env_time._again1,boldness1,agg1,prey_reaction1,prey_attack1,voracity1,n_viable_eggsacs,lab_longevity)
data.r1$trial <- rep(1,nrow(data.r1)) # add a column for trial number
data.r1<-dplyr::rename(data.r1, NE_start = Novel_env_time._to.stop1,NE_Resume=Novel_env_time._again1, Bold=boldness1, Agg=agg1,prey_reaction=prey_reaction1,prey_attack=prey_attack1,voracity=voracity1)

data.r2<-select(data, ID, sex, Novel_env_time._to.stop2,Novel_env_time._again2,boldness2,agg2,prey_reaction2,prey_attack2,voracity2,n_viable_eggsacs,lab_longevity)
data.r2$trial<-rep(2, nrow(data.r2))
data.r2<-dplyr::rename(data.r2, NE_start = Novel_env_time._to.stop2,NE_Resume=Novel_env_time._again2, Bold=boldness2, Agg=agg2,prey_reaction=prey_reaction2,prey_attack=prey_attack2,voracity=voracity2)

data.full<-bind_rows(data.r1,data.r2) #combine data from both trials to run mixed effect models
data.f<-subset(data.full, sex==1) #create separate data frames for males and females
data.m<-subset(data.full, sex==2)

hist(data.f$n_viable_eggsacs) #poisson
data.f$eggs_t<-(data.f$n_viable_eggsacs)^2
hist(data.f$eggs_t)#more poisson
hist(data.f$lab_longevity) #assume normal?
hist(data.f$NE_start)
data.f$NE_start_t<-log((data.f$NE_start)+0.5)

# 4) models 1, 11-13 female repro
m1<-glmer(n_viable_eggsacs~(NE_start)+(1|ID), data=data.f, family="poisson")
##can't get model to converge!

m11<-glmer(eggs_t~(NE_Resume)+(1|ID), data=data.f, family="poisson")
r.squaredGLMM(m11) #R2=1.880345e-05; neg (switch sign to reflect exploration)
summary(m11) #obs=64, ind=33
rpt(NE_Resume~(1|ID), grname = "ID", data=data.f) #r=0.816

m12<-glmer(eggs_t~Bold+(1|ID), data=data.f, family="poisson")
r.squaredGLMM(m12) #R2=5.754206e-06, pos, obs=64, ind=33
summary(m12)

m13<-glmer(eggs_t~Agg+(1|ID), data=data.f, family="poisson")
r.squaredGLMM(m13) #R2=8.799589e-05, pos, obs=66, ind=33
summary(m13)

# 5) models 2-6 female longevity
m2<- glmer(lab_longevity~scale(NE_start)+(1|ID), data=data.f, family="poisson")
r.squaredGLMM(m2) #R2=8.658333e-05 positive, obs=62, ind=32
summary(m2)
rpt(NE_start~(1|ID), grname = "ID", data=data.f) #r=0.23

#check model estimates obtained from Moiron et al
prior1 <-list(R = list(V = diag(c(1,1),2,2), nu = 3, fix = 2),
              G=list(G1 =list(V = diag(2), nu=3, alpha.mu = c(0,0), alpha.V = diag(25^2,2))))


moiron<- MCMCglmm(cbind(NE_start, lab_longevity) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=data.f, prior = prior1, verbose = FALSE,nitt=3200000,thin=2000,burnin=300000)
plot(moiron)

c1 <- posterior.cor(moiron$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(moiron$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)


m3<-glmer(lab_longevity~Agg+(1|ID), data=data.f, family="poisson")
r.squaredGLMM(m3) #R2=0.0002561427 negative, obs=64, ind=32
summary(m3)
rpt(Agg~(1|ID), grname = "ID", data=data.f) #r=0.673

m4<-lmer(lab_longevity~NE_Resume+(1|ID), data=data.f) #model does not converge

m5<-glmer(lab_longevity~Bold+(1|ID), data=data.f, family="poisson")
r.squaredGLMM(m5) #R2=0.0002046803 pos, obs=62, ind=32
summary(m5)
rpt(Bold~(1|ID), grname = "ID", data=data.f)#r=0.381

m6<-glmer(lab_longevity~scale(voracity)+(1|ID), data=data.f, family="poisson")
r.squaredGLMM(m6) #R2=2.289710e-05, pos, obs=42, ind=22
summary(m6)
rpt(voracity~(1|ID), grname = "ID", data=data.f) #r=0

# 6) models 7-10 male longevity
m7<-glmer(lab_longevity~scale(NE_start)+(1|ID), data=data.m, family="poisson")
r.squaredGLMM(m7) #R2=0.0005792962 neg, obs=39, ind=20
summary(m7) 
rpt(NE_start~(1|ID), grname = "ID", data=data.m) #r=0.4

m8<-glmer(lab_longevity~scale(NE_Resume)+(1|ID), data=data.m, family="poisson")
r.squaredGLMM(m8) #R2=1.697611e-05; neg (swith sign to POS to refelct exploration), obs=39, ind=20
summary(m8)
rpt(NE_Resume~(1|ID), grname = "ID", data=data.m) #r=0.105

m9<-glmer(lab_longevity~scale(Agg)+(1|ID), data=data.m, family="poisson")
r.squaredGLMM(m9) #R2=0.0001161308, neg, obs=40, ind=20
summary(m9)
rpt(Agg~(1|ID), grname = "ID", datatype = c("Poisson"), data=data.m) #r=0.777

m10<-glmer(lab_longevity~Bold+(1|ID), data=data.m, family="poisson")
r.squaredGLMM(m10) #R2=2.513219e-05; neg, obs=39, ind=20
summary(m10)
rpt(Bold~(1|ID), grname = "ID", data=data.m) #r=0.493