# Marshall et al 2016, Rec 323

# Load libraries====
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(here)

#function to convert Zr to r        
Zr.to.r<-function(Zr){
  r<-(exp(2*Zr)-1)/(exp(2*Zr)+1)}

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

hist(maleBS$bs.freq) #poisson
maleBS$bs.freq_t<-log(maleBS$bs.freq)
hist(maleBS$bs.freq_t) #still poisson
maleBS$bs.freq_t2<-sqrt(maleBS$bs.freq)
hist(maleBS$bs.freq_t2) #still poisson
maleBS$bs.freq_t3<-sqrt(maleBS$bs.freq+0.5)
hist(maleBS$bs.freq_t3) #still poisson
maleBS$bs.freq_t4<-1/(maleBS$bs.freq)
hist(maleBS$bs.freq_t4) #reverse poisson
maleBS$bs.freq_t5<-(maleBS$bs.freq)^1/3
hist(maleBS$bs.freq_t5) #still poisson...
maleBS$bs.freq_t6<-sqrt(maleBS$bs.freq+3/8)
hist(maleBS$bs.freq_t6)
  
hist(maleBS$bs.sess) #right skewed, but pretty normal
hist(maleMG$mg.freq) #poisson
hist(maleMG$mg.sess)#poisson

hist(femaleBS$bs.freq) #poisson

femaleBS$bs.freq_t<-sqrt(femaleBS$bs.freq+0.5)
hist(femaleBS$bs.freq_t)

hist(femaleBS$bs.sess) #almost normal

# partition variance for male babysitting sessions
m1<-MCMCglmm(cbind(bs.sess, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=maleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m1)
#plot(resid(m1))
#Moiron Zr estimate= -0.187 (bolded to indicate it was flipped)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m1$VCV[,1:4]) #-0.33
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m1$VCV[,5:8]) #0
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

MaleBS.sess.Zr<-Zr.to.r(-0.33) #-0.318

#male babysitting frequency
m2<-MCMCglmm(cbind(bs.freq, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=maleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m2)
# Zr estimate from Moiron 0.204 (bolded to indicate sign was flipped)

c3 <- posterior.cor(m2$VCV[,1:4]) #-0.32
round(apply(c3,2,mean),2)
round(apply(c3,2, quantile, c(0.025, 0.975)),2)
c4 <- posterior.cor(m2$VCV[,5:8])
round(apply(c4,2,mean),2)
round(apply(c4,2, quantile, c(0.025, 0.975)),2)

MaleBS.freq.Zr<-Zr.to.r(-0.32) #-0.3095

#female babysitting session
m3<-MCMCglmm(cbind(bs.sess, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=femaleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m3)
# Zr estimate from Moiron -0.365 (bolded to indicate sign was flipped)

c5 <- posterior.cor(m3$VCV[,1:4]) #0.33
round(apply(c5,2,mean),2)
round(apply(c5,2, quantile, c(0.025, 0.975)),2)
c6 <- posterior.cor(m3$VCV[,5:8])
round(apply(c6,2,mean),2)
round(apply(c6,2, quantile, c(0.025, 0.975)),2)

#female babysiting frequency
m4<-MCMCglmm(cbind(bs.freq_t, survival) ~ (trait-1), random = ~us(trait):indiv ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=femaleBS, prior = prior4, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m4)
# Zr estimate from Moiron -0.031 (bolded to indicate sign flipped)

c7 <- posterior.cor(m4$VCV[,1:4]) #0.4 (without transfomration) # 0.36 with transformation
round(apply(c7,2,mean),2)
round(apply(c7,2, quantile, c(0.025, 0.975)),2)
c8 <- posterior.cor(m4$VCV[,5:8])
round(apply(c8,2,mean),2)
round(apply(c8,2, quantile, c(0.025, 0.975)),2)