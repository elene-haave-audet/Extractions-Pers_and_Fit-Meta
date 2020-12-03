# DiRienzo 2018 Rec 72
# Phenotypic and among-ind extractions

# Load libraries====
library(here)
library(MCMCglmm)
library(rptR)
library(lme4)

# Set wd====
dir<-here()

# Load priors====

# Load data====
data.f<-read.csv("Data/DiRienzoF.csv")
cor.test(data.f$attack_number, data.f$success) #r=0.0731635

data.m<-read.csv("Data/DiRienzoM.csv")
hist(data.m$prop_decon)
data.m$decon_t<-((data.m$prop_decon)+0.5)^(1/3)
hist(data.m$decon_t) #somewhat more normal

m.r<-rpt(decon_t~(1|male_ID), grname = "male_ID", data=data.m)
m.r #0.161

# phenotypic cor====
m2<-glmer(success~decon_t+ (1|male_ID),family="binomial", data=data.m)
summary(m2)
r.squaredGLMM(m2)

plot<-ggplot(data.m, aes(x=decon_t, y=success))+
  geom_smooth(method = "glm")
plot

# Amg-ind====
m1<-MCMCglmm(cbind(decon_t, success) ~ (trait-1), random = ~us(trait):male_ID ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=data.m, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m1)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m1$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m1$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)