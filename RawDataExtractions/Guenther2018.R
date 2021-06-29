# Guenther 2018 rec 99
# Variance partitioning

# Load libraries====
require(MCMCglmm)
require(lme4)
require(VGAM)
require(pscl)
require(rptR)
library(here)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/data99.csv")
names(data)

# Load priors====
prior.miw<-list(R=list(V=diag(2), nu=2.002), G=list(G1=list(V=diag(2), nu=2.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000)))
prior.invgamma<-list(R=list(V=diag(2), nu=1.02), G=list(G1=list(V=diag(2), nu=1.02)))


##################### NO_LATENCY BEHAVIOUR DO NOT USE THAT BEHAVIOUR
hist(data$No_latency)
data$No_latencyln<-log(data$No_latency)
hist(data$No_latencyln)
#not good transfo
#don't know which family to use to model bimodal data??????
m1<-MCMCglmm(cbind(No_latency, littersize) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("geometric", "gaussian"), data=data, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m1)

# NO_CONTACT BEHAVIOUR====
hist(data$No_contacts)
hist(data$littersize)
#poisson, do not transform, keep as

###NO_contact~litter size --> switch sign to reflect neophpbia, not neophilia
m2<-MCMCglmm(cbind(No_contacts, littersize) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("poisson", "poisson"), data=data, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m2)

###posterior correlation matrix - 1 through 16 is among individual, 10-18 is within-individual
c1 <- posterior.cor(m2$VCV[,1:4]) #-0.06
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m2$VCV[,5:8]) #0.04
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

# phenotypic correlation
m3<-glmer(littersize~No_contacts+(1|ID), family="poisson", data=data)
summary(m3)

rpt(No_contacts~(1|ID), grname = "ID", data=data, datatype = c("Poisson"))


# HAND behaviour====
hist(data$Hand) #assume poisson
data$Hand_t<-(data$Hand)^1/2
hist(data$Hand_t)

# phenotypic correlation
m4<-glmer(littersize~Hand + (1|ID),data=data, family = "poisson")
summary(m4)

# Hand~litter size
m5<-MCMCglmm(cbind(Hand, littersize) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("gaussian", "poisson"), data=data, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m5)

###posterior correlation matrix - 1 through 16 is among individual, 10-18 is within-individual
c1 <- posterior.cor(m5$VCV[,1:4]) #-0.11
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m5$VCV[,5:8]) #-0.16
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

rpt(Hand~(1|ID), grname = "ID", data=data)

# STRUGGLE behaviour====
hist(data$struggle)
data$struggleln<-log(data$struggle)
hist(data$struggleln)
#not good transfo
#let's try square root transfo
data$StruggleS<-sqrt(data$struggle)
hist(data$StruggleS)
#nope
#let's try cube root
data$StruggleC<-sign(data$Struggle)*abs(data$Struggle)^(1/3)
hist(data$StruggleC)
#nope, so will have to run MCMCglmm with poisson or other family types

# Struggle~litter size
#doesn't work with poisson or binomial family (error data must be 
#non-negative integers), so used "gaussian", "exponential")
m6<-MCMCglmm(cbind(struggle, littersize) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("exponential", "gaussian"), data=data, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m6)
#okish

###posterior correlation matrix - 1 through 16 is among individual, 10-18 is within-individual
c1 <- posterior.cor(m6VCV[,1:4]) #-0.51
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m6$VCV[,5:8]) #0.04
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

# phenotypic correlation
m7<-glmer(littersize~struggle+(1|ID), data=data, family="poisson", na.action = "na.omit")
summary(m7)

rpt(struggle~(1|ID), grname="ID",data=data)

# OF_DISTANCE behaviour====
hist(data$OF_distance)
data$OF_distanceln<-log(data$OF_distance)
hist(data$OF_distanceln)
#not good transfo
#let's try square root transfo
data$OF_distanceS<-sqrt(data$OF_distance)
hist(data$OF_distanceS)
#nope
#let's try cube root
data$OF_distanceC<-sign(data$OF_distance)*abs(data$OF_distance)^(1/3)
hist(data$OF_distanceC)
#nope, so will have to run MCMCglmm with exponential (same as poisson but for continuous data
#or hurdle, =0 and >0)
m8<-MCMCglmm(cbind(OF_distance, littersize) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("exponential", "gaussian"), data=data, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m8)

###posterior correlation matrix - 1 through 16 is among individual, 10-18 is within-individual
c1 <- posterior.cor(m8$VCV[,1:4]) #0.65
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m8$VCV[,5:8]) #-0.15
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

# phenotypic correlation
m9<-glmer(littersize~OF_distanceC+ (1|ID), family="poisson", data=data, na.action = "na.omit")
summary(m9)

rpt(OF_distance~(1|ID), grname = "ID", data=data)

# LF_TRIPS behaviour====
hist(data$LF_trips)
data$LF_tripsln<-log(data$LF_trips)
hist(data$LF_tripsln)
#ish but not great

#poisson distribution or neg binomial?
model.pois=glmer(LF_trips~littersize +(1|ID), na.action=na.exclude, data=data, family=poisson)
summary(model.pois)
#overdisped, deviance >df residuals
model.nb=glmer.nb(LF_trips~littersize +(1|ID), na.action=na.exclude, data=data)
summary(model.nb)
#still overdispersed,let's try zero-inflated poisson
model.zip=zeroinfl(LF_trips~littersize |ID, na.action=na.exclude, data=data)
summary(model.zip)

#let's try without transforming the data and using zero inflated poisson distribution (count data)
m10<-MCMCglmm(cbind(LF_trips, littersize) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("poisson", "poisson"), data=data, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m10)

c1 <- posterior.cor(m10$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m10$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)
