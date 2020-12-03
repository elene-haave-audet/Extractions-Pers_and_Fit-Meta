# Le Coeur 2015
# Variance partitioning

# Load libraries====
library(here)
library(MCMCglmm)
library(lme4)
library(rptR)

# Set wd====
dir<-here()

# Load data====
data<- read.csv("Data/LeCoeur2015.csv")

hist(data$trappy) #Poisson dist
hist(data$ars) #Poisson dist

# Models====
m1<-MCMCglmm(cbind(trappy, ars) ~ (trait-1), random = ~us(trait):ID ,rcov = ~us(trait):units, family = c("poisson", "poisson"), data=data, prior = Prior2, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m1) #again, priors look better!

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m1$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m1$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

cor.test(data$trappy, data$ars, method = c("spearman")) #r=-0.01

r1<-rpt(trappy~(1|ID), grname=c("ID"), data=data, datatype=c("Poisson"))
r1
hist(data$trappy)

r2<-lmer(trappy_t~ -1+(1|ID), data=data)
summary(r2)
plot(r2)