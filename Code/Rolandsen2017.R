# Rolandsen 2017 Rec 2
# Variance partitioning

# Load libraries====
library(here)
library(MCMCglmm)

# Load priors====

# Load data====
data<-read.csv("Data/Rol2017surv.csv")
data2<-read.csv("Data/Rol2017twin.csv")

data$DeadAlive<-as.factor(data$DeadAlive)
data$tactic<-as.factor(data$tactic)
data$twins<-as.factor(data$twins)
data$year<-as.factor(data$year)
data$mother<-as.factor(data$mother)

data2$tactic<-as.factor(data2$tactic)
data2$twins<-as.factor(data2$twins)
data2$year<-as.factor(data2$year)
data2$mother<-as.factor(data2$mother)

# Models====
m4<-MCMCglmm(cbind(tactic,DeadAlive) ~ (trait-1), random = ~us(trait):mother ,rcov = ~us(trait):units, family = c("categorical", "categorical"), data=data, prior = Prior2, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m4) #none of the 6 priors we have had worked for this model

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m4$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m4$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

m5<-MCMCglmm(cbind(tactic,twins) ~ (trait-1), random = ~us(trait):mother ,rcov = ~us(trait):units, family = c("categorical", "categorical"), data=data2, prior = Prior2, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m5)