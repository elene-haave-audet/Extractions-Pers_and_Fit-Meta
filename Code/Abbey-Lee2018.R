# Abbey-Lee 2018 Rec 5
# Variance partitioning

# Load libraries====
library(here)
library(tidyverse)
library(MCMCglmm)
library(rptR)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/AbbeyLee2018.csv")

ind=data %>% 
  group_by(MaleID) %>% 
  summarize(Male=n_distinct(MaleID))

# Models====
m1<-rpt(Exploration~(1|MaleID), data=data, grname="MaleID")
m1

hist(data$Exploration) #can assume normality

# prior
#use when behaviour is binomial and survival is gaussian
Prior3 <- list(R=list(V=diag(2),nu=3, fix = 2),G=list(G1 =list(V = diag(2), nu=3, alpha.mu = c(0,0), alpha.V = diag(c(1000, 25^2)))))

m3<-MCMCglmm(cbind(Exploration, EPY) ~ (trait-1), random = ~us(trait):MaleID ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=data, prior =Prior3, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m3) #Maria's prior seems to work okay for this model

c1 <- posterior.cor(m3$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m3$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)