# Santicchia 2018 rec 219
# Phenotypic and amg-ind extractions

# Load libraries====
library(here)
library(lme4)
library(MCMCglmm)
library(rptR)

# Load priors====

# Survival data====
San2018Surv<-read.csv("Data/Santicchia2018_Survival.csv")
table(San2018Surv$Sex) #63 F, 78 M

San2018M<-subset(San2018Surv,Sex=="Male")
San2018F1<-subset(San2018Surv,Sex=="Female")

# Behaviour data====
San2018Behav<-read.csv("Data/Santicchia2018_Repeat.csv")
table(San2018Behav$Sex) #19 F, 30 M

# Females====
San2018Litter<-read.csv("Data/Santicchia2018_Fem.csv")
San2018F<-merge(San2018Litter, San2018Behav, by="id")
San2018Merge<-merge(San2018Behav, San2018Surv, by="id")
San2018F.B<-subset(San2018Merge, Sex.x=="female")

San2018R1<-subset(San2018F, Repeat==1)
San2018R1$Trap<-San2018R1$Square.root.traps.1
San2018R1$Cap<-San2018R1$ln.captures.1

San2018R2<-subset(San2018F, Repeat==2)
San2018R2$Trap<-San2018R2$square.root.traps2
San2018R2$Cap<-San2018R2$ln.captures.2

San2018FComb<-rbind(San2018R1, San2018R2)
San2018FComb<-left_join(San2018FComb, San2018Surv, by="id")
hist(San2018FComb$Trap) #normal enough (exploration)
hist(San2018FComb$Cap) #looks normal (boldness)

m1<-glmer(Litter.produced~ scale(Trap)+ (1|id), family="binomial", data=San2018FComb)
summary(m1) #pos
r.squaredGLMM(m1) #R2=0.1271735

m2<-glmer(Litter.produced~scale(Cap)+ (1|id), family = "binomial", data= San2018FComb)
summary(m2) #pos
r.squaredGLMM(m2) #R2=0.04025647

rpt(Trap~(1|id), grname = "id", data=San2018FComb) #r=0.258
rpt(Cap~(1|id), grname="id", data=San2018FComb) #r=0.481

m3<-MCMCglmm(cbind(Trap, Litter.produced) ~ (trait-1), random = ~us(trait):id ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=San2018FComb, prior = prior4, verbose = FALSE,nitt=1003000,thin=1000,burnin=30000)
plot(m3)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m3$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m3$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

m4<-MCMCglmm(cbind(Cap, Litter.produced) ~ (trait-1), random = ~us(trait):id ,rcov = ~us(trait):units, family = c("gaussian", "categorical"), data=San2018FComb, prior = prior.miw, verbose = FALSE,nitt=103000,thin=100,burnin=3000)
plot(m4)

###posterior correlation matrix - 1 through 4 is among individual, 5-8 is within-individual
c1 <- posterior.cor(m4$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
c2 <- posterior.cor(m4$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

m5<-lmer(Local.survival..months.~Trap+(1|id), data=San2018FComb)
summary(m5) #pos
r.squaredGLMM(m5) #R2=8.627807e-30

m6<-glmer(Local.survival..months.~Cap+(1|id), data=San2018FComb, family="poisson")
summary(m6)v#pos
r.squaredGLMM(m6) #R2=0.01551943

# Males====
San2018M1<- select(San2018Behav,id, Sex, ln.captures.1, Square.root.traps.1)
San2018M1<-filter(San2018M1,Sex=="male")
San2018M1$Cap<-San2018M1$ln.captures.1
San2018M1$Trap<-San2018M1$Square.root.traps.1
San2018M1$Repeat<-rep(1, nrow(San2018M1))


San2018M2<-select(San2018Behav, id, Sex, ln.captures.2, square.root.traps2)
San2018M2<-filter(San2018M2, Sex=="male")      
San2018M2$Cap<-San2018M2$ln.captures.2
San2018M2$Trap<-San2018M2$square.root.traps2
San2018M2$Repeat<-rep(2,nrow(San2018M2))

San2018MComb<-bind_rows(San2018M1, San2018M2) #traps=exploration
San2018MComb<-left_join(San2018MComb, San2018Surv, by="id") #caps=boldness

m7<-glmer(Local.survival..months.~Cap+(1|id), family="poisson", data=San2018MComb)
summary(m7)#pos
r.squaredGLMM(m7) #0.005139763 
rpt(Cap~(1|id), grname = "id", data=San2018MComb) #r=0.676

m8<-glmer(Local.survival..months.~Trap+(1|id), family="poisson", data=San2018MComb)
summary(m8) #pos
r.squaredGLMM(m8) #0.003512940
rpt(Trap~(1|id), grname = "id", data=San2018MComb) #r=0.786