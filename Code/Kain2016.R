# Kain & McCoy 2016 rec 324
# Phenotypic estimates

# Load libraries====
library(here)
library(lme4)
library(MuMIn)
library(tidyverse)

# Set wd====
dir<-(here)

# Load data====
Kain<-read.csv("Data/Kain2016.csv")

# Clean data====
#create a variable for survival at time of observation
Kain=mutate(Kain, survival=if_else(activity=="D", 0, 1))

#subset for dead individuals to extract time of death
dead<-subset(Kain, survival==0)

dead.time=dead %>% 
  group_by(individual) %>% 
  dplyr::summarize(time=first(time)) #now I have time of death for all the individuals that died

#keep the last observation from each individual to get survival df

count(Kain, individual)#60 individuals
table(survival$survival) #51 dead, 9 live

survival=Kain %>% 
  select(individual, survival, time) %>% 
  group_by(individual) %>% 
  slice(c(n())) %>% 
  ungroup() #end w/ 60 observations

#all individuals have time 2820 at end of experiment; remove dead animals and replace with time of death
live<-subset(survival, survival==1)
live<-select(live, individual,time)

#now combine live and dead animals to have time for all individuals
longevity<-as.data.frame(rbind(live, dead.time))

behav=left_join(longevity,Kain, by="individual")
behav$time_death<-(behav$time.x+120) #make all values positive by 
#adding 120 minnutes, which was length of observation period prior 
#adding predator cue

stoch<-subset(behav, treatment=="stochastic")
const<-subset(behav, treatment=="constant")
ctrl<-subset(behav, treatment=="control")
high<-subset(behav, treatment=="high")

hist(behav$location)
hist(behav$time_death)

# stochastic treatment####
m1<-glmer(time_death~scale(location)+(1|individual),data=stoch, family="poisson")
r.squaredGLMM(m1) #R2= 4.864433e-10
rpt(location~(1|individual), grname = "individual", data=stoch) #r=0.255
summary(m1)

# constant treatment####
m2<-glmer(time_death~scale(location)+ (1|individual), family="poisson", data=const)
r.squaredGLMM(m2) #R2=5.923116e-11
summary(m2)
rpt(location~(1|individual), grname = "individual", data=const) #r=0.305

# control treatment####
m3<-glmer(time_death~scale(location)+(1|individual), family="poisson", data=ctrl)
r.squaredGLMM(m3) #R2= 1.069043e-10
summary(m3)
rpt(location~(1|individual), grname="individual", data=ctrl)
#r=0.198

# high treatment####
m4<-glmer(time_death~scale(location)+(1|individual), family="poisson", data=high)
r.squaredGLMM(m4) #R2=2.028931e-09
summary(m4)
rpt(location~(1|individual), grname = "individual", data=high) #r=0.149
