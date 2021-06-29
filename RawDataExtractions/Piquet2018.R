# Piquet 2018 rec 197
# Phenotypic extractions

# Load libraries====
library(here)
library(lme4)
library(rptR)
library(MuMIn)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/Piquet2018.csv")
data.f<-subset(data, sex==2)
data.m<-subset(data, sex==1)

#female- exploration (time spent in NE)
female.exp<-data.f %>% 
  filter(!is.na(time_novel))

hist(data.f$time_novel)
m1<-glmer(survival_2017~scale(time_novel)+(1|id), data=female.exp, family="binomial")
r.squaredGLMM(m1) #fails to converge, but did not receive that warning when I ran the model
summary(m1)
rpt(time_novel~(1|id), grname = "id", data=data.f) #r=0.465

#female- escape speed (antipredator behaviour) #REMOVED FROM EXTRACTIONS B/C IT'S A PERFORMANCE TRAIT
hist(data.f$speed)
m2<-glmer(survival_2017~speed+(1|id), data=data.f, family="binomial")
r.squaredGLMM(m2) #model failed to converge, but ran fine above
summary(m2)
rpt(speed~(1|id), grname = "id", data=data.f) #r=0.51

#female- exploration (yes/no entered NE)
m3<-glmer(survival_2017~crossornot+(1|id), data=data.f, family="binomial")
r.squaredGLMM(m3) #same as above
summary(m3)
rpt(crossornot~(1|id), grname = "id", data=data.f,datatype = "Binary") #r=

#female- breath rate (stress)
m4<-glmer(survival_2017~breath_rate+(1|id), family="binomial", data=data.f)
r.squaredGLMM(m4) #same as above
summary(m4)
rpt(breath_rate~(1|id), grname = "id", data=data.f) #r=0.168

#males
data.m<-filter(data.m, !is.na(survival_2017))

#male- exploration (time in NE)
male.explore<-data.m %>% 
  filter(!is.na(time_novel))

m5<-glmer(survival_2017~time_novel+(1|id), data=male.explore, family="binomial")
r.squaredGLMM(m5) #converged
summary(m5)
rpt(time_novel~(1|id), grname = "id", data=data.m) #r=0.558

#male- escape speed: PERFORMANCE TRAIT- REMOVED FROM TABLE
m6<-glmer(survival_2017~speed+(1|id), data=data.m, family="binomial")
r.squaredGLMM(m6)
summary(m6)
rpt(speed~(1|id), grname = "id", data = data.m) #r=0.457

#male- breath rate
male.breath<-data.m %>% 
  filter(!is.na(breath_rate))

m7<-glmer(survival_2017~breath_rate+(1|id), data=male.breath, family="binomial")
r.squaredGLMM(m7)
summary(m7)
rpt(breath_rate~(1|id), grname = "id", data=data.m) #r=0.453

#male- exploratiom (entered NE yes/no)
male.cross<-data.m %>% 
  filter(!is.na(crossornot))

m8<-glmer(survival_2017~crossornot+(1|id), data=male.cross, family="binomial")
r.squaredGLMM(m8)
summary(m8)
rpt(crossornot~(1|id), grname = "id", data=data.m, datatype = c("Binary")) #r=
