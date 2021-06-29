# Monceau 2017 Rec 166
# Phenotypic extractions

# Load libraries====
library(here)
library(lme4)
library(MuMIn)

# Set wd====
dir<-here()

# Load data====
Monc<-read.csv("Data/Monceau2017.csv")
MoncF<-subset(Monc, sex=="f")
MoncM<-subset(Monc, sex=="m")

# Female survival====
hist(MoncF$survival) #not particularly normal
MoncF$survival_t<-(MoncF$survival)^2
hist(MoncF$survival_t)
hist(MoncM$survival) #also not particulrly normal
hist(MoncF$activity)
hist(Monc$larvae) #not too bad for normality

m1<-glmer(survival_t~ scale(activity)+ (1|id), family="poisson", data=MoncF)
plot(m1)
summary(m1) #beta=0.0001763 SE= 0.0050638 obs=82, ind=41 

hist(MoncF$exploration)
MoncF$exploration_t<-sqrt(MoncF$exploration)
hist(MoncF$exploration_t)

m2<-glmer(survival_t~exploration_t+ (1|id), family="poisson", data=MoncF)
plot(m2)
summary(m2) #beta=0.0000328 SE= 0.002795

hist(MoncF$food_neophobia)
MoncF$food_neophobia_t<-log(MoncF$food_neophobia)
hist(MoncF$food_neophobia_t)

m3<-glmer(survival_t~food_neophobia_t+ (1|id), family="poisson", data=MoncF)
plot(m3)
summary(m3) #beta=0.00003873 SE= 0.002869 switch sign to make low neophobia reference

hist(MoncF$gregariousness)
MoncF$gregariousness_t<-log(MoncF$gregariousness)
hist(MoncF$gregariousness_t)

m4<-glmer(survival_t~gregariousness_t+ (1|id), family="poisson", data=MoncF)
plot(m4)
summary(m4) #beta= -0.00005551 SE= 0.002279 reverse sign since small number indicates 'highly gregarious'

# Male survival====
hist(MoncM$survival)
MoncM$survival_t<-(MoncM$survival)^2
hist(MoncM$survival_t)

hist(MoncM$activity)
MoncM$activity_t<-(MoncM$activity)^1/2
hist(MoncM$activity_t)

m5<-glmer(survival_t~scale(activity_t)+ (1|id), family="poisson", data=MoncM)
plot(m5)
summary(m5) #beta= -0.00004055 SE= 0.003508 obs=82 ind=41

hist(MoncM$exploration)

m6<-glmer(survival_t~scale(exploration)+ (1|id), family="poisson", data=MoncM)
plot(m6)
summary(m6) #beta= 0.00001448 SE= 0.002366
r.squaredGLMM(m6)


m7<-glmer(survival_t~scale(food_neophobia)+ (1|id), family="poisson", data=MoncM)
plot(m7)
summary(m7) #beta= -0.00001944 SE= 0.002642 switch sign to make low neophobia a higher value

m8<-glmer(survival_t~scale(gregariousness)+ (1|id), family="poisson", data=MoncM)
plot(m8)
summary(m8) #beta= -0.0001381 SE= 0.0041662 switch sign because high value=less sociable

# Female repro====
m9<-lmer(larvae~activity+(1|id), data=MoncF)
summary(m9) #beta=0.00 se=0.00

m11<-lmer(larvae~exploration+(1|id), data=MoncF)
summary(m11) #beta=0.00 se=0.00

m13<-glmer(larvae~scale(gregariousness_t)+(1|id), data=MoncF, family="poisson")
summary(m13) #beta=-0.002478 SE= 0.022997 switch sign to reflect sociability
hist(MoncF$gregariousness)
MoncF$gregariousness_t<-log((MoncF$gregariousness)+0.5)
hist(MoncF$gregariousness_t)

m15<-glmer(larvae~food_neophobia_t+(1|id), data=MoncF, family="poisson")
summary(m15) #beta=0.0005 se=0.028
hist(MoncF$food_neophobia)
MoncF$food_neophobia_t<-log((MoncF$food_neophobia)+0.5)
hist(MoncF$food_neophobia_t)

# Male repro====
m10<-glmer(larvae~activity_t+(1|id), data=MoncM, family="poisson")
summary(m10) #beta=-0.0003 se=0.0016
hist(MoncM$activity)
MoncM$activity_t<-(MoncM$activity)^1/2
hist(MoncM$activity_t)

m12<-lmer(larvae~exploration+(1|id), data=MoncM)
summary(m12) #beta=0.00 se=0.00

m14<-glmer(larvae~gregariousness_t+(1|id), data=MoncM, family="poisson")
summary(m14) #beta=0.0033737 SE=0.023401 switch sign to reflect sociability
MoncM$gregariousness_t<-log((MoncM$gregariousness)+0.5)
hist(MoncM$gregariousness_t)

m16<-glmer(larvae~food_neophobia_t+(1|id), data=MoncM, family="poisson")
summary(m16) #beta=0.004 se=0.0257
MoncM$food_neophobia_t<-log((MoncM$food_neophobia)+0.5)
hist(MoncM$food_neophobia_t)