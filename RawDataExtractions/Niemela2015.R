# Niemela 2015 rec 183
# Phenotypic extractions

# Load libraries====
library(here)
library(lme4)
library(rptR)

# Set wd====
dir<-here()

# Load data====
NimF2015<-read.csv("Data/Niemela2015_females.csv")
NimM2015<-read.csv("Data/Niemela2015_males.csv")

hist(NimF2015$Longevity) #try poisson?
hist(NimM2015$Longevity) #closer to poisson than females

hist(NimF2015$FID)
hist(NimF2015$BUTT_OUT)
hist(NimF2015$BODY_OUT)

# Models====
m1<-glmer(Longevity~ log(ACTIVITY)+ (1|INDIVIDUAL), family="poisson",data=NimF2015)
plot(m1)
summary(m1)
rpt(log(ACTIVITY)~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimF2015)
r.squaredGLMM(m1)

m2<-glmer(Longevity~log(ACTIVITY)+(1|INDIVIDUAL), family = "poisson", data=NimM2015)
summary(m2)
rpt(ACTIVITY~(1|INDIVIDUAL), grname="INDIVIDUAL", data=NimM2015, datatype = "Poisson")

hist(NimF2015$HOME_RANGE)
m3<-glmer(Longevity~HOME_RANGE+(1|INDIVIDUAL), family="poisson", data=NimF2015)
summary(m3)
rpt(HOME_RANGE~(1|INDIVIDUAL), grname = "INDIVIDUAL", data = NimF2015)
#r=0.369

m4<-glmer(Longevity~HOME_RANGE+(1|INDIVIDUAL), family="poisson", data=NimM2015)
summary(m4)
rpt(HOME_RANGE~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimM2015)
#r=0.234

hist(NimF2015$FID)
m5<-glmer(Longevity~scale(FID)+(1|INDIVIDUAL), family="poisson", data=NimF2015)
summary(m5) #switch sign since low FID indicates low boldness
rpt(FID~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimF2015)
#R  = 0.126

m6<-glmer(Longevity~scale(FID)+(1|INDIVIDUAL), family="poisson", data=NimM2015)
summary(m6) #switch sign as above
rpt(FID~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimM2015)
#r=0.163

m7<-glmer(Longevity~scale(BUTT_OUT)+(1|INDIVIDUAL), family="poisson", data=NimF2015)
summary(m7) #time to emerge after disturbance: switch sign to reflect boldness
rpt(BUTT_OUT~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimF2015)
#r=0.051

m8<-glmer(Longevity~scale(BUTT_OUT)+(1|INDIVIDUAL), family="poisson", data=NimM2015)
summary(m8)#switch sign as above
rpt(BUTT_OUT~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimM2015)
#r=0.016

m9<-glmer(Longevity~scale(BODY_OUT)+(1|INDIVIDUAL), family="poisson", data=NimF2015)
summary(m9) #swtich sign as above
rpt(BODY_OUT~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimF2015)
#r= 0.06

m10<-glmer(Longevity~scale(BODY_OUT)+(1|INDIVIDUAL), family="poisson", data=NimM2015)
summary(m10) #switch sign as above
rpt(BODY_OUT~(1|INDIVIDUAL), grname = "INDIVIDUAL", data=NimM2015)
#r=0.017