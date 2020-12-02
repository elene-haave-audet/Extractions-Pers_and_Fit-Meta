# Le Gailliard 2015 Rec 143
# Phenotypic extractions

# Load libraries====
library(here)

# Set wd====
dir<-here()

LGSurv<-read.csv("Data/LG2015Behav.csv")
LGRep<-read.csv("Data/LG2015Repro.csv")

LGAct<-read.csv("Data/ActivityNewborns.csv")
LGBold<-read.csv("Data/BoldnessNewborns.csv")
LGSoc<-read.csv("Data/SociabilityNewborns.csv")
LGRecap<-read.csv("Data/Recaptures.csv")
LGRepro<-read.csv("Data/Reproduction.csv")

LGsurvF<-subset(LGSurv, Sex=="F")
LGsurvM<-subset(LGSurv, Sex=="M")

cor.test(LGsurvF$Act, LGsurvF$Surv) #r=0.05179593
cor.test(LGsurvF$Bold, LGsurvF$Surv) #r=-0.006183899 reverse sign to reflect boldness
cor.test(LGsurvF$Soc, LGsurvF$Surv) #r=0.008438001

cor.test(LGsurvM$Act, LGsurvM$Surv) #r=0.03620088
cor.test(LGsurvM$Bold, LGsurvM$Surv) #r=-0.05183443 reverse sign to reflect boldness
cor.test(LGsurvM$Soc, LGsurvM$Surv) #r=0.0749004

cor.test(LGRep$Act, LGRep$LitSize) #r=0.1089149
cor.test(LGRep$Bold, LGRep$LitSize) #r=-0.4965431 reverse sign to reflect boldness
cor.test(LGRep$Soc, LGRep$LitSize) #r=-0.1176825