# Debecker 2016 Rec 70
# Phenotypic extractions

# Load libraries====
library(here)

# Set wd====
dir<-here()

# D. elegans====
data.e<-read.csv("Data/DebeckerElegansCont.csv")
data.ef<-subset(data.e, Sex=="Female")
data.em<-subset(data.e, Sex=="Male")

cor.test(data.ef$Avg_Projection, data.ef$Adult_lifespan) #r=-0.1638117 
cor.test(data.ef$Avg_Latency, data.ef$Adult_lifespan) #r=-0.1349558 (switch sign)
cor.test(data.em$Avg_Projection, data.em$Adult_lifespan) #r=0.4162892
cor.test(data.em$Avg_Latency, data.em$Adult_lifespan) #r=-0.1686977 (switch sign)

# D. genei====
data.g<-read.csv("DebeckerGeneiCont.csv")
data.gf<-subset(data.g, Sex=="Female")
data.gm<-subset(data.g, Sex=="Male")

cor.test(data.gf$Avg_Projection, data.gf$Adult_lifespan) #r=-0.4442666
cor.test(data.gf$Avg_Latency, data.gf$Adult_lifespan) #r=0.482926 (switch sign)
cor.test(data.gm$Avg_Projection, data.gm$Adult_lifespan) #r=-0.004833125
cor.test(data.gm$Avg_Latency, data.gm$Adult_lifespan) #r=0.0625705 (switch sign)

data.gr<-read.csv("DebeckerGraellsiiCont.csv")
data.grf<-subset(data.gr, Sex=="Female")
data.grm<-subset(data.gr, Sex=="Male")

cor.test(data.grf$Avg_Projection, data.grf$Adult_lifespan) #r=-0.126299
cor.test(data.grf$Avg_Latency, data.grf$Adult_lifespan) #r=-0.3290612 (switch sign)
cor.test(data.grm$Avg_Projection, data.grm$Adult_lifespan) #r=0.08668952
cor.test(data.grm$Avg_Latency, data.grm$Adult_lifespan) #r=-0.09584818 (switch sign)

# D. pumillio====
data.p<-read.csv("DebeckerPumillioCont.csv")
data.pf<-subset(data.p, Sex=="Female")
data.pm<-subset(data.p, Sex=="Male")

cor.test(data.pf$Avg_Projection, data.pf$Adult_lifespan) #r=-0.5240291
cor.test(data.pf$Avg_Latency, data.pf$Adult_lifespan) #r=-0.1716046 (switch sign)
cor.test(data.pm$Avg_Projection, data.pm$Adult_lifespan) #r=-0.3574302
cor.test(data.pm$Avg_Latency, data.pm$Adult_lifespan) #r=-0.5183914 (switch sign)
