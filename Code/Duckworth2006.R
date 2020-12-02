#===================
#Duckworth 2006
#==================

setwd("G:/Shared drives/Personality & Fitness Meta-Analysis/R_Pers&FitExtractions/Data")

data<-read.csv("Duckworth2006.csv")

m1<-glm(Repro~Agg, family="poisson", data=data)

cor.test(data$Agg, data$Repro)
#r=-0.7347