# Zavorka 2016

# Libraries====
library(here)

# Set wd====
dir<-here()

# data====
zavorka<-read.csv("data_Sten2014.csv")
str(zavorka)

#Logistic regression
glm.fit<-glm(final_survival~pca_dim1_activ, data=zavorka, family="binomial")
summary(glm.fit)

#the estimates given are log odd ratio, to get odd ratio only
exp(coef(glm.fit))