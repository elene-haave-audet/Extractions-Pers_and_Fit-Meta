# Phenotypic correlation for Reale & Festa-Bianchet 2003
# April 28, 2020

# Load libraries====
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(GeneNet)
library(here)

# Set wd
dir<-here()

# Load data====
data<-read.csv("Data/Reale2003.csv")

# Clean and separate data====
## survival
## remove "T"--> translocated individuals for which surv is unknown
# get last observation for each individual (i.e. dead or alive in '98)
surv<-data %>% 
  filter(type_of_data=="Survival") %>% 
  filter(Survival!="T") %>% 
  group_by(ID) %>% 
  slice(n()) %>% 
  select(ID, Survival)

## boldness (number of times captured in a season)
bold<-data %>% 
  filter(type_of_data=="Boldness") %>% 
  select(ID, Year, NB_captures)
## some females died before 97, so do not have info on survival
## join survival and boldness data
bold_join<-left_join(bold, surv, by="ID") #remove NAs for survival
bold_join<-filter(bold_join, !is.na(Survival))  

## docility (female measured several times in a given year)
doc<-data %>% 
  filter(type_of_data=="docility") %>% 
  select(ID, Year, Docility)
## join survival and dociliy data
doc_join<-left_join(doc, surv, by="ID")

# Run glms to find cor b/w fitness and behav====
## survival~boldness+(1|ID)
m1<-glmer(Survival~scale(NB_captures)+(1|ID), data=bold_join, family="binomial")
summary(m1) #NInd=67, Nobs=288 sign=pos
r.squaredGLMM(m1) #R2=9.06e-18
rpt(NB_captures~(1|ID), grname = "ID", data=bold_join, datatype = c("Poisson"))
#r=0

## survival~docility+(1|ID)
m2<-glmer(Survival~Docility+(1|ID), data=doc_join, family="binomial")
summary(m2) #NInd=56, Nobs=304 sign=pos, but will be reversed to reflect stress
r.squaredGLMM(m2) #R2=1.0587e-17
rpt(Docility~(1|ID), grname = "ID", data=doc_join, datatype = c("Poisson"))
hist(doc_join$Docility)
#r=0.696
