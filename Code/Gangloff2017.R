# Gangloff et al 2017, Rec 90

# Load libraries====
library(tidyverse)
library(lme4)
library(MuMIn)
library(rptR)
library(here)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/Gangloff2017.csv")

# Prep data====
#calculate number of surviving offspring after 12 months for each mother
NOffspring<-data %>% 
  group_by(MomID) %>% 
  summarise(NSurv=sum(Survive_12Months))

#take average of 4 trials of tongue flicks and escape latency
data$AvgFlicks<-rowMeans(data[,c("Flicks_1", "Flicks_2", "Flicks_3", "Flicks_4")])

data$AvgEscLatency<-rowMeans(data[,c("EscLatency_1", "EscLatency_2", "EscLatency_3", "EscLatency_4")])

# take 1 row per mom and join with offspring data
data.clean<-data %>% 
  group_by(MomID) %>% 
  slice(c(n())) %>% 
  ungroup()

data.clean<-full_join(data.clean, NOffspring, by="MomID")

# Extract estimates====
m1<-glm(NSurv~AvgFlicks, data=data.clean, family="poisson")
