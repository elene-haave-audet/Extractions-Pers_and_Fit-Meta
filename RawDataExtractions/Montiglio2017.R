# Montiglio 2017 Rec 170
# Phenotypic extractions

# Load libraries====
library(here)
library(tidyverse)

# Set wd====
dir<-here()

# Load data====
data<-read.csv("Data/data170_mating_rate.csv")

# Prep data====
proportion= data %>% 
  group_by(id_exp, mating) %>% 
  summarise(n=n()) %>% 
  mutate(rel.freq=n/sum(n)) %>% 
  mutate(p_mate= if_else(mating==0 & rel.freq==1, 0, rel.freq)) %>% 
  subset(mating==1| p_mate==0)


ind=data %>% 
  group_by(id_exp) %>% 
  summarize(Male=n_distinct(id_exp)) #n=295

Activity=data %>% 
  group_by(id_exp) %>% 
  summarize(Act=max(coef.activity))

Agg=data %>% 
  group_by(id_exp) %>% 
  summarize(Agg=max(coef.agg))

data2<-merge(Activity, Agg, by="id_exp")  
data3<-merge(data2, proportion, by="id_exp")

# correlations====
hist(data3$p_mate)
cor.test(data3$Act, data3$p_mate) #r=0.1663807
cor.test(data3$Agg, data3$p_mate) #r=0.1976249