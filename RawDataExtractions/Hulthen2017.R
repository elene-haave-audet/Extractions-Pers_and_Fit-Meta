# Hulthen 2017 Rec 113
# Phenotypic extractions

# Load libraries====
library(here)
library(tidyverse)

# Set wd
dir<-here()

# Load data====
Hulthen<-read.csv("Data/Hulthen2017.csv")

hulthen_ind<-Hulthen %>% 
  group_by(PIT_ID) %>% 
  summarize(n_distinct(Year))

cor.test(Hulthen$Predated.by.cormorant..Yes.No., Hulthen$Boldness.score)
#r=0.101724