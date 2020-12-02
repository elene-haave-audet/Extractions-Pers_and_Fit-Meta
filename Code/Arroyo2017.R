#===============================================
# Variance partitioning for Arroyo 2017, Rec 18
# March 25, 2020, Pandemic year :P...
#===============================================

# 1) Set working directory to data location
setwd("G:/Shared drives/Personality & Fitness Meta-Analysis/R_Pers&FitExtractions/Data")

# 2) Load libraries
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(GeneNet)

# 3) Load data
data<-read.table("Arroyo2017_PersFit.txt", header=TRUE)

#does not appear to be a way to link individuals
#with known identity to repeated repro output
#nest ID repeated across years, but could be different individuals