#================================================
# Variance Partitioning for Pruitt 2012, Rec 204
# Elene Haave Audet, Mar 23, 2020
# Day 12 pandemic
#================================================

# 1) Set working directory to data location
setwd("G:/Shared drives/Personality & Fitness Meta-Analysis/R_Pers&FitExtractions/Data")

# 2) Load libraries
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)

# 3) Load Data
#no individual IDs given, 
#therefore cannot partition variation at among level

