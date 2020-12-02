# Variance partitioning for Arroyo 2017, Rec 18
# March 25, 2020, Pandemic year :P...

# Set-up----
# Load libraries====
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(GeneNet)
library(here)

# Set wd====
dir<-here()

# 3) Load data
data<-read.table("Data/Arroyo2017_PersFit.txt", header=TRUE)

#does not appear to be a way to link individuals
#with known identity to repeated repro output
#nest ID repeated across years, but could be different individuals