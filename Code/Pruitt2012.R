# Variance Partitioning for Pruitt 2012, Rec 204
# Elene Haave Audet, Mar 23, 2020
# Day 12 pandemic

# Load libraries====
library(tidyverse)
library(lme4)
library(MuMIn)
library(MCMCglmm)
library(rptR)
library(here)

# Set wd
dir<-here()

# Load Data
# no individual IDs given, 
# therefore cannot partition variation at among level
