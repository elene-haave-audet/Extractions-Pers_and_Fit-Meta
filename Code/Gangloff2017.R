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
data<-read.csv"Data/Gangloff2017.csv"
