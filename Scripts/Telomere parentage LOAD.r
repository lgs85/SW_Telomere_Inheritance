rm(list=ls())

dd0 <- read.csv('Data/REL_TL_EB_EAF_25 Feb FULL DATA.csv') #main dataset
terr <- read.csv('Data/SW TL all territories.csv')
status <- read.csv('Data/SW TL all breedstatus.csv')
insects <- read.csv('Data/Insects.csv')
pars <- read.csv('Data/Hannah pedigree with sex and layyear.csv')
statusdate <- read.csv('Data/SW TL breedstatus with date and Territory.csv')
allcatches <- read.csv('Data/SW_all_juv_catches.csv')

# Load relevant libraries -------------------------------------------------
library(arm)
library(plyr)
library(ggplot2)
library(gvlma)
library(survival)
library(flexsurv)
library(car)