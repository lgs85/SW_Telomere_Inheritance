rm(list=ls())

dd <- read.csv('Data/SW TL main data for analysis 2.csv') #main dataset
terr <- read.csv('Data/SW TL all territories.csv')
status <- read.csv('Data/SW TL all breedstatus.csv')
insects <- read.csv('Data/Insects.csv')
pars <- read.csv('Data/Hannah parentage no EPP with sex.csv')
statusdate <- read.csv('Data/SW TL breedstatus with date and Territory.csv')
allcatches <- read.csv('Data/SW_all_juv_catches.csv')

# Load relevant libraries -------------------------------------------------
library(MuMIn)
library(arm)
library(plyr)
library(ggplot2)
library(gvlma)
library(survival)
library(flexsurv)
library(car)