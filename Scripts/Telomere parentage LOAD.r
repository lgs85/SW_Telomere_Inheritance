rm(list=ls())
dd <- read.csv('Data/SW TL main data for analysis.csv')
pars <- read.csv('Data/Hannah parentage no EPP.csv')
TLmaster <- read.csv('Data/Master telomeres Emma and Ellie.csv')
insects <- read.csv('Data/Insects.csv')
terr <- read.csv('Data/SW TL all territories.csv')

library(knitr)
library(ggplot2)
