#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################*


# Telomere data -----------------------------------------------------------

dd$TLKB <- dd$TL/1000


# Weird variable names ----------------------------------------------------

colnames(dd)[colnames(dd) == 'OccasionDate'] <- 'CatchDate'
colnames(dd)[colnames(dd) == 'MaxOfMaxOfSeenDate'] <- 'DeathDate'
colnames(dd)[colnames(dd) == 'MinOfFieldPeriodID'] <- 'FieldPeriodID'


# Catch Year, Catch month and death year ----------------------------------

dd$CatchYear <- as.numeric(substr(dd$CatchDate,7,10))
dd$DeathYear <- as.numeric(substr(dd$DeathDate,7,10))





# Age data ----------------------------------------------------------------

dd$Age <- dd$CatchYear-dd$LayYear

#Sort out and order ageclass levels
dd <- droplevels(dd[dd$Ageclass != '',])
levels(dd$Ageclass) <- c('A','CH','FL','FL','FL','SA')
dd$Ageclass <- factor(dd$Ageclass,levels = c('CH','FL','FL','FL','SA','A'))

dd$Age[dd$Ageclass != 'A'] <- 0






# Survival and lifespan ---------------------------------------------------

dd$RemainingLife <- dd$DeathYear-dd$CatchYear
dd$SurvivedNext <- ifelse(dd$RemainingLife>0,1,0)
dd$Lifespan <- dd$DeathYear-dd$LayYear
dd$Died <- ifelse(dd$DeathYear<2013,1,0)





# Sex ---------------------------------------------------------------------

dd$Sex <- ifelse(dd$SexEstimate == 1,'Males','Females')




# Tarsus ------------------------------------------------------------------

dd$Tarsus <- NA
for(i in 1:nrow(dd))
{
  ifelse(is.na(dd$RightTarsus[i]),
         ifelse(is.na(dd$LeftTarsus[i]),
                dd$Tarsus[i] <- NA,
                dd$Tarsus[i] <- dd$LeftTarsus[i]),
         dd$Tarsus[i] <- dd$RightTarsus[i])
}

dd$RightTarsus <- NULL
dd$LeftTarsus <- NULL




# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows
terr <- subset(terr,SummerIndex>0.9) #Get rid of winter seasons

insects <- subset(insects,FieldPeriodID != 26)

# take average for year
yearmean <- tapply(terr$TQcorrected,terr$Year,mean)
terrmean <- tapply(terr$TQcorrected,terr$TerritoryID,mean)

dd$TQspace <- NA
dd$TQtime <- NA
dd$Insect <- NA
dd$InsectPrev <- NA

for(i in 1:nrow(dd))
{
  if(dd$LayYear[i] %in% names(yearmean))
  {
    dd$TQtime[i] <- yearmean[names(yearmean) == dd$LayYear[i]]
  }
  
  if(dd$TerritoryID[i] %in% names(terrmean))
  {
    dd$TQspace[i] <- terrmean[names(terrmean) == dd$TerritoryID[i]]
  }
  
  if(dd$FieldPeriodID[i] %in% insects$FieldPeriodID)
  {
    dd$Insect[i] <- insects$MeanInsects[insects$FieldPeriodID == dd$FieldPeriodID[i]]
    
    dd$InsectPrev[i] <- insects$MeanInsects[which(insects$FieldPeriodID == dd$FieldPeriodID[i])-1]
  } 
}

dd$InsectF <- ifelse(dd$Insect > 14, 'High','Low')
dd$TQspace <- log(dd$TQspace)
dd$TQtime <- log(dd$TQtime)


# Remove unwanted data ----------------------------------------------------


dd <- droplevels(subset(subset(dd,TL>1000),TL<15000))
#dd <- subset(dd,BodyMass>11)
dd <- subset(dd,LayYear <2009)


# Parentage ---------------------------------------------------------------

AvgTL <- tapply(dd$TL,dd$BirdID,mean)
nTL <- tapply(dd$TL,dd$BirdID,length)
nTL <- nTL[nTL>1]

ddlong <- subset(dd,BirdID %in% names(nTL))
temp <- unlist(lapply(ddlong$BirdID,findTL))
ddlong$cenTL <- ddlong$TL-temp

pars <- subset(pars,offspring %in% dd$BirdID)
ddpar <- subset(dd,BirdID %in% pars$offspring)

ddpar$mother <- unlist(lapply(ddpar$BirdID,findpar))
ddpar$father <- unlist(lapply(ddpar$BirdID,findpar,parent = 'father'))
ddpar$EPP <- ifelse(is.na(ddpar$father),'Extra pair','Within pair')

ddpar <- subset(ddpar,mother %in% names(AvgTL))
ddpar$mumTL <- unlist(lapply(ddpar$mother,findTL))
dads <- subset(ddpar,!(is.na(ddpar$father)))

dads <- subset(dads,father %in% names(AvgTL))
dads$dadTL <- unlist(lapply(dads$father,findTL))

chicks <- subset(ddpar,Ageclass == 'CH')
fl <- subset(ddpar,Ageclass %in% c('FL','SA'))




