#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################*


#Average repeats of blood samples
av <- ave(dd$TL,c(dd$BloodID,dd$Status))
dd$TL <- av
dd <- unique(dd)


# Weird variable names ----------------------------------------------------

colnames(dd)[colnames(dd) == 'OccasionDate'] <- 'CatchDate'
colnames(dd)[colnames(dd) == 'MaxOfMaxOfSeenDate'] <- 'DeathDate'
colnames(dd)[colnames(dd) == 'MinOfFieldPeriodID'] <- 'FieldPeriodID'


# Catch Year, Catch date and death year ----------------------------------


dd$CatchYear <- as.numeric(substr(dd$CatchDate,7,10))
dd$DeathYear <- as.numeric(substr(dd$DeathDate,7,10))
dd$CatchDate <- as.Date(dd$CatchDate,"%d/%m/%Y")

dd$Season <- ifelse(as.numeric(format(dd$CatchDate,'%m')) %in% c(4:10),
                    'Summer','Winter')

head(dd)
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



# Remove unwanted data/outliers ----------------------------------------------------


dd <- droplevels(subset(subset(dd,TL>1000),TL<15000))
dd <- subset(dd,BodyMass>5)
dd <- subset(dd,Tarsus>17)


# Telomere data -----------------------------------------------------------

dd$TLKB <- dd$TL/1000
dd$LogTL <- log(dd$TL)

dd$cenTL <- NA

#Centre telomere length by Field period
for(i in 1:nrow(dd))
{
  currentdata <- subset(dd,FieldPeriodID == FieldPeriodID[i])
  dd$cenTL[i] <- (dd$TLKB[i] - mean(currentdata$TLKB))/sd(currentdata$TLKB)
}


mymed <- median(dd$TL)
dd$TLF <- ifelse(dd$TL>mymed,'Long telomeres','Short telomeres')

# Age data ----------------------------------------------------------------

dd$Age <- dd$CatchYear-dd$LayYear


#Sort out and order ageclass levels
dd <- droplevels(dd[dd$Ageclass != '',])
levels(dd$Ageclass) <- c('A','CH','FL','FL','FL','SA')
dd$Ageclass <- factor(dd$Ageclass,levels = c('CH','FL','FL','FL','SA','A'))
dd$Fledged <- ifelse(dd$Ageclass == 'CH','Nestlings',
                     ifelse(dd$Ageclass == 'A','Adults',
                            'Fledglings'))
dd$Fledged <- factor(dd$Fledged,levels = c('Nestlings','Fledglings','Adults'))


dd$Agemonths <- ifelse(dd$Ageclass == 'CH',1,
                       ifelse(dd$Ageclass == 'FL',3,
                              ifelse(dd$Ageclass == 'SA',9,dd$Age*12)))
dd$AgemonthF <- ifelse(dd$Ageclass == 'CH','<1',
                       ifelse(dd$Ageclass == 'FL','1-9',
                              ifelse(dd$Ageclass == 'SA','9-12','>12')))
dd$AgemonthF <- factor(dd$AgemonthF,levels = c('<1','1-9','9-12','>12'))

dd <- subset(dd,Agemonths>0)

#Delta age
for(i in 1:nrow(dd))
{
  currentdata = subset(dd,BirdID == dd$BirdID[i])
  dd$MeanAge[i] <- mean(currentdata$Age)
  dd$DeltaAge[i] <- dd$Age[i]-dd$MeanAge[i]
}



# Survival and lifespan ---------------------------------------------------

dd$RemainingLife <- dd$DeathYear-dd$CatchYear
dd$SurvivedNext <- ifelse(dd$RemainingLife>0,1,0)
dd$Lifespan <- dd$DeathYear-dd$LayYear
dd$Died <- ifelse(dd$DeathYear<2013,1,0)





# Sex ---------------------------------------------------------------------

dd$Sex <- ifelse(dd$SexEstimate == 1,'Males','Females')


# Helpers and social group size -------------------------------------------

status <- subset(status,BreedGroupID %in% dd$BreedGroupID)

for(i in 1:nrow(dd))
{
  currentBG <- dd$BreedGroupID[i]
  currentdata <- subset(subset(status,BreedGroupID == currentBG),BirdID !=dd$BirdID[i])
  dd$GroupSize[i] <- nrow(currentdata)
  dd$Helper[i] <- nrow(subset(currentdata,Status == 'H'))
  dd$Otherdds[i] <- nrow(subset(currentdata,Status %in% c('CH','FL','OFL')))
  dd$NonHelper[i] <- nrow(subset(currentdata,Status %in% c('AB','ABX')))
}

dd <- subset(dd,GroupSize>0)
dd$HelperF <- ifelse(dd$Helper>0,'Helpers present','Helpers absent')


# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows


insects <- subset(insects,FieldPeriodID != 26)
insects$Insectcen <- (insects$MeanInsects-mean(insects$MeanInsects))/sd(insects$MeanInsects)

# take average for year
terrmean <- tapply(terr$TQcorrected,terr$TerritoryID,mean)

dd$TQ <- NA
dd$TQI <- NA
dd$cenTQ <- NA
dd$Insect <- NA

for(i in 1:nrow(dd))
{
  if(dd$TerritoryID[i] %in% names(terrmean))
  {
    dd$TQ[i] <- terrmean[names(terrmean) == dd$TerritoryID[i]]
    dd$TQI[i] <- dd$TQ[i]/dd$GroupSize[i]
  }
  
  if(dd$FieldPeriodID[i] %in% insects$FieldPeriodID)
  {
    dd$Insect[i] <- insects$Insectcen[insects$FieldPeriodID == dd$FieldPeriodID[i]] 
  } 
}






# Parentage ---------------------------------------------------------------

AvgTL <- tapply(dd$TL,dd$BirdID,max)
nTL <- tapply(dd$TL,dd$BirdID,length)
nTL <- nTL[nTL>1]

ddlong <- subset(dd,BirdID %in% names(nTL))
temp <- unlist(lapply(ddlong$BirdID,findTL))
ddlong$cenTL <- ddlong$TL-temp

ddpar <- subset(dd,BirdID %in% pars$offspring)

ddpar$mother <- unlist(lapply(ddpar$BirdID,findpar))
ddpar$father <- unlist(lapply(ddpar$BirdID,findpar,parent = 'father'))
ddpar$EPP <- ifelse(is.na(ddpar$father),'Extra pair','Within pair')

ddpar <- subset(ddpar,mother %in% names(AvgTL))
ddpar$mumTL <- unlist(lapply(ddpar$mother,findTL))
dadTL <- unlist(lapply(ddpar$father,findTL))
ddpar$dadTL <- NA

for(i in 1:nrow(ddpar))
{
  if(ddpar$father[i] %in% names(dadTL))
  {
    ddpar$dadTL[i] <- dadTL[names(dadTL) == ddpar$father[i]]
  }
  mumdata <- subset(dd,BirdID == ddpar$mother[i])[1,]
  daddata <- subset(dd,BirdID == ddpar$father[i])[1,]
  ddpar$mumage[i] <- ddpar$CatchYear[i] - mumdata$LayYear
  ddpar$dadage[i] <- ddpar$CatchYear[i] - daddata$LayYear
  ddpar$daddelta[i] <- 
}




# Subset juveniles --------------------------------------------------------------


juv <- droplevels(subset(ddpar,Ageclass %in% c('CH','FL','OFL','SA')))
adults <- droplevels(subset(dd,Ageclass == 'A'))

mymed <- mean(juv$cenTL,na.rm=T)
juv$TLF <- ifelse(juv$cenTL > mymed,'Long telomeres','Short telomeres')

mymed <- median(juv$dadage,na.rm=T)
juv$dadAgeF <- ifelse(juv$dadage>mymed,'Old','Young')

mymed <- median(juv$mumage,na.rm=T)
juv$mumAgeF <- ifelse(juv$mumage>mymed,'Old','Young')




# Subset Fledglings and subadults -----------------------------------------------------------

xf <- subset(juv,Status == 'XF')

juv <- subset(juv,!is.na(Tarsus))
juv <- subset(juv,!is.na(BodyMass))
juv <- subset(juv,Status!='XF')


for(i in 1:nrow(juv))
{
  currentdata <- subset(juv,FieldPeriodID == juv$FieldPeriodID[i])
  juv$cenTQ[i] <- (juv$TQ[i] - mean(currentdata$TQ))/sd(currentdata$TQ)
  juv$cenTarsus[i] <- (juv$Tarsus[i] -  mean(currentdata$Tarsus))/sd(currentdata$Tarsus)
  juv$cenHelper[i] <- (juv$Helper[i] -  mean(currentdata$Helper))/sd(currentdata$Helper)
  juv$cenNonHelper[i] <- (juv$NonHelper[i] -  mean(currentdata$NonHelper))/sd(currentdata$NonHelper)
}

juv <- subset(juv,cenTQ<4)


FlSA <- subset(juv,Ageclass!='CH')
chicks <- subset(juv,Ageclass == 'CH')

juvall <- droplevels(juv[complete.cases(juv),])
FlSAall <- droplevels(FlSA[complete.cases(FlSA),])
chickall <- droplevels(chicks[complete.cases(chicks),])



# Look at telomere loss ---------------------------------------------------

x1 <- aggregate(CatchDate~BirdID,adults,min)
x2 <- merge(adults,x1)
x2 <- x2[!(duplicated(x2$BirdID)),]


x3 <- aggregate(TLKB~BirdID,juv,max)
x4 <- merge(x3,juv)
x4 <- x4[!(duplicated(x4$BirdID)),]

x2 <- x2[x2$BirdID %in% x4$BirdID,]
x4 <- x4[x4$BirdID %in% x2$BirdID,]

x3 <- x3[order(x3$BirdID),]
x4 <- x4[order(x4$BirdID),]

xx1 <- x4$TLKB-mean(x4$TLKB)
xx2 <- x2$TLKB-mean(x2$TLKB)          

rho <- cor(x4$TLKB,x2$TLKB)

(rho*xx1)-xx2

Loss <- data.frame(x4,
                   Loss = x4$TL-x2$TL,
                   TimeDiff = as.numeric(x2$CatchDate-x4$CatchDate),
                   D=(rho*xx1)-xx2)
Loss$TROC <- with(Loss,D/TimeDiff)
for(i in 1:nrow(Loss))
{
  currentdata <- Loss[Loss$FieldPeriodID == Loss$FieldPeriodID[i],]
  Loss$cenTROC[i] <- (Loss$TROC[i]-mean(currentdata$TROC)/sd(currentdata$TROC))
}

Loss <- Loss[!(is.na(Loss$cenTROC)),]
Loss <- subset(subset(Loss,TROC > -0.01),TROC < 0.02)
Loss <- subset(subset(Loss,TimeDiff<1500),TimeDiff>365)


# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate,x1,x2,x3,x4)




# Field period average data -----------------------------------------------

juvseason <- ddply(juv,
                   .(FieldPeriodID,Season),
                   summarize,
                   TLKBmean = mean(TLKB),
                   TLKBse = se(TLKB),
                   Tarsus = mean(Tarsus),
                   Insect = mean(Insect),
                   Lifespan = mean(RemainingLife),
                   Lifespanse = se(RemainingLife),
                   CatchYear = mean(CatchYear),
                   n = length(TLKB))
juvseason$cenTL <- juvseason$TLKBmean-mean(juvseason$TLKBmean)
juvseason <- subset(juvseason,n>5)

chickseason <- ddply(chicks,
                     .(FieldPeriodID,Season),
                     summarize,
                     TLKBmean = mean(TLKB),
                     TLKBse = se(TLKB),
                     Tarsus = mean(Tarsus),
                     Insect = mean(Insect),
                     Lifespan = mean(RemainingLife),
                     Lifespanse = se(RemainingLife),
                     CatchYear = mean(CatchYear),
                     n = length(TLKB))
chickseason$cenTL <- chickseason$TLKBmean-mean(chickseason$TLKBmean)
chickseason <- subset(chickseason,n>5)


flseason <- ddply(FlSA,
                  .(FieldPeriodID,Season),
                  summarize,
                  TLKBmean = mean(TLKB),
                  TLKBse = se(TLKB),
                  Tarsus = mean(Tarsus),
                  Insect = mean(Insect),
                  Lifespan = mean(RemainingLife),
                  Lifespanse = se(RemainingLife),
                  CatchYear = mean(CatchYear),
                  n = length(TLKB))
flseason$cenTL <- flseason$TLKBmean-mean(flseason$TLKBmean)
flseason <- subset(flseason,n>5)

loss.season <- ddply(Loss,
                     .(FieldPeriodID,Season),
                     summarize,
                     TROCmean = mean(TROC),
                     TROCse = se(TROC),
                     Insect = mean(Insect),
                     Lifespan = mean(RemainingLife),
                     Lifespanse = se(RemainingLife),
                     n = length(TROC),
                     CatchYear = mean(CatchYear))
loss.season <- subset(loss.season,TROCse<0.005)



######################################################
ddpar2 <- dd
allpars <- c(pars$mother,pars$father)
allpars <- unique(allpars[!(is.na(allpars))])


ddpar2 <- subset(dd,BirdID %in% allpars)
ddpar2$numoffspring <- NA

for(i in 1:nrow(ddpar2))
{
  if(ddpar2$Sex[i] == 'Males')
  {
  currentdata <- subset(pars,father == ddpar2$BirdID[i])
  ddpar2$numoffspring[i] <- nrow(currentdata)
  } else
  {
    if(ddpar2$Sex[i] == 'Females')
    {
      currentdata <- subset(pars,mother == ddpar2$BirdID[i])
      ddpar2$numoffspring[i] <- nrow(currentdata)
    }
  }
  ddpar2$offpropmale[i] <- mean(currentdata$SexEstimate)
  
}

ddpar2 <- ddpar2[order(ddpar2$BirdID,ddpar2$CatchDate),]
ddpar2 <- ddpar2[!(duplicated(ddpar2$BirdID)),]
ddpar2$nmales <- with(ddpar2,offpropmale*numoffspring)
ddpar2$nfemales <- with(ddpar2,numoffspring-nmales)
ddpar2$sexratio <- with(ddpar2,log((nmales+1)/(nfemales+1)))

ddpar2$sexbias <- with(ddpar2,(nmales-nfemales)*numoffspring)
