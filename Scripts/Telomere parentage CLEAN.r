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



# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows


insects <- subset(insects,FieldPeriodID != 26)
insects$Insectcen <- (insects$MeanInsects-mean(insects$MeanInsects))/sd(insects$MeanInsects)

# take average for year
terrmean <- tapply(terr$TQcorrected,terr$TerritoryID,mean)

dd$TQ <- NA
dd$Insect <- NA

for(i in 1:nrow(dd))
{
  if(dd$TerritoryID[i] %in% names(terrmean))
  {
    dd$TQ[i] <- terrmean[names(terrmean) == dd$TerritoryID[i]]
  }
  
  if(dd$FieldPeriodID[i] %in% insects$FieldPeriodID)
  {
    dd$Insect[i] <- insects$Insectcen[insects$FieldPeriodID == dd$FieldPeriodID[i]] 
  } 
}




# Parentage ---------------------------------------------------------------


pars$MumAge <- with(pars,LayYear-MumLayYear)
pars$DadAge <- with(pars,LayYear-DadLayYear)



#Get parents
ddpar <- subset(dd,BirdID %in% pars$offspring)
ddpar$mother <- unlist(lapply(ddpar$BirdID,findpar))
ddpar$father <- unlist(lapply(ddpar$BirdID,findpar,parent = 'father'))
ddpar$EPP <- ifelse(is.na(ddpar$father),'Extra pair','Within pair')

#Get parental telomere lngth from both early and later life
ddDate <- dd[order(dd$BirdID,dd$CatchDate),]
earlies <- subset(ddDate,Ageclass!='A')
lates <- subset(ddDate,Ageclass == 'A')

#Early life parental TL
ParTL <- with(earlies,aggregate(TL,list(BirdID),mean))
colnames(ParTL) <- c('BirdID','TL')
ddpar$EmumTL <- unlist(lapply(ddpar$mother,findTL))
ddpar$EdadTL <- unlist(lapply(ddpar$father,findTL))

#Adult parental TL
ParTL <- with(lates,aggregate(TL,list(BirdID),mean))
colnames(ParTL) <- c('BirdID','TL')
ddpar$LmumTL <- unlist(lapply(ddpar$mother,findTL))
ddpar$LdadTL <- unlist(lapply(ddpar$father,findTL))


for(i in 1:nrow(ddpar))
{
  mumdata <- subset(dd,BirdID == ddpar$mother[i])[1,]
  daddata <- subset(dd,BirdID == ddpar$father[i])[1,]
  ddpar$mumage[i] <- ddpar$LayYear[i] - mumdata$LayYear
  ddpar$dadage[i] <- ddpar$LayYear[i] - daddata$LayYear
  ddpar$mumlife[i] <- mumdata$Lifespan
  ddpar$dadlife[i] <- daddata$Lifespan
}

# Subset juveniles --------------------------------------------------------------


#juv <- droplevels(subset(ddpar,Ageclass %in% c('CH','FL','OFL','SA')))
juv <- subset(ddpar,Age<2)
adults <- droplevels(subset(dd,Ageclass == 'A'))

mymed <- mean(juv$cenTL,na.rm=T)
juv$TLF <- ifelse(juv$cenTL > mymed,'Long telomeres','Short telomeres')

mymed <- median(juv$dadage,na.rm=T)
juv$dadAgeF <- ifelse(juv$dadage>mymed,'Old','Young')

mymed <- median(juv$mumage,na.rm=T)
juv$mumAgeF <- ifelse(juv$mumage>mymed,'Old','Young')




# Subset Fledglings and subadults -----------------------------------------------------------

xf <- subset(juv,Status == 'XF')
juv <- subset(juv,Status!='XF')

juv$mother <- factor(juv$mother)
juv$father <- factor(juv$father)



######################################################
pars <- pars[order(pars$mother),]
pars$MumOffOrd[1] <- 1

x = 1

for(i in 2:nrow(pars))
{
  if(pars$mother[i] == pars$mother[i-1])
  {
    x = x + 1
  } else
  {
    x = 1
  }
  
  pars$MumOffOrd[i] <- x
  
}



pars <- pars[order(pars$father,pars$DadAge),]
pars$DadOffOrd <- NA
pars$Numoffspring <- NA
pars$DadOffOrd[1] <- 1

x = 1

for(i in 2:nrow(pars))
{
  if(is.na(pars$father[i])) break
  if(pars$father[i] == pars$father[i-1])
  {
    x = x + 1
  } else
  {
    pars$Numoffspring[pars$father == pars$father[i-1]] <- x
    x = 1
  }
  
  pars$DadOffOrd[i] <- x
  
}

cummean <- function(x) cumsum(x)/seq_along(x)

pars$CumSex <- unlist(tapply(pars$SexEstimate,pars$father,cummean))

pars2 <- subset(pars,Numoffspring >6)

ggplot(pars2,aes(x = DadAge,y = SexEstimate,col=factor(father))) +
  geom_line()+
  theme_lgs()+
  stat_sum(geom='line')


myn <- tapply(pars$SexEstimate,pars$MumAge,length)
myn2 <- tapply(pars$SexEstimate,pars$DadAge,length)
mymean <- tapply(pars$SexEstimate,pars$MumAge,mean,na.rm=T)
mymean2 <- tapply(pars$SexEstimate,pars$DadAge,mean,na.rm=T)


mymean <- mymean[myn>5]
mymean2 <- mymean2[myn2>5]

plot(mymean,type='l',ylim = c(0,1))
points(mymean2,col='red',type='l')
abline(0.5,0,lty=2)

subset(pars,DadAge == 11)
