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
c

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
dd <- subset(dd,LayYear>1992)

dd$Condition <- lm(BodyMass~Tarsus, data = dd,na.action = na.omit)$resid

# Telomere data -----------------------------------------------------------

mymed <- median(dd$TL)
dd$TLF <- ifelse(dd$TL>mymed,'Long telomeres','Short telomeres')
dd$TLKB <- dd$TL/1000

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


# Survival and lifespan ---------------------------------------------------

dd$RemainingLife <- dd$DeathYear-dd$CatchYear
dd$SurvivedNext <- ifelse(dd$RemainingLife>0,1,0)
dd$Lifespan <- dd$DeathYear-dd$LayYear
dd$Died <- ifelse(dd$DeathYear<2015,1,0)





# Sex ---------------------------------------------------------------------

dd$Sex <- ifelse(dd$SexEstimate == 1,'Males','Females')



# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows
terr$cenTQ <- NA
for(i in 1:nrow(terr))
{
  currentTQ <- subset(terr,FieldPeriodID == terr$FieldPeriodID[i])$TQcorrected
  terr$cenTQ[i] <- (terr$TQcorrected[i] - mean(currentTQ))/sd(currentTQ)
}

insects$Insectcen <- (insects$MeanInsects-mean(insects$MeanInsects))/sd(insects$MeanInsects)

# take average for year
terrmean <- tapply(terr$TQ,terr$TerritoryID,mean)


iyear <- tapply(insects$Insectcen,insects$Year,mean)

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

pars <- subset(pars,prob>0.79)

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
lates <- subset(ddDate,Ageclass =='A')

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

ParTL<- with(lates,aggregate(Condition,list(BirdID),mean))
colnames(ParTL) <- c('BirdID','TL')
ddpar$mumCon <- unlist(lapply(ddpar$mother,findTL))
ddpar$dadCon <- unlist(lapply(ddpar$father,findTL))


for(i in 1:nrow(ddpar))
{
  mumdata <- subset(dd,BirdID == ddpar$mother[i])[1,]
  daddata <- subset(dd,BirdID == ddpar$father[i])[1,]
  ddpar$mumage[i] <- ddpar$LayYear[i] - mumdata$LayYear
  ddpar$dadage[i] <- ddpar$LayYear[i] - daddata$LayYear
  ddpar$mumfood[i] <- mumdata$TQ
  ddpar$mumlife[i] <- ifelse(mumdata$Died == 1, mumdata$Lifespan, NA)
  ddpar$dadlife[i] <- ifelse(daddata$Died == 1,daddata$Lifespan, NA)
}

ddpar$LmumTLKB <- ddpar$LmumTL/1000

# Subset juveniles --------------------------------------------------------------


#juv <- droplevels(subset(ddpar,Ageclass %in% c('CH','FL','OFL','SA')))
juv <- subset(ddpar,Age<1)
adults <- droplevels(subset(dd,Ageclass == 'A'))

mymed <- mean(juv$TL,na.rm=T)
juv$TLF <- ifelse(juv$TL > mymed,'Long telomeres','Short telomeres')




# Subset Fledglings and subadults -----------------------------------------------------------

xf <- subset(juv,Status == 'XF')
juv <- subset(juv,Status!='XF')

# Helpers and social group size -------------------------------------------

status <- subset(status,BreedGroupID %in% juv$BreedGroupID)

for(i in 1:nrow(juv))
{
  currentBG <- juv$BreedGroupID[i]
  currentdata <- subset(subset(status,BreedGroupID == currentBG),BirdID !=juv$BirdID[i])
  juv$GroupSize[i] <- nrow(currentdata)
  juv$Helper[i] <- nrow(subset(currentdata,Status == 'H'))
  juv$OtherJuvs[i] <- nrow(subset(currentdata,Status %in% c('CH','FL','OFL')))
  juv$NonHelper[i] <- nrow(subset(currentdata,Status %in% c('AB','ABX')))
}

juv$EbothTL <- apply(juv[,c('EmumTL','EdadTL')],1,mean)
juv$LbothTL <- apply(juv[,c('LmumTL','LdadTL')],1,mean)

juv$TQI <- juv$TQ/(juv$Helper+1)


juv$InsectF <- ifelse(juv$Insect>0, 'Good Years','Bad Years')





# Look at telomere loss ---------------------------------------------------

x1 <- aggregate(CatchDate~BirdID,adults,min)
x2 <- merge(adults,x1)
x2 <- x2[!(duplicated(x2$BirdID)),]


x3 <- aggregate(TL~BirdID,juv,max)
x4 <- merge(x3,juv)
x4 <- x4[!(duplicated(x4$BirdID)),]

x2 <- x2[x2$BirdID %in% x4$BirdID,]
x4 <- x4[x4$BirdID %in% x2$BirdID,]

x3 <- x3[order(x3$BirdID),]
x4 <- x4[order(x4$BirdID),]

xx1 <- x4$TL-mean(x4$TL)
xx2 <- x2$TL-mean(x2$TL)          

rho <- cor(x4$TL,x2$TL)

(rho*xx1)-xx2

Loss <- data.frame(x4,
                   Loss = x4$TL-x2$TL,
                   TimeDiff = as.numeric(x2$CatchDate-x4$CatchDate),
                   D=(rho*xx1)-xx2)
Loss$TROC <- with(Loss,D/TimeDiff)



Loss <- subset(subset(Loss,TimeDiff<1460),TimeDiff>365)
 



# Sex ratio by year -------------------------------------------------------

sr <- tapply(allcatches$SexEstimate,allcatches$LayYear,mean,na.rm=T)
srn <- tapply(allcatches$SexEstimate,allcatches$LayYear,length)

sr <- sr[srn>20]
sr <- sr[names(sr)>1989]


sims <- rep(NA,5000)
upperCI <- rep(NA,length(sr))
lowerCI <- rep(NA,length(sr))

set.seed(111)

for(i in 1:length(sr))
{
  currentyear <- names(sr[i])
  n <- srn[names(srn) == currentyear]
  for(j in 1:5000) sims[j] <- mean(sample(c(0,1),n,replace = T))
  upperCI[i] <- quantile(sims,0.95)
  lowerCI[i] <- quantile(sims,0.05)
}

ddFig1 <- data.frame(Year = as.numeric(names(sr)),sr,upperCI,lowerCI,row.names = NULL)

juv2 <- subset(juv,EPP == 'Within pair')

