#Hypothesis 1 - Early-life telomere length is related to parental telomere length

#Check for outliers - cleveland plots
plot(c(1:length(juv$TL))~juv$TL)
plot(c(1:length(juv$LmumTL))~juv$LmumTL)
plot(c(1:length(juv$LdadTL))~juv$LdadTL)
plot(c(1:length(juv$EmumTL))~juv$EmumTL)
plot(c(1:length(juv$EdadTL))~juv$EdadTL)

#One outlier for LmumTL - have a look
subset(juv,LmumTL>12000)$mother #932

#Take a look to see if mother has multiple catches
subset(dd,BirdID == 932) #Only caught once - can't verify so exclude

juv <- subset(juv,LmumTL<12000)

layout(matrix(1:4,2,2))
lm1 <- lm(TL~LmumTL+LdadTL+Ageclass,data=juv,na.action=na.exclude)
plot(lm1)

#Right-hand skey - try logging variables
juv$LogTL <- log10(juv$TL)
juv$LogLmumTL <- log10(juv$LmumTL)
juv$LogLdadTL <- log10(juv$LdadTL)

#Run again
lm_juvTL_LparTL <- lmer(LogTL ~ LogLmumTL * Ageclass +
                        (LogLmumTL|mother),
                        data=juv,
                        na.action=na.exclude)

standardize(lm_juvTL_LparTL)
summary(lm_juvTL_LparTL)
