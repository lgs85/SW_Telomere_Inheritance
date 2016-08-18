#TELOMERE PARENTAGE MODELS


# Maternal and paternal TL ------------------------------------------------
cor_mumdad <- with(juv,cor.test(LmumTL,LdadTL))

with(FL,cor.test(LmumTL,EmumTL))


# Sex ratio and parental TL -----------------------------------------
sr_FLmumTL <- glmer(SexEstimate ~  LmumTL  + mumage +(1|mother) + (RTL|LayYear),data = juv,family = 'binomial',na.action=na.exclude)
sr_FLmumTQ <- glmer(SexEstimate ~ mumage + TQ + (1|mother) + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)

summary(sr_FLmumTL)

sr_FLdadTL <- glmer(SexEstimate ~ LdadTL + dadage + (1|father),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLdadTQ <- glmer(SexEstimate ~ dadage + TQ + Helper + (1|father),data = FL,family = 'binomial',na.action=na.exclude)

summary(sr_FLdadTL)

sr_FLparTL <- glmer(SexEstimate ~ parTL + mumage + dadage + Helper + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLparTQ <- glmer(SexEstimate ~ TQ + mumage + dadage + Helper + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)


lm_parTLF <- lm(RTL~parTL,data=subset(FL, Sex == 'Females'),na.action = na.exclude)
lm_parTLM <- lm(RTL~parTL,data=subset(FL, Sex == 'Males'),na.action = na.exclude)



# Juvenile TL and parental adult TL ---------------------------------------
FLTL_LdadTL <- lmer(RTL~LdadTL + (1|father) + (1|LayYear),data = FL,na.action=na.exclude)
FLTL_LmumTL <- lmer(RTL~LmumTL + TQ + (1|mother) + (1|LayYear),data = FL,na.action=na.exclude)
FLTL_parTL <- lmer(RTL~parTL + Sex + Helper + TQ + (1|mother) + (1|LayYear),data = FL,na.action=na.exclude)

summary(FLTL_LdadTL)

# Telomere length and survival --------------------------------------------
FLsurv <- glmer(SurvivedNext~RTL*Sex+(1|LayYear),data=FL,family = binomial,na.action=na.exclude)
FLsurvsex <- glmer(SurvivedNext~Sex+(1|LayYear),data=FL,family = binomial,na.action=na.exclude)




# Parental body condition and sex ratio -----------------------------------
sr_FLmumcon <- glmer(SexEstimate ~  mumcon + mumage + (1|mother) + (RTL|LayYear),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLdadcon <- glmer(SexEstimate ~ dadcon + (1|father),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLparcon <- glmer(SexEstimate ~ parcon + mumage + dadage + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)

summary(sr_FLparcon)

# Maternal telomere length and territory quality --------------------------
mumTL_TQ <- lm(LmumTL~TQ,data=FL,na.action=na.exclude)

FLearly <- subset(FL,LayYear<2001)
mumTL_TQ_early <- lm(LmumTL~TQ,data=FLearly,na.action=na.exclude)

summary(mumTL_TQ_early)



# Nestling models ---------------------------------------------------------

sr_NLparTL <- glm(SexEstimate~parTL,data=NL,na.action=na.exclude)
NLTL_parTL <- lm(RTL~parTL,data=NL,na.action=na.exclude)
NLsurv <- glm(SurvivedNext~RTL,data=NL,na.action=na.exclude)
