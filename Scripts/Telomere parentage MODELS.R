#TELOMERE PARENTAGE MODELS


# Maternal and paternal TL ------------------------------------------------
cor_mumdad <- with(FL,cor.test(LmumTL,LdadTL))



# Sex ratio and parental TL -----------------------------------------
sr_FLmumTL <- glmer(SexEstimate ~  LmumTLKB + mumage + Helper + (1|mother) + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLmumTQ <- glmer(SexEstimate ~ mumage + Helper + TQ + (1|mother) + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)

sr_FLdadTL <- glmer(SexEstimate ~ LdadTLKB + Helper + (1|father),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLdadTQ <- glmer(SexEstimate ~ dadage + TQ + Helper + (1|father),data = FL,family = 'binomial',na.action=na.exclude)

sr_FLparTL <- glmer(SexEstimate ~ parTL + mumage + dadage + Helper + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLparTQ <- glmer(SexEstimate ~ TQ + mumage + dadage + Helper + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)


lm_parTLF <- lm(TLKB~parTL,data=subset(FL, Sex == 'Females'),na.action = na.exclude)
lm_parTLM <- lm(TLKB~parTL,data=subset(FL, Sex == 'Males'),na.action = na.exclude)



# Juvenile TL and parental adult TL ---------------------------------------
FLTL_LdadTL <- lmer(TLKB~LdadTLKB * Sex + Helper + TQ + (1|father) + (1|LayYear),data = FL,na.action=na.exclude)
FLTL_LmumTL <- lmer(TLKB~LmumTLKB * Sex + Helper + TQ + (1|mother) + (1|LayYear),data = FL,na.action=na.exclude)
FLTL_parTL <- lmer(TLKB~parTL * Sex + Helper + TQ + (1|mother) + (1|LayYear),data = FL,na.action=na.exclude)

# Telomere length and survival --------------------------------------------
FLsurv <- glmer(SurvivedNext~TLKB*Sex+(1|LayYear),data=FL,family = binomial,na.action=na.exclude)
FLsurvsex <- glmer(SurvivedNext~Sex+(1|LayYear),data=FL,family = binomial,na.action=na.exclude)




# Parental body condition and sex ratio -----------------------------------
sr_FLmumcon <- glmer(SexEstimate ~  mumcon + mumage + Helper + (1|mother) + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLdadcon <- glmer(SexEstimate ~ dadcon + Helper + (1|father),data = FL,family = 'binomial',na.action=na.exclude)
sr_FLparcon <- glmer(SexEstimate ~ parcon + mumage + dadage + Helper + (1|LayYear),data = FL,family = 'binomial',na.action=na.exclude)


# Maternal telomere length and territory quality --------------------------
mumTL_TQ <- lm(LmumTLKB~TQ,data=FL,na.action=na.exclude)

FLearly <- subset(FL,LayYear<2001)
mumTL_TQ_early <- lm(LmumTLKB~TQ,data=FLearly,na.action=na.exclude)

summary(mumTL_TQ_early)



# Nestling models ---------------------------------------------------------

sr_NLparTL <- glm(SexEstimate~parTL,data=NL,na.action=na.exclude)
NLTL_parTL <- lm(TLKB~parTL,data=NL,na.action=na.exclude)
NLsurv <- glm(SurvivedNext~TLKB,data=NL,na.action=na.exclude)
