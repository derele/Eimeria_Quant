### Code to analyse
## 1) Correlation among DNA quantification methods
## 2) Course of infection by DNA quantification method 

library("ggplot2")
library("dplyr")
library("vegan")
library("gridExtra")
library("tscount")
library(gplots)
library(lme4)
library(MASS)
library("stargazer")
library("plm")
library("AER")
library("nparLD")

mytab <- read.csv("tmp/temporarydataset.csv")

# making little cheating
mytab$OPG[mytab$dpi==2] <- 0
mytab$OPG[mytab$dpi==1] <- 0


################################ Data  analyses start here #################
###########
###########
###### Question 1: Can DNA and oocyst predict weightloss & and is there an interaction between the 2?
# 1.1 Weightloss: as a rate
# 1.2 weightloss as the peak --> ALICE

#### classical glmm with ID as random effect
weightglmm <- lmer(weightloss~Genome_copies_mean*OPG + (1|EH_ID), data= mytab)
summary(weightglmm)

#### linear regression, does not include time or ID.
myweight_lm <- lm(weightloss~Genome_copies_mean*OPG, data=mytab)
coeftest(myweight_lm, vcov=vcovHC, type="HC1")

#### Linear regression with ID as fixed effect #######################
# elimninates the risk of biases due to ommited factors that vary across animals but not over time
# linear regression with ID as a fixed effect: reports dummy coefficients for each ID. which is annoying
myweight_lmf <- lm(weightloss~Genome_copies_mean*OPG + EH_ID - 1, data=mytab)
#coeftest(myweight_lm, vcov=vcovHC, type="HC1")
summary(myweight_lmf)

#OLS to the demeaned data
##obtain demeaned data
mytab_demeaned <- with(mytab, data.frame(weightloss=weightloss- ave(weightloss, EH_ID),
                                         Genome_copies_mean=Genome_copies_mean-ave(Genome_copies_mean, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                                         OPG=OPG-ave(OPG, EH_ID, FUN=function(x) mean(x, na.rm=T))))
# estimate the regression
summary(lm(weightloss~Genome_copies_mean*OPG - 1, data=mytab_demeaned))


# alternative, use the plm package, plm() uses the entity-demeaned OLS algorithm and thus does not report dummy coefficients. 
myweigh_plm_ID <- plm(weightloss~Genome_copies_mean * OPG,
               data=mytab,
               index=c("EH_ID"),
               model="within",
               effect = "individual")
coeftest(myweigh_plm_dpi, vcov. = vcovHC, type = "HC1")

# for time only
myweigh_plm_dpi <- plm(weightloss~Genome_copies_mean * OPG,
               data=mytab,
               index=c("dpi"),
               model="within",
               effect = "time")
coeftest(myweigh_plm_dpi, vcov. = vcovHC, type = "HC1")


#### Regression with time and ID fixed effects###################################
# Controlling for variables that are constant across entities but vary over time
# can be done by including time fixed effects.
# The combined model (time and ID fixed effects) allows to eliminate bias from
#unobservables that change over time but are constant over entities and it controls
#for factors that differ across entities but are constant over time. Such models
#can be estimated using the OLS algorithm 

class(mytab$EH_ID)
class(mytab$dpi)                              

mytab$dpi <- as.factor(mytab$dpi)

# with lm we get the dummy variables, which is annoying
myweigh_lmft <- lm(weightloss~Genome_copies_mean * OPG + EH_ID + dpi -1,
               data=mytab)
summary(myweigh_lmft)

# plm does not report any dummy variable which is cool
myweigh_plm <- plm(weightloss~Genome_copies_mean * OPG,
               data=mytab,
               index=c("EH_ID", "dpi"),
               model="within",
               effect = "twoways")
# obtain a summary based on clusterd standard errors
# (adjustment for autocorrelation + heteroskedasticity)
coeftest(myweigh_plm, vcov. = vcovHC, type = "HC1")

# interestingly when we include time fixed effects (controls for effects that change
#over time and not because of ID) we don't get a significant effect for OPG and the
#interaction between Genome_copies and OPG.

### non parametric analysis of longitudinal data in factorial experiments
### Brunner et al. 2002



############ Alice 1.2 question



########################################
#######################################
########### Question 2: do DNA and OOcyst predict infection stage?

mytab$inf <- NA

for (i in 1:nrow(mytab)){
    if (mytab$dpi[i]== 1|| mytab$dpi[i]==2||mytab$dpi[i]==3||mytab$dpi[i]==4) {
        mytab$inf[i] <- "early"
    } else if (mytab$dpi[i]==0){
        mytab$inf[i] <- "NI"
    } else if (mytab$dpi[i]==5||mytab$dpi[i]==6||mytab$dpi[i]==7) {
        mytab$inf[i] <- "peak"
    } else if (mytab$dpi[i]==8||mytab$dpi[i]==9||mytab$dpi[i]==10) {
        mytab$inf[i] <- "late"
    } else {mytab$inf[i] <- NA}
}


##########################################
#########################################
######### Question 3: Are DNA and Oocyst correlated?

cor.test(mytab$OPG, mytab$Genome_copies_mean, method="spearman")

plot(log10(1+ mytab$Genome_copies_mean), log(mytab$OPG))

lm(weightloss~Genome_copies_mean * OPG + EH_ID + dpi -1,
               data=mytab)

summary(myweigh_lmft)


# granger causality
# hypothesis: DNA does not cause OPG for all individuals

pgrangertest(OPG~Genome_copies_mean, data=mytab, index=c("EH_ID", "dpi"))

pgrangertest(Genome_copies_mean~OPG, data=mytab, index=c("EH_ID", "dpi"))


# plm does not report any dummy variable which is cool
myweigh_plm <- plm(weightloss~Genome_copies_mean * OPG,
               data=mytab,
               index=c("EH_ID", "dpi"),
               model="within",
               effect = "twoways")


## testing with time series: create a mean for each individualö
meantab <- sdt[,c("Genome_copies_mean", "OPG", "dpi", "weightloss")]
meantab <- unique(meantab)

head(meantab)

meantab$dpi <- as.factor(meantab$dpi)

meantab  <- aggregate(meantab[, c(1, 2, 4)], list(meantab$dpi), mean, na.action=na.exclude)

tab <- meantab %>% group_by(dpi)

mytab

mtab <- with(mytab, data.frame(weightloss=ave(weightloss, dpi),
                                         Genome_copies_mean=ave(Genome_copies_mean, dpi, FUN=function(x) mean(x, na.rm=T)),
                               OPG=ave(OPG, dpi, FUN=function(x) mean(x, na.rm=T))))

mtab <- unique(mtab)
mtab$dpi <- seq(0,10,1)

# let's just pretend that OPG at day 1 and 2 are 0
mtab[2, 3] <- 0
mtab[3, 3] <- 0


acf(na.omit(mtab$Genome_copies_mean))

acf(na.omit(mtab$OPG))

acf(na.omit(mtab$weightloss))

