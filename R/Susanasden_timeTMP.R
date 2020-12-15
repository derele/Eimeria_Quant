### Code to analyse
## 1) Correlation among DNA quantification methods
## 2) Course of infection by DNA quantification method

library("ggplot2")
#library("dplyr")
library("vegan")
library("gridExtra")
#library("tscount")
#library(gplots)
library(lme4)
library(MASS)
#library("stargazer")
library("plm")
library("AER")
library("nparLD")

#mytab <- read.csv("tmp/temporarydataset.csv")

source("R/1_Data_preparation.R")
source("R/2_qPCR_data_preparation.R")

ls()

str(sdt)

sdt$dpi=as.factor(sdt$dpi)
###### Question 1: Can DNA and oocyst predict weightloss & and is there an interaction between the 2?
# 1.1 Weightloss: as a rate
# 1.2 weightloss as the peak --> ALICE

#### classical glmm with ID as random effect
weightglmm <- lmer(weightloss~Genome_copies_gFaeces*OPG + (1|EH_ID), data= sdt)
summary(weightglmm)

#### linear regression, does not include time or ID.
myweight_lm <- lm(weightloss~Genome_copies_gFaeces*OPG, data=sdt)
coeftest(myweight_lm, vcov=vcovHC, type="HC1")

#### Linear regression with ID as fixed effect #######################
# elimninates the risk of biases due to ommited factors that vary across animals but not over time
# linear regression with ID as a fixed effect: reports dummy coefficients for each ID. which is annoying
myweight_lmf <- lm(weightloss~Genome_copies_gFaeces*OPG + EH_ID - 1, data=sdt)
#coeftest(myweight_lm, vcov=vcovHC, type="HC1")
summary(myweight_lmf)

#OLS to the demeaned data
##obtain subject demeaned data
sdt_Sdemeaned <- with(sdt, data.frame(weightloss=weightloss- ave(weightloss, EH_ID),
                      Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                      OPG=OPG-ave(OPG, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                      dpi=dpi,
                      EH_ID=EH_ID))
# estimate the regression
summary(lm(weightloss~Genome_copies_gFaeces*OPG, data=sdt_Sdemeaned))

# time demeaned data
sdt_Tdemeaned <- with(sdt, data.frame(weightloss=weightloss- ave(weightloss, dpi),
                                         Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, dpi, FUN=function(x) mean(x, na.rm=T)),
                                         OPG=OPG-ave(OPG, dpi, FUN=function(x) mean(x, na.rm=T))))
# estimate the regression
summary(lm(weightloss~Genome_copies_gFaeces*OPG, data=sdt_Tdemeaned))


#### Regression with time and ID fixed effects###################################
# Controlling for variables that are constant across entities but vary over time
# can be done by including time fixed effects.
# The combined model (time and ID fixed effects) allows to eliminate bias from
#unobservables that change over time but are constant over entities and it controls
#for factors that differ across entities but are constant over time. Such models
#can be estimated using the OLS algorithm

sdt_STdemeaned <- with(sdt_Sdemeaned, data.frame(weightloss=weightloss- ave(weightloss, dpi),
                      Genome_copies_gFaeces=Genome_copies_gFaeces-ave(Genome_copies_gFaeces, dpi, FUN=function(x) mean(x, na.rm=T)),
                      OPG=OPG-ave(OPG, dpi, FUN=function(x) mean(x, na.rm=T)),
                      dpi=dpi,
                      EH_ID=EH_ID))


summary(sdt_STdemeaned)
summary(sdt_Sdemeaned)

cor.test(sdt_STdemeaned$Genome_copies_gFaeces, sdt_STdemeaned$OPG, method="spearman")

cor.test (sdt$Genome_copies_gFaeces, sdt$OPG, method="spearman")
nrow((na.omit(sdt[,c("Genome_copies_gFaeces", "OPG")])))

nrow((na.omit(sdt[,c("Genome_copies_gFaeces", "OPG", "weightloss")])))


# plotting correlations
jpeg("fig/GG_OPG_cor.jpeg",
     width = 5, height = 6, units = "in", pointsize = 10,
     res = 500)
ggplot(sdt, aes(x=log(1+Genome_copies_gFaeces), y=log(1+OPG), colour=dpi))+
    geom_point(size=2, alpha=0.8)+
    annotate("text", x=10, y=15, label="rho=0.72, p<0.001", hjust = "left")+
    labs(x="Genome copies (log+1)", y="OPG (log+1)")+
    theme_classic()
dev.off()

log10(sdt_STdemeaned$Genome_copies_gFaeces+max(na.omit(sdt_STdemeaned$Genome_copies_gFaeces)))

sdt_STdemeaned$Genome_copies_gFaeces

ggplot(sdt_STdemeaned, aes(x=(Genome_copies_gFaeces+ max(na.omit(Genome_copies_gFaeces))), y=(OPG+max(na.omit(OPG)))))+
  geom_point(size=5, alpha=0.5)+
  scale_y_log10(name = "log10 (Oocyst per gram faeces + 1) \n (Flotation)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(name = "log10 (Genome copies per gram faeces + 1) \n (qPCR)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  labs(tag= "C)")+
  theme_bw()+
  stat_cor(label.x = 9.99,  label.y = 6.1, method = "spearman",
           aes(label= paste(..r.., ..p.label.., sep= "~`,`~")))+
  theme(text = element_text(size=16))+
  annotation_logticks()+
  geom_text (x = 9.9, y = 6.1, show.legend = F,
             label = paste ("Spearman's rho ="))-> C

C

####OPGs modeled by Genome copies
sdt%>%
  ggplot(aes(Genome_copies_gFaeces+1, OPG+1))+
  geom_smooth(method = lm, col= "black")+
  scale_y_log10(name = "log10 (Oocyst per gram faeces + 1) \n (Flotation)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(name = "log10 (Genome copies per gram faeces + 1) \n (qPCR)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  stat_cor(label.x = 4.5,  label.y = 6, method = "spearman",
           aes(label= paste(..r.., ..p.label.., sep= "~`,`~")))+
  theme(text = element_text(size=16))+
  annotation_logticks()+
  geom_text (x = 3.75, y = 6.05, show.legend = F,
             label = paste ("Spearman's rho ="))-> A

sdt%>%
  ggplot(aes(Genome_copies_gFaeces+1, OPG+1, fill=dpi))+
  geom_smooth(method = lm, se=FALSE, aes(Genome_copies_gFaeces+1, OPG+0.1, color=dpi))+
  scale_y_log10(name = "log10 (Oocyst per gram faeces + max) \n (Flotation)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(name = "log10 (Genome copies per gram faeces + max) \n (qPCR)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_point(shape=21, size=5) +
  annotation_logticks()+
  theme_bw()+
  labs(tag= "B)")+
  theme(text = element_text(size=16)) -> B

##Figure 4# Spearman's Correlation between genome copies and OPG overall and by dpi
pdf(file = "fig/Figure_4abc.pdf", width = 10, height = 20)
grid.arrange(A,B,C)
dev.off()



# estimate the regression
sdtST=na.omit(sdt[,c("Genome_copies_gFaeces", "OPG", "weightloss")])
ST.lm=lm(weightloss~Genome_copies_gFaeces*OPG, data=sdtST)
int.lm=lm(weightloss~1, data=sdtST)

OPG.lm=lm(weightloss~Genome_copies_gFaeces, data=sdtST)
GC.lm=lm(weightloss~OPG, data=sdtST)
I.lm=lm(weightloss~Genome_copies_gFaeces+OPG, data=sdtST)

anova(ST.lm, GC.lm)
anova(ST.lm, OPG.lm)
anova(ST.lm, I.lm)

summary(ST.lm)
anova(ST.lm, int.lm)

#plot(ST.lm)

class(sdt$EH_ID)
class(sdt$dpi)

# plot weight loss and genome copies
jpeg("fig/GG_weightloss.jpeg",
     width = 5, height = 6, units = "in", pointsize = 10,
     res = 500)
ggplot(sdtST, aes(x=log(1+Genome_copies_gFaeces), y=weightloss))+
    geom_point(size=2, alpha=0.8)+
    annotate("text", x=10, y=15, label="F=25.1, p<0.001", hjust = "left")+
    labs(x="Genome copies (log+1)", y="Weight loss")+
    theme_classic()
dev.off()

jpeg("fig/OPG_weightloss.jpeg",
     width = 5, height = 6, units = "in", pointsize = 10,
     res = 500)
ggplot(sdtST, aes(x=log(1+OPG), y=weightloss))+
    geom_point(size=2, alpha=0.8)+
    annotate("text", x=10, y=15, label="F=3.9, p=0.05", hjust = "left")+
    labs(x="OPG (log 1+)", y="Weight loss")+
    theme_classic()
dev.off()


# plot residuals
sdtST$predicted <- predict(ST.lm)
sdtST$residuals <- residuals(ST.lm)
sdtST$logGC <- log(1+sdtST$Genome_copies_gFaeces)
library(tidyr)

sdtST %>%
    gather(key="iv", value = "x", -logGC, -predicted, -residuals) %>%
    ggplot(aes(x=x, y=logGC))+
    geom_segment(aes(xend=x, yend=predicted), alpha=0.2)+
    geom_point(aes(color=residuals))+
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    guides(color = FALSE) +
    geom_point(aes(y = predicted), shape = 1) +
     facet_grid(~ iv, scales = "free_x") +
     theme_bw()

ggplot(sdtST, aes(x = logGC, y = weightloss)) +  # Set up canvas with outcome variable on y-axis
    geom_segment(aes(xend = logGC, yend = predicted), alpha = .2) +  # alpha to fade lines
    geom_point() +
    geom_point(aes(y = predicted), shape = 1) +
    theme_classic()  # Add theme for cleaner look

# with lm we get the dummy variables, which is annoying
myweigh_lmft <- lm(weightloss~Genome_copies_mean * OPG + EH_ID + dpi -1,
               data=mytab)
summary(myweigh_lmft)

# interestingly when we include time fixed effects (controls for effects that change
#over time and not because of ID) we don't get a significant effect for OPG and the
#interaction between Genome_copies and OPG.

### non parametric analysis of longitudinal data in factorial experiments
### Brunner et al. 2002
# LD-F1 design fers to the experimental design with one sub-plot factor
#(longitudinal data forone homogeneous group of subjects).
#1 hypothesis: no time effect

ex.f1 <- ld.f1(y=mytab$weightloss, time=mytab$dpi, subject=mytab$EH_ID, time.order=c(0,1,2,3,4,5,6,7,8,9,10))

summary(ex.f1)
plot(ex.f1)
print(ex.f1)



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
# spearman correlation
cor.test(mytab$OPG, mytab$Genome_copies_mean, method="spearman")

summary(mytab$Genome_copies_mean)

# time demeaned data

# estimate the regression
summary(lm(OPG~Genome_copies_gFaeces - 1, data=sdtST))

#corr <- ggplot(mytab, aes(x= log(1+Genome_copies_mean), y= log(1+OPG))) +
#    geom_point(aes(colour=factor(dpi)), size = 2, alpha=0.8)+
#    xlab("qpcr DNA, log(+ 1)") +
#    ylab("Oocysts, log(+ 1)") +
#    theme_classic()

#corr_t <- ggplot(mytab_td, aes(x= log(max(Genome_copies_mean, na.rm=T) + Genome_copies_mean), y= log(max(OPG) + OPG))) +
#    geom_point(size = 2, alpha=0.8)+
#    xlab("qpcr DNA, log(+ max)") +
#    ylab("Oocysts, log(+ max)") +
#    theme_classic()

#corr_t

#jpeg("fig/Figure_cor.jpeg",
#     width = 6, height = 5, units = "in", pointsize = 10,
#     res = 500)
#ggarrange(corr, corr_t, labels = c("a", "b"))
#dev.off()

library(ggpubr)

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

##### to do:
## stage of infection
## residual plots for LM's for supplementary material
