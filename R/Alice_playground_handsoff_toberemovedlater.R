## Alice Balard 16 October 2020 
#  -----------------
# DNA based and classical coprological techniques provide different but
# complementary quantification – an example from rodent Coccidia (Eimeria)
#  -----------------

## We should be in the main project folder "Eimeria_Quant"
getwd()
#setwd("../")

# color blind palette: https://i.stack.imgur.com/zX6EV.png 
colorBlindPal =  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Import data by sourcing previous codes
source("R/1_Data_preparation.R")
source("R/2_qPCR_data_preparation.R") 
source("R/3_Data_analysis_qPCR_flotation.R") 
# sdt contains our data; sdt$weightloss is in % or original weight

#  -----------------
### Preparation of dataframe with all max values:
myData = sdt
# NB. Let's not consider which parent is which, but make A_B mouse = B_A mouse
# we don't have enough individuals to test this effect, and we are not interested in it anyway!
x <- strsplit(myData$Strain, "_")
y <- lapply(x, sort)
z <- unlist(lapply(y, FUN = function(x){paste(x, collapse="-")}))
myData$Mouse_genotype <- z
rm(x, y, z)
## Order the levels to be more clear in later plots (parents will be low and down 
## on the legend, hybrids in between...)
myData$Mouse_genotype <- factor(
  myData$Mouse_genotype,
  levels = c("NMRI", "WSB", "WP", "PWD1", "SCHUNT-SCHUNT", 
             "STRA-STRA", "SCHUNT-STRA", "BUSNA-STRA","PWD-SCHUNT",
             "BUSNA-PWD", "BUSNA-BUSNA", "PWD-PWD"),
  labels = c("NMRI", "MMd_F0 (Ws-Ws)", "Mmm-Mmd_Hybrid (WP)", "MMm_F0 (Pw1-Pw1)", "MMd_F0 (Sc-Sc)", 
             "MMd_F0 (St-St)", "MMd_F1 (Sc-St)", "Mmm-Mmd_F1 (Bu-St)", "Mmm-Mmd_F1 (Pw-Sc)",
             "MMm_F1 (Bu-Pw)", "MMm_F0 (Bu-Bu)", "MMm_F0 (Pw-Pw)"))
# dpi as numbers
myData$dpi = as.numeric(as.character(myData$dpi))

#  -----------------

# 1. Calculate maximum value & at which dpi it happens
# 1.1 OPG
datMaxOPG = myData %>% dplyr::group_by(EH_ID) %>%
  filter(OPG == max(OPG, na.rm = T)) %>%
  dplyr::select(EH_ID, dpi, OPG, Mouse_genotype)%>% data.frame()
table(datMaxOPG$dpi)
# 1.2 fecal DNA
datMaxDNA = myData %>% group_by(EH_ID) %>%
  dplyr::filter(Genome_copies_gFaeces == max(Genome_copies_gFaeces, na.rm = T)) %>%
  dplyr::select(EH_ID, dpi, Genome_copies_gFaeces, Mouse_genotype)%>% data.frame()
table(datMaxDNA$dpi)
# 1.3 relative weight loss
datMaxWL = myData %>% group_by(EH_ID) %>%
  filter(weightloss == max(weightloss, na.rm = T)) %>%
  dplyr::select(EH_ID, dpi, weightloss, Mouse_genotype)%>% data.frame()
table(datMaxWL$dpi)

datMaxALL = merge(merge(datMaxOPG[c("EH_ID", "OPG", "Mouse_genotype")], 
                        datMaxDNA[c("EH_ID", "Genome_copies_gFaeces", "Mouse_genotype")]),
                  datMaxWL[c("EH_ID", "weightloss", "Mouse_genotype")])

# ----------------------------
# Q2(a): Does DNA predict a health effect on the host “overall” (e.g. maximum)?
# ----------------------------

# Extract metrics: # what is the peak day? # what is the peak value?
ggplot(datMaxALL, aes(x = OPG, y = Genome_copies_gFaeces)) +
  geom_point() + theme_bw()
cor.test(datMaxALL$OPG, datMaxALL$Genome_copies_gFaeces, method = "spearman")

ggplot(datMaxALL, aes(x = OPG, y = weightloss)) +
  geom_point() + theme_bw()
cor.test(datMaxALL$OPG, datMaxALL$weightloss, method = "spearman")

ggplot(datMaxALL, aes(x = Genome_copies_gFaeces, y = weightloss)) +
  geom_point() + theme_bw()
cor.test(datMaxALL$Genome_copies_gFaeces, datMaxALL$weightloss, method = "spearman")

# Can we predict weight loss by a combination of OPG and fecDNA?
library(lmtest)
modFull = lm(weightloss ~ OPG * Genome_copies_gFaeces, data = datMaxALL)
modminusOPG = lm(weightloss ~ Genome_copies_gFaeces, data = datMaxALL)
modminusDNA = lm(weightloss ~ OPG, data = datMaxALL)
modminusInter = lm(weightloss ~ OPG + Genome_copies_gFaeces, data = datMaxALL)
list(signifOG = lrtest(modFull, modminusOPG),
     signifDNA = lrtest(modFull, modminusDNA),
     signifInter = lrtest(modFull, modminusInter))

summary(modFull)
summary(modminusInter)

par(mfrow= c(2,2))
plot(modFull) # saved as FigS_modelFit_alice.pdf
par(mfrow= c(1,1))


library(ggeffects)
ggpredict(modFull)


# Results:
# We tested if the maximum weight loss during the experiment was influenced by OPG, 
# Genome_copies_gFaeces and their interaction. To test the significance of the marginal 
# contribution of each parameter to the full model, each parameter was removed from 
# the full model, and the difference between full and reduced model was assessed using 
# likelihood ratio tests (G). 


  
# We found that the maximum weight loss during infection was influenced by
# both OPG and fecal Eimeria DNA; Their interaction was not found significant
# (LRT: OPG: G=6.6, df=2, P=0.036; 
# Genome_copies_gFaeces: G=18.2, df=2, P<0.001;
# interaction: G=2.5, df=1, P=0.11)

# plots
par(mfrow = c(2, 2))
plot(modFull) # looks OK
par(mfrow = c(1, 1)) 

# test same analyses than in paper 2 BUT with genomic data
testSignifWithinParas_updated <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- MASS::glm.nb(OPG ~ Mouse_genotype, data = dataframe)
    mod0 <- MASS::glm.nb(OPG ~ 1, data = dataframe)
  } else if (which == "RES_DNA"){
    modFULL <- MASS::glm.nb(Genome_copies_gFaeces ~ Mouse_genotype, data = dataframe)
    mod0 <- MASS::glm.nb(Genome_copies_gFaeces ~ 1, data = dataframe)
  } else if (which == "IMP"){
    modFULL <- lm(weightloss ~ Mouse_genotype, data = dataframe)
    mod0 <- lm(weightloss ~ 1, data = dataframe)
  } else if (which == "TOL"){
    modFULL <- lm(weightloss ~ 0 + OPG : Mouse_genotype, data = dataframe)
    mod0 <- lm(weightloss ~ 0 + OPG, data = dataframe)
  } else if (which == "TOL_DNA"){
    modFULL <- lm(weightloss ~ 0 + Genome_copies_gFaeces : Mouse_genotype, data = dataframe)
    mod0 <- lm(weightloss ~ 0 + Genome_copies_gFaeces, data = dataframe)
  }
  G <- lrtest(modFULL, mod0)
  return(list(modfull = modFULL, LRT = G))
}

# RES
testSignifWithinParas_updated(datMaxALL, "RES")
# RES_DNA
testSignifWithinParas_updated(datMaxALL, "RES_DNA")
# IMP
testSignifWithinParas_updated(datMaxALL, "IMP")
# TOL
testSignifWithinParas_updated(datMaxALL, "TOL")
# TOL_DNA
testSignifWithinParas_updated(datMaxALL, "TOL_DNA")

library(ggeffects)
library(cowplot)
DFpredicted = ggpredict(testSignifWithinParas_updated(datMaxALL, "RES")$modfull, terms = "Mouse_genotype")
p1=ggplot(DFpredicted, aes(x, predicted)) + geom_point() + theme_bw() +
  geom_errorbar(aes(ymin= conf.low, ymax=conf.high)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

DFpredicted2 = ggpredict(testSignifWithinParas_updated(datMaxALL, "RES_DNA")$modfull, terms = "Mouse_genotype")
p2=ggplot(DFpredicted2, aes(x, predicted)) + geom_point() + theme_bw() +
  geom_errorbar(aes(ymin= conf.low, ymax=conf.high)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot_grid(p1, p2, labels = c("A: OPG", "B: DNA"))

  # ------------------
# APPENDIX (Junk...)
# ------------------

# A simple method for distinguishing within- versus between-subjecteffects using 
# mixed modelsMartijn van de Pol*, Jonathan Wright

# https://royalsocietypublishing.org/doi/10.1098/rspb.2015.2151?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
# about trajectory analyses; very similar data to us
# Even a bit less because they had to measure within a mouse and had to use glowing bacteria.
# In terms of tolerance and resistance we might have more... DNA (intensity) and the actual surviving (tolerated) parasites
# the trajectories could help for the question 3 (later on...)
# Would you span the trajectories into 3 dimensinal space between oocysts, DNA and WL?
# yes probably, that would be very interesting to see i the trajectories correlate
# The article is well cited... especially now, five years after it appeared.
# I wonder how the hamming distance between the two dimensional vectors are interpreted... have to read on...
# The problem is (AND WHY THIS WON'T WORK), that the directions seem to be different on each day... the mouse could (lose weight&become more infected), (lose weight&become less infected), (gain weight&become more infected) or (gain weight & become less infected), or die
# This changes relatively freely in their model.
# I guess it's an up and down of both variables over the time of infection...
# In ours it would be boring as we have this clear peak dynamic.... the relative height of peak  is more important then differences in the trajectory in our case!
# I think our trajectories will likely be quite uniform (all mice first gain DNA, gain oocysts, lose weight, then all lose DNA, lose oocysts, gain weight).
# The question would be more HOW MUCH? and whether differences in the how much influence each other (or rather oocysts and/or/combined with DNA influence WL) (edited) 
# The ups and downs create trajectories in "tolerance resistance space"... those can be different in these bacteria. I don't think they will be very different in our model. I think what we are after is quantitative differences between qualitatively similar trajectories.
# Oh but it applies just for TS data
# Yep maybe it's indeed not applicable to our question. I will have a closer look when my head won't hurt ^^
# It's also only a short TS but it's more TS like in that it has less of a phasing structure than our problem.
# with phasing I mean clear UP DOWN pattern as in our data.