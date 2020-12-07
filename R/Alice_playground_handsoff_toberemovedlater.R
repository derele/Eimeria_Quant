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

datMaxALL1 = merge(merge(datMaxOPG[c("EH_ID", "OPG", "Mouse_genotype")], 
                        datMaxDNA[c("EH_ID", "Genome_copies_gFaeces", "Mouse_genotype")]),
                  datMaxWL[c("EH_ID", "weightloss", "Mouse_genotype")])

# ----------------------------
# NB. NOW WITH TIME WINDOW
myData0 = myData # keep full for comparison
# time window, dpi 3 to dpi 10, were we have non null data
myData <- myData[myData$dpi %in% 3:10,]
# ----------------------------

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

datMaxALL2 = merge(merge(datMaxOPG[c("EH_ID", "OPG", "Mouse_genotype")], 
                         datMaxDNA[c("EH_ID", "Genome_copies_gFaeces", "Mouse_genotype")]),
                   datMaxWL[c("EH_ID", "weightloss", "Mouse_genotype")])

# Comparison with/without window
comp = merge(datMaxALL1, datMaxALL2, by = c("EH_ID", "Mouse_genotype"))

table(comp$OPG.x != comp$OPG.y) # does not change OPG results
table(comp$Genome_copies_gFaeces.x != comp$Genome_copies_gFaeces.y) # does not change DNA results
table(comp$weightloss.x != comp$weightloss.y) # change for 2 animals weight loss

comp[comp$weightloss.x != comp$weightloss.y,]

check = myData0[myData0$EH_ID %in% c("LM0210", "LM0219"),]
ggplot(check, aes(x = dpi, y = weight, col =EH_ID)) + geom_line() + geom_point()

# ----------------------------
### Q1: Is the DNA coming from the counted oocyst?
# Q1a: Correlation of the two (Susie's code)
# Q1b: Prediction of one b the other with a lag?
# Granger causality: DROPPED. For details, see "Appendix" part at the end of the code
# ----------------------------

# Q1b
# Granger causality

# https://towardsdatascience.com/fun-with-arma-var-and-granger-causality-6fdd29d8391c
# test autocorrelation
inds = unique(myData$EH_ID)
par(mfrow=c(1,2))
for (i in inds){
  series = sdtWin[myData$EH_ID %in% i,]
  series = series[order(series$dpi),"OPG"]
  # replace NA per 0 
  series[is.na(series)] = 0
  #Store acf values
  acf = acf(series, plot = F)
  #plot the series using title, x & ylabels, use lines and a marker at yaxis=0
  plot(series, main="Series-1", xlab="Time", ylab="Values")
  lines(series)
  abline(h=0, col="red")
  #Plot acf values
  plot(acf, main = as.character(i))
}
par(mfrow=c(1,1))

series = myData[sdtWin$EH_ID %in% "LM0224", c("OPG", "dpi")]
series = series[myData(series$dpi),]
# Check if the ACF plot has no bands crossing the threshold boundary except 0th lag
# (which is the correlation of the observation with itself).

# install.packages("tseries")
# library(tseries)
# adf.test computes the Augmented Dickey-Fuller test for the null that x has a unit root.
# Unit root tests are tests for stationarity in a time series. A time series has 
# stationarity if a shift in time doesn't cause a change in the shape of the distribution;
# unit roots are one cause for non-stationarity. These tests are known for having low 
# statistical power
# inds = unique(myData$EH_ID)
# for (i in inds){
#   series = myData[myData$EH_ID %in% i,]
#   series = series[order(series$dpi),"OPG"]
#   # replace NA per 0 
#   series[is.na(series)] = 0
#   # perform test
#   print(c(i, round(adf.test(series)$p.value, 4)))
# }
## rejection of stationary TS in 3/22 individuals

# for (i in inds){
#   series = sdt[sdt$EH_ID %in% i,]
#   series = series[order(series$dpi),"Genome_copies_mean"]
#   # replace NA per 0 
#   series[is.na(series)] = 0
#   # perform test
#   print(c(i, round(adf.test(series)$p.value, 4)))
# }
## no rejection of stationary TS

# https://en.wikipedia.org/wiki/Granger_causality
# Better terms: "precedence", or, as Granger himself later claimed in 1977, "temporally related".
# Rather than testing whether Y causes X, the Granger causality tests whether Y forecasts X.

# https://www.statisticshowto.com/granger-causality/#:~:text=Granger%20causality%20is%20a%20%E2%80%9Cbottom,see%20if%20they%20are%20correlated.
# Make sure your time series is stationary before proceeding. Data should be transformed to eliminate the possibility of autocorrelation. You should also make sure your model doesn’t have any unit roots, as these will skew the test results.
# The basic steps for running the test are:
# State the null hypothesis and alternate hypothesis. For example, y(t) does not Granger-cause x(t).
# Choose the lags. This mostly depends on how much data you have available. 
## One way to choose lags i and j is to run a model order test (i.e. use a model order selection method). It might be easier just to pick several values and run the Granger test several times to see if the results are the same for different lag levels. The results should not be sensitive to lags.
# Find the f-value. Two equations can be used to find if βj = 0 for all lags j

# If you have a large number of variables and lags, your F-test can lose power. An alternative would be to run a chi-square test, constructed with likelihood ratio or Wald tests. Although both versions give practically the same result, the F-test is much easier to run.

library(plm)
pgrangertest(OPG~Genome_copies_gFaeces, data=myData, index=c("EH_ID", "dpi"))

pgrangertest(OPG~Genome_copies_gFaeces, data=myData, index=c("EH_ID", "dpi"))
# Error in pgrangertest(OPG ~ Genome_copies_gFaeces, data = myData, index = c("EH_ID",  : 
# Condition for test = "Ztilde" not met for all individuals: length of time series must be larger than 5+3*order (>5+3*1=8)

# ----------------------------
# Q2b: Does DNA predict a health effect on the host “overall” (e.g. maximum)?
# ----------------------------

# Extract metrics: # what is the peak day? # what is the peak value?
ggplot(datMaxALL1, aes(x = OPG, y = Genome_copies_gFaeces)) +
  geom_point() + theme_bw()
cor.test(datMaxALL1$OPG, datMaxALL1$Genome_copies_gFaeces, method = "spearman")
# on window: does not change results
ggplot(datMaxALL2, aes(x = OPG, y = Genome_copies_gFaeces)) +
  geom_point() + theme_bw()
cor.test(datMaxALL2$OPG, datMaxALL2$Genome_copies_gFaeces, method = "spearman")

ggplot(datMaxALL1, aes(x = OPG, y = weightloss)) +
  geom_point() + theme_bw()
cor.test(datMaxALL1$OPG, datMaxALL1$weightloss, method = "spearman")
# on window: change results but not conclusion
ggplot(datMaxALL2, aes(x = OPG, y = weightloss)) +
  geom_point() + theme_bw()
cor.test(datMaxALL2$OPG, datMaxALL2$weightloss, method = "spearman")

ggplot(datMaxALL1, aes(x = Genome_copies_gFaeces, y = weightloss)) +
  geom_point() + theme_bw()
cor.test(datMaxALL1$Genome_copies_gFaeces, datMaxALL1$weightloss, method = "spearman")
# on window: change results but not conclusion
ggplot(datMaxALL2, aes(x = Genome_copies_gFaeces, y = weightloss)) +
  geom_point() + theme_bw()
cor.test(datMaxALL2$Genome_copies_gFaeces, datMaxALL2$weightloss, method = "spearman")

# Can we predict weight loss by a combination of OPG and fecDNA?
library(lmtest)
modFull = lm(weightloss ~ OPG * Genome_copies_gFaeces, data = datMaxALL1)
modminusOPG = lm(weightloss ~ Genome_copies_gFaeces, data = datMaxALL1)
modminusDNA = lm(weightloss ~ OPG, data = datMaxALL1)
modminusInter = lm(weightloss ~ OPG + Genome_copies_gFaeces, data = datMaxALL1)
list(signifOG = lrtest(modFull, modminusOPG),
     signifDNA = lrtest(modFull, modminusDNA),
     signifInter = lrtest(modFull, modminusInter))
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