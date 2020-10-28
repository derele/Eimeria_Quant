## Alice Balard 16 October 2020 
#  -----------------
# DNA based and classical coprological techniques provide different but
# complementary quantification – an example from rodent Coccidia (Eimeria)
#  -----------------
### Q1: Is the DNA coming from the counted oocyst?
# Q1a: Correlation of the two
# Q1b: (How?) does this change over the time of infection
### Q2: What are DNA measurements good for, if anything?
# Q2a: Does DNA predict a health effect on the host (better than oocyst counts) on 
# a daily basis?
# Q2b: Does DNA predict a health effect on the host “overall” (e.g. maximum)?
### Q3: Is a combination of the two measurements helpful?
# Q3a: to predict at which time of an infection we are (imagine we don’t know)
# Q3b: to predict health effect
#  -----------------

## We should be in the main project folder "Eimeria_Quant"
getwd()
#setwd("../")

## Import data by sourcing previous codes
source("R/1_Data_preparation.R")
source("R/2_qPCR_data_preparation.R") # 17 oct 2020 -> error
source("R/3_Data_analysis_qPCR_flotation.R") # 17 oct 2020 -> error
# sdt contains our data
# sdt$weightloss is in % or original weight
#  -----------------

# bits for Victor to check in Data preparation:
par(mfrow = c(4,6))
for (i in 1:22) {
d <- sdt[sdt$EH_ID %in% unique(sdt$EH_ID)[i], ]
plot(as.numeric(as.character(d$dpi)), d$OPG)
}
ggplot(sdt, aes(x=dpi, y=OPG)) + geom_point() + geom_line(aes(group=EH_ID))
table(sdt$EH_ID, sdt$dpi)
## Data Article2 Alice
dataAl <- read.csv("../Article_RelatedParasitesResTol/data/ExpeDF_005_Alice.csv")
dataAl[dataAl$EH_ID %in% "LM0226", c("dpi", "OPG")]
sdt[sdt$EH_ID %in% "LM0226", c("dpi", "OPG")]
  # NB: how to have data for all days... We cannot compare
#  -----------------

# TEMPORARY Load backed-up data 
sdt <- read.csv("temp/sdt_bak.csv")
  ## Still an error, 228 only

### Q1: Is the DNA coming from the counted oocyst?
# Q1a: Correlation of the two
ggplot(sdt, aes(x=dpi, y=OPG)) + geom_point() + geom_line(aes(group=EH_ID))
ggplot(sdt, aes(x=dpi, y=Genome_copies_mean)) + geom_point() + geom_line(aes(group=EH_ID))
ggplot(sdt, aes(x=dpi, y=weightloss)) + geom_point() + geom_line(aes(group=EH_ID))

ggplot(sdt, aes(x=OPG, y=Genome_copies_mean)) + geom_point()
cor(sdt$OPG, sdt$Genome_copies_mean, use = "complete.obs", method = "pearson")
# no direct obvious correlation

# Granger causality

# https://towardsdatascience.com/fun-with-arma-var-and-granger-causality-6fdd29d8391c
# test autocorrelation
inds = unique(sdt$EH_ID)
par(mfrow=c(1,2))
for (i in inds){
  series = sdt[sdt$EH_ID %in% i,]
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

series = sdt[sdt$EH_ID %in% "LM0224", c("OPG", "dpi")]
series = series[order(series$dpi),]

# Check if the ACF plot has no bands crossing the threshold boundary except 0th lag
# (which is the correlation of the observation with itself).





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
pgrangertest(OPG~Genome_copies_mean, data=sdt, index=c("EH_ID", "dpi"))

# hypothesis: DNA does not cause OPG for all individuals
# Susie: Granger causality
# hypothesis: DNA does not cause OPG for all individuals
# pgrangertest(OPG~Genome_copies_mean, data=mytab, index=c("EH_ID", "dpi"))
# 
# Panel Granger (Non-)Causality Test (Dumitrescu/Hurlin (2012))                                                                                                                                     
# 
# data:  OPG ~ Genome_copies_mean                                                                                                                                                                           
# Ztilde = 93.114, p-value < 2.2e-16                                                                                                                                                                        
# alternative hypothesis: Granger causality for at least one individual 

# ----------------------------
# Q2b: Does DNA predict a health effect on the host “overall” (e.g. maximum)?
# ----------------------------

# Tips from stats meeting 14 oct 2020
# fit within individual, extract metric
# extract summary statistics, meaningful, then relate that to maxWL
# Cons: neglect the uncertainty; Pros:simple and biologically meaningful
# -> cf 

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
