## Alice Balard 16 October 2020 
#  -----------------
# DNA based and classical coprological techniques provide different but
# complementary quantification – an example from rodent Coccidia (Eimeria)
#  -----------------

## We should be in the main project folder "Eimeria_Quant"
getwd()
#setwd("../")

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

modNull = lm(weightloss ~ 1, data = datMaxALL)
modFull = lm(weightloss ~ OPG * Genome_copies_gFaeces, data = datMaxALL)
modminusOPG = lm(weightloss ~ Genome_copies_gFaeces, data = datMaxALL)
modminusDNA = lm(weightloss ~ OPG, data = datMaxALL)
modminusInter = lm(weightloss ~ OPG + Genome_copies_gFaeces, data = datMaxALL)

# test difference LRT or anova
anova(modFull, modNull)
lrtest(modFull, modNull)
anova(modFull, modNull, test ="LRT")
anova(modFull, modNull, test ="Chisq")
# Homebrew log-likelihood test
like.diff = logLik(modFull) - logLik(modNull)
df.diff = modNull$df.residual - modFull$df.residual
pchisq(as.numeric(like.diff) * 2, df=df.diff, lower.tail=F)
  
summary(modFull)
summary(modminusInter)

## Conclusion
# anova(mod1, mod2) -> Wald test, we report a F-statistics
# anova(mod1, mod2, test="LRT" or "Chisq") -> LRT but computed with a certain method (sum of square whatever)
# lrtest(mod1, mod2) -> LRT computed with a second method
# we go for lrtest and report a Chisq

# Results:
list(signifFull = lrtest(modFull, modNull),
     signifOG = lrtest(modFull, modminusOPG),
     signifDNA = lrtest(modFull, modminusDNA),
     signifInter = lrtest(modFull, modminusInter))

par(mfrow= c(2,2))
#plot(modFull) # saved as FigS_modelFit_alice.pdf
par(mfrow= c(1,1))

# plot prediction
d <- datMaxALL[c("OPG", "Genome_copies_gFaeces", "weightloss")]

d$predicted <- predict(modFull)   # Save the predicted values
d$residuals <- residuals(modFull) # Save the residual values

#https://drsimonj.svbtle.com/visualising-residuals
#pdf(file = "fig/plotResidAlice_temp.pdf", width = 8, height = 5)
plotResidAlice_temp <- d %>% 
  gather(key = "iv", value = "x", -weightloss, -predicted, -residuals) %>%  # Get data into shape
  ggplot(aes(x = x, y = weightloss)) +  # Note use of `x` here and next line
  geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
  geom_point(aes(color = residuals)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  guides(color = FALSE) +
  geom_point(aes(y = predicted), shape = 1) +
  facet_grid(~ iv, scales = "free_x") +  # Split panels here by `iv`
  theme_bw()
#dev.off()

# save figure in a temp directory
saveRDS(plotResidAlice_temp, "fig/plotResidAlice_temp.RDS")

# relative importance of predictors
library(relaimpo)

calc.relimp(modminusInter)
calc.relimp(modminusInter, rela=TRUE)

# The total proportion of variance explained by the model with all two predictors is 59.53%. 
# OPG contributed to 13%, and DNA for 46%. After normalisation, we found that the proportion of contribution of each predictor to the overall 
# R2 is 78% for DNA, and 22% for OPG. 

#For two predictors, after we get their relative importance measured by R2
#, we might want to test whether one predictor is significantly more important than the other. However, unlike t-test, it is rather difficult to find an analytical test statistic for a test. Instead, bootstrap can be used. The package relaimpo includes two functions -- boot.relimp() and booteval.relimp() -- for the task.  The first function conducts the bootstrap and the second one gets the confidence intervals.

bootresults<-boot.relimp(modminusInter, b=1000) 
ci<-booteval.relimp(bootresults, norank=T)
ci
plot(ci)

# Ulrike Grömping (2006). Relative Importance for Linear Regression in R: The
# Package relaimpo. Journal of Statistical Software, 17(1), 1--27.
# Your version of package relaimpo: R package version 2.2-3

# for standardisation:
datMaxALL_standard = data.frame(scale(datMaxALL[c("Genome_copies_gFaeces", "OPG", "weightloss")]))

modFull_std = lm(weightloss ~ OPG * Genome_copies_gFaeces, data = datMaxALL_standard)
modminusOPG_std = lm(weightloss ~ Genome_copies_gFaeces, data = datMaxALL_standard)
modminusDNA_std = lm(weightloss ~ OPG, data = datMaxALL_standard)
modminusInter_std = lm(weightloss ~ OPG + Genome_copies_gFaeces, data = datMaxALL_standard)
list(signifOG = lrtest(modFull_std, modminusOPG_std),
     signifDNA = lrtest(modFull_std, modminusDNA_std),
     signifInter = lrtest(modFull_std, modminusInter_std))
# all good, does not change

summary(modFull)
summary(modFull_std)
summary(modminusInter_std)

# http://dmcglinn.github.io/quant_methods/lessons/standardized_beta_coefficients.html#:~:text=Standardized%20%CE%B2%20coefficients,a%20standard%20deviation%20of%201.
# weird, should not change if standardised!!


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