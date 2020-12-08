### Code to analyse
## 1) Correlation among qPCR and oocyst flotation quantification
### library(ggsci)

if(!exists("sample.data")){
  source("R/1_Data_preparation.R")
}

if(!exists("sdt")){
    source("R/2_qPCR_data_preparation.R")
}

###Let's start plotting and analysing the data!
### 1) Correlation among Eimeria quantification methods
####Genome copies modeled by OPGs 
sdt%>%
  ggplot(aes(OPG+1, Genome_copies_gFaeces))+
  geom_smooth(method = lm, col= "black")+
  scale_x_log10(name = "log10 (Oocyst per gram faeces + 1) (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies per gram faeces) (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks()

##Model 1: Genome copies/g faeces modeled by OPG
DNAbyOPG <- lm(log10(Genome_copies_gFaeces)~log10(OPG+1),
               data = sdt, na.action = na.exclude)
summary(DNAbyOPG)
##Model 2: Genome copies/g faeces modeled by OPG with DPI interaction
DNAbyOPG_dpi <- lm(log10(Genome_copies_gFaeces)~log10(OPG+1)*dpi,
                   data = sdt, na.action = na.exclude)
summary(DNAbyOPG_dpi)
##Comparison of models
anova(DNAbyOPG, DNAbyOPG_dpi)

####OPGs modeled by Genome copies 
sdt%>%
  ggplot(aes(Genome_copies_gFaeces, OPG+1))+
  geom_smooth(method = lm, col= "black")+
  scale_y_log10(name = "log10 (Oocyst per gram faeces + 1) (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(name = "log10 (Genome copies per gram faeces) (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks()

##Model 3: OPG modeled by Genome copies/g faeces
OPGbyDNA <- lm(log10(OPG+1)~log10(Genome_copies_gFaeces),
               data = sdt, na.action = na.exclude)
summary(OPGbyDNA)
##Model 4: Genome copies/g faeces modeled by OPG with DPI interaction
OPGbyDNA_dpi <- lm(log10(OPG+1)~log10(Genome_copies_gFaeces)*dpi,
                   data = sdt, na.action = na.exclude)
summary(OPGbyDNA_dpi)

##Comparison of models
anova(OPGbyDNA, OPGbyDNA_dpi)

##Linear models by DPI
sdt%>%
    ggplot(aes(Genome_copies_gFaeces, OPG+1, fill=dpi))+
    geom_smooth(method = lm, se=FALSE, aes(Genome_copies_gFaeces, OPG+0.1, color=dpi))+
    scale_y_log10(name = "log10 (Oocyst per gram faeces + 1) (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_x_log10(name = "log10 (Genome copies per gram faeces) (qPCR)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    geom_point(shape=21, size=5) +
    theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks()

##Weak correlation between measurments by DPI :S 

###Course of infection 
##Genome copies/g of faeces
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  wilcox_test(Genome_copies_gFaeces ~ dpi)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Genome_copies_gFaeces_DPI_Comparison.csv")

##Select just comparison against DPI 0
stats.test%>%
  filter(group1%in%c("0"))-> stats.test 

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= Genome_copies_gFaeces))+
  scale_y_log10("log10 Genome copies/g Faeces (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  scale_color_brewer(palette = "Paired")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> A

##Significant mean difference from day 3 and on... Basically DPI 0, 1 and 2 DNA measurments are the same!

## Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
#sdt%>%
#  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
#  dplyr::select(EH_ID, dpi,OPG)%>%
#  dplyr::arrange(EH_ID)%>%
#  dplyr::arrange(dpi)%>% ##for comparison 
#  dplyr::mutate(OPG= OPG+1)%>%
#  rstatix::wilcox_test(OPG ~ dpi, ref.group = "0")%>%
#  adjust_pvalue(method = "bonferroni") %>%
# add_significance()%>%
#  add_xy_position(x = "dpi")#-> stats.test

##Save statistical analysis
#x <- stats.test
#x$groups<- NULL
#write.csv(x, "Tables/OPG_DPI_Comparison.csv")

##Select just comparison against DPI 0
#stats.test%>%
#  filter(group1%in%c("0"))-> stats.test 

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= OPG+1))+
  scale_y_log10("log10 (Oocyst per gram faeces + 1) (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  scale_color_brewer(palette = "Paired")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> B

##Significant mean difference from day 4 and on... Basically DPI 0 to 3 No OPG and DPI 4 equal to 10!

##Weight loss 
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi, weightloss)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  wilcox_test(weightloss ~ dpi)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Weightloss_DPI_Comparison.csv")

##Select just comparison against DPI 0
stats.test%>%
  filter(group1%in%c("0"))-> stats.test 

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi,weightloss)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= weightloss))+
  geom_boxplot()+
  geom_point(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  xlab("Day post infection")+
  ylab("Weight loss (%)")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  scale_color_brewer(palette = "Paired")+
  labs(tag= "C)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> C

##Figure # Course of Eimeria Infection in genome copies, OPG, and weight loss
#pdf(file = "fig/Figure_3.pdf", width = 10, height = 15)
grid.arrange(A,B,C)
#dev.off()
rm(A,B, C)

###################################### Extra code ###########################################
## DNA as a predictor of weightloss
sdt%>%
   ggplot(aes(Genome_copies_gFaeces, weightloss))+
   geom_smooth(method = lm, color= "black")+
   scale_y_continuous(name = "Weight loss to 0 dpi")+
   scale_x_log10(name = "log10 Genome copies/g Faeces (qPCR)", 
                 breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
   geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
   labs(tag= "B)")+
   theme_bw()+
   theme(text = element_text(size=16))#+
   #stat_cor(label.x = 2.0, label.y = 84, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
   #stat_regline_equation(label.x = 2.0, label.y = 88)

##Model 5: Genome copies/g of feaces as predictor of weight loss 
WlbyDNA <- lm(weightloss~log10(Genome_copies_gFaeces),
               data = sdt, na.action = na.exclude)
summary(WlbyDNA)
##Model 6: Genome copies/g of feaces as predictor of weight loss with DPI interaction
WlbyDNA_dpi <- lm(weightloss~log10(Genome_copies_gFaeces)*dpi,
                   data = sdt, na.action = na.exclude)
summary(WlbyDNA_dpi)

##Comparison of models
anova(WlbyDNA, WlbyDNA_dpi)

## OPG as a predictor of weightloss
sdt%>%
  ggplot(aes(OPG, weightloss))+
  geom_smooth(method = lm, color= "black")+
  scale_y_continuous(name = "Weight loss to 0 dpi")+
  scale_x_log10(name = "log10 Oocysts per gram of faeces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

## sdt$RelWeight<- (sdt$weight/sdt$weight_dpi0)*100 ##Check with Alice 
## summary(lm(RelWeight~Genome_copies_mean, sdt))

############## Create Time-Series (Alice and Susi analysis)####################
#require("reshape")
#sdt%>%
#  dplyr::select(EH_ID, dpi, Genome_copies_mean)%>%
#  dplyr::arrange(EH_ID)%>%
#  dplyr::arrange(dpi)%>%
#  na.exclude()-> dna
#dna<- reshape(dna, idvar = "EH_ID", timevar = "dpi", direction = "wide")

#sdt%>%
#  dplyr::select(EH_ID, dpi, OPG)%>%
#  dplyr::arrange(EH_ID)%>%
#  dplyr::arrange(dpi)%>%
#  na.exclude()-> oocysts
#oocysts<- reshape(oocysts, idvar = "EH_ID", timevar = "dpi", direction = "wide")

#ts.data<- inner_join(oocysts, dna, by= "EH_ID")

#ts.data%>%
#  dplyr::rowwise()%>%
#  dplyr::mutate(Sum_Oocysts= sum(c(OPG.0,OPG.3,OPG.4,OPG.5,OPG.6,OPG.7,OPG.8,OPG.9,OPG.10)))-> ts.data
#ts.data<- na.omit(ts.data)

#ts.data%>%
#  dplyr::rowwise()%>%
#  dplyr::mutate(Max_Oocysts= max(c(OPG.0,OPG.3,OPG.4,OPG.5,OPG.6,OPG.7,OPG.8,OPG.9,OPG.10)))-> ts.data

##Total OPGs during infection are predicted by DNA at different dpi? Genome copies per dpi as individual predictors

## sum.opg <- glm.nb(formula = Sum_Oocysts~ Genome_copies_mean.0+
##       Genome_copies_mean.1+
##       Genome_copies_mean.2+
##       Genome_copies_mean.3+
##       Genome_copies_mean.4+
##       Genome_copies_mean.5+
##       Genome_copies_mean.6+
##       Genome_copies_mean.7+
##       Genome_copies_mean.8+
##       Genome_copies_mean.9+
##       Genome_copies_mean.10, data = ts.data, na.action = na.exclude)

## summary(sum.opg)
## plot(sum.opg)

## drop1(sum.opg, test= "LRT")

## library(sjPlot)
## library(sjmisc)
## library(sjlabelled)
## tab_model(sum.opg, show.intercept = T, show.est = T, show.stat = T, show.fstat = T, show.aic = T, show.obs = T, show.loglik = T)

## ##extract p values for bonferroni correction
## #p.sum.opg<- as.data.frame(coef(summary(sum.opg))[,'Pr(>|z|)'])
## #colnames(p.sum.opg)<- "P_unadjusted"
## #p.sum.opg$P_adjusted<-p.adjust(p.sum.opg$`P_unadjusted`, method = "bonferroni")

## ##Max OPG during infection are predicted by DNA at different dpi? Genome copies per dpi as individual predictors 
## max.opg <- glm.nb(formula = Max_Oocysts~ Genome_copies_mean.0+
##                  Genome_copies_mean.1+
##                  Genome_copies_mean.2+
##                  Genome_copies_mean.3+
##                  Genome_copies_mean.4+
##                  Genome_copies_mean.5+
##                  Genome_copies_mean.6+
##                  Genome_copies_mean.7+
##                  Genome_copies_mean.8+
##                  Genome_copies_mean.9+
##                  Genome_copies_mean.10, data = ts.data, na.action = na.exclude)

## summary(max.opg)
## plot(max.opg)
## tab_model(max.opg, show.intercept = T, show.est = T, show.stat = T, show.fstat = T, show.aic = T, show.obs = T, show.loglik = T)

## ##extract p values for bonferroni correction
## #p.max.opg<- as.data.frame(coef(summary(max.opg))[,'Pr(>|z|)'])
## #colnames(p.max.opg)<- "P_unadjusted"
## #p.max.opg$P_adjusted<-p.adjust(p.max.opg$`P_unadjusted`, method = "bonferroni")

## ##Oocysts at pick infection (dpi6) are predicted by DNA at different dpi? Genome copies per dpi as individual predictors 
## dpi6.opg <- glm.nb(formula = OPG.6~ Genome_copies_mean.0+
##                     Genome_copies_mean.1+
##                     Genome_copies_mean.2+
##                     Genome_copies_mean.3+
##                     Genome_copies_mean.4+
##                     Genome_copies_mean.5+
##                     Genome_copies_mean.6+
##                     Genome_copies_mean.7+
##                     Genome_copies_mean.8+
##                     Genome_copies_mean.9+
##                     Genome_copies_mean.10, data = ts.data, na.action = na.exclude)

## summary(dpi6.opg)
## plot(dpi6.opg)
## tab_model(dpi6.opg, show.intercept = T, show.est = T, show.stat = T, show.fstat = T, show.aic = T, show.obs = T, show.loglik = T)

## ##extract p values for bonferroni correction
## #p.dpi6.opg<- as.data.frame(coef(summary(dpi6.opg))[,'Pr(>|z|)'])
## #colnames(p.dpi6.opg)<- "P_unadjusted"
## #p.dpi6.opg$P_adjusted<-p.adjust(p.dpi6.opg$`P_unadjusted`, method = "bonferroni")

## set.seed(2020)
## ts.data%>%
##   dplyr::select(EH_ID,Genome_copies_mean.4, Sum_Oocysts, Max_Oocysts, OPG.6)%>%
##   ggplot(aes(Genome_copies_mean.4, Sum_Oocysts))+
##   geom_smooth(method = lm, col="black")+
##   scale_x_log10(name = "log10 Genome copies/µL gDNA dpi 4 (qPCR)", 
##                 breaks = scales::trans_breaks("log10", function(x) 10^x),
##                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##   scale_y_log10(name = "log10 Sum Oocyst per gram feces (Flotation)", 
##                 breaks = scales::trans_breaks("log10", function(x) 10^x),
##                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##   geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
##   labs(tag= "A)")+
##   theme_bw()+
##   theme(text = element_text(size=16))+
##   stat_cor(label.x = 4.25, label.y = 5.0, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
##   stat_regline_equation(label.x = 4.25, label.y = 5.25)+
##   stat_cor(label.x = 4.25,  label.y = 4.75,method = "spearman")+
##   annotation_logticks()+
##   coord_cartesian(ylim = c(10000, 10000000))-> tssum

## set.seed(2020)
## ts.data%>%
##   dplyr::select(EH_ID,Genome_copies_mean.4, Sum_Oocysts, Max_Oocysts, OPG.6)%>%
##   ggplot(aes(Genome_copies_mean.4, Max_Oocysts))+
##   geom_smooth(method = lm, col="black")+
##   scale_x_log10(name = "log10 Genome copies/µL gDNA dpi 4 (qPCR)", 
##                 breaks = scales::trans_breaks("log10", function(x) 10^x),
##                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##   scale_y_log10(name = "log10 Max Oocyst per gram feces (Flotation)", 
##                 breaks = scales::trans_breaks("log10", function(x) 10^x),
##                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##   geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
##   labs(tag= "B)")+
##   theme_bw()+
##   theme(text = element_text(size=16))+
##   stat_cor(label.x = 4.25, label.y = 4.75, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
##   stat_regline_equation(label.x = 4.25, label.y = 5.0)+
##   stat_cor(label.x = 4.25,  label.y = 4.5,method = "spearman")+
##   annotation_logticks()+
##   coord_cartesian(ylim = c(10000, 10000000))-> tsmax

## set.seed(2020)
## ts.data%>%
##   dplyr::select(EH_ID,Genome_copies_mean.4, Sum_Oocysts, Max_Oocysts, OPG.6)%>%
##   ggplot(aes(Genome_copies_mean.4, OPG.6))+
##   geom_smooth(method = lm, col="black")+
##   scale_x_log10(name = "log10 Genome copies/µL gDNA dpi 4 (qPCR)", 
##                 breaks = scales::trans_breaks("log10", function(x) 10^x),
##                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##   scale_y_log10(name = "log10 Oocyst per gram feces dpi 6 (Flotation)", 
##                 breaks = scales::trans_breaks("log10", function(x) 10^x),
##                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##   geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
##   labs(tag= "C)")+
##   theme_bw()+
##   theme(text = element_text(size=16))+
##   stat_cor(label.x = 4.0, label.y = 4.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
##   stat_regline_equation(label.x = 4.0, label.y = 4.75)+
##   stat_cor(label.x = 4.0,  label.y = 4.25,method = "spearman")+
##   annotation_logticks()+
##   coord_cartesian(ylim = c(10000, 10000000))-> ts6

## ts.data%>%
##  ggplot(aes(Sum_Oocysts, Max_Oocysts))+
##  geom_smooth(method = lm, col= "black")+
##  scale_x_log10(name = "log10 Sum Oocyst per gram feces", 
##                breaks = scales::trans_breaks("log10", function(x) 10^x),
##                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##  scale_y_log10(name = "log10 Max Oocyst per gram feces", 
##                breaks = scales::trans_breaks("log10", function(x) 10^x),
##                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
##  labs(tag= "C)")+
##  theme_bw()+
##  theme(text = element_text(size=16))+
##  stat_cor(label.x = 5.5, label.y = 1.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
##  stat_regline_equation(label.x = 5.5, label.y = 2)+
##  stat_cor(label.x = 5.5,  label.y = 1,method = "spearman")+
##  annotation_logticks()-> opgmaxsum

## ##Save plots
##Define later which plots from here will be included 