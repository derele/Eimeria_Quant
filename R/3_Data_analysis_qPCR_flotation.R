### Code to analyses
## 1) Correlation among qPCR and oocyst flotation quantification
### library(ggsci)

if(!exists("sample.data")){
  source("R/1_Data_preparation.R")
}

if(!exists("sdt")){
    source("R/2_qPCR_data_preparation.R")
}

###Let's start plotting and analysing the data!
### 1) Course of infection 
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
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "none")+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> A

##Significant mean difference from day 3 and on... Basically DPI 0, 1 and 2 DNA measurments are the same!

## Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
sdt%>%
  filter(dpi%in%c("0","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi,OPG)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  dplyr::mutate(OPG= OPG+1)#%>% ##To check
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
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "none")+
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
  theme(text = element_text(size=16), legend.position = "none")+
  stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "0", paired = F, na.rm = TRUE)-> C

##Figure 2: Course of Eimeria Infection in genome copies, OPG, and weight loss
#pdf(file = "fig/Figure_2.pdf", width = 10, height = 15)
grid.arrange(A,B,C)
#dev.off()
rm(A,B, C)

### 2) Correlation among Eimeria quantification methods
####Genome copies modeled by OPGs 
sdt%>%
  ggplot(aes(OPG, Genome_copies_gFaeces))+
  geom_smooth(method = lm, col= "black")+
  scale_x_log10(name = "log10 (Oocyst per gram faeces) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies per gram faeces)  \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  stat_cor(label.x = 4.5,  label.y = 2,method = "spearman",
           aes(label= paste(..r.., ..p.label.., sep= "~`,`~")))+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks()+
  geom_text (x = 3.75, y = 2.05, show.legend = F,
             label = paste ("Spearman's rho ="))-> A

##Model 1: Genome copies/g faeces modeled by OPG
DNAbyOPG <- lm(log10(Genome_copies_gFaeces)~log10(OPG+1),
               data = sdt, na.action = na.exclude)
summary(DNAbyOPG)
##Model 2: Genome copies/g faeces modeled by OPG with DPI interaction
DNAbyOPG_dpi <- lm(log10(Genome_copies_gFaeces)~log10(OPG+1)*dpi,
                   data = sdt, na.action = na.exclude)
summary(DNAbyOPG_dpi)
##Comparison of models
anova(DNAbyOPG, DNAbyOPG_dpi,  test="LRT")

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
             label = paste ("Spearman's rho ="))

corOPGbyDNA <-cor.test(log10(sdt$Genome_copies_gFaeces+1), log10(sdt$OPG+1),  method = "spearman")

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
colores<- c("4"="#00BD5C", "5"= "#00C1A7", "6"= "#00BADE", "7"= "#00A6FF", 
         "8" = "#B385FF", "9"= "#EF67EB", "10" = "#FF63B6")

sdt%>%
  mutate(dpi = fct_relevel(dpi, "0","1", "2", "3", "4", "5", 
                                   "6", "7", "8", "9", "10"))%>%
  ggplot(aes(OPG, Genome_copies_gFaeces, fill=dpi))+
  geom_point(shape=21, size=5) +
  geom_smooth(method = lm, se=FALSE, aes(OPG, Genome_copies_gFaeces, color=dpi))+
  scale_color_manual(values = colores, guide= "none")+
  scale_x_log10(name = "log10 (Oocyst per gram faeces) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies per gram faeces) \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  labs(tag= "B)")+
  theme(text = element_text(size=16))+
  annotation_logticks() -> B

AB<- ggarrange(A, B, common.legend = TRUE, ncol = 1, nrow = 2)


sdt%>% 
  nest(-dpi)%>% 
  mutate(cor=map(data,~cor.test(log10(.x$Genome_copies_gFaeces+1), log10(.x$OPG+1), method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()-> x

x$data<- NULL
x$cor<- NULL
x%>%
  arrange(dpi)-> corOPGbyDNA_DPI

write.csv(corOPGbyDNA_DPI, "Tables/Q1_OPG_DNA_Correlation_DPI.csv",  row.names = F)
##Non significant correlation between measurements by DPI

##Figure 4# Spearman's Correlation between genome copies and OPG overall and by dpi
#pdf(file = "fig/Figure_4.pdf", width = 10, height = 15)
grid.arrange(A,B)
#dev.off()
rm(A,B)

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
   theme(text = element_text(size=16))

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