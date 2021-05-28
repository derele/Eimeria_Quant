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
  filter(!is.na(Genome_copies_gFaeces))%>%
  wilcox_test(Genome_copies_gFaeces ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
#write.csv(x, "Tables/Genome_copies_gFaeces_DPI_Comparison.csv")
stats.test%>%
  dplyr::mutate(y.position = log10(y.position))%>%
  dplyr::mutate(dpi = c("1","2","3","4", "5","6", "7", "8", "9", "10"))-> stats.test

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces, Infection)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  filter(!is.na(Genome_copies_gFaeces))%>% 
  #filter(Genome_copies_gFaeces!=0)%>% 
  ggplot(aes(x= dpi, y= Genome_copies_gFaeces+1))+
  scale_y_log10("log10 (Genome copies/g Faeces + 1) (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot()+
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Infection, fill= dpi), color= "black")+
  scale_shape_manual(values = c(21, 24))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray", alpha= 0.5)+
  labs(tag= "A)", shape= "qPCR status (Melting curve)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top")+
  annotation_logticks(sides = "l")+
  guides(fill=FALSE)+
  stat_pvalue_manual(stats.test, label= "p.adj.signif", x= "dpi", y.position = 100000000000)-> A

##Significant mean difference from DPI 2!!!!! not all samples but some sowed signal!!!!

## Oocysts
##Wilcoxon test (Compare mean per DPI with DPI 0 as reference)
# sdt%>%
#   filter(dpi%in%c("0","1", "2", "3","4", "5","6", "7", "8", "9", "10"))%>%
#   dplyr::select(EH_ID, dpi,OPG)%>%
#   dplyr::arrange(EH_ID)%>%
#   dplyr::arrange(dpi)%>% ##for comparison 
#   dplyr::mutate(OPG= OPG+1)%>% ##To check
#   filter(!is.na(OPG))->#%>% 
  #wilcox_test(OPG ~ dpi, alternative = "two.sided", ref.group = "0")%>%
#  adjust_pvalue(method = "bonferroni") %>%
# add_significance()%>%
#  add_xy_position(x = "dpi")#-> stats.test

##Save statistical analysis
#x <- stats.test
#x$groups<- NULL
#write.csv(x, "Tables/OPG_DPI_Comparison.csv")

sdt%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_gFaeces)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= OPG+1))+
  scale_y_log10("log10 (Oocyst/g Faeces + 1) (Flotation)", 
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
  wilcox_test(weightloss ~ dpi, alternative = "two.sided", ref.group = "0")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "dpi")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
#write.csv(x, "Tables/Weightloss_DPI_Comparison.csv")

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
rm(A,B,C, x, stats.test)

### 2) Correlation among Eimeria quantification methods
##Model 1: Genome copies/g faeces modeled by OPG
DNAbyOPG <- lm(log10(Genome_copies_gFaeces+1)~log10(OPG+1),
               data = sdt, na.action = na.exclude)
summary(DNAbyOPG)

sdt$predictedM1 <- predict(DNAbyOPG)   # Save the predicted values
sdt$residualsM1 <- residuals(DNAbyOPG) # Save the residual values

##Plot model
####Genome copies modeled by OPGs 
sdt%>%
  ggplot(aes(OPG+1, Genome_copies_gFaeces+1))+
  geom_smooth(method = lm, col= "black")+
  scale_x_log10(name = "log10 (Oocyst/g Faeces + 1) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces +1)  \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  annotation_logticks()-> A

##Plot residuals
##Mean residuals for plot 
sdt%>%
  group_by(dpi) %>% 
  summarise(residualsM1_mean = mean(na.omit(residualsM1)))%>%
  inner_join(sdt, by= "dpi")%>%
  filter(dpi%in%c("0","1","2","3","4", "5","6", "7", "8", "9", "10"))%>%
  dplyr::select(dpi, residualsM1, residualsM1_mean)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  mutate(residualim= 0)%>%
  ggplot(aes(x= dpi, y= residualsM1))+
  geom_jitter(width = 0.5, shape=21, size=2.5, aes(fill= dpi), alpha= 0.75, color= "black")+
  geom_segment(aes(y= residualsM1_mean, yend= residualim, xend= dpi), color= "black", size= 1) +
  geom_point(aes(x = dpi, y = residualsM1_mean), size=4)+
  geom_rect(aes(xmin=-0.1,xmax=4.5,ymin=-Inf,ymax=Inf),alpha=0.01,fill="grey")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  xlab("Day post infection")+
  scale_y_continuous(name = "Residuals\n (Genome copies/g Faeces)")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")-> B

##Model 2: Genome copies/g faeces modeled by OPG without DPI interaction
DNAbyOPG_dpi <- lm(log10(Genome_copies_gFaeces+1)~log10(OPG+1)+dpi,
                   data = sdt, na.action = na.exclude)
summary(DNAbyOPG_dpi)

##Model 3: Genome copies/g faeces modeled by DPI 
DNAbydpi <- lm(log10(Genome_copies_gFaeces+1)~dpi,
                   data = sdt, na.action = na.exclude)
summary(DNAbydpi)

##Model 4: Genome copies/g faeces modeled by OPG with DPI interaction
DNAbyOPGxdpi <- lm(log10(Genome_copies_gFaeces+1)~log10(OPG+1)*dpi,
                   data = sdt, na.action = na.exclude)
summary(DNAbyOPGxdpi)

##Comparison of models
# test difference LRT or anova
DNANULL<- lm(log10(Genome_copies_gFaeces+1)~1,
             data = sdt, na.action = na.exclude)

lrtest(DNANULL, DNAbyOPG) 
lrtest(DNANULL, DNAbyOPG_dpi) 
lrtest(DNANULL, DNAbyOPGxdpi) 
lrtest(DNAbyOPG, DNAbyOPG_dpi) #--> Report this table in the results 
lrtest(DNAbyOPG, DNAbyOPGxdpi) #--> Report this table in the results 
lrtest(DNAbyOPG_dpi, DNAbyOPGxdpi) #--> Report this table in the results 

### GLMM
require(lme4)
require(sjPlot)
DNAbyOPG_dpi_glmm <- lmer(log10(Genome_copies_gFaeces+1)~log10(OPG+1) + dpi + (1|EH_ID),
                          data = sdt, na.action = na.exclude, REML=TRUE)

summary(DNAbyOPG_dpi_glmm)

plot_models(DNAbyOPG_dpi, DNAbyOPG_dpi_glmm)
dev.off()

##Plot model by DPI
colores<- c("4"="#00BD5C", "5"= "#00C1A7", "6"= "#00BADE", "7"= "#00A6FF", 
         "8" = "#B385FF", "9"= "#EF67EB", "10" = "#FF63B6")

sdt%>%
  mutate(dpi = fct_relevel(dpi, "0","1", "2", "3", "4", "5", 
                                   "6", "7", "8", "9", "10"))%>%
  ggplot(aes(OPG+1, Genome_copies_gFaeces+1, fill=dpi))+
  geom_point(shape=21, size=5) +
  geom_smooth(method = lm, se=FALSE, aes(OPG, Genome_copies_gFaeces, color=dpi))+
  scale_color_manual(values = colores, guide= "none")+
  scale_x_log10(name = "log10 (Oocyst/g Faeces + 1) \n (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 (Genome copies/g Faeces) \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  labs(tag= "C)")+
  theme(text = element_text(size=16), legend.position = "none")+
  annotation_logticks()-> C

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

#write.csv(corOPGbyDNA_DPI, "Tables/Q1_OPG_DNA_Correlation_DPI.csv",  row.names = F)
##Non significant correlation between measurements by DPI
###################################### Extra code ###########################################
## DNA as a predictor of weightloss
#sdt%>%
#   ggplot(aes(Genome_copies_gFaeces, weightloss))+
#   geom_smooth(method = lm, color= "black")+
#   scale_y_continuous(name = "Weight loss to 0 dpi")+
#   scale_x_log10(name = "log10 Genome copies/g Faeces (qPCR)", 
#                 breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#   geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
#   labs(tag= "B)")+
#   theme_bw()+
#   theme(text = element_text(size=16))

##Model 5: Genome copies/g of feaces as predictor of weight loss 
#WlbyDNA <- lm(weightloss~log10(Genome_copies_gFaeces),
#               data = sdt, na.action = na.exclude)
#summary(WlbyDNA)
##Model 6: Genome copies/g of feaces as predictor of weight loss with DPI interaction
#WlbyDNA_dpi <- lm(weightloss~log10(Genome_copies_gFaeces)*dpi,
#                   data = sdt, na.action = na.exclude)
#summary(WlbyDNA_dpi)

##Comparison of models
#anova(WlbyDNA, WlbyDNA_dpi)

## OPG as a predictor of weightloss
#sdt%>%
#  ggplot(aes(OPG, weightloss))+
#  geom_smooth(method = lm, color= "black")+
#  scale_y_continuous(name = "Weight loss to 0 dpi")+
#  scale_x_log10(name = "log10 Oocysts per gram of faeces (Flotation)", 
#                breaks = scales::trans_breaks("log10", function(x) 10^x),
#                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
#  labs(tag= "B)")+
#  theme_bw()+
#  theme(text = element_text(size=16))