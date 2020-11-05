## Code for 
## 1) Standard curve of qPCR Eimeria

### We should have ONE model able to predict for different species and
### cyclers, IF they are significantly different

## 2) Determine Eimeria amount for infection experiment samples

### then make the predictions based on the data (including cycler and
### species) and the models (with those factors if significant)

## only load packages that are needed!!!
library(ggpubr)
library(rcompanion)
library(dplyr)
library(gridExtra)

##Load data
if(!exists("sample.data")){
    source("R/1_Data_preparation.R")
}
##Standard curves
data.std<- read.csv("data/Eimeria_quantification_Std_Curve_data.csv")
data.std%>%
    dplyr::mutate(Genome_copies= Oocyst_count*8)-> data.std

##Define numeric and factor variables 
num.vars <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", "Genome_copies")
fac.vars <- c("Well", "Sample.Name", "Detector", "Task",  "Std_series","Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction")  

## as.numeric alone will likely fail if stringsAsfactors is TRUE! 
data.std[, num.vars] <- apply(data.std[, num.vars], 2,
                                 function (x) as.numeric(as.character(x)))
data.std[, fac.vars] <- apply(data.std[, fac.vars], 2, as.factor)

##Correct zero in NTC with not-detected results 
data.std$Ct[data.std$Ct == 0] <- NA
data.std$Ct_mean[data.std$Ct_mean == 0 & data.std$Task == "NTC"] <- NA
data.std$Sd_Ct[data.std$Sd_Ct == 0] <- NA

##Correct labels to have them homogeneous 
data.std$Sample.Name<- gsub(pattern = " ", replacement = "_", x = data.std$Sample.Name)

## Select just standards data
## Estimate the number of genome copies per ng of gDNA
data.std.lm<- subset(data.std, Task== "Standard") ## Select just data from standards 
data.std.lm %>% 
  select(Sample.Name, Task, Ct, Cycler, Oocyst_count, Parasite, Genome_copies)%>%
  dplyr::mutate(Oocyst_DNA= Oocyst_count*(3.8E-4))%>% ##Estimation of DNA (ng) derived from Oocyst
  dplyr::mutate(DNA_PCR= Oocyst_DNA/30)%>% ##DNA (ng) in PCR considering 1uL from a stock of 30uL
  ##Considering that 1 ng of Eimeria gDNA is equivalent to 2.11E4 genome copies
  dplyr::mutate(Genome_copies_ngDNA= (2.11E4)*DNA_PCR)-> data.std.lm  

##Inter-sample variation and spiked samples
data.unk<-read.csv("data/Eimeria_quantification_Sample_data.csv")
  
##Define numeric and factor variables 
num.vars2 <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", 
                "Oocyst_1", "Oocyst_2", "Oocyst_3", "Oocyst_4", "Oocyst_5", "Oocyst_6", "Oocyst_7", "Oocyst_8", "Dilution_factor", "Volume", "Sporulated")
  fac.vars2 <- c("Well", "Sample.Name", "Detector", "Task", "Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction", "Strain")  

  data.unk[, num.vars2] <- apply(data.unk[, num.vars2], 2,
                               function (x) as.numeric(as.character(x)))
  data.unk[, fac.vars2] <- apply(data.unk[, fac.vars2], 2, as.factor)

##Information from intersample variation experiment 
data.unk%>%
  dplyr::select(Sample.Name, Task, Ct, Cycler, Parasite, Sample_type, Extraction, Tm, 
                Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8, 
                Sporulated, Dilution_factor, Volume, Strain)%>%
  filter(Sample_type=="Oocysts" & Task=="Unknown")%>%
  dplyr::group_by(Sample.Name)%>%
  dplyr::mutate(N= n())%>%
  dplyr::mutate(Total_oocysts= (sum(Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6,
                                    Oocyst_7, Oocyst_8))/N)%>%
  ##Concentration of Oocyst in the solution
  dplyr::mutate(Oocyst_count= (((Total_oocysts*10000)/8)*Dilution_factor))%>%
  ##Concentration of sporulated oocyst in the solution
  dplyr::mutate(Sporulated_count= (((Sporulated*10000)/8)*Dilution_factor))%>%
  dplyr::mutate(Sporulation_rate= (Sporulated_count/Oocyst_count)*100)%>%
  dplyr::mutate(Sporulation_rate= as.numeric(Sporulation_rate))-> data.unk.lm

##Spiked data 
data.unk%>%
  select(Sample.Name, Task, Ct, Cycler,Parasite, Sample_type, Feces_weight, Extraction, Oocyst_count, Qubit, NanoDrop)%>%
  filter(Sample_type=="Feces" & Task=="Unknown")-> data.spk
  
##Infection experiment
data.inf.exp<-read.csv("data/Eimeria_quantification_Inf_exp_data.csv")
data.inf.exp%>%
    select(Content, Sample, Plate_number, Cq, Melt_Temperature)%>%
    dplyr::rename(Ct= Cq, labels= Sample, Task= Content, Tm= Melt_Temperature)%>%
    dplyr::mutate(Cycler= "BioRad")-> data.inf.exp
  
##Define numeric and factor variables 
num.vars3 <- c("Ct", "Tm")
fac.vars3 <- c("labels", "Task", "Plate_number", "Cycler")  
data.inf.exp[, num.vars3] <- apply(data.inf.exp[, num.vars3], 2,
                               function (x) as.numeric(as.character(x)))
data.inf.exp[, fac.vars3] <- apply(data.inf.exp[, fac.vars3], 2, as.factor)
  

rm(fac.vars, num.vars, fac.vars2, num.vars2, fac.vars3, num.vars3)

####### Standard curves #######
set.seed(2020)
## Comput simple linear models from standards
## "Genome copies modeled by Ct"

##Ct modeled by Oocyst counts; data from different Cyclers
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Ct, color= Cycler)) +
  geom_smooth(method = "lm", se = T) +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), 
              aes(size= 20,fill= Cycler, shape= Parasite), color= "black", alpha= 0.5)+
  #stat_cor(label.x = 5, label.y = c(35,30,25), 
  #         aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+ # Add correlation coefficient
  #stat_regline_equation(label.x = 5, label.y = c(36.5,31.5,26.5))+ # Add Regression equation lm Ct~log10(Oocyst_count)
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "b")

##Ct modeled by Oocyst_counts and extra predictors to be considered 
##Model 1: Ct modeled by oocyst count simple without other predictor
##considering all data 
lm.Ct<- lm(Ct~log10(Oocyst_count), data.std.lm)

##Model 2: Ct modeled by oocyst count and parasite as predictors
lm.CtPar<- lm(Ct~log10(Oocyst_count)+Parasite, data.std.lm)

##Model 3: Ct modeled by oocyst count and cycler as predictors
lm.CtCyc<- lm(Ct~log10(Oocyst_count)+Cycler, data.std.lm)

##Model 4: Ct modeled by oocysts counts, parasite and cycler used as predictors
lm.CtAll<- lm(Ct~log10(Oocyst_count)+Parasite+Cycler, data.std.lm)

##Model 5: Ct modeled by oocysts counts and parasite/cycle interaction (Check with Alice and Susi)
lm.CtInt<- lm(Ct~log10(Oocyst_count)+Parasite*Cycler, data.std.lm)

##Comparison of models 
compareLM(lm.CtAll, lm.CtPar, lm.CtCyc, lm.Ct, lm.CtInt)

##Model 3 fit better the data... Cycler has major impact (confirm somehow our expectations)!

##Real standard curve##
##Genome copies modeled by Ct and extra predictors to be considered 
##Model 6: Genome copies modeled by Ct simple without other predictor
##considering all data 
lm.SC<- lm(log10(Genome_copies)~Ct, data.std.lm)

##Model 7: Genome copies modeled by Ct and parasite as predictors
lm.SCPar<- lm(log10(Genome_copies)~Ct+Parasite, data.std.lm)

##Model 8: Genome copies modeled by Ct and cycler as predictors
lm.SCCyc<- lm(log10(Genome_copies)~Ct+Cycler, data.std.lm)

##Model 9: Genome copies modeled by Ct, parasite and cycler used as predictors
lm.SCAll<- lm(log10(Genome_copies)~Ct+Parasite+Cycler, data.std.lm)

##Model 10: Genome copies modeled by Ct and parasite/cycle interaction (Check with Alice and Susi)
lm.SCInt<- lm(log10(Genome_copies)~Ct+Parasite*Cycler, data.std.lm)

##Comparison of models 
compareLM(lm.SC, lm.SCPar, lm.SCCyc, lm.SCAll, lm.SCInt)

##Model 8 fit better the data... Cycler has major impact (again confirm expectations)!
##Linear model (Standard curve for the rest of experiments)
data.std.lm%>%
  ggplot(aes(x = Ct, y = Genome_copies)) +
  geom_smooth(method = "lm", se = T, color="black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_y_log10("log 10 Eimeria genome copies", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "l")-> A

##Linear model Genome copies per ng modeled by Oocyst count 
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Genome_copies_ngDNA)) +
  geom_smooth(method = "lm", se = F, color= "black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_x_log10("log 10 Eimeria Oocysts Count (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10("log 10 Eimeria genome copies/ng gDNA \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "B)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "bl")-> B

##Model 11: Genome copies modeled by Oocyst count, parasite and cycle 
lm.SCOoc<- lm(log10(Genome_copies_ngDNA)~log10(Oocyst_count)+Parasite+Cycler, data.std.lm)

### Using MODEL 8 to predict using different levels of the factor cycler
data.std.lm$predicted<- 10^predict(lm.SCCyc)
data.std.lm$residuals<- 10^residuals(lm.SCCyc)

##Predicted genome copies
##modeled by Ct
ggplot(data.std.lm, aes(x = Ct, y = predicted)) +
  geom_smooth(method = "lm", se = T, color= "black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_y_log10("log 10 Predicted Eimeria genome copies", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "l")

##by Oocyst count
ggplot(data.std.lm, aes(x = Oocyst_count, y = predicted)) +
  geom_smooth(method = "lm", se = T, color= "black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10("log 10 Predicted Eimeria genome copies", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "bl")

## ### Figure 1 Final Standard curves 
## pdf(file = "fig/Figure_1.pdf", width = 8, height = 10)
## grid.arrange(A, B)
## dev.off()

## If it is necessary some of the previous figures could be included as supplementary  

##Mean comparison standars against NTC
set.seed(2020)
data.std%>%
  select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Date)%>%
  filter(Task%in%c("Standard", "NTC"))%>%
  ggplot(aes(x = Sample.Name, y = Ct)) +
  scale_x_discrete(name = "Standard",
                   labels= c("Eimeria_10_0"= "Oocysts 10⁰", "Eimeria_10_1"= "Oocysts 10¹",
                             "Eimeria_10_2"= "Oocysts 10²", "Eimeria_10_3"= "Oocysts 10³",
                             "Eimeria_10_4"= "Oocysts 10⁴", "Eimeria_10_5"= "Oocysts 10⁵",
                             "Eimeria_10_6"= "Oocysts 10⁶", "H2O"= "NTC")) +
  scale_y_continuous(name = "Ct")+ 
  geom_jitter(shape=21, position=position_jitter(0.2), color= "black", alpha= 0.5,
              aes(size= 25, fill= Cycler))+
  theme_bw() +
  theme(text = element_text(size=16),legend.position = "none")+
  theme(axis.text.x = element_text(angle=90))+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange",
               shape=16, size=0.5, color="black")+
  labs(tag = "A)")+
  geom_hline(yintercept = 30, linetype = 2)+
  stat_compare_means(method = "anova",
                     aes(label = paste0(..method.., ",\n","p=",..p.format..)),
                     label.y= 33, label.x = 7)+
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "H2O", 
                     label.y = c(40, 39, 33, 30.5, 24, 20, 16, 0))
###Determine that 10^0 and 10^1 meassurments are basically like NTC when all the information is taken into account

###### Intersample variation experiment #####
##Predict genome copies using model 8
data.unk.lm$Genome_copies<- 10^predict(lm.SCCyc, data.unk.lm)
##Adjust genome copies per ng of DNA
data.unk.lm%>%
dplyr::mutate(Oocyst_DNA= Oocyst_count*(3.8E-4))%>% ##Estimation of DNA (ng) derived from Oocyst
  dplyr::mutate(DNA_PCR= Oocyst_DNA/30)%>% ##DNA (ng) in PCR considering 1uL from a stock of 30uL
  dplyr::mutate(Genome_copies_ngDNA= Genome_copies*DNA_PCR)-> data.unk.lm  

data.unk.lm%>%
  bind_rows(data.std.lm)-> data.unk.lm

data.unk.lm%>%
  filter(Task%in%c("Unknown"))%>%
  ggplot(aes(Oocyst_count, Genome_copies_ngDNA), geom=c("point", "smooth"))+
  scale_x_log10(name = "log10 Eimeria Oocysts count (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log 10 Eimeria genome copies/ng gDNA \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  geom_point(shape=21, position=position_jitter(0.1), aes(size= 25, fill= Task), color= "black", alpha= 0.5)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(color= "black", method = "lm")+            
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)")+
  annotation_logticks(sides = "bl")

set.seed(2020)
data.unk.lm%>%
  dplyr::select(Sample.Name,Task, Oocyst_count, Genome_copies, Genome_copies_ngDNA)%>%  
  filter(Task%in%c("Standard", "Unknown"))%>%
  ggplot(aes(x = Oocyst_count, y = Genome_copies_ngDNA), geom=c("point", "smooth")) +
  scale_x_log10(name = "log10 Eimeria Oocyst count (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Eimeria genome copies/ng gDNA \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  geom_jitter(shape=21, position=position_jitter(0.2),
              aes(fill= Task), size= 5, color= "black", alpha= 0.5)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(aes(color= Task, fill= Task), method = "lm")+            
  #stat_cor(label.x = 0.75, label.y = c(5.5, 6.5),
  #         aes(label = paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+
  ## Add correlation coefficient
  #stat_regline_equation(aes(color = Task), label.x = 0.75, label.y = c(5,6))+
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "B)")+
  annotation_logticks(sides = "bl")-> C

##Model 12: Intersample variation considering Parasite, strain, cycler and sporulation rate as predictors
lm.ISV<- lm(formula = log10(Genome_copies_ngDNA)~log10(Oocyst_count)+Parasite+Strain+Cycler+Sporulation_rate, data = data.unk.lm)

##Compair model 11 (perfect fit) vs model 12 
##compareLM(lm.SCOoc, lm.ISV)
##Predict genome copies from oocyst count and compair against previous prediction (Not finished, check!)
##data.unk.lm$Genome_copies_ng_pred<- 10^predict(lm.ISV, data.unk.lm)

##Figure 2 Intersample variation
##pdf(file = "fig/Figure_2.pdf", width = 10, height = 8)
##grid.arrange(C)
##dev.off()

########## Spiked samples Experiment #########

##Compair Spiked samples qPCR estimation with real oocyst count by two extraction methods
set.seed(2020)
##Predict genome copies using model 8
data.spk$Genome_copies<- 10^predict(lm.SCCyc, data.spk)
data.spk%>%
  dplyr::mutate(Genome_copies_gFaeces= Genome_copies/Feces_weight)-> data.spk

##Difference between 2 hatching strategies 
data.spk%>%
    ggplot(aes(x = Oocyst_count, y = Genome_copies_gFaeces, color=Extraction), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocysts (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies (qPCR)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2), color= "black", aes(size= 25, fill= Extraction), alpha= 0.5)+
    theme_bw() +
    geom_smooth(aes(color= Extraction, fill= Extraction), method = "lm")+
    facet_grid(cols = vars(Extraction))+
    stat_cor(label.x = log10(100), label.y = log10(50000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Extraction))+        # Add correlation coefficient
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    stat_regline_equation(aes(color = Extraction), label.x = log10(100), label.y = log10(75000))+
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "A)")+
    annotation_logticks(sides = "bl")

##All data together
data.spk%>%
  ggplot(aes(x = Oocyst_count, y = Genome_copies_gFaeces), geom=c("point", "smooth")) +
  scale_x_log10(name = "log10 Eimeria Oocysts (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Eimeria genome copies (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  geom_jitter(shape=21, position=position_jitter(0.2), color= "black", aes(size= 25, fill= Extraction), alpha= 0.5)+
  theme_bw() +
  geom_smooth(color= "black", method = "lm")+
  stat_cor(label.x = log10(100), label.y = log10(50000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+ # Add correlation coefficient
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  stat_regline_equation(label.x = log10(100), label.y = log10(70000))+
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)")+
  annotation_logticks(sides = "bl")

##Comparison between standard curve and spiked samples from ceramic bead extracted data
data.spk%>%
    filter(Sample_type=="Feces" & Task=="Unknown" & Extraction!="Glass_beads")%>%
    dplyr::mutate(Genome_copies_ngDNA= Genome_copies/50, ## copies by ng of fecal DNA considering 1uL from 50 ng/uL DNA
                  DNA_sample= Qubit*40, ## Estimate total gDNA of sample
                  DNA_g_feces= DNA_sample/Feces_weight,
                  ## Transform it to ng fecal DNA by g of faeces
                  Genome_copies_gFaeces= Genome_copies_ngDNA*DNA_g_feces, ## Estimate genome copies by g of faeces
                  OPG=Oocyst_count/Feces_weight) -> data.spk.lm ## Estimate oocyst per g of faeces for spiked samples

data.spk.lm%>%
    bind_rows(data.std.lm)-> data.spk.lm

##
set.seed(2020)
data.spk.lm%>%
    dplyr::select(Sample.Name,Task, Oocyst_count, Genome_copies, Genome_copies_ngDNA)%>%  
    filter(Task%in%c("Standard", "Unknown"))%>%
    ggplot(aes(x = Oocyst_count, y = Genome_copies_ngDNA), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocyst count (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies/ng gDNA \n (qPCR)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2),
                aes(fill= Task), size= 5, color= "black", alpha= 0.5)+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(aes(color= Task, fill= Task), method = "lm")+            
    #stat_cor(label.x = 0.75, label.y = c(5.5, 6.5),
    #         aes(label = paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+
    ## Add correlation coefficient
    #stat_regline_equation(aes(color = Task), label.x = 0.75, label.y = c(5,6))+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "A)")+
    annotation_logticks(sides = "bl")-> D

##Figure 3 comparison between Eimeria genome copies from oocyst DNA and from fecal DNA 
#pdf(file = "fig/Figure_3.pdf", width = 10, height = 8)
#grid.arrange(D)
#dev.off()

##Model 13: Genome copies/ng gDNA modeled by Oocyst count, cycler and parasite as predictors
#lm.spk<- lm(formula = log10(Genome_copies_ngDNA)~ log10(Oocyst_count), data = subset(data.spk.lm, Task== "Unknown"))

######### Infection experiment data############
## Define real positive and negatives based on Tm 
data.inf.exp %>% 
    dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                        Tm >= 80   ~ "Negative",
                                        Tm < 80 ~ "Positive")) -> data.inf.exp 

##Estimate number of genome copies with qPCR Ct value (Model 8)
data.inf.exp$Genome_copies<- 10^predict(lm.SCCyc, data.inf.exp)

##Summarise genome copies by sample  
data.inf.exp %>%
    select(Genome_copies,labels) %>% # select variables to summarise
    na.omit()%>%
    dplyr::group_by(labels)%>%
    dplyr::summarise_each(funs(Genome_copies_min = min, Genome_copies_q25 = quantile(., 0.25),
                               Genome_copies_median = median, Genome_copies_q75 = quantile(., 0.75), 
                               Genome_copies_max = max, Genome_copies_mean = mean, Genome_copies_sd = sd)) -> Sum.inf
##qPCR data from 236 samples
## Tm values were not sumarised to avoid problems in the function 

##Join summirised data
data.inf.exp<- inner_join(data.inf.exp, Sum.inf, by= "labels")

##Eliminate an unprocessed sample and controls
data.inf.exp%>%
    select(labels, Genome_copies_mean, Infection)%>%
    filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl","FML"))%>% ## Replace NAs in real negative samples to 0 
    dplyr::mutate(Genome_copies_mean= replace_na(Genome_copies_mean, 0))%>%
    ##Get unique labels from qPCR data
    distinct(labels, .keep_all = TRUE)-> data.inf.exp

##Merging Infection experiment oocyst and weight loss data with qPCR data
##Check differences between two dataframes
setdiff(sample.data$labels, data.inf.exp$labels)

##Samples MTU, FJL, EHM, DRT, CEY were not taken 
##Sample CPY was taken but not extracted (Faeces not found in boxes)
##Sample FLM was collected and extracted but not processed for qPCR (DNA not found)
##We end with 235 samples out of the original 242

###Join all the data in the same dataframe
sdt<- inner_join(sample.data, data.inf.exp, by="labels") ## Add qPCR data
###Tiny adjustment  
sdt$dpi<- as.factor(sdt$dpi)

sdt%>%
  dplyr::mutate(Genome_copies_ngDNA= Genome_copies_mean/50, ## copies by ng of fecal DNA considering 1uL from 50 ng/uL DNA
                DNA_sample= Conc_DNA*30, ## Estimate total gDNA of sample
                DNA_g_feces= DNA_sample/fecweight_DNA,
                ## Transform it to ng fecal DNA by g of faeces
                Genome_copies_gFaeces= Genome_copies_ngDNA*DNA_g_feces) -> sdt ## Estimate genome copies by g of faeces

##Transform to zero OPGs for DPI 1 and 2 
sdt$OPG[sdt$dpi==1] <- 0
sdt$OPG[sdt$dpi==2] <- 0                 
##Remove dataframes with data not related to the infection experiment data that won't be used in the following scripts
rm(data.std, data.std.lm, data.unk, data.unk.lm, data.spk, data.spk.lm, Sum.inf)
