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

## Select just standards data
## Estimate the number of genome copies per ng of gDNA
data.std.lm<- subset(data.std, Task== "Standard") ## Select just data from standards 
data.std.lm %>% 
  select(Sample.Name, Task, Ct, Cycler, Oocyst_count, Parasite, Genome_copies)%>%
  dplyr::mutate(Oocyst_DNA= Oocyst_count*(3.8E-4))%>% ##Estimation of DNA (ng) derived from Oocyst
  dplyr::mutate(DNA_PCR= Oocyst_DNA/30)%>% ##DNA (ng) in PCR considering 1uL from a stock of 30uL
  ##Considering that 1 ng of Eimeria gDNA is equivalent to 2.11E4 genome copies
  dplyr::mutate(Genome_copies_ng= (2.11E4)*DNA_PCR)-> data.std.lm  

##Inter-sample variation
data.unk.lm<-read.csv("data/Eimeria_quantification_Sample_data.csv")
  
##Define numeric and factor variables 
num.vars2 <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", 
                "Oocyst_1", "Oocyst_2", "Oocyst_3", "Oocyst_4", "Oocyst_5", "Oocyst_6", "Oocyst_7", "Oocyst_8", "Dilution_factor", "Volume", "Sporulated")
  fac.vars2 <- c("Well", "Sample.Name", "Detector", "Task", "Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction", "Strain")  

  data.unk.lm[, num.vars2] <- apply(data.unk.lm[, num.vars2], 2,
                               function (x) as.numeric(as.character(x)))
data.unk.lm[, fac.vars2] <- apply(data.unk.lm[, fac.vars2], 2, as.factor)

data.unk.lm%>%
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

##Infection experiment
data.inf<-read.csv("data/Eimeria_quantification_Inf_exp_data.csv")
data.inf%>%
    select(Content, Sample, Plate_number, Cq, Melt_Temperature)%>%
    dplyr::rename(Ct= Cq, labels= Sample, Task= Content, Tm= Melt_Temperature)-> data.inf
  
##Define numeric and factor variables 
num.vars3 <- c("Ct", "Tm")
fac.vars3 <- c("labels", "Task", "Plate_number")  
data.inf[, num.vars3] <- apply(data.inf[, num.vars3], 2,
                               function (x) as.numeric(as.character(x)))
data.inf[, fac.vars3] <- apply(data.inf[, fac.vars3], 2, as.factor)
  

rm(fac.vars, num.vars, fac.vars2, num.vars2, fac.vars3, num.vars3)

####### Standard curves #######
set.seed(2020)
## Comput simple linear models from standards
## "Genome copies modeled by Ct"

##Ct modeled by Oocyst counts; data from different Cyclers
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Ct, color= Cycler, shape= Parasite)) +
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
  annotation_logticks(sides = "b")-> A

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
lm.SC<- lm(log10(Genome_copies_ng)~Ct, data.std.lm)

##Model 7: Genome copies modeled by Ct and parasite as predictors
lm.SCPar<- lm(log10(Genome_copies_ng)~Ct+Parasite, data.std.lm)

##Model 8: Genome copies modeled by Ct and cycler as predictors
lm.SCCyc<- lm(log10(Genome_copies_ng)~Ct+Cycler, data.std.lm)

##Model 9: Genome copies modeled by Ct, parasite and cycler used as predictors
lm.SCAll<- lm(log10(Genome_copies_ng)~Ct+Parasite+Cycler, data.std.lm)

##Model 10: Genome copies modeled by Ct and parasite/cycle interaction (Check with Alice and Susi)
lm.SCInt<- lm(log10(Genome_copies_ng)~Ct+Parasite*Cycler, data.std.lm)

##Comparison of models 
compareLM(lm.SC, lm.SCPar, lm.SCCyc, lm.SCAll, lm.SCInt)

##Model 8 fit better the data... Cycler has major impact (again confirm expectations)!
##Linear model (Standard curve for the rest of experiments)
data.std.lm%>%
  ggplot(aes(x = Ct, y = Genome_copies_ng)) +
  geom_smooth(method = "lm", se = T, color="black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_y_log10("log 10 Eimeria genome copies/ng gDNA", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "l")->A

##Linear model Genome copies per ng modeled by Oocyst count 
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Genome_copies_ng)) +
  geom_smooth(method = "lm", se = F, color= "black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10("log 10 Eimeria genome copies/ng gDNA", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  labs(tag = "B)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "bl")-> B

##Model 11: Genome copies modeled by Oocyst count, parasite and cycle 
lm.SCOoc<- lm(log10(Genome_copies_ng)~log10(Oocyst_count)+Parasite+Cycler, data.std.lm)

### Using MODEL 8 to predict using different levels of the factor cycler
data.std.lm$predicted<- 10^predict(lm.SCCyc)
data.std.lm$residuals<- 10^residuals(lm.SCCyc)

##Predicted genome copies per ng of DNA
##modeled by Ct
ggplot(data.std.lm, aes(x = Ct, y = predicted)) +
  geom_smooth(method = "lm", se = T, color= "black") +
  guides(color = "none", size = "none") +  # Size legend also removed
  scale_y_log10("log 10 Predicted \n Eimeria genome copies/ng gDNA", 
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
  scale_y_log10("log 10 Predicted \n Eimeria genome copies/ng gDNA", 
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

###### Intersample variation experiment #####
##Predict genome copies per ng of gDNA using model 8
data.unk.lm$Genome_copies_ng<- 10^predict(lm.SCCyc, data.unk.lm)

data.unk.lm%>%
  ggplot(aes(Oocyst_count, Genome_copies_ng), geom=c("point", "smooth"))+
  scale_x_log10(name = "log10 Eimeria Oocysts count (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log 10 Eimeria genome copies/ng gDNA \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  #geom_jitter(shape=21, position=position_jitter(0.1), aes(size= 25, fill= Task), color= "black", alpha= 0.5)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(color= "black", method = "lm")+            
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)")+
  annotation_logticks(sides = "bl")->C

##Model 12: Intersample variation considering Parasite, strain, cycler and sporulation rate as predictors
lm.ISV<- lm(formula = log10(Genome_copies_ng)~log10(Oocyst_count)+Parasite+Strain+Cycler+Sporulation_rate, data = data.unk.lm)

##Compair model 11 (perfect fit) vs model 12 
##compareLM(lm.SCOoc, lm.ISV)
##Predict genome copies from oocyst count and compair against previous prediction (Not finished, check!)
##data.unk.lm$Genome_copies_ng_pred<- 10^predict(lm.ISV, data.unk.lm)

##Figure 2 Intersample variation
##pdf(file = "fig/Figure_2.pdf", width = 10, height = 8)
##grid.arrange(C)
##dev.off()

########## Mock samples Experiment #########
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
    geom_jitter(shape=21, position=position_jitter(0.2), color= "black",
                aes(size= 25, fill= Std_series))+
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
                       label.y = c(34, 31, 27, 24, 21, 17, 15, 0))
###Determine that 10^0 meassurments are basically like NTC

##Compair mock samples qPCR estimation with real oocyst count by two extraction methods
set.seed(2020)
data.unk%>%
    select(Sample.Name, Task, Ct,Qty,Cycler,Parasite, Sample_type, Feces_weight, Extraction, Oocyst_count)%>%
    filter(Sample_type=="Feces" & Task=="Unknown")%>%
    dplyr::mutate(Qty= 10^((Ct-36)/-3.1))%>%
    ggplot(aes(x = Oocyst_count, y = Qty, color=Extraction), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocysts (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies/µL gDNA", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2), color= "black", aes(size= 25, fill= Extraction), alpha= 0.5)+
    theme_bw() +
    geom_smooth(aes(color= Extraction, fill= Extraction), method = "lm")+
    facet_grid(cols = vars(Extraction))+
    stat_cor(aes(color = Extraction), label.x = log10(100),  label.y = log10(100000),method = "spearman")+
    stat_cor(label.x = log10(100), label.y = log10(50000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Extraction))+        # Add correlation coefficient
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
                                        #stat_regline_equation(aes(color = Std_series), label.x = 3, label.y = c(35, 31))
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "A)")+
    annotation_logticks(sides = "bl")

##Standard curve and ceramic beads data
data.std%>%
    dplyr::select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,
                  Parasite,Tm,Sample_type,Feces_weight,Extraction,Date)%>%
    filter(Task%in%c("Standard", "NTC") & Cycler=="ABI" & Std_series%in%c("A","B"))%>%
    dplyr::mutate(Qty= Qty*8)%>% ##Transform to Genome copies per uL gDNA qPCR 
    dplyr::group_by(Parasite)-> Std.mock

data.unk%>%
    dplyr::select(Sample.Name,Task,Ct,Qty,Cycler,Oocyst_count,
                  Parasite,Tm,Sample_type,Feces_weight,Extraction,Date, NanoDrop, Strain)%>%
    filter(Sample_type=="Feces" & Task=="Unknown" & Extraction!="Glass_beads")%>%
     ## Transform Ct to Genome copies per uL gDNA qPCR
    dplyr::mutate(Qty= 10^((Ct-36)/-3.1),
                  ## GC_ngDNA= Genome_copies/NanoDrop, Estimate Genome
                  ## copies by ng of fecal DNA
                  GC_ngDNA= Qty/50, ## Estimate Genome copies by ng of fecal DNA
                  DNA_sample= NanoDrop*40, ## Estimate total gDNA of sample
                  DNA_g_feces= DNA_sample/Feces_weight,
                  ## Transform it to ng fecal DNA by g of feces
                  GC_gfeces= GC_ngDNA*DNA_g_feces, ## Estimate genome copies by g of feces
                  OPG=Oocyst_count/Feces_weight) ->
    data.mock ## Estimate oocyst per g of feces for mock samples

data.mock$predicted.Gc<- 10^predict(lm.GC1, data.mock)
## data.mock$residuals.Gc<- 10^residuals(lm.GC1, data.mock) 
## ## Error in match.arg(type) : 'arg' must be NULL or a character vector

data.mock%>%
    bind_rows(data.std.lm)-> data.mock

## rm(Std.mock)

##
set.seed(2020)
data.mock%>%
    dplyr::select(Sample.Name,Task,Qty,Ct,Oocyst_count, predicted.Gc)%>%  
    filter(Task%in%c("Standard", "Unknown"))%>%
    ggplot(aes(x = Oocyst_count, y = predicted.Gc), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocyst count (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies/µL gDNA \n (qPCR)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2),
                aes(fill= Task), size= 5, color= "black", alpha= 0.5)+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(aes(color= Task, fill= Task), method = "lm")+            
    ## stat_cor(aes(color = Task), label.x = 2,  label.y = c(7, 6),method = "spearman")+
    stat_cor(label.x = 0.75, label.y =  5.5,
             aes(label = paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+
    ## Add correlation coefficient
    stat_regline_equation(aes(color = Task), label.x = 0.75, label.y = 6)+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "B)")+
    annotation_logticks(sides = "bl")->E#+
                                        #theme(legend.position = c(0.85, 0.25), legend.direction = "vertical",
                                        # Change legend key size and key width
                                        #legend.key.size = unit(0.25, "cm"),
                                        #legend.key.width = unit(0.15,"cm"))

pdf(file = "fig/Figure_2.2.pdf", width = 10, height = 8)
grid.arrange(E)
dev.off()

## summary(lm(formula = log10(Qty)~log10(Oocyst_count),
##            data = subset(data.mock, Task== "Standard")))

## # Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
## #  0 (non-NA) cases

## modelstd<- lm(formula = log10(Qty)~log10(Oocyst_count),
##               data = subset(data.mock, Task== "Standard"))

## Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
##   0 (non-NA) cases
 
summary(lm(formula = log10(Qty)~log10(Oocyst_count),
           data = subset(data.mock, Task== "Unknown" & Oocyst_count >0)))

## # Produces an error... but is also nowhere used after creating it??!!
## modelmock<- lm(formula = log10(Qty)~log10(Oocyst_count),
##                data = subset(data.mock, Task== "Unknown" & Oocyst_count >0))

## ## Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 

data.mock%>%
    dplyr::select(Sample.Name, Qty, Oocyst_count, Task)%>%
    filter(Task== "Unknown"& Oocyst_count >0)%>%
    dplyr::mutate(Qty_estimated= 10^(0.9+log10(Oocyst_count)), Percent_error= ((Qty_estimated- Qty)/Qty_estimated)*100)->Error_mock

mean_ci(Error_mock$Percent_error)
mean_sd(Error_mock$Percent_error)

data.mock%>%
    select(Sample.Name,Task,Qty,Ct,Oocyst_count, GC_gfeces, OPG)%>%  
    filter(Task== "Unknown")%>%
    ggplot(aes(x = OPG, y = Qty), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocysts per gram of feces (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies/µL gDNA \n (qPCR)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Task), color= "black", alpha= 0.5)+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(color= "black", method = "lm")+            
    stat_cor(label.x = 2,  label.y = 6,method = "spearman")+
    stat_cor(label.x = 2, label.y = 5.75,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+        # Add correlation coefficient
    stat_regline_equation(label.x = 2, label.y = 6.25)+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "C)")+
    annotation_logticks(sides = "bl")-> G

## pdf(file = "fig/Figure_2.pdf", width = 20, height = 15)
## grid.arrange(D,E,G, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 3)))
## dev.off()

## # Error in grob$wrapvp <- vp : object of type 'closure' is not subsettable
## ## In addition: There were 21 warnings (use warnings() to see them)

### why are those created in the first place if they are removed here?
### D was acutally never created!!!
## rm(data.mock, data.unk, data.std, D, E, G)

## Infection experiment --- well is this now what it is all about?

## Considering the standard curve generated with the data from the BioRad Cycler
## Ct = 42x -4(log10Number of genome number per uL gDNA) Figure 1.1B
## Number of genome copies = 10^((Ct-42)/-4)

## Estimate number of genome copies with qPCR Ct value 

## Define real positive and negatives based on Tm 
data.inf %>% 
    dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                        Tm >= 80   ~ "Negative",
                                        Tm < 80 ~ "Positive")) -> data.inf #%>%
    #dplyr::mutate(Qty= 10^((Ct-42)/-4), Genome_copies= 10^((Ct-42)/-4)) -> data.inf

data.inf$Genome_copies<- 10^predict(lm.GC1, data.inf)
## data.inf$residuals<- 10^residuals(lm.GC1, data.inf) ## breaks

data.inf %>%
    select(Genome_copies,labels) %>% # select variables to summarise
    na.omit()%>%
    dplyr::group_by(labels)%>%
    dplyr::summarise_each(funs(Genome_copies_min = min, Genome_copies_q25 = quantile(., 0.25),
                               Genome_copies_median = median, Genome_copies_q75 = quantile(., 0.75), 
                               Genome_copies_max = max, Genome_copies_mean = mean, Genome_copies_sd = sd)) -> Sum.inf

## Tm values were not sumarised to avoid problems in the function 

data.inf<- inner_join(data.inf, Sum.inf, by= "labels")

data.inf%>%
    select(labels, Genome_copies_mean, Infection)%>%
    filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl","FML"))%>% ## Replace NAs in real negative samples to 0 
    dplyr::mutate(Genome_copies_mean= replace_na(Genome_copies_mean, 0))-> data.inf.exp

### why on earth creat something to delete it then? What is really needed here?
rm(data.inf, data.std, Sum.inf)


### There were 11 warnings (use warnings() to see them)
### I didn't take care of these for now!! TODO!!!
