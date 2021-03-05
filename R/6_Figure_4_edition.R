### Code to generate figure 4 based on Alice and Susi analysis

if(!exists("plotResidAlice_temp.RDS")){
  A<- readRDS("fig/plotResidAlice_temp.RDS") ## This RDS file is generated in code 4_Alice_playground_handsoff_toberemovelater.R
}

if(!exists("ResDemSusana.rds")){
  B<- readRDS("fig/ResDemSusana.rds") ## This RDS file is generated in code 5_Susanasden_timeTMP.R
}



##Make tiny adjustment
# New facet label names for supp variable
A$data%>%
  mutate(iv = case_when(iv == "OPG" ~ "Maximum Oocysts per g of faeces",
                        iv == "Genome_copies_gFaeces"  ~ "Maximum Genome copies per g of faeces"))-> A$data

A+
  facet_grid(~ iv, scales = "free_x")+
  xlab("Measurments")+
  ylab("Maximum weight loss (%)")+
  labs(tag= "A)")+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "none")-> A

B$data%>%
  mutate(iv = case_when(iv == "OPG" ~ "Oocysts per g of faeces",
                        iv == "Genome_copies_gFaeces"  ~ "Genome copies per g of faeces"))-> B$data
B+
  facet_grid(~ iv, scales = "free_x")+
  xlab("Measurments")+
  ylab("Weight loss relative to DPI 0 (%)")+
  labs(tag= "B)")+
  theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "none")-> B

##Figure 4: Prediction of impact in health (weight loss) by Eimeria genome copies and OPG
#pdf(file = "fig/Figure_4.pdf", width = 10, height = 15)
grid.arrange(A,B)
#dev.off()
rm(A,B)