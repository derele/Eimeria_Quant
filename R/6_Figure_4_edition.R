### Code to generate figure 4 based on Alice and Susi analysis

if(!exists("plotResidAlice_temp.RDS")){
  A<- readRDS("fig/plotResidAlice_temp.RDS") ## This RDS file is generated in code 4_Prediction_health_by_Max_OPG_fecDNA.R
}

if(!exists("ResDemSusana.rds")){
  B<- readRDS("fig/ResDemSusana.rds") ## This RDS file is generated in code 5_Prediction_health_by_OPG_fecDNA.R
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
#ggarrange(A,B, ncol = 1, nrow = 2)-> tmp.fig

##To visualize it externally 
#ggsave(filename = "Rplots.pdf", tmp.fig, width = 8, height = 10)

#ggsave(file = "fig/Figure_4.pdf", tmp.fig, width = 8, height = 10)
rm(A,B)