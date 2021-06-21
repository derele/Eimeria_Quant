### Code to generate figure 4 based on Alice and Susi analysis

if(!exists("plotResidAlice_temp.RDS")){
  A<- readRDS("fig/plotResidAlice_temp.RDS") ## This RDS file is generated in code 4_Prediction_health_by_Max_OPG_fecDNA.R
}

if(!exists("ResDemSusana.rds")){
  B<- readRDS("fig/ResDemSusana.rds") ## This RDS file is generated in code 5_Prediction_health_by_OPG_fecDNA.R
}

##Figure 4: Prediction of impact in health (weight loss) by Eimeria genome copies and OPG
ggarrange(A, B, ncol = 1, nrow = 2, legend = "none")-> tmp.fig

tmp.fig

##To visualize it externally 
#ggsave(filename = "Rplots.pdf", tmp.fig, width = 8, height = 10)

ggsave(file = "fig/Figure_4_beforeLeg.pdf", tmp.fig,  width = 8, height = 8)
rm(A,B)
