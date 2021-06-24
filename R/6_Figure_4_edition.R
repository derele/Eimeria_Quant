### Code to generate figure 4 based on Alice and Susi analysis

if(!exists("plotResidAlice_temp")){
  A<- readRDS("fig/plotResidAlice_temp.RDS") ## This RDS file is generated in code 4_Prediction_health_by_Max_OPG_fecDNA.R
}

if(!exists("ResDemSusana")){
  B<- readRDS("fig/ResDemSusana.rds") ## This RDS file is generated in code 5_Prediction_health_by_OPG_fecDNA.R
}

if(!exists("legend")){
  C<- readRDS("fig/legend.rds") ## This RDS file is generated in code 5_Prediction_health_by_OPG_fecDNA.R
}

##Figure 4: Prediction of impact in health (weight loss) by Eimeria genome copies and OPG
require(grid)
require(gridExtra)
require(cowplot)

plot_grid(C, A, B, align = "v", nrow = 3, rel_heights = c(1, 5, 5))-> tmp.fig

#tmp.fig

##To visualize it externally 
#ggsave(filename = "Rplots.pdf", tmp.fig, width = 8, height = 10)

ggplot2::ggsave(file = "fig/Figure_4.pdf", tmp.fig,  width = 10, height = 10)
rm(A,B,C, tmp.fig)