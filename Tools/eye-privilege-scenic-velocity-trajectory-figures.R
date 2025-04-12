# This script contains code to reproduce figures reported in Ocular privilege manuscript.
# The data corresponding to the objects plotted here are provided in the "scenic_monocle_trajectory.RData" file.
# The detailed code for standard data processing carried out to generate these figure objects are provided in the main workflow file.
# The detailed script to generate scenic outputs is also provided as For Fig.6 pySCENIC.py.
# Load libraries
library(scCustomize)
library(ggplot2)
library(ComplexHeatmap)

load("scenic_monocle_trajectory.RData")

#########################################
# Fig. 6B
# Choose regulons to plot
Reguloneye = read.delim("regulonActivity_byGroup_Scaled_eyenaive0928.txt",header = T,sep="\t",row.names = 1)
Reguloneye <- Reguloneye[,c("naiveT","anergic_quiet","anergic_IL17","Treg","Proliferative")]

pdf("fig6b-regulon-heatmap.pdf")
Heatmap(
  Reguloneye,
  name                         = "Regulon activity",
  col                          = colorRampPalette(c( "steelblue3","white", "tomato3"))(200),
  show_row_names               = T,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 10),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = T,
  cluster_columns              = F)
dev.off()

#########################################
# Fig. 6C
# Choose regulons to plot
regulonsToPlot = c('Foxp3(+)', 'Prdm1(+)', 'Nfatc1(+)', 'Nfatc2(+)', 'Nr4a2(+)', 'Egr2(+)', 'Rara(+)','Rxra(+)')
# Extract corresponding data
regulon_plots<-FeaturePlot_scCustom(eyenaive,features = regulonsToPlot,colors_use = viridis_inferno_light_high)
# Save the plot
ggsave("fig-6c-regulon-activity.pdf", width = 10, height = 8)

#########################################
# Fig. 7C
# Trajectory plot colored by cluster
# Plot data object presented in the manuscript is in plot3, to create the plot with other parameters, run the below code with monocle2
#plot3 <- plot_cell_trajectory(HSMM_myo4, color_by = "clusters")+scale_color_manual(values = c("#bfa708","#22c981", "#23acde", "#de6de8"))
plot3
ggsave("fig7c-trajectory-clusters.pdf")

#########################################
# Fig. 7D
# Trajectory plot colored by pseudotime
# Plot data object presented in the manuscript is in plot1, to create the plot with other parameters, run the below code with monocle2
#plot1 <- plot_cell_trajectory(HSMM_myo4, color_by = "Pseudotime")+scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
plot1
ggsave("fig7d-trajectory-pseudotime.pdf")

#########################################
# Fig. 7E
# Heatmap showing bifurcation of gene expression dynamics along pseudotime
# The plot data object is shared here as 'p30'. Detailed code is presented in the main workflow script.
p30$ph
ggsave("fig7e-trajectory-branch.pdf", p30$ph)

#########################################
# Fig. 7F
# Detailed code that was used to generate this plot is provided 
# NOTE: The below code was tested to work with monocle2 in windows
pdf("fig7f-kineticsplot.pdf", width = 4, height = 6)
kinetic_genes <- row.names(subset(fData(HSMM_myo4),gene_short_name %in% c("Foxp3","Tcf7", "Tox2", "Nr4a3")))
plot_genes_branched_pseudotime(HSMM_myo4[kinetic_genes,],branch_point = 1,color_by = "clusters",ncol = 1, cell_size=0.2,branch_labels = c("angery state","Treg state" )) +scale_color_manual(values = c("#23acde","#bfa708","#de6de8", "#22c981"))
dev.off()
