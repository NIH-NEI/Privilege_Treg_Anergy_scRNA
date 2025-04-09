# This R script contains data and code used to generate all Figures for Ocular immune privilege in action manuscript.
library(Seurat)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(patchwork)
library(EnhancedVolcano)
library(data.table)
library(pheatmap)

load("eyenaive.RData")
Idents(eyenaive) <- "RNA_snn_res.0.2"

#########################################
# Fig. 1D
# Show sorted naive before injection and retrieved from eye after on week
DimPlot(object = eyenaive, pt.size=0.2, reduction = "umap", group.by = "type",cols = c("#2dc2a9","salmon"))
ggsave("fig1d-umap-eyenaive-type.png")

#########################################
# Fig. 1E
# Primed phenotype and Treg and proliferation marker
FeaturePlot(eyenaive, features=c("Cd44", "Sell","Foxp3", "Mki67"), ncol = 2)
ggsave("fig1e-featureplot-eyenaive-cd44.pdf", height = 5, width = 7)

#########################################
# Fig. 1F
# Cluster identification
eyenaive$eyenaive_annotation=recode(as.character(eyenaive$RNA_snn_res.0.2),
                                    "0"="nfc1",
                                    "1"="nfc2",
                                    "2"="naiveT",
                                    "3"="Treg",
                                    "4"="proliferative")

Idents(eyenaive) <- "eyenaive_annotation"
DimPlot(object = eyenaive, pt.size=0.2, reduction = "umap",label = T, cols = c("#23acde","#bfa708","#de6de8", "#22c981","salmon"))
ggsave("fid1f-umap-eyenaive-annotated.pdf", height = 4, width = 7)

#########################################
# Fig. 2A
# Cluster identification
Idents(eyenaive) <- "RNA_snn_res.0.2"
Idents(eyenaive) <- factor(Idents(eyenaive),levels=c(2,0,1,3,4))
Vlngene <- c('Foxp3', 'Tbx21','Rorc', 'Gata3', 'Bcl6', 'Il2', 'Ifng', 'Il17a', 'Csf2', 'Il10', 'Tgfb1')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("fig2a.pdf", height = 20, width = 3)

#########################################
# Fig. 2B
# Th lineage-specific markers
# Th1 related
FeaturePlot(eyenaive, features=c("Tbx21", "Ccr5", "Cxcr3"), ncol = 3)+NoLegend()
ggsave("fig2b-th1-related.pdf", height = 3, width = 12)
# Th17 related
FeaturePlot(eyenaive, features=c("Rorc", "Ccr6", "Il17a"), ncol = 3)+NoLegend()
ggsave("fig2b-th17-related.pdf", height = 3, width = 12)
# Treg related
FeaturePlot(eyenaive, features=c("Il2ra", "Il10", "Tgfb1"), ncol = 3)+NoLegend()
ggsave("fig2b-treg-related.pdf", height = 3, width = 12)

#########################################
# Fig. 2C
# naive related markers
FeaturePlot(eyenaive, features=c("Ccr7", "Lef1", "Tcf7"), ncol = 3)+NoLegend()
ggsave("fig2c-naive-related.pdf", height = 3, width = 12)

#########################################
# Fig. 3A
# Volcano plot
# nfc1 vs naive
Vnfc1vsnaive<-read.delim("DEG_nfc1vsnaive.txt",sep = "\t")
pval_0 <- subset(Vnfc1vsnaive, Vnfc1vsnaive$p_val_adj == 0)
# creating a random list of numbers around 1e-300. 
pval_jitter <- round(rnorm(nrow(pval_0), mean =150, sd = 35),0) 
pval_jitter <- 1*10^-pval_jitter
# replace old p-values by new values
pval_0$p_val_adj <- pval_jitter
# creating a new data set for visualization 
Vnfc1vsnaive <- subset(Vnfc1vsnaive,Vnfc1vsnaive$p_val_adj!=0)
Vnfc1vsnaive <- rbind(Vnfc1vsnaive, pval_0)
# Set custom colors and labels
keyvals <- ifelse(
  Vnfc1vsnaive$avg_log2FC < -1, 'dodgerblue3',
  ifelse(Vnfc1vsnaive$avg_log2FC > 1, 'brown2',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'brown2'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey'] <- 'Non-DEG'
names(keyvals)[keyvals == 'dodgerblue3'] <- 'Downregulated'
p1 <- EnhancedVolcano(Vnfc1vsnaive,rownames(Vnfc1vsnaive),
                      x ="avg_log2FC", y ="p_val_adj",title = 'nfc1 vs naive',
                      selectLab = c("Nrgn","Cblb","Nt5e","Nr4a1","Nr4a2","Nr4a3","Ctla4","Stat3","Tox",
                                    "Tgfb1","Maf","Lag3","Ikzf2","Foxp1","Sell","Pik3r5","Satb1","Pag1","S100a4","Cd200",
                                    "Nrp1","Tnfrsf4","Pdcd1","Ccl5","Tnfrsf1b", "Tcf7","Itgav","Itga4","Ccr7","Aff3","Il7r","Tox2"),
                      pCutoff = 10e-10,FCcutoff = 1,cutoffLineCol = 'darkgrey',vlineCol = 'darkgrey', pointSize = 1.0,labSize = 3,vline= c(-0.25,0.25),
                      colCustom = keyvals,drawConnectors = TRUE,arrowheads = FALSE,
                      lengthConnectors = unit(3, 'npc'),legendIconSize = 0, widthConnectors = 0.5,max.overlaps=50)
p2 <- p1 + ggplot2::coord_cartesian(xlim=c(-3, 3.5)) +
  ggplot2::scale_x_continuous(breaks=seq(-3, 3.5,1))+
  ggplot2::theme_bw()
p2
ggsave("fig3a-nfc1.pdf", plot = p2, width = 6, height = 4)

# nfc2 vs naive
Vnfc2vsnaive<-read.delim("DEG_nfc2vsnaive.txt",sep = "\t")
pval_0 <- subset(Vnfc2vsnaive, Vnfc2vsnaive$p_val_adj == 0)
# creating a random list of numbers around 1e-300. 
pval_jitter <- round(rnorm(nrow(pval_0), mean =150, sd = 35),0) 
pval_jitter <- 1*10^-pval_jitter
# replace old p-values by new values
pval_0$p_val_adj <- pval_jitter
# creating a new data set for visualization 
Vnfc2vsnaive <- subset(Vnfc2vsnaive,Vnfc2vsnaive$p_val_adj!=0)
Vnfc2vsnaive <- rbind(Vnfc2vsnaive, pval_0)
# Set custom colors and labels
keyvals <- ifelse(Vnfc2vsnaive$avg_log2FC < -1, 'dodgerblue3',
                  ifelse(Vnfc2vsnaive$avg_log2FC > 1, 'brown2','grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'brown2'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey'] <- 'Non-DEG'
names(keyvals)[keyvals == 'dodgerblue3'] <- 'Downregulated'
p1 <- EnhancedVolcano(Vnfc2vsnaive,rownames(Vnfc2vsnaive),
                      x ="avg_log2FC", y ="p_val_adj",title = 'nfc2 versus naive',
                      selectLab = c("Nrgn","Cblb","Nt5e","Nr4a1","Nr4a2","Nr4a3","Ctla4","Il17a","Lgals3","Prkcq","Itgav","Il2rb",
                                    "Tgfb1","Maf","Lag3","Foxp1","Sell","Pik3r5","Satb1","Pag1","S100a4",
                                    "Nrp1","Tnfrsf4","Pdcd1","Tnfrsf1b", "Tcf7","Rgs1","Ccr7","Aff3","Tox2"),
                      pCutoff = 10e-10,FCcutoff = 1,cutoffLineCol = 'darkgrey',vlineCol = 'darkgrey', pointSize = 1.0,labSize = 3,vline= c(-0.25,0.25),
                      colCustom = keyvals,drawConnectors = TRUE,arrowheads = FALSE,
                      lengthConnectors = unit(3, 'npc'),legendIconSize = 0, widthConnectors = 0.5,max.overlaps=100)
p2 <- p1 + ggplot2::coord_cartesian(xlim=c(-3, 4)) +
  ggplot2::scale_x_continuous(breaks=seq(-3, 4,1))+
  ggplot2::theme_bw()
p2          
ggsave("fig3a-nfc2.pdf", plot = p2, width = 6, height = 4)

# Treg vs naive
VTregvsnaive<-read.delim("DEG_Tregvsnaive.txt",sep = "\t")
pval_0 <- subset(VTregvsnaive, VTregvsnaive$p_val_adj == 0)
# creating a random list of numbers around 1e-300. 
pval_jitter <- round(rnorm(nrow(pval_0), mean =150, sd = 35),0) 
pval_jitter <- 1*10^-pval_jitter
# replace old p-values by new values
pval_0$p_val_adj <- pval_jitter
# creating a new data set for visualization 
VTregvsnaive <- subset(VTregvsnaive,VTregvsnaive$p_val_adj!=0)
VTregvsnaive <- rbind(VTregvsnaive, pval_0)
# Set custom colors and labels
keyvals <- ifelse(VTregvsnaive$avg_log2FC < -1, 'dodgerblue3',
                  ifelse(VTregvsnaive$avg_log2FC > 1, 'brown2','grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'brown2'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey'] <- 'Non-DEG'
names(keyvals)[keyvals == 'dodgerblue3'] <- 'Downregulated'
p1 <- EnhancedVolcano(VTregvsnaive,rownames(VTregvsnaive),
                      x ="avg_log2FC", y ="p_val_adj",title = 'Treg versus naive',
                      selectLab = c("Nt5e","Nr4a1","Nr4a2","Nr4a3","Ctla4","Tnfrsf4","Lgals3","Prkcq","Itgav","Il2rb","Bach2",
                                    "Tgfb1","Maf","Lag3","Foxp1","Sell","Pik3r5","Gzmb","Pag1","S100a4","S100a6","Ccl5",
                                    "Nrp1","Tnfrsf4","Pdcd1","Tnfrsf1b", "Tcf7","Rgs1","Itga4","Ccr7","Aff3","Foxp3","Tbx21","Lef1"),
                      pCutoff = 10e-10,FCcutoff = 1,cutoffLineCol = 'darkgrey',vlineCol = 'darkgrey', pointSize = 1.0,labSize = 3,vline= c(-0.25,0.25),
                      colCustom = keyvals,drawConnectors = TRUE,arrowheads = FALSE,
                      lengthConnectors = unit(3, 'npc'),legendIconSize = 0, widthConnectors = 0.5,max.overlaps=100)
p2 <- p1 + ggplot2::coord_cartesian(xlim=c(-4, 5)) +
  ggplot2::scale_x_continuous(breaks=seq(-4, 5,1))+
  ggplot2::theme_bw()
p2          
ggsave("fig3a-treg.pdf", plot = p2, width = 6, height = 4)

#########################################
# Fig. 3B,C,D,E
# GSEA Pathway Profiles
# Define the GSEA pathway that need to merge
pathway_names <- c("SAFFORD_T_LYMPHOCYTE_ANERGY", 
                   "GSE20366_TREG_VS_NAIVE_CD4_TCELL_HOMEOSTATIC_CONVERSION_UP", 
                   "GSE7852_TREG_VS_TCONV_UP",
                   "GSE14308_TH17_VS_NAIVE_CD4_TCELL_UP")

# the input path
nfc1_path <- "GSEA_Fig3_input_20250407/nfc1.GseaPreranked.1694187881280/"
nfc2_path <- "GSEA_Fig3_input_20250407/nfc2.GseaPreranked.1694134932965/"
treg_path <- "GSEA_Fig3_input_20250407/treg.GseaPreranked.1694353624854/" 
prolifer_path <- "GSEA_Fig3_input_20250407/prolifermannualy.GseaPreranked.1695817747363/" 
group_paths <- c(nfc1_path, nfc2_path,treg_path,prolifer_path)

#pathway_name <-pathway_names[1]
for (pathway_name in pathway_names) {
  file_name <- paste0(pathway_name,".tsv")
  merged_data <- data.table()
  
  for (group_path in group_paths) {
    file <- file.path(group_path, file_name)
    if (grepl("treg", group_path)) {
      group_name <- "treg"
    } else if (grepl("nfc1", group_path)) {
      group_name <- "nfc1"
    } else if (grepl("nfc2", group_path)) {
      group_name <- "nfc2"
    } else if (grepl("prolifer", group_path)) {
      group_name <- "prolifer"
    } else {
      group_name <- "unknown"
    }
    
    if (file.exists(file)) {
      data <- read.delim(file, sep = "\t", header = TRUE)
      data$pathway <- pathway_name
      data$group <- group_name
      merged_data <- rbind(merged_data, data)
    } else {
      print(paste(group_name, "file", file_name, "does not exist"))
    }
  }
  
  # plotting
  data<-merged_data
  t_name<- pathway_name
  
  color_palette=c("nfc1" = "#bfa708", "nfc2" = "#6fb878", "treg" = "#3d59a0","prolifer"= "#a519e6")
  p1<-ggplot(data) +
    aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES) + 
    scale_color_manual(values = color_palette) + 
    geom_line(aes(color = group),size = 2) +
    ylim(0,0.6)+
    labs( y = "Enrichment score (ES)",  title = t_name,x="Ranked in ordered dataset") +
    theme_bw(base_size = 26)+ 
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          #legend.position = c(0.9, 0.8), 
          legend.title = element_blank(), 
          legend.background = element_blank(), 
          plot.title=element_text(hjust=0.5))+ 
    scale_x_continuous(expand = c(0, 0)) + 
    geom_hline(yintercept = 0, linetype = "dashed")  # add ES as baseline 0
  p1
  p2 <- ggplot(data,aes(x = RANK.IN.GENE.LIST,y = group)) +
    scale_color_manual(values = color_palette) + 
    geom_vline(aes(xintercept = RANK.IN.GENE.LIST,color = group),show.legend = F,size = 1) +
    theme_classic(base_size = 12) +
    facet_wrap(~group,ncol = 1) +
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.ticks.length = unit(0.25,'cm'),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  p2
  rect.bm.col = c("#CC3333", "white", "#003366")
  rect.bm.label = c("Up regulated","Down regulated")
  rect_b <- data.frame(x = 1:max(data$RANK.IN.GENE.LIST), val = seq(-2, 2,length.out = length(1:max(data$RANK.IN.GENE.LIST))))
  label_b <- data.frame(x = c(0, max(data$RANK.IN.GENE.LIST)), y = 1, label = rect.bm.label,hjust = c(0, 1))
  
  rect_bp <- ggplot(rect_b) + 
    geom_tile(aes(x = x,y = 1, fill = val), show.legend = F) + 
    scale_fill_gradient2(low = rect.bm.col[1],mid = rect.bm.col[2], high = rect.bm.col[3], midpoint = 0) + 
    #theme(plot.background = element_rect(colour = "black")) + 
    geom_label(data = label_b, aes(x = x,y = y, label = label, hjust = hjust), fontface = "bold.italic") + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,0)) + 
    theme_bw(base_size = 12) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title.y = element_blank(),
          panel.grid = element_blank(), axis.line.x = element_blank(),
          strip.background = element_blank(), strip.text = element_blank(),
          panel.spacing = unit(0.1, "cm"), 
          plot.margin = margin(t = 0,r = 0.2, b = 0.5, l = 0.2, unit = "cm")) + 
    xlab("Rank in Ordered Dataset")
  rect_bp
  out_p <- p1/p2/rect_bp+plot_layout(heights = c(0.7,0.2,0.05))
  out_p 
  # export fig
  ggsave(filename = paste0("fig3a_",t_name, ".png"), plot = out_p , width = 12, height = 10, dpi = 300)
}

#########################################
# Fig. 5A
pathway = read.delim("ipa-immune.txt",header = T,sep="\t",row.names = 1)
ipa_pathway<-pheatmap(pathway, color = colorRampPalette(c("steelblue3", "white", "tomato2"))(500),
         cluster_rows = FALSE, cluster_cols = FALSE)
ggsave("fig5a-ipa.pdf", ipa_pathway, height = 5, width = 6)

#########################################
# Fig. 5B
Vlngene <- c('Ctla4','Lag3','Pdcd1','Fasl', 'Cd200')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("fig5b.pdf", height = 8, width = 3)

#########################################
# Fig. 5C
Vlngene <- c('Tgfbr1','Tgfbr2','Tgfbr3')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("fig5c.pdf", height = 5, width = 3)

#########################################
# Fig. 5D
Vlngene <- c('Cblb','Dgkz',"Nfatc1","Nfatc2","Nr4a1","Nr4a2","Nr4a3")
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("fig5d.pdf", height = 10, width = 3)

#########################################
# Fig. 5E
# Dotplot
dot_markers = c("Grap2","Prkcq","Prkcb","Card11","Pik3cd","Pik3r5","Pdk1","Pag1")
DotPlot(eyenaive, features = dot_markers)+coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(),  
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))  
labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
ggsave("fig5e-dotplot1.pdf",width = 4, height = 3)

#########################################
# Fig. 6D
FeaturePlot(eyenaive, features=c("Foxp3","Nrgn"), ncol = 2)+NoLegend() 
ggsave("fig6d-foxp3-nrgn.pdf", height = 3, width = 7)

#########################################
# Supplementary Fig. 1A
# Top markers heatmap
Idents(eyenaiveheatmap) <- factor(Idents(eyenaive),levels=c(2,0,1,3,4))
markersall_eyenaive<-read.delim("markersall_eyenaive.txt",sep = "\t")

# Extract the top 10 marker for each cluster and plot heatmap 
markersall_eyenaive %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
top10

DoHeatmap(eyenaiveheatmap, features = top10$gene, group.by = "RNA_snn_res.0.2") + scale_fill_gradientn(colors = c("steelblue4", "white", "firebrick3"))+  NoLegend()
ggsave("sfig1a-top10-heatmap.png", height = 6, width = 5)

#########################################
# Supplementary Fig. 1B
# Th17 associated
Vlngene <- c('Il17f','Il22',"Il23r")
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-th17-associated.pdf", height = 5, width = 3)
# Th2 associated
Vlngene <- c('Ccr3','Ccr8',"Il4", 'Il5')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-th2-associated.pdf", height = 7, width = 3)
# Tfh associated
Vlngene <- c('Cxcr5','Cxcl13',"Il21")
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-tfh-associated.pdf", height = 5, width = 3)
# Th22 associated
Vlngene <- c('Ahr','Ccr4',"Ccl7",'Il13','Il22')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-th22-associated.pdf", height = 8, width = 3)
# Th9 associated
Vlngene <- c('Spi1','Il4ra',"Il9")
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-th9-associated.pdf", height = 5, width = 3)

#########################################
# Supplementary Fig. 1C
# Other Tgf-b genes
Vlngene <- c('Tgfb2','Tgfb3')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-other-tgfb-associated.pdf", height = 4, width = 3)
# Tr1 associated gene
Vlngene <- c('Itga2')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-other-tr1-associated.pdf", height = 2, width = 3)
# Il35 genes
Vlngene <- c('Il12a','Ebi3')
VlnPlot(eyenaive,features =  Vlngene,pt.size=0,ncol=1)+NoLegend()
ggsave("sfig1b-il35-associated.pdf", height = 4, width = 3)

#########################################
sessionInfo()

