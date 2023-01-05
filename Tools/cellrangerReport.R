# This script Reads all the cellranger summary metric files for all the samples
# and plots the aggregated data for easy comparison
# Load the below libraries, install them if not already installed
#install.packages("ggdendro")
library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(ggplot2)
library(ggdendro)
library(plotly)
library(htmlwidgets)
library(data.table)

##################################### Setup
# Set the input, output paths
directorypath="/Users/nagarajanv/OneDrive - National Institutes of Health/Zixuan/Results/cellranger/"
outputpath="/Users/nagarajanv/OneDrive - National Institutes of Health/Zixuan/Results/"
sampleinfopath="/Users/nagarajanv/OneDrive - National Institutes of Health/Zixuan/SampleInformation.tab"

# Read the results folders 
file_list <- list.files(path=directorypath)
# List result folders
#file_list = c("S1_1")
# Set the column headers
filehead=c("SampleName","Cells","Confidently mapped antisense","Confidently mapped to exonic regions","Confidently mapped to genome","Confidently mapped to intergenic regions","Confidently mapped to intronic regions","Confidently mapped to transcriptome","Mapped to genome","Median UMI counts per cell","Median genes per cell","Median reads per cell","Number of reads from cells called from this sample","Total genes detected","Cell-associated barcodes identified as multiplets","Cell-associated barcodes not assigned any CMOs","Cells assigned to other samples","Cells assigned to this sample","Number of reads","Number of short reads skipped","Q30 RNA read","Q30 UMI","Q30 barcodes","Cell-associated barcodes identified as multiplets","Cell-associated barcodes not assigned any CMOs","Cells assigned to a sample","Confidently mapped antisense","Confidently mapped reads in cells","Confidently mapped to exonic regions","Confidently mapped to genome","Confidently mapped to intergenic regions","Confidently mapped to intronic regions","Confidently mapped to transcriptome","Estimated number of cells","Mapped to genome","Mean reads per cell","Number of reads","Number of reads in the library","Sequencing saturation","Valid UMIs","Valid barcodes","Cell-associated barcodes identified as multiplets","Cells assigned to a sample","Estimated number of cell-associated barcodes","Median CMO UMIs per cell","Number of samples assigned at least one cell","Singlet capture ratio","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","Number of reads","Number of short reads skipped","Q30 RNA read","Q30 UMI","Q30 barcodes","Cell-associated barcodes identified as multiplets","Cell-associated barcodes not assigned any CMOs","Cells assigned to a sample","Estimated number of cell-associated barcodes","Fraction CMO reads","Fraction CMO reads usable","Fraction reads from multiplets","Fraction reads in cell-associated barcodes","Fraction unrecognized CMO","Mean reads per cell-associated barcode","Median CMO UMIs per cell-associated barcode","Number of reads","Samples assigned at least one cell","Sequencing saturation","Valid UMIs","Valid barcodes")
filehead2=c("Cells","Confidently mapped antisense","Confidently mapped to exonic regions","Confidently mapped to genome","Confidently mapped to intergenic regions","Confidently mapped to intronic regions","Confidently mapped to transcriptome","Mapped to genome","Median UMI counts per cell","Median genes per cell","Median reads per cell","Number of reads from cells called from this sample","Total genes detected","Cell-associated barcodes identified as multiplets","Cell-associated barcodes not assigned any CMOs","Cells assigned to other samples","Cells assigned to this sample","Number of reads","Number of short reads skipped","Q30 RNA read","Q30 UMI","Q30 barcodes","Cell-associated barcodes identified as multiplets","Cell-associated barcodes not assigned any CMOs","Cells assigned to a sample","Confidently mapped antisense","Confidently mapped reads in cells","Confidently mapped to exonic regions","Confidently mapped to genome","Confidently mapped to intergenic regions","Confidently mapped to intronic regions","Confidently mapped to transcriptome","Estimated number of cells","Mapped to genome","Mean reads per cell","Number of reads","Number of reads in the library","Sequencing saturation","Valid UMIs","Valid barcodes","Cell-associated barcodes identified as multiplets","Cells assigned to a sample","Estimated number of cell-associated barcodes","Median CMO UMIs per cell","Number of samples assigned at least one cell","Singlet capture ratio","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","CMO signal-to-noise ratio","Cells assigned to CMO","Fraction reads in cell-associated barcodes","Number of reads","Number of short reads skipped","Q30 RNA read","Q30 UMI","Q30 barcodes","Cell-associated barcodes identified as multiplets","Cell-associated barcodes not assigned any CMOs","Cells assigned to a sample","Estimated number of cell-associated barcodes","Fraction CMO reads","Fraction CMO reads usable","Fraction reads from multiplets","Fraction reads in cell-associated barcodes","Fraction unrecognized CMO","Mean reads per cell-associated barcode","Median CMO UMIs per cell-associated barcode","Number of reads","Samples assigned at least one cell","Sequencing saturation","Valid UMIs","Valid barcodes")

# Initiate an empty dataframe
filealldata = data.frame()

##################################### Extract and map short names
# Read sample information file
sampleinformation=read.csv(sampleinfopath,sep="\t")
# Extract shortname columns
samplenametype=unique(paste(sampleinformation$submitted_name,sampleinformation$ShortNameType,sep=","))
# Convert to dataframe
samplenametype=as.data.frame(samplenametype)
# Initiate emply lists to save hc and uv sample names
sampleS1list=list()
sampleS2list=list()
sampleS3list=list()

##################################### Extract and aggregate summary metrics
# Read individual summary metric file
for(i in 1:length(file_list))
{
  # File path
  filename=paste(directorypath,file_list[i],"/metrics_summary.csv",sep = "")
  # Store sample name
  filesample=unlist(file_list[i])
  # Identify and store short sample name
  filesampleoneline=samplenametype[grep(filesample,samplenametype[,1]), ]
  filesampleoneline=unlist(as.vector(strsplit(filesampleoneline, ',')))
  filesampleoneline=filesampleoneline[2]
  print(filesampleoneline)
  filesample=filesampleoneline
  # Create S1 name list
  if(grepl("S1",filesample))
  {
    sampleS1list<-append(sampleS1list,filesample)
  }
  # Create S2 name list
  if(grepl("S2",filesample))
  {
    sampleS2list<-append(sampleS2list,filesample)
  }
  # Read summary metric data
  #filedata <- read_csv(filename, col_types = cols(`Valid Barcodes` = col_number(), `Sequencing Saturation` = col_number(), `Q30 Bases in Barcode` = col_number(), `Q30 Bases in RNA Read` = col_number(), `Q30 Bases in UMI` = col_number(), `Reads Mapped to Genome` = col_number(), `Reads Mapped Confidently to Genome` = col_number(), `Reads Mapped Confidently to Intergenic Regions` = col_number(), `Reads Mapped Confidently to Intronic Regions` = col_number(), `Reads Mapped Confidently to Exonic Regions` = col_number(), `Reads Mapped Confidently to Transcriptome` = col_number(), `Reads Mapped Antisense to Gene` = col_number(), `Fraction Reads in Cells` = col_number(), `Antibody: Valid Barcodes` = col_number(), `Antibody: Sequencing Saturation` = col_number(), `Antibody: Q30 Bases in Barcode` = col_number(), `Antibody: Q30 Bases in Antibody Read` = col_number(), `Antibody: Q30 Bases in UMI` = col_number(), `Antibody: Fraction Antibody Reads` = col_number(), `Antibody: Fraction Antibody Reads Usable` = col_number(), `Antibody: Fraction Antibody Reads in Aggregate Barcodes` = col_number(), `Antibody: Fraction Unrecognized Antibody` = col_number(), `Antibody: Antibody Reads in Cells` = col_number()))
  filedata=transpose(fread(filename, select=c(6)))
  filedata <- data.frame(lapply(filedata, function(x) { gsub("( .*)", "", x) }))
  filedata <- data.frame(lapply(filedata, function(x) { gsub("%|,", "", x) }))
  
  # Aggregate summary metric data in to filealldata variable
  if(length(filedata[2,])>=79)
  {
    filedataframe=as.data.frame(filedata[1,])
    filedatanumbers=as.numeric(filedataframe[1,])
    fileonerow=c(filesample,filedatanumbers)
    filealldata=rbind(filealldata,c(filesample,filedatanumbers))
  }
}

##################################### Plot and save interactive multipanel image as html file
# Create list for ordering hc and uv samples
sampleS1list=unlist(sampleS1list)
sampleS2list=unlist(sampleS2list)
sampleS3list=unlist(sampleS3list)
sampletypelist=append(sampleS1list,sampleS2list)
print(sampletypelist)

# Add column headers
filealldata
forplot=filealldata
colnames(forplot)=filehead

# Exclude samplename column and include all other columns for plotting
dfm <- pivot_longer(forplot, -SampleName, names_to="variable", values_to="value")
intplot=ggplot(dfm, aes(x = factor(SampleName, levels=sampletypelist),y = as.numeric(value))) + 
  geom_bar(aes(fill = variable), position="stack", stat="identity")+
  facet_wrap(~variable, scales="free_y")+
  theme(legend.position = "none",
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90),
    strip.text = element_text(size=8),
    axis.title.x = element_blank())
# Generate the plot and save it as html in the "Reports" folder
l <- plotly::ggplotly(intplot)
htmlwidgets::saveWidget(as_widget(l), paste(outputpath,"cellrangerReport.html",sep=""))

# Cluster sample metrics and generate a dendrogram
nohead=filealldata[,-1]
temp2=filealldata[,1]
forheatmap <- sapply( nohead, as.numeric )
forheatmap
# Add column headers
rownames(forheatmap)=temp2
colnames(forheatmap)=filehead2
# Convert data to matrix
data=as.matrix(forheatmap)
# Calculate the vector similarity using euclidean distance method
dd <- dist(scale(data), method = "euclidean")
# Draw the dendrogram
hc <- hclust(dd, method = "ward.D2")
# Save the dendrogram as an image in the "Reports" folder
png(paste(outputpath,"cellrangerClusterReport.png",sep=""),res=300)
ggdendrogram(hc)
dev.off()

