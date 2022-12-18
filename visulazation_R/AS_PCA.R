library(splitstackshape )
library(matrixStats)

setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/rMATS_files")
Nova1_AS <- read.delim("Nova1_AS.csv", header = TRUE, row.names = 1, sep = ",")
selected = c ("Acsl3", " Nrbp2", "Trem2", "Nrxn1", "Nrxn3", "Nrxn2", "Syngap1", "Nova1", "Nova2", "Emx1", "Vip", "Pvr", "Prox1", "Prox2", "Sst")
EnhancedVolcano(Nova1_AS, x= "IncLevelDifference", y = "PValue", lab = Nova1_AS$geneSymbol, selectLab = selected)
Nova1_AS=Nova1_AS[which(duplicated(Nova1_AS$GeneID)==FALSE),]
dim(Nova1_AS)
rowCounts= data.frame(Nova1_AS$IC_SAMPLE_1, Nova1_AS$IC_SAMPLE_2)
rowCounts=cSplit(rowCounts, "Nova1_AS.IC_SAMPLE_1", sep=",")
rowCounts=cSplit(rowCounts, "Nova1_AS.IC_SAMPLE_2", sep=",")
rowCounts= cbind(rowCounts, Nova1_AS$geneSymbol)
rownames(rowCounts)= Nova1_AS$GeneID
#rowCounts= na.omit(rowCounts)
#rowCounts[is.na(rowCounts)] <- 0
genelist= rowCounts$V2
rowCounts=rowCounts[,-9]
condition <- factor(c("Nova1", "Nova1", "Nova1", "Nova1","WT_sham", "WT_sham", "WT_sham", "WT_sham"  ))
coldata <- data.frame(row.names= colnames(rowCounts), condition )
mat= data.frame(rowCounts)
mat[is.na(mat)] <- 0
#Calculate Z score 
mat.z=(mat-rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
Heatmap(mat.z, cluster_rows = T, cluster_columns = T,  column_labels = colnames(mat.z), name = "Z-score", row_labels = genelist )

#if I want to select the significant genes 
sigs.df2 <- Nova1_AS[(Nova1_AS$PValue <0.05 ) & (abs(Nova1_AS$IncLevelDifference) > 0.1),]
#SynapticGenes=  c("Nova1", "Macf1", "Ercc5", "Setx", "Dnajc5", "Hspa8" )
#sigs.df2=Nova1_AS[Nova1_AS$symbol %in% SynapticGenes ,]
mat= data.frame(rowCounts)
rownames(mat)=rownames(rowCounts)
mat[is.na(mat)] <- 0
mat<- mat[rownames(sigs.df2),]
mat= na.omit(mat)
mat.z=(mat-rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
Heatmap(mat.z, cluster_rows = T, cluster_columns = T,  column_labels = colnames(mat.z), name = "Z-score" )


#Nova2
setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/rMATS_files")
Nova1_AS <- read.delim("Nova2_AS.csv", header = TRUE, row.names = 1, sep = ",")
selected = c ("Acsl3", " Nrbp2", "Trem2", "Nrxn1", "Nrxn3", "Nrxn2", "Syngap1", "Nova1", "Nova2", "Emx1", "Vip", "Pvr", "Prox1", "Prox2", "Sst")
EnhancedVolcano(Nova1_AS, x= "IncLevelDifference", y = "PValue", lab = Nova1_AS$geneSymbol, selectLab = selected)
Nova1_AS=Nova1_AS[which(duplicated(Nova1_AS$GeneID)==FALSE),]
dim(Nova1_AS)
rowCounts= data.frame(Nova1_AS$IC_SAMPLE_1, Nova1_AS$IC_SAMPLE_2)
rowCounts=cSplit(rowCounts, "Nova1_AS.IC_SAMPLE_1", sep=",")
rowCounts=cSplit(rowCounts, "Nova1_AS.IC_SAMPLE_2", sep=",")
rowCounts= cbind(rowCounts, Nova1_AS$geneSymbol)
rownames(rowCounts)= Nova1_AS$GeneID
#rowCounts= na.omit(rowCounts)
#rowCounts[is.na(rowCounts)] <- 0
genelist= rowCounts$V2
rowCounts=rowCounts[,-9]
condition <- factor(c("Nova2", "Nova2", "Nova2", "Nova2","WT_sham", "WT_sham", "WT_sham", "WT_sham"  ))
coldata <- data.frame(row.names= colnames(rowCounts), condition )
mat= data.frame(rowCounts)
mat[is.na(mat)] <- 0
rownames(mat)= rownames(rowCounts)
#Calculate Z score 
mat.z=(mat-rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
Heatmap(mat.z, cluster_rows = T, cluster_columns = T,  column_labels = colnames(mat.z), name = "Z-score", row_labels = genelist )

#if I want to select the significant genes 
sigs.df2 <- Nova1_AS[(Nova1_AS$PValue <0.005 ) & (abs(Nova1_AS$IncLevelDifference) > 0.4),]#SynapticGenes=  c("Nova1", "Macf1", "Ercc5", "Setx", "Dnajc5", "Hspa8" )
#sigs.df2=Nova1_AS[Nova1_AS$symbol %in% SynapticGenes ,]
mat<- mat[sigs.df2$GeneID,]
mat.z=(mat-rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
colnames(mat.z) = c("Nova2", "Nova2", "Nova2", "Nova2","WT_sham", "WT_sham", "WT_sham", "WT_sham"  )
Heatmap(mat.z, cluster_rows = T, cluster_columns = T,  column_labels = colnames(mat.z), name = "Z-score", row_labels = sigs.df2$geneSymbol )
pheatmap(mat.z,labels_row =  sigs.df2$geneSymbol)


#GO Analysis
#gene_to_test <- rownames(sigs.df[sigs.df$log2FoldChange>0.5,])
gene_to_test <- rownames(mat)
GO_results <- enrichGO(gene = gene_to_test, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "ALL")
View(as.data.frame(GO_results))
plot(barplot(GO_results, showCategory = 30 ))


gene_to_test <- sigs.df2$geneSymbol
GO_results <- enrichGO(gene = gene_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL")
plot(barplot(GO_results, showCategory = 30 ))

#ECS
setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/rMATS_files")
Nova1_AS <- read.delim("ECS_AS.csv", header = TRUE, row.names = 1, sep = ",")

#NovaDB
setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/rMATS_files")
Nova1_AS <- read.delim("NovaDB_AS.csv", header = TRUE, row.names = 1, sep = ",")

#ECS_DBL
setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/rMATS_files")
Nova1_AS <- read.delim("ECS_DB_AS.csv", header = TRUE, row.names = 1, sep = ",")
#####################################################################3
#Compare the GE (sigs.df) vs AS (Nova1_AS)
Nova1_AS=Nova1_AS[which(duplicated(Nova1_AS$GeneID)==FALSE),]

sigs.df <- sigs.df[which(sigs.df$symbol %in% Nova1_AS$geneSymbol),]
Nova1_AS <- Nova1_AS[which(sigs.df$symbol %in% Nova1_AS$geneSymbol),]

Combined_table<- data.frame(sigs.df$symbol,sigs.df$log2FoldChange, sigs.df$pvalue, Nova1_AS$IncLevelDifference, Nova1_AS$PValue)
rownames(Combined_table) <- rownames(sigs.df)

#write.csv(Combined_table,"GEvsAS_Nova2.csv", row.names = TRUE)

#View(Combined_table[(Combined_table$sigs.df.log2FoldChange>0.2 ) & (abs(Combined_table$Nova1_AS.IncLevelDifference) > 0.09),])
High<-Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) > 0.2 ) & (abs(Combined_table$Nova1_AS.IncLevelDifference) > 0.04),]
Low<- Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) < 0.2 ) & (abs(Combined_table$Nova1_AS.IncLevelDifference) < 0.04),]
Diff_GE <- Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) > 0.2 ) & (abs(Combined_table$Nova1_AS.IncLevelDifference) < 0.04),]
Diff_AS <- Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) < 0.2) & (abs(Combined_table$Nova1_AS.IncLevelDifference) > 0.04),]


##Scalling 
Scalled_table= Combined_table
Scalled_table$sigs.df.log2FoldChange= scale(abs(Scalled_table$sigs.df.log2FoldChange))
Scalled_table$Nova1_AS.IncLevelDifference= scale(abs(Scalled_table$Nova1_AS.IncLevelDifference))



aa= Scalled_table[, c(2,4)]
aa$ID= rownames(aa)
colnames(aa)= c("FC_GE", "FC_AS", "ID")
library(reshape2)
reshape  <- melt(aa, id.var = "ID")
library(ggplot2)
ggplot(reshape, aes(x = ID, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") 


#I need to remove abs()
Scalled_table= Combined_table
Scalled_table$sigs.df.log2FoldChange= scale(Scalled_table$sigs.df.log2FoldChange)
Scalled_table$Nova1_AS.IncLevelDifference= scale(Scalled_table$Nova1_AS.IncLevelDifference)

aa= Scalled_table[, c(2,4)]
aa$GeneSymbol= Scalled_table$sigs.df.symbol
colnames(aa)= c("FC_GE", "FC_AS", "GeneSymbol")
library(reshape2)
reshape  <- melt(aa, id.var = "GeneSymbol")
library(ggplot2)
ggplot(reshape, aes(x = GeneSymbol, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") 

#To get dipper
Scalled_table= Combined_table
Scalled_table=Scalled_table[1:35,]

#Cleanup 
#All
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library("org.Mm.eg.db")
library(EnhancedVolcano)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)

setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/CountMatrix_STAR")
Counts <- read.delim("CountReads_csv.csv", header = TRUE, row.names = 1, sep = ",")
#Reduce #replocates 
Counts <- Counts[,-1]
library(tidyverse)
New_counts <-matrix(0, nrow= nrow(Counts), ncol=1)
for (i in 1:70) {
  for (u in 1:35){
    
    
    
    select_counts=rowMeans(Counts %>% select(i,i+1))
    
    
  }
  New_counts <- cbind(New_counts, select_counts)
  
}

New_counts= New_counts[,-1]
index= which(1:ncol(New_counts)%%2!=1)
Counts2=New_counts[,index]
Counts2= cbind(Counts2, Counts[,71:74])
condition <- factor(c("Nova1", "Nova1", "Nova1","Nova2", "Nova2", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO" ,"ECS_2KO" ,"ECS_2KO" , "N1N2", "N1N2", "N1N2", "N1N2", "Nova1", "Nova1",  "Nova2", "Nova2","N1N2", "N1N2", "N1N2", "N1N2", "Nova1", "Nova1", "Nova2", "Nova2", "Nova2", "Nova2", "ECS_WT", "ECS_WT", "ECS_WT", "ECS_WT", "WT_sham", "WT_sham", "WT_sham", "WT_sham"  ))
colnames(Counts2)= condition

#Nova1
Counts= Counts2
Counts= Counts[,c(1:3,18:19, 26:27 ,36:39)]
condition <- factor(c("Nova1", "Nova1", "Nova1", "Nova1", "Nova1", "Nova1", "Nova1","WT_sham", "WT_sham", "WT_sham", "WT_sham"  ))
coldata <- data.frame(row.names= colnames(Counts), condition )

adds <- DESeqDataSetFromMatrix(countData = round(Counts), colData = coldata, design = ~condition)
adds<- DESeq(adds)
vsdata <- vst(adds, blind = FALSE)
res <- as.matrix(results(adds))
sigs <- na.omit(res)
sigs<- as.data.frame(sigs)
#sigs <- sigs[sigs$padj < 0.05, ]

#Plotting Volcano 
sigs.df <- as.data.frame(sigs)
sigs.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")
library(splitstackshape )
library(matrixStats)

setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/rMATS_files")
Nova1_AS <- read.delim("Nova1_AS.csv", header = TRUE, row.names = 1, sep = ",")
Nova1_AS=Nova1_AS[which(duplicated(Nova1_AS$GeneID)==FALSE),]

sigs.df <- sigs.df[which(sigs.df$symbol %in% Nova1_AS$geneSymbol),]
Nova1_AS <- Nova1_AS[which(sigs.df$symbol %in% Nova1_AS$geneSymbol),]

Combined_table<- data.frame(sigs.df$symbol,sigs.df$log2FoldChange, sigs.df$pvalue, Nova1_AS$IncLevelDifference, Nova1_AS$PValue)
rownames(Combined_table) <- rownames(sigs.df)

#I need to remove abs()
Scalled_table= Combined_table
Scalled_table$sigs.df.log2FoldChange= scale(Scalled_table$sigs.df.log2FoldChange)
Scalled_table$Nova1_AS.IncLevelDifference= scale(Scalled_table$Nova1_AS.IncLevelDifference)

aa= Scalled_table[, c(2,4)]
aa$GeneSymbol= Scalled_table$sigs.df.symbol
colnames(aa)= c("FC_GE", "FC_AS", "GeneSymbol")
library(reshape2)
reshape  <- melt(aa, id.var = "GeneSymbol")
library(ggplot2)
ggplot(reshape, aes(x = GeneSymbol, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") 

#Sort in high to low 
ggplot(reshape, 
       aes(x=reorder(GeneSymbol, value),y= value), 
       fill = variable) + 
  geom_bar(stat = "identity", position = "dodge")  +
  ggtitle("Total ")

reshape %>%
  mutate(class = fct_reorder(GeneSymbol, value)) %>%
  ggplot( aes(x=reorder(GeneSymbol, value), y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge")

reshape %>%
  mutate(class = fct_reorder(GeneSymbol, value)) %>%
  ggplot( aes(x=reorder(GeneSymbol, value), y=value, fill=variable)) + 
  geom_col()

ggplot(reshape, aes(x = reorder(GeneSymbol, -value), y = value, fill=variable)) + geom_bar(stat = "identity")

#To get dipper
Scalled_table= Combined_table
Scalled_table=Scalled_table[1:35,]

GE_Genes= Combined_table$sigs.df.symbol[which(Combined_table$sigs.df.pvalue < 0.008)]

AS_Genes=Combined_table$sigs.df.symbol[which(Combined_table$Nova1_AS.PValue < 0.000003)]

#write.csv(GE_Genes,"GE_Genes.csv", row.names = FALSE)
#write.csv(AS_Genes,"AS_Genes.csv", row.names = FALSE)
require(dplyr)
require(forcats)

reshape %>% 
  mutate(ordering = -as.numeric(variable) + value,
         GeneSymbol = fct_reorder(GeneSymbol, ordering, .desc = T)) %>% 
  ggplot(aes(GeneSymbol, value, fill = variable)) + geom_col()


p= reshape %>% 
  mutate(ordering = -as.numeric(variable) + value,
         GeneSymbol = fct_reorder(GeneSymbol, ordering, .desc = T)) %>% 
  ggplot(aes(GeneSymbol, value, fill = variable)) + geom_col()

require(scales) # for removing scientific notation
p <- p + scale_y_continuous(labels = comma)
# manually generate breaks/labels
labels <- seq(1901, 2000, length.out=10)
# and set breaks and labels
p <- p + scale_x_discrete(breaks=labels, labels=as.character(labels))
p


library(tidyverse)

ggplot(reshape %>% arrange(GeneSymbol, desc(value)) %>%
         mutate(GeneSymbol=factor(GeneSymbol, levels=unique(GeneSymbol))), 
       aes(x=GeneSymbol,y=value, fill = variable)) +
  geom_bar(stat="identity", show.legend=TRUE) +
  
  facet_grid(. ~ variable, scales="free_x", space="free_x") +
  scale_y_continuous(limits=c(-0.005, 1.05*max(reshape$value)), expand=c(0,0)) +
  theme_classic() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA))


p= ggplot(reshape %>% arrange(GeneSymbol, desc(value)) %>%
            mutate(GeneSymbol=factor(GeneSymbol, levels=unique(GeneSymbol))), 
          aes(x=GeneSymbol,y=value, fill = variable)) +
  geom_bar(stat="identity", show.legend=TRUE) +
  
  facet_grid(. ~ variable, scales="free_x", space="free_x") +
  scale_y_continuous(limits=c(-0.005, 1.05*max(reshape$value)), expand=c(0,0)) +
  theme_classic() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA))

require(scales) # for removing scientific notation
p <- p + scale_y_continuous(labels = comma)
# manually generate breaks/labels
labels <- seq(1901, 2000, length.out=10)
# and set breaks and labels
p <- p + scale_x_discrete(breaks=labels, labels=as.character(labels))
p

#Zoom in positive area
#ADD presynapthic genes 

SynapticGenes= c("Nrxn1", "Macf1", "Ercc5", "Setx" )
AS_Genes=Combined_table$sigs.df.symbol[which(Combined_table$Nova1_AS.IncLevelDifference >0.1)]
AS_Genes= c(AS_Genes,SynapticGenes )
rr=reshape[which(reshape$GeneSymbol%in%AS_Genes),]
rr %>% 
  mutate(ordering = -as.numeric(variable) + value,
         GeneSymbol = fct_reorder(GeneSymbol, ordering, .desc = T)) %>% 
  ggplot(aes(GeneSymbol, value, fill = variable)) + geom_col()


#Barplot in zoom genes for NOVA2
GE= rr[which(rr$variable=="FC_GE"),]
AS= rr[which(rr$variable=="FC_AS"),]
rr2=rbind(GE[1,], AS[1,])
dd= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[2,], AS[2,])
dd2= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[3,], AS[3,])
dd3= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))


rr2=rbind(GE[101,], AS[101,])
dd4= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[30,], AS[30,])
dd5= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[67,], AS[67,])
dd6= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[18,], AS[18,])
dd7= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[19,], AS[19,])
dd8= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[54,], AS[54,])
dd9= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))


library("ggpubr")

ggarrange(dd, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
          ncol = 3, nrow = 3)
###

#Corrlation Pvalue
Nova1_coor= sigs.df
Nova1_coor= Nova1_coor[1:22529,5 ]
Nova2_coor= sigs.df
Nova2_coor= Nova2_coor[1:22529,5 ]
N1N2_coor= sigs.df
N1N2_coor=N1N2_coor[1:22529, 5]
ECS_coor= sigs.df
ECS_coor= ECS_coor[1:22529,5 ]
ECSWT_coor= sigs.df
ECSWT_coor= ECSWT_coor[1:22529,5 ]

Cor_Events= cbind(Nova1_coor, Nova2_coor, N1N2_coor, ECS_coor, ECSWT_coor)
colnames(Cor_Events)= c("N1", "N2", "N1N2", "ECS2ko", "ECSWT")
Cor_Events= as.data.frame(Cor_Events)


corrplot(cor(Cor_Events[1:12,]))
corrplot(cor(Cor_Events),
         method = "number",
         type = "upper" # show only upper side
)


#Corr one by one 

corrplot(cor(as.data.frame(Nova1_coor[1:12]), as.data.frame(Nova2_coor[1:12])),
         method = "number",
         type = "upper" # show only upper side
)
corrplot(cor(Cor_Events))

#Corrlation LogFC
Nova1_coor= sigs.df
Nova1_coor= Nova1_coor[1:22529,2 ]
Nova2_coor= sigs.df
Nova2_coor= Nova2_coor[1:22529,2 ]
N1N2_coor= sigs.df
N1N2_coor=N1N2_coor[1:22529, 2]
ECS_coor= sigs.df
ECS_coor= ECS_coor[1:22529,2 ]
ECSWT_coor= sigs.df
ECSWT_coor= ECSWT_coor[1:22529,2 ]

Cor_Events= cbind(Nova1_coor, Nova2_coor, N1N2_coor, ECS_coor, ECSWT_coor)
colnames(Cor_Events)= c("N1", "N2", "N1N2", "ECS2ko", "ECSWT")
Cor_Events= as.data.frame(Cor_Events)

corrplot(cor(Cor_Events[1:12,]))

#GO Analysis
#gene_to_test <- rownames(sigs.df[sigs.df$log2FoldChange>0.5,])
gene_to_test <- rownames(sigs.df)
GO_results <- enrichGO(gene = gene_to_test, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "ALL")
GO_terms=as.data.frame(GO_results)
#view(as.data.frame(grep("synap", GO_terms$Description, value = TRUE)))
length(grep("synap", GO_terms$Description, value = TRUE))
CT= c(grep("synap", GO_terms$Description, value = TRUE))

All_Synap_CT= GO_terms[ which(GO_terms$Description%in%CT),]
presynaptic_CT= All_Synap_CT[ which(All_Synap_CT$Description%in% c(grep("presynap", GO_terms$Description, value = TRUE))),]
postsynaptic_CT= All_Synap_CT[ which(All_Synap_CT$Description%in% c(grep("postsynap", GO_terms$Description, value = TRUE))),]
synaptic_CT= All_Synap_CT[  !(All_Synap_CT$Description %in% c(grep( "presynap", GO_terms$Description, value = TRUE))),]
synaptic_CT= synaptic_CT[  !(synaptic_CT$Description %in% c(grep( "postsynap", GO_terms$Description, value = TRUE))),]

#PRESYNAPTIC
library(stringr)
str_extract(presynaptic_CT$geneID[1], "[A-Z]+\\d+")
str_extract_all(presynaptic_CT$geneID[1], "[A-Z]+\\d+")
presynap_list= str_extract_all(presynaptic_CT$geneID, "[A-Z]+\\d+")
presynap_list=unique(unlist(presynap_list))
#run condition
sigs.df2= sigs.df[which(rownames(sigs.df) %in% presynap_list),]

#POSTSYNAPTIC
library(stringr)
str_extract(postsynaptic_CT$geneID[1], "[A-Z]+\\d+")
str_extract_all(postsynaptic_CT$geneID[1], "[A-Z]+\\d+")
presynap_list= str_extract_all(postsynaptic_CT$geneID, "[A-Z]+\\d+")
presynap_list=unique(unlist(presynap_list))
#run condition
sigs.df2= sigs.df[which(rownames(sigs.df) %in% presynap_list),]


#SYNAPTIC
library(stringr)
str_extract(synaptic_CT$geneID[1], "[A-Z]+\\d+")
str_extract_all(synaptic_CT$geneID[1], "[A-Z]+\\d+")
presynap_list= str_extract_all(synaptic_CT$geneID, "[A-Z]+\\d+")
presynap_list=unique(unlist(presynap_list))
#run condition
sigs.df2= sigs.df[which(rownames(sigs.df) %in% presynap_list),]



#RUN FOR ALL
setwd("C:/Users/ALGHN0A/Desktop/P1_083122_NovaProtein/Nextflow_otput_QC/rMATS_files")
Nova1_AS <- read.delim("ECS_DB_AS.csv", header = TRUE, row.names = 1, sep = ",")
Nova1_AS=Nova1_AS[which(duplicated(Nova1_AS$GeneID)==FALSE),]

sigs.df2 <- sigs.df2[which(sigs.df2$symbol %in% Nova1_AS$geneSymbol),]
Nova1_AS <- Nova1_AS[which(sigs.df2$symbol %in% Nova1_AS$geneSymbol),]

Combined_table<- data.frame(sigs.df2$symbol,sigs.df2$log2FoldChange, sigs.df2$pvalue, Nova1_AS$IncLevelDifference, Nova1_AS$PValue)
rownames(Combined_table) <- rownames(sigs.df2)

#I need to remove abs()
Scalled_table= Combined_table
Scalled_table$sigs.df.log2FoldChange= scale(Scalled_table$sigs.df.log2FoldChange)
Scalled_table$Nova1_AS.IncLevelDifference= scale(Scalled_table$Nova1_AS.IncLevelDifference)

aa= Scalled_table[, c(2,4)]
aa$GeneSymbol= Scalled_table$sigs.df2.symbol
colnames(aa)= c("FC_GE", "FC_AS", "GeneSymbol")
library(reshape2)
reshape  <- melt(aa, id.var = "GeneSymbol")
library(ggplot2)

reshape %>%
  mutate(class = fct_reorder(GeneSymbol, value)) %>%
  ggplot( aes(x=reorder(GeneSymbol, value), y=value, fill=variable )) +  labs(x ="Genes") +
  geom_bar(stat = "identity", position = "dodge") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

reshape %>%
  mutate(class = fct_reorder(GeneSymbol, value)) %>%
  ggplot( aes(x=reorder(GeneSymbol, value), y=value, fill=variable)) +  labs(x ="Genes") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_col() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


library(tidyverse)
require(scales) # for removing scientific notation
p <- p + scale_y_continuous(labels = comma)
# manually generate breaks/labels
labels <- seq(1901, 2000, length.out=10)
# and set breaks and labels
p <- p + scale_x_discrete(breaks=labels, labels=as.character(labels))
p


#fouces on specific genes

rr=reshape
rr %>% 
  mutate(ordering = -as.numeric(variable) + value,
         GeneSymbol = fct_reorder(GeneSymbol, ordering, .desc = T)) %>% 
  ggplot(aes(GeneSymbol, value, fill = variable)) + geom_col() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


#Barplot in zoom genes for NOVA2
GE= rr[which(rr$variable=="FC_GE"),]
AS= rr[which(rr$variable=="FC_AS"),]
rr2=rbind(GE[47,], AS[47,])
dd= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())

rr2=rbind(GE[36,], AS[36,])
dd2= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())

rr2=rbind(GE[27,], AS[27,])
dd3= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())


rr2=rbind(GE[34,], AS[34,])
dd4= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())

rr2=rbind(GE[103,], AS[103,])
dd5= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())

rr2=rbind(GE[38,], AS[38,])
dd6= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())

rr2=rbind(GE[32,], AS[32,])
dd7= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())

rr2=rbind(GE[109,], AS[109,])
dd8= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())


rr2=rbind(GE[43,], AS[43,])
dd9= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title=element_blank())


library("ggpubr")

ggarrange(dd, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
          ncol = 3, nrow = 3)

#SD
rr2= reshape

my_se <- rr2 %>%
  group_by(GeneSymbol) %>%
  summarise(
    sd=sd(value))

my_se <- rr2 %>%
  group_by(GeneSymbol) %>%
  summarise(n=n(),
            sd=sd(value),
            se=sd/sqrt(n))


rr2$sd= my_se$sd[which(my_se$GeneSymbol %in% rr2$GeneSymbol)]
rr2$se= my_se$se[which(my_se$GeneSymbol %in% rr2$GeneSymbol)]


ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.2,
                position=position_dodge(.9)) 

ggplot(data=rr3, aes(x=GeneSymbol, y=value, fill=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("red" , "lightblue4")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.2,
                position=position_dodge(.9)) 