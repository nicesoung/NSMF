rm(list=ls())

library(tidyr)
library(dplyr)
library(biomaRt)
library(genefu)
library(data.table)
library(jsonlite)

rsem.results = Sys.glob("D:/TCGA/GDCdata/TCGA-COAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/*/*.tsv")
rsem.table <- fread(rsem.results[1], header = T, sep = "\t")
merge.df = data.frame(rsem.table[, 1])
colnames(merge.df) <- "gene_id"


for(rsem.result in rsem.results){
  
  
  sample.id <- basename(rsem.result)
  rsem.df <- fread(rsem.result, header = T)
  #rsem.df <- rsem.df[, -"transcript_id(s)"]
  sub.df <- rsem.df[, c("gene_id", "tpm_unstranded")]
  colnames(sub.df) <- c("gene_id", sample.id)
  merge.df <- cbind(merge.df, sub.df[, 2])
}
dim(merge.df)


merge.df[,1] <- data.frame(rsem.table[, 2])
merge.df <- unique(merge.df)
merge.df <- merge.df[-c(1),]
merge.df <- t(merge.df)
colnames(merge.df) <- merge.df[1,]
merge.df <- merge.df[-c(1),]
merge.df <- cbind(rownames(merge.df), data.frame(merge.df, row.names=NULL))



x.dir <- "D:/TCGA/GDCdata/TCGA-COAD/harmonized/Transcriptome_Profiling"

x.files <- Sys.glob(paste0(x.dir,"/*.json"))

for (x.file in x.files){
  temp.dir <- dirname(x.file)
  temp.id <- basename(temp.dir)
  
  temp.json <- fromJSON(x.file)
  temp.df <- do.call(rbind, temp.json$associated_entities)
  temp.id.df <- data.frame(file_id = temp.json$file_id, 
                           file_name = temp.json$file_name, 
                           entity_submitter_id = temp.df$entity_submitter_id)
  
  write.table(x = temp.id.df, file = paste0(temp.dir,"/", temp.id, ".sample.id"), quote = F, sep = "\t", row.names = F)
  
  temp <- apply(temp.id.df, 1, function(x) {
    temp.split <- strsplit(x[2], split = "\\.")
    temp.raw <- paste0(temp.dir, "/", x[1], "/", temp.split[[1]][1], ".", temp.split[[1]][2], ".", temp.split[[1]][3])
    temp.new <- paste0(temp.dir, "/", x[1], "/", x[3], ".", temp.split[[1]][2], ".", temp.split[[1]][3])
    
    system(paste("ln -s", temp.raw, temp.new, sep = " "))
    
    return(temp.new)
  })
}

TCGA <- merge(merge.df,temp.id.df, by.x = "rownames(merge.df)", by.y = "file_name")
TCGA <- data.frame(TCGA)
L <- dim(TCGA)[2]
TCGA <- TCGA[,c(L:(L-1),1:(L-2))]
TCGA$Type <- TCGA$entity_submitter_id
a <- dim(TCGA)[2]
TCGA <- TCGA[,c(1,a,2:(a-1))]
s <- TCGA$Type
TCGA$Type <- substr(s, 14, 16)
unique(TCGA$Type)

Tumor <- subset(TCGA,Type != "11A" & Type != "11B")
Normal <- subset(TCGA,Type == "11A"| Type == "11B")



unique(Tumor$Type)
dim(Tumor)

TPM <- Tumor[,-c(2,3,4)]
rownames(TPM) <- TPM$entity_submitter_id
TPM <- TPM[,-c(1)]


#TPM <- TCGA[,-c(2,3,4)]
#rownames(TPM) <- TPM$entity_submitter_id
#TPM <- TPM[,-c(1)]
write.csv(TPM,"C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_TPM_INPUT_RSS_SCORE.csv")

TPM <- read.csv("C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_TPM_INPUT_RSS_SCORE.csv")
TPM <- t(TPM)
colnames(TPM) <- TPM[1,]

TPM_COAD_TUMOR <- TPM[-c(1),]
write.csv(TPM_COAD_TUMOR,"C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_TPM_INPUT_RSS_SCORE_final.csv")
##########################################################################################################################
##########################################################################################################################
##############Replication stress 계산산
pathGeneTable = read.csv("C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/CCLE_Pathway_Gene_List.csv", header = T)
pathGeneTable

sampleRNA <- read.csv("C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_TPM_INPUT_RSS_SCORE_final.csv", header = T)
head(sampleRNA)
names(sampleRNA)[1] <- "Gene"
sampleRNA[,-c(1)]
martID <- useMart("ensembl",dataset="hsapiens_gene_ensembl") ## human
# martID <- useMart("ensembl",dataset="mmusculus_gene_ensembl") ## mouse

## annotation list
annotIDs <- getBM(attributes=c('hgnc_symbol', 'entrezgene_id', 'ensembl_gene_id', 'gene_biotype'),
                  filters = 'hgnc_symbol',
                  value = sampleRNA$Gene,
                  mart = martID) ## human
head(annotIDs)
# annotIDs <- getBM(attributes=c('hgnc_symbol', 'entrezgene_id', 'ensembl_gene_id', 'gene_biotype'),
#                   filters = 'hgnc_symbol',
#                   value = unique(sampleRNA[,1]),
#                   mart = martID) ## human

# annotIDs <- getBM(attributes=c('mgi_symbol', 'entrezgene_id', 'ensembl_gene_id', 'gene_biotype'), 
#                   filters = 'mgi_symbol', 
#                   value = unique(sampleRNA[,1]), 
#                   mart = martID) ## mouse

colnames(annotIDs)[1] <- "probe"
colnames(annotIDs)[2] <- "EntrezGene.ID"
colnames(annotIDs)[3] <- "Gene"

annotIDs <- annotIDs[!duplicated(annotIDs$probe),]
head(annotIDs)
# annotIDs <- annotIDs[!duplicated(annotIDs$EntrezGene.ID),]
# annotIDs <- annotIDs[annotIDs$probe %in% unique(annotIDs$probe), ]

sampleRNA <- merge(sampleRNA, annotIDs, by.x = "Gene", by.y = "probe")
head(sampleRNA)
# head(temp[,c(1, 178:180)])
L <- dim(sampleRNA)[2]
L
head(sampleRNA[,c(1, (L-2):L)])


sampleRNA$Gene.y <- sampleRNA$Gene
sampleRNA$Gene <- sampleRNA$EntrezGene.ID
head(sampleRNA)
dim(sampleRNA)

sampleRNA[,c((L-1))]
names(sampleRNA)[L-1] <- "probe"
sampleRNA <- sampleRNA[!duplicated(sampleRNA$probe),]
#sampleRNA <- sampleRNA[!duplicated(sampleRNA$EntrezGene.ID),]
dim(sampleRNA)
geneAnnotation <- sampleRNA$probe
# sampleRNA <- sampleRNA[,-c(178:180)]
L <- dim(sampleRNA)[2]
sampleRNA <- sampleRNA[,-c((L-2):L)]
sampleRNA[,2:dim(sampleRNA)[2]]
sampleRNA <- sampleRNA[,2:dim(sampleRNA)[2]]
# dataGroups <- c(groupAdeno, groupDuctal) ### for Adeno and Ductal
rownames(sampleRNA) <- geneAnnotation
expData <- sampleRNA
expData <- log2(expData+1) # for Chae or CCLE data 
####################ex##########################################################################################


###scoring

pN <- dim(pathGeneTable)[1]
pN
for(i in 1:pN){
  print(i)
  pathGeneTable[i,pathGeneTable[i,] != ""]
  l <- length(pathGeneTable[i,pathGeneTable[i,]!=""])
  pathName <- pathGeneTable[i,1:2]
  pathGeneSym <- as.character(pathGeneTable[i,3:l])
  
  symIDs <- getBM(attributes=c('hgnc_symbol', 'entrezgene_id', 'ensembl_gene_id', 'gene_biotype'), 
                  filters = 'hgnc_symbol', 
                  value = pathGeneSym, 
                  mart = martID)
  
  symIDs <- cbind(symIDs, symIDs[,1])
  colnames(symIDs)[1] <- "probe"
  colnames(symIDs)[2] <- "EntrezGene.ID"
  colnames(symIDs)[5] <- "coefficient"
  symIDs[,5] <- 1
  
  tmp1 <- unique(symIDs[,c(1,2,5)])
  tmp1 <- tmp1[!duplicated(tmp1$probe),]
  rownames(tmp1) <- tmp1$probe
  rownames(annotIDs) <- annotIDs$probe
  unique(annotIDs$probe)
  length(annotIDs$probe)
  length(unique(annotIDs$probe))
  
  print("sig")
  tmp <- sig.score(x=tmp1, data=t(expData), annot=annotIDs,
                   do.mapping=TRUE, signed=TRUE, verbose=TRUE)
  
  if(i==1){
    scoreTable <- tmp$score
  }else{
    scoreTable <- cbind(scoreTable, tmp$score)
  }
}


scoreTable
# expData <- sampleRNA[sampleRNA$Symbol %in% unique(annotIDs$probe), ]
colnames(scoreTable) <- pathGeneTable$V1
# write.csv(pathGeneTable[, c(1,2)], file = "pathName_chae_TCGA.csv")
# write.csv(scoreTable, file = "scoreTable_chae_TCGA.csv")
#write.csv(pathGeneTable[, c(1,2)], file = "/BiO/Live/nicesoung/NSMF/text/pathName_CCLE_Colorectal_Adenocarcinoma.csv")
write.csv(scoreTable, file = "C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_ScoreTable.csv")
write.csv(expData, file = "C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_expData.csv")



scoreTable

##################################################################################################################################
#################################################################################################################################
############# Heatmap final
library(tidyr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(viridis)

TCGA_COAD_TUMOR_scoreTable <- read.csv("C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_ScoreTable.csv", header = T)

TCGA_COAD_TUMOR_scoreTable$Type <- TCGA_COAD_TUMOR_scoreTable$X
s <- TCGA_COAD_TUMOR_scoreTable$Type
TCGA_COAD_TUMOR_scoreTable$Type <- substr(s, 14, 16)
unique(TCGA_COAD_TUMOR_scoreTable$Type)

Tumor <- subset(TCGA_COAD_TUMOR_scoreTable,Type != "11A" & Type != "11B")

Tumor <- Tumor[1:22]
rownames(Tumor) <- Tumor$X
Tumor <- Tumor[,-c(1)]
Tumor <- t(Tumor)
Tumor <- Tumor[,order(colSums(Tumor))]

tmp_r <- cbind.data.frame(Tumor)
#col_ann <- cbind.data.frame(type = rep("Tumor", ncol(Tumor))))
#rownames(col_ann) <- colnames(tmp_r)
#breaksList = seq(-8, 8, by = 0.5)
xx <- pheatmap(tmp_r, 
               #annotation_col = col_ann, 
               scale = "row", 
               cluster_rows = F, 
               cluster_cols =F, 
               fontsize_col = 10, 
               show_rownames = T, 
               show_colnames = F , 
               main = "TCGA-COAD" ,
               cellheight=10, 
               cellwidth = 1,
               color = colorRampPalette(c("blue", "white", "red"))(50),
               silent = TRUE
               
)


save_pheatmap_pdf <- function(x, filename, width=20, height=20) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


save_pheatmap_pdf(xx, "C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/pheatmap_TCGA_TUMOR.pdf")
write.csv(tmp_r, file = "C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_RSSSCORE_INPUT_Heatmap.csv")

########################################################################################################################################
########################################################################################################################################
#######correlation_score

expData <- read.csv("C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_expData.csv")
scoreTable <- read.csv("C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_ScoreTable.csv")
rownames(scoreTable) <- scoreTable$X
scoreTable
scoreTable <- scoreTable[,-c(1)]

rownames(expData) <- expData$X
expData <- expData[,-c(1)]


#expData[12353,]
#plot(as.numeric(expData[12353,]),rowSums(scoreTable))
#cor.test(as.numeric(expData[12353,]),rowSums(scoreTable))
#rowSums(scoreTable)
extmp <- expData
a <- dim(extmp)[1]
a
for(i in 1:a){
  Gene <- rownames(extmp)[i] 
  Gene_exp <- t(extmp[Gene,])
  correlation_score <- cor(Gene_exp, rowSums(scoreTable))
  P_value <- cor.test(Gene_exp, rowSums(scoreTable))$p.value  
  colnames(correlation_score)[1] <- "correlation_score" 
  correlation_table1 <- cbind(correlation_score, P_value)  
  
  if(i==1){
    correlation_table <- correlation_table1
  }else{
    correlation_table <- rbind(correlation_table, correlation_table1)
  }
}

rownames(correlation_table)
head(correlation_table)

write.csv(correlation_table,"C:/Users/user/OneDrive/바탕 화면/RSS_SCORE/TCGA/TCGA_COAD_TUMOR_Correlation.csv")

correlation_table
