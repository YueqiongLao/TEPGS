library(tidyverse)
library(ISOpureR) 
library(impute) 
library(pRRophetic)  
library(SimDesign) 
library(ggplot2) 
library(cowplot) 

##############################################CTRP2.0##########################
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
#
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  
ctrp.ccl.anno <- read.table("./datafiles/CTRP_ccl_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) 
ctrp.cpd.anno <- read.delim("./datafiles/CTRP_cpd_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) 

############inputs

# 1. expression data
CRC <- read.table("./datafiles/CRC.data.txt",sep = "\t") 
group <- read.table("./datafiles/CRC_TEPGS_group.txt",sep = "\t")
testExpr <- CRC[rownames(CRC) %in% rownames(group),]
expr <- data.frame(t(testExpr))

# 2.Loading drug susceptibility AUC matrix and data processing
ctrp.auc <- read.table("./datafiles/CTRP_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

## a. Excluding compounds with > 20% missing data
ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.auc)]

## b. Removing data from haematopoietic_and_lymphoid_tissue
rmccl <- paste0("CCL",na.omit(ctrp.ccl.anno[which(ctrp.ccl.anno$ccle_primary_site == "haematopoietic_and_lymphoid_tissue"),"master_ccl_id"]))
rownames(ctrp.auc) <- paste0("CCL",rownames(ctrp.auc))
ctrp.auc <- ctrp.auc[setdiff(rownames(ctrp.auc),rmccl),]

## c. imputing the missing AUC values using the k-nearest neighbour (k-NN) imputation method
ctrp.auc.knn <- impute.knn(as.matrix(ctrp.auc))$data

## d. 
ctrp.auc.knn <- ctrp.auc.knn/ceiling(max(ctrp.auc.knn)) 


# Chemotherapeutic response prediction

# Load the expression profile of the CCLE cell line as the training set
ccl.expr <- read.table("./datafiles/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 

# Load gene annotation files for gene ID conversion
Ginfo <- read.table("./datafiles/overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 

# convert ensembl to gene symbol
ccl.expr <- ccl.expr[,-1]; rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)
comgene <- intersect(rownames(ccl.expr),rownames(Ginfo))
ccl.expr <- ccl.expr[comgene,]
ccl.expr$gene <- Ginfo[comgene,"genename"]; ccl.expr <- ccl.expr[!duplicated(ccl.expr$gene),]; rownames(ccl.expr) <- ccl.expr$gene; ccl.expr <- ccl.expr[,-ncol(ccl.expr)]

## CTRP
keepgene <- apply(ccl.expr, 1, mad) > 0.5 # 
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) 
trainPtype <- as.data.frame(ctrp.auc.knn)
ccl.name <- ccl.miss <- c() 
for (i in rownames(trainPtype)) {
  if(!is.element(gsub("CCL","",i),ctrp.ccl.anno$master_ccl_id)) {
    cat(i,"\n")
    ccl.miss <- c(ccl.miss, i) 
    ccl.name <- c(ccl.name, i) 
  } else {
    ccl.name <- c(ccl.name,  ctrp.ccl.anno[which(ctrp.ccl.anno$master_ccl_id == gsub("CCL","",i)),"ccl_name"]) 
  }
}

cpd.name <- cpd.miss <- c() 
for (i in colnames(trainPtype)) {
  if(!is.element(i,ctrp.cpd.anno$master_cpd_id)) {
    cat(i,"\n")
    cpd.miss <- c(cpd.miss, i) 
    cpd.name <- c(cpd.name, i) 
  } else {
    cpd.name <- c(cpd.name,  ctrp.cpd.anno[which(ctrp.cpd.anno$master_cpd_id == i),"cpd_name"]) 
  }
}

rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),] 
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)] 
comccl <- intersect(rownames(trainPtype),colnames(trainExpr)) 
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

# Test set
keepgene <- apply(expr, 1, mad) > 0.5
testExpr <- log2(expr[keepgene,] + 1)
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
all(rownames(trainExpr)==rownames(testExpr))
testExpr <- as.matrix(testExpr)
all(colnames(trainExpr)==rownames(trainPtype))

outTab <- NULL

# 
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) 
  ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = testExpr,
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.00001 
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab
write.table(ctrp.pred.auc,"./results/CTRP/ctrp.pred.auc.txt",sep = "\t",quote = F)






