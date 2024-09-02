rm(list = ls())
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics) 
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

# Set the paths
work.path <- "./OptimalML"; setwd(work.path) 
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")

#########################Step1 Loads scripts for model training and validation

source(file.path(code.path, "Machine learning script.R"))

#########################Step2 Datasets preparation 

##Training Cohort -------------------------------------------------------------

#expression
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_expr <- Train_expr[rowSums(Train_expr > 0) > ncol(Train_expr) * 0.1, ]

#survival
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t")
Train_surv <- Train_surv[Train_surv$OS.time > 0, c("OS", "OS.time")] 
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,]

## Validation Cohort -------------------------------------------------------

#expression
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

#survival
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t",check.names = F,stringsAsFactors = F)
Test_surv <- Test_surv[Test_surv$OS.time > 0, c("Cohort","OS", "OS.time")] 
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,]

# Extract the same genes in training and validation datasets
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,]) 

#########################Step3. Model training and validation

## method list --------------------------------------------------------

methods <- read.xlsx(file.path(code.path, "method.list.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
head(methods)

########################Step3.1 Train the model -----------------------
model <- list()
set.seed(seed = 123)

for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method 
  method <- strsplit(method, "\\+")[[1]] 
  
  Variable = colnames(Train_expr)
  for (i in 1:length(method)){
    if (i < length(method)){
      selected.var <- RunML(method = method[i], 
                            Train_expr = Train_expr, 
                            Train_surv = Train_surv, 
                            mode = "Variable",       
                            timeVar = "OS.time", statusVar = "OS")
      if (length(selected.var) > 5) Variable <- intersect(Variable, selected.var)
    } else {
      model[[method_name]] <- RunML(method = method[i],
                                    Train_expr = Train_expr[, Variable],
                                    Train_surv = Train_surv,
                                    mode = "Model",
                                    timeVar = "OS.time", statusVar = "OS")
    }
  }
}
saveRDS(model, file.path(res.path, "model.rds"))

########################Step3.2 Evaluate the model -----------------------------

model <- readRDS(file.path(res.path, "model.rds"))

Cindexlist <- list()
for (method in methods){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], 
                                  Test_expr = Test_expr, 
                                  Test_surv = Test_surv, 
                                  Train_expr = Train_expr, 
                                  Train_surv = Train_surv, 
                                  Train_name = "TCGA", 
                                  cohortVar = "Cohort", 
                                  timeVar = "OS.time",
                                  statusVar = "OS") # 
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "Cindex.mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

########################Step4 Plot ---------------------------------------------

Cindex_mat <- read.delim2("./Cindex.mat.txt",sep = "\t",header = T) ###Test cohort
rownames(Cindex_mat) <- Cindex_mat$Model
avg_Cindex <- apply(Cindex_mat, 1, mean)          
avg_Cindex <- sort(avg_Cindex, decreasing = T)   
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]     
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) 
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "#89619D", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") 
names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)

cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              #col = c("#1CB8B2", "#FFFFFF", "#EEB849"), 
              col = c("#8EA4C0","#FFFFFF", "#AC5048"), 
              rect_gp = gpar(col = "black", lwd = 1), 
              cluster_columns = FALSE, cluster_rows = FALSE, 
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path(fig.path, "Cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 4, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())






