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

## Set the working path
work.path <- "/home/data/t170345/project/comPrgML_2.0/"; setwd(work.path) 
code.path <- file.path(work.path, "Codes") 
data.path <- file.path(work.path, "InputData") 
res.path <- file.path(work.path, "Results") 
fig.path <- file.path(work.path, "Figures") 

## Load the script for model comparison
source(file.path(code.path, "Compare.R"))

## Read the training and test sets ---------------------------------------------

## Training Cohort -------------------------------------------------------------
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
rownames(Train_surv) <- Train_surv$ID
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]
colnames(Train_surv) <- c("sample","OS","OS.time")
Train_surv <- Train_surv[,-1]

## Validation Cohort -----------------------------------------------------------
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
rownames(Test_surv) <- Test_surv$ID
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]
Test_surv <- Test_surv[,c("ID","Cohort","OS","OS.time")]
colnames(Test_surv) <- c("sample","Cohort","OS","OS.time")
Test_surv <- Test_surv[,-1]

# Extract the same genes
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,]) 

# Data standardization
Train_set <-  scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) 
Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)

## read Public Signature ------------------------------------------------------------

## pubSIG1
pubSIG <- read.xlsx(file.path(data.path, "Table S4.xlsx"), startRow = 2) 
pubSIG <- split(pubSIG[, c("Symbol", "Coef")], pubSIG$Author)

## My Signature ----------------------------------------------------------------
mySIGname = "TEPGS" 
myAlgorithm = "CoxBoost+StepCox[forward]"

## mySIG
mySIG <- read.table(file.path(res.path, "TEPGS.RS_mat.txt"), header = T, check.names = F)
mySIG <- setNames(object = mySIG[[myAlgorithm]], nm = rownames(mySIG))

## Integrate signatures --------------------------------------------------------
signatures <- pubSIG
signatures[[mySIGname]] <- mySIG

## Calculate the C-index -------------------------------------------------------
model <- list(); cinfo <- list() 
log.file <- file.path(res.path, "makeCox.log") 
if (file.exists(log.file)) file.remove(log.file) 
log.file <- file(log.file, open = "a")
sink(log.file, append = TRUE, type = "message")
for (i in names(signatures)){
  if (class(signatures[[i]]) == "data.frame"){
    model[[i]] <- makeCox(Features = signatures[[i]]$Symbol, 
                          coefs = signatures[[i]]$Coef,      
                          SIGname = i,                       
                          unmatchR = 0.2,                    
                          Train_expr = Train_set,           
                          Train_surv = Train_surv,           
                          statusVar = "OS",                  
                          timeVar = "OS.time")               
  }else{
    model[[i]] = signatures[[i]]
  }
  
  cinfo[[i]] <- calCindex(model = model[[i]],                
                          name = i,                          
                          Test_expr = Test_set,              
                          Test_surv = Test_surv,             
                          Train_expr = Train_set,           
                          Train_surv = Train_surv,           
                          Train_name = "TCGA",               
                          #Train_expr = NULL,                
                          #Train_surv = NULL,                
                          CohortVar = "Cohort",              
                          metaCohort = TRUE,                 
                          statusVar = "OS",                  
                          timeVar = "OS.time")               
  message("")
}
closeAllConnections()

cinfo <- do.call(rbind, cinfo)
write.table(cinfo[,1:5], file = file.path(res.path,"cinfo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F) 
cinfo <- split(cinfo, cinfo$Cohort)

# Drawing -------------------------------------------------------------------------
CohortCol <- brewer.pal(n = length(cinfo), name = "Paired") 
names(CohortCol) <- names(cinfo)

# 
plots <- lapply(cinfo, function(plot.data){
  plot.data$method <- 
    factor(plot.data$method,
           levels = plot.data$method[order(plot.data$C, decreasing = F)])
  # compares two concordance indices: the statistical test is a two-sided Student t test for dependent samples.
  C.compare <- plot.data$C[plot.data$method == mySIGname]
  se.compare <- plot.data$se[plot.data$method == mySIGname]
  n.compare <- plot.data$n[plot.data$method == mySIGname]
  RS.compare <- plot.data$RS[plot.data$method == mySIGname][[1]]
  r.combined <- unlist(lapply(plot.data$RS, function(x) cor(x, RS.compare)))
  var.combined <- plot.data$se^2 + se.compare^2 - 2*r.combined*plot.data$se*se.compare
  p <- pt(abs((plot.data$C-C.compare))/(sqrt(var.combined)), n.compare - 1, lower.tail = F) * 2
  plot.data$label <- cut(p, breaks = c(0, 0.05, 0.01, 0.001, 0.0001))
  plot.data$label <- plyr::mapvalues(x = plot.data$label,
                                     from = c("(0,0.0001]", "(0.0001,0.001]", "(0.001,0.01]", "(0.01,0.05]"), 
                                     to = c("****", "***", "**", "*"))
  return(ggplot(plot.data, aes(x = method, y = C, fill = Cohort)) +
           geom_errorbar(aes(ymin = C - 1.96 * se, ymax = C + 1.96 * se), width = .1) +
           geom_point(color = CohortCol[unique(plot.data$Cohort)], size = 2.5) +
           geom_text(aes(x = method, y = max(plot.data$C + 1.96 * plot.data$se - 0.05), label = label)) +
           geom_hline(yintercept = 0.6, linetype = "dashed") +
           ggtitle(label = unique(plot.data$Cohort)) +
           coord_flip() + 
           theme_classic() +
           theme(panel.border = element_rect(fill = NA, size = 1),
                 axis.title = element_blank(),
                 legend.position = "none"))
})
# Merge and save forest diagrams
plot_grid(plotlist = plots, nrow = 1)
ggsave(file.path(fig.path, "comparison.pdf"), width = 16, height = 26)


