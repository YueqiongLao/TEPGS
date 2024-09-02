library(dplyr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(futile.logger)

############################################GSE68086 Volcano map################

DEGs <- read.delim2("./GSE68086.CRC.DEGs.txt",sep = "\t",header = T)
DEGs[,c(6:13)] <- as.numeric(unlist(DEGs[,c(6:13)]))
df <- DEGs[,c(1,2,10:13)]
colnames(df) <- c("Id","Symbol","log2FoldChange","logCPM","pvalue","FDR")
##Determine the screening threshold：p＜0.05，|log2FC|＞0.5
pvalue = 0.05
log2FC = 0.5
##
df$group <- case_when(
  df[,3] > log2FC & df$pvalue < pvalue ~ "up",
  df[,3] < -log2FC & df$pvalue < pvalue ~ "down",
  TRUE ~ 'none'
)
df$'-log10(pvalue)' <- -log10(df$pvalue) 
df$group <- factor(df$group, levels = c("up","down","none"))

####ggplot2
mycol <- c("#D2352C","#497EB3","#d8d8d8")
mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "top",legend.title = element_blank()
  )

pdf(paste("./","TEP.GSE68086.Volcano.map",".pdf",sep=""), width = 6, height = 6)
p1 <- ggplot(data = df,
             aes(x = log2FoldChange, y = -log10(pvalue), color = group)) + 
  geom_point(size = 2.2) 
p2 <- p1 +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  mytheme
p3 <- p2 +
  geom_hline(yintercept = min(df$`-log10(pvalue)`),size = 0.7,color = "black",lty = "dashed") + 
  geom_vline(xintercept = c(-log2FC, log2FC),size = 0.7,color = "black",lty = "dashed") 
p3
dev.off()

##################################TCGA-CRC Volcano map #########################

DEG.CRC <- read.table("./DEseq2_sig.TCGA_CRC_TvsN.txt",sep = "\t") 
df <- tibble::rownames_to_column(DEG.CRC,"ID")
colnames(df) <- c("Symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","FDR")

##Determine the screening threshold：p＜0.05，|log2FC|＞0.5
pvalue = 0.05
log2FC = 0.5

#
df$group <- case_when(
  df[,3] > log2FC & df$pvalue < pvalue ~ "up",
  df[,3] < -log2FC & df$pvalue < pvalue ~ "down",
  TRUE ~ 'none'
)
df$'-log10(pvalue)' <- -log10(df$pvalue) 
df$group <- factor(df$group, levels = c("up","down","none"))

##ggplot2
mycol <- c("#D2352C","#497EB3","#d8d8d8")
mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "top",legend.title = element_blank()
  )
df <- arrange(df,df$`-log10(pvalue)`)

pdf(paste("./","TEP.TCGA_CRC.Volcano.map",".pdf",sep=""), width = 6, height = 6)
p1 <- ggplot(data = df,
             aes(x = log2FoldChange, y = -log10(pvalue), color = group)) + 
  geom_point(size = 2.2) 
p2 <- p1 +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  mytheme
p3 <- p2 +
  geom_hline(yintercept = c(-log10(pvalue)),size = 0.7,color = "black",lty = "dashed") + 
  geom_vline(xintercept = c(-log2FC, log2FC),size = 0.7,color = "black",lty = "dashed") 
p3
dev.off()

##################################Venn diagram##################################

DEGs <- read.delim2("./GSE68086.CRC.DEGs.txt",sep = "\t",header = T)
DEGs[,c(6:13)] <- as.numeric(unlist(DEGs[,c(6:13)]))
DEseq <- DEGs[with(DEGs,abs(logFC) > 0.5 & FDR < 0.05),]
up <- DEseq[DEseq$logFC >0,] 
down <- DEseq[DEseq$logFC <0,]

DEG.CRC <- read.table("./DEseq2_sig.TCGA_CRC_TvsN.0.5.txt",sep = "\t") 
DEG.CRC <- DEG.CRC[with(DEG.CRC,abs(log2FoldChange) > 0.5 & padj < 0.05),]
DEG.CRC.up <- DEG.CRC[DEG.CRC$log2FoldChange>0,]
DEG.CRC.down <- DEG.CRC[DEG.CRC$log2FoldChange<0,]

two.up <- DEG.CRC.up[rownames(DEG.CRC.up) %in% up$HGNC.symbol,] 
two.down <- DEG.CRC.down[rownames(DEG.CRC.down) %in% down$HGNC.symbol,]

#######################venn
venn <- function(..., 
                 labels = NULL,
                 Rplot = TRUE, 
                 PDF = NULL, 
                 PPT = NULL, 
                 cex.label = 10, 
                 cex.num = 8, 
                 alpha.fill = 0.7, 
                 col.fill = c("#D2352C","#497EB3", "#F1BB72", "#F3B1A0"), 
                 alpha.stroke = 1, 
                 size.stroke = 1,
                 col.stroke = "black",
                 show_percentage = FALSE
){
  
  getName <- function(...){
    as.character(substitute(list(...)))[-1]
  }
  
  data_list = list(...)
  if(is.null(labels)){
    names(data_list) = getName(...)
  } else {
    names(data_list) = labels
  }
  
  # build venn table
  Len = length(data_list)
  vennM = matrix("", Len, Len, dimnames = list(getName(...), getName(...)))
  for(i in 1:Len){
    for(j in i:Len){
      vennM[i, j] = length(intersect(data_list[[i]], data_list[[j]]))
    }
  }
  resM = apply(vennM, 1, as.numeric)
  dimnames(resM) = dimnames(vennM)
  
  if(Len > 4){
    message("More than 4 objects, will not plot venn.")
    return(resM)
  }
  
  # plot venn
  require(ggvenn)
  p <- suppressWarnings(ggvenn(
    data = data_list,
    show_percentage = show_percentage,
    set_name_size = cex.label,
    text_size = cex.num,
    stroke_alpha = alpha.stroke,
    stroke_size = size.stroke,
    fill_alpha = alpha.fill,
    fill_color = col.fill,
    stroke_linetype = "solid",
  ))
  
  if(Rplot){
    print(p)
  }
  
  if(!is.null(PDF)){
    message(paste0("File save to '", getwd(), "/", PDF, "'"))
    pdf(PDF, 8, 8)
    print(p)
    dev.off()
  }
  
  if(!is.null(PPT)){
    message(paste0("File save to '", getwd(), "/", PPT, "'"))
    require(export)
    graph2ppt(p, file = PPT, width = 8, height = 8)
  }
  
  return(resM)
}

s1  <- up$HGNC.symbol
s2 <- rownames(DEG.CRC.up)
venn(s1, s2, PDF = "./venn_up.pdf")

s1 <- down$HGNC.symbol
s2 <- rownames(DEG.CRC.down)
venn(s1, s2, PDF = "./venn_down.pdf")







