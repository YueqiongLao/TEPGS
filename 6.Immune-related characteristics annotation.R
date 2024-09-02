library(GSVA)

##Load the CRC Expression Matrix (TPM Normalized) and Immune Cell Marker Files
CRC <- read.table("./CRC.data.txt",sep = "\t") 
geneset <- read.table("./ssGSEA.cell.txt",sep = "\t") 
geneset <- geneset %>% column_to_rownames("V1")%>%t()

a <- geneset
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
View(set)
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
View(l)
class(l)

####
dat <- as.matrix(t(CRC)) 
ssgsea<- gsva(dat,l, method='ssgsea', kcdf='Gaussian',abs.ranking=TRUE) #kcdf='Gaussian'

#Min-Max standardization
resm <- ssgsea
for (i in colnames(ssgsea)) {
  resm[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
}

save(resm,file = c("TCGA-CRC_immune_cell_result.resm.Rdata")) 

####
load("TCGA-CRC_immune_cell_result.resm.Rdata")
resm <- data.frame(t(ssgsea))
resm$Barcode <- rownames(resm)
group_list <- read.table("./TCGA_CRC_TEPGS_group.txt",sep = "\t")
group_list$Barcode <- rownames(group_list)

##
data <- left_join(group_list[,c(13:14)],resm,by=c("Barcode"))
data <- na.omit(data)
TME_New <-  melt(data,id.vars = c("riskScore","Barcode"))
colnames(TME_New) <- c("Group","Sample","Celltype","Composition")  
head(TME_New)

###plot

pdf(paste("./results/","ssGSEA_TEPGS",".pdf",sep=""), width = 16, height = 4)
P <- ggplot(TME_New,aes(x=Celltype,y=Composition,fill=Group))+
  geom_boxplot(color='black',width = 1,outlier.size=0.3,outlier.fill="white",outlier.color="white")+
  scale_fill_manual(values=c("#AC5048", "#8EA4C0"))+
  scale_color_manual(values=c("black","black"))+
  stat_compare_means(aes(group=Group),label = "p.signif",method = 'wilcox.test',size =3,bracket.size=0.5)+
  theme(panel.background  = element_rect(fill = 'white'),
        axis.ticks = element_line(colour =  'black'),
        axis.line =element_line(colour =  'black'),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab("cell proportion")
P
dev.off()








