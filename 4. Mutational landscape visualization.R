library(circlize)
library(ComplexHeatmap)

#########mutational landscape in TEPGS-high group###############################
mat <- read.table('./TEP_high/onco_matrix.txt',sep = '\t')
#########Set the color
col = c("Missense_Mutation" = "#4A79A9", "Splice_Site" = "#ECBA7E", 
        "Nonsense_Mutation" = "#C23529","Frame_Shift_Del" = "#6AA153",
        "Frame_Shift_Ins" = "#895295","In_Frame_Del"="#a6d854","Multi_Hit" = "#945330")
#
alter_fun = list(
  background = alter_graphic("rect", fill = "#e9e9e9"),   
  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
  Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
  Nonsense_Mutation = alter_graphic("rect",fill = col["Nonsense_Mutation"]),
  Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
  Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
  In_Frame_Del=alter_graphic("rect", fill = col["In_Frame_Del"]),
  Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"])
)

# title
column_title = "OncoPrint for TCGA-CRC.TEP.high"

# Modify the tag name
heatmap_legend_param = list(title = "Alternations", 
                            at = c('Missense_Mutation', 
                                   'Multi_Hit',                       
                                   'Nonsense_Mutation',
                                   'Frame_Shift_Del', 
                                   'Frame_Shift_Ins', 'Splice_Site','In_Frame_Del'), 
                            labels = c('Missense_Mutation', 'Multi_Hit', 
                                       'Nonsense_Mutation',
                                       'Frame_Shift_Del', 
                                       'Frame_Shift_Ins', 'Splice_Site','In_Frame_Del'))

# Waterfall plot
pdf('./TEP_high/oncoPrint.pdf',width = 6,height = 4.5)
oncoPrint(mat,
          alter_fun = alter_fun, 
          alter_fun_is_vectorized = FALSE,
          pct_gp = gpar(fontsize = 6), 
          col = col,
          column_title = column_title)
dev.off()

# Add clinical information
pdata <- read.table("./TEP_high/TEP.high.cli.txt",sep = "\t")
pdata$Status <- ifelse(pdata$Status==0,"Alive","Dead")
pdata$Tumor_Sample_Barcode <- pdata$Tumor_Sample_Barcode %>% gsub("-",".",.)
rownames(pdata) <- pdata$Tumor_Sample_Barcode
#
comsample <- intersect(colnames(mat),rownames(pdata))
mat <- mat[,comsample]
pdata <- pdata[comsample,]
all(colnames(mat)==rownames(pdata))

#Define annotation information
ha <-HeatmapAnnotation(Age=pdata$Age,
                       Gender=pdata$Gender,
                       Stage=pdata$TNM,
                       Status=pdata$Status,
                       col = list(Gender = c("male" = "#c6eae4", "female" = "#80cdc1"),
                                  Status=c("Alive"="#dfdfdf", "Dead"="#b9b9b9"),
                                  Stage=c("I_II"="#8BC0EB","III_IV"="#DD80A0"),
                                  Age= colorRamp2(c(min(pdata$Age),median(pdata$Age), max(pdata$Age)), c("#c4d3e9","#f1a169","#b43a29"))),
                       TMB=anno_barplot(pdata$total_perMB_log,
                                        border=FALSE, gp = gpar(fill="grey2")),
                       show_annotation_name = TRUE,
                       annotation_name_gp = gpar(fontsize = 10))

#draw(ha)
pdf('./TEP_high/com_oncoPrint.high.pdf',width = 8,height = 6)
oncoplot_anno <- oncoPrint(mat,
                           alter_fun = alter_fun, 
                           alter_fun_is_vectorized = TRUE,
                           col = col,
                           remove_empty_columns = TRUE, remove_empty_rows = TRUE, 
                           row_order = 1:nrow(mat),
                           column_title = column_title,
                           top_annotation = ha)
oncoplot_anno
dev.off()



######################mutational landscape in TEPGS-low group###################

mat <- read.table('./TEP_low/onco_matrix_low.txt',sep = '\t',header = T)

#########Set the color
col = c("Missense_Mutation" = "#4A79A9", "Splice_Site" = "#ECBA7E", 
        "Nonsense_Mutation" = "#C23529","Frame_Shift_Del" = "#6AA153",
        "Frame_Shift_Ins" = "#895295","In_Frame_Del"="#a6d854","Multi_Hit" = "#945330",
        "In_Frame_Ins"="#BB9CC2","Nonstop_Mutation"="#E9B5B1","Translation_Start_Site"="#E6DBE9")
#
alter_fun = list(
  background = alter_graphic("rect", fill = "#e9e9e9"),   
  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
  Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
  Nonsense_Mutation = alter_graphic("rect",fill = col["Nonsense_Mutation"]),
  Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
  Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
  In_Frame_Del=alter_graphic("rect", fill = col["In_Frame_Del"]),
  In_Frame_Ins=alter_graphic("rect", fill = col["In_Frame_Ins"]),
  Nonstop_Mutation=alter_graphic("rect", fill = col["Nonstop_Mutation"]),
  Translation_Start_Site=alter_graphic("rect", fill = col["Translation_Start_Site"]),
  Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"])
)
# title
column_title = "OncoPrint for TCGA-CRC-TEPGS.low"

# Modify the tag name
heatmap_legend_param = list(title = "Alternations", 
                            at = c('Missense_Mutation', 
                                   'Multi_Hit',                       
                                   'Nonsense_Mutation',
                                   'Frame_Shift_Del', 
                                   'Frame_Shift_Ins', 'Splice_Site','In_Frame_Del',
                                   'In_Frame_Ins','Nonstop_Mutation','Translation_Start_Site'), 
                            labels = c('Missense_Mutation', 'Multi_Hit', 
                                       'Nonsense_Mutation',
                                       'Frame_Shift_Del', 
                                       'Frame_Shift_Ins', 'Splice_Site','In_Frame_Del',
                                       'In_Frame_Ins','Nonstop_Mutation','Translation_Start_Site'))

# Waterfall plot
pdf('./TEP_low/oncoPrint.low.pdf',width = 6,height = 4.5)
oncoPrint(mat,
          alter_fun = alter_fun, 
          alter_fun_is_vectorized = FALSE,
          pct_gp = gpar(fontsize = 6), 
          col = col,
          column_title = column_title)
dev.off()

# Add clinical information
pdata <- read.table("./TEP_low/TEP.low.cli.txt",sep = "\t")
pdata$Status <- ifelse(pdata$Status==0,"Alive","Dead")
pdata$Tumor_Sample_Barcode <- pdata$Tumor_Sample_Barcode %>% gsub("-",".",.)
rownames(pdata) <- pdata$Tumor_Sample_Barcode
###
comsample <- intersect(colnames(mat),rownames(pdata))
mat <- mat[,comsample]
pdata <- pdata[comsample,]
all(colnames(mat)==rownames(pdata))

#Define annotation information
ha <-HeatmapAnnotation(Age=pdata$Age,
                       Gender=pdata$Gender,
                       Stage=pdata$TNM,
                       Status=pdata$Status,
                       col = list(Gender = c("male" = "#c6eae4", "female" = "#80cdc1"),
                                  Status=c("Alive"="#dfdfdf", "Dead"="#b9b9b9"),
                                  Stage=c("I_II"="#8BC0EB","III_IV"="#DD80A0"),
                                  Age= colorRamp2(c(min(pdata$Age),median(pdata$Age), max(pdata$Age)), c("#c4d3e9","#f1a169","#b43a29"))),
                       TMB=anno_barplot(pdata$total_perMB_log,
                                        border=FALSE, gp = gpar(fill="grey2")),
                       show_annotation_name = TRUE,
                       annotation_name_gp = gpar(fontsize = 10))

#draw
pdf('./TEP_low/com_oncoPrint.low.pdf',width = 8,height = 6)
oncoplot_anno <- oncoPrint(mat,
                           alter_fun = alter_fun, 
                           alter_fun_is_vectorized = TRUE,
                           col = col,
                           remove_empty_columns = TRUE, remove_empty_rows = TRUE,
                           row_order = 1:nrow(mat),
                           column_title = column_title,
                           top_annotation = ha)
oncoplot_anno
dev.off()




