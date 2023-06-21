load("All Data.RData")

library(ggplot2)
library(survcomp)

#tumor vs normal analysis
m1_tumor_normal <- data.frame(metabolite = rownames(m1),
                              pathway = metanno_m1$SUB_PATHWAY[match(rownames(m1),metanno_m1$COMP_IDstr)],
                              p.value = sapply(1:nrow(m1),function(i){wilcox.test(as.numeric(m1[i,])~MetaData_M1$TISSUE_TYPE)$p.value}),
                              dm = sapply(1:nrow(m1),function(i){mean(as.numeric(m1[i,which(MetaData_M1$TISSUE_TYPE=="TUMOR")]))-
                                  mean(as.numeric(m1[i,which(MetaData_M1$TISSUE_TYPE=="NORMAL")]))}))
m1_tumor_normal <- m1_tumor_normal[-grep("X -",metanno_m1$BIOCHEMICAL[match(rownames(m1),metanno_m1$COMP_IDstr)]),] #remove unmatched metabolites
m1_tumor_normal$p.adj <- p.adjust(m1_tumor_normal$p.value, method = "BH")
m1_tumor_normal$sig <- ifelse(m1_tumor_normal$p.adj < 0.05 & abs(m1_tumor_normal$dm) > 0.5, "FDR < 0.05","Not Significant")

m4_tn_factor <- MetaData_M4$TISSUE_STATUS[-which(colnames(m4) %in% m4_outliers_id)]
m4 <- m4[,-which(colnames(m4) %in% m4_outliers_id)]
m4_tumor_normal <- data.frame(metabolite = rownames(m4),
                              pathway = metanno_m4$SUB_PATHWAY[match(rownames(m4),metanno_m4$COMP_IDstr)],
                              p.value = sapply(1:nrow(m4),function(i){wilcox.test(as.numeric(m4[i,])~m4_tn_factor)$p.value}),
                              dm = sapply(1:nrow(m4),function(i){mean(as.numeric(m4[i,which(m4_tn_factor == "TISSUE_TUMOR")]))-
                                  mean(as.numeric(m4[i,which(m4_tn_factor == "TISSUE_NORMAL")]))}))
m4_tumor_normal$p.adj <- p.adjust(m4_tumor_normal$p.value, method = "BH")
m4_tumor_normal$sig <- ifelse(m4_tumor_normal$p.adj < 0.05 & abs(m4_tumor_normal$dm) > 0.5, "FDR < 0.05","Not Significant")

m5_tn_factor <- MetaData_M5$SAMPLE_DESCRIPTION[-which(colnames(m5) %in% m5_outliers_id)]
m5 <- m5[,-which(colnames(m5) %in% m5_outliers_id)]
m5_tumor_normal <- data.frame(metabolite = rownames(m5),
                              pathway = metanno_m5$SUB_PATHWAY[match(rownames(m5),metanno_m5$COMP_IDstr)],
                              p.value = sapply(1:nrow(m5),function(i){wilcox.test(as.numeric(m5[i,])~m5_tn_factor)$p.value}),
                              dm = sapply(1:nrow(m5),function(i){mean(as.numeric(m5[i,which(m5_tn_factor == "TISSUE_TUMOR")]))-
                                  mean(as.numeric(m5[i,which(m5_tn_factor == "TISSUE_NORMAL")]))}))
m5_tumor_normal$p.adj <- p.adjust(m5_tumor_normal$p.value, method = "BH")
m5_tumor_normal$sig <- ifelse(m5_tumor_normal$p.adj < 0.05 & abs(m5_tumor_normal$dm) > 0.5, "FDR < 0.05","Not Significant")

#PCA of tumor vs normal
m4_pca <- prcomp(t(m4))
m4_pca <- m4_pca$x
m4_pca_table <- data.frame(sample = rownames(m4_pca),
                           PC1 = m4_pca[,1],
                           PC2 = m4_pca[,2],
                           type = m4_tn_factor)
ggplot(m4_pca_table, aes(x=PC1,y=PC2,col =type)) + geom_point(size=1) + theme_Publication(base_size = 7) + scale_color_brewer(palette="Pastel1") + ggtitle("M4")

m5_pca <- prcomp(t(m5))
m5_pca <- m5_pca$x
m5_pca_table <- data.frame(sample = rownames(m5_pca),
                           PC1 = m5_pca[,1],
                           PC2 = m5_pca[,2],
                           type = m5_tn_factor)
ggplot(m5_pca_table, aes(x=PC1,y=PC2,col =type)) + geom_point(size=1) + theme_Publication(base_size = 7) + scale_color_brewer(palette="Pastel1") + ggtitle("M5")



#differential abundance plot code
calculate_da <- function(pathways,data){
  da_scores <- as.data.frame(matrix(NA, nrow = length(pathways),ncol = 5))
  colnames(da_scores) <- c("Metabolite_Count","Up","Down","DA_Score","Pathway")
  da_scores$Metabolite_Count <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i],])[1]})
  da_scores$Up <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i] & data$dm > 0,])[1]}) 
  da_scores$Down <- sapply(1:length(pathways), function(i){(da_scores$Metabolite_Count[i]-da_scores$Up[i])})
  da_scores$DA_Score <- sapply(1:length(pathways), function(i){(da_scores$Up[i]-da_scores$Down[i])/da_scores$Metabolite_Count[i]})
  da_scores$Pathway <-  pathways 
  da_scores$Pathway <- factor(da_scores$Pathway, levels = da_scores[order(da_scores$DA_Score),5])
  return(da_scores)
}

plot_da <- function(data,title){
  ggplot(data, aes(x=Pathway, y=DA_Score)) +
    geom_segment(aes(x=Pathway, xend=Pathway, y=0, yend=DA_Score), color=ifelse(data$DA_Score > 0,"red","blue")) +
    geom_point(color=ifelse(data$DA_Score > 0,"red","blue")) + geom_hline(yintercept = 0, col = "grey") + coord_flip() +
    theme_Publication(base_size = 7) +
    xlab("") +
    ylab("Differential Abundance Score") +
    ggtitle(paste(title))
}



#compare fold-changes in m1 vs m4/m5
m4_intersect <- intersect(m1_tumor_normal$metabolite,m4_tumor_normal$metabolite)
m5_intersect <- intersect(m1_tumor_normal$metabolite,m5_tumor_normal$metabolite)

rownames(m1_tumor_normal) <- m1_tumor_normal$metabolite
rownames(m4_tumor_normal) <- m4_tumor_normal$metabolite
rownames(m5_tumor_normal) <- m5_tumor_normal$metabolite

cor.test(m1_tumor_normal[m4_intersect,4],m4_tumor_normal[m4_intersect,4],method="spearman")
qplot(m1_tumor_normal[m4_intersect,4],m4_tumor_normal[m4_intersect,4]) + theme_Publication() + xlab("M1") + ylab("M4") +
  geom_smooth(method=lm)

cor.test(m1_tumor_normal[m5_intersect,4],m5_tumor_normal[m5_intersect,4],method="spearman")
qplot(m1_tumor_normal[m5_intersect,4],m5_tumor_normal[m5_intersect,4]) + theme_Publication() + xlab("M1") + ylab("M5") +
  geom_smooth(method=lm)


m4vsm1 <- data.frame(m1 = m1_tumor_normal[m4_intersect,4],
                     m4 = m4_tumor_normal[m4_intersect,4])
write.table(m4vsm1,file = "/Users/tangc1/Documents/RCC Final/Source/FS1b1.txt", sep = "\t", quote = FALSE)

m5vsm1 <- data.frame(m1 = m1_tumor_normal[m5_intersect,4],
                     m5 = m5_tumor_normal[m5_intersect,4])
write.table(m5vsm1,file = "/Users/tangc1/Documents/RCC Final/Source/FS1b2.txt", sep = "\t", quote = FALSE)


#treated vs untreated
m4_cc_tumors <- MetaData_M4$SAMPLE_NAME[match(master_histo[master_histo$METABOLON_COHORT=="M4" & master_histo$HISTOLOGY == "Clear cell" & master_histo$TUMOR_NORMAL=="TUMOR",3],MetaData_M4$CLIENT_IDENTIFIER)]
m4_tumors <- m4[,which(colnames(m4) %in% m4_cc_tumors)] #less columns bc it removes outliers
m4_treatment_status <- MetaData_M4$TREATMENT[match(colnames(m4_tumors),MetaData_M4$SAMPLE_NAME)]
m4_treatment <- data.frame(p.value = sapply(intersecting_metabolites,function(i){wilcox.test(as.numeric(m4_tumors[i,])~m4_treatment_status)$p.value}))
m4_treatment$p.adj <- p.adjust(m4_treatment$p.value, method = "BH")
m4_treatment$metabolite <- metanno_m4$BIOCHEMICAL[match(intersecting_metabolites,metanno_m4$COMP_IDstr)]
m4_treatment$dm <- sapply(intersecting_metabolites,function(i){mean(as.numeric(m4[i,which(m4_treatment_status=="ANTI-PD1")]))-mean(as.numeric(m4[i,which(m4_treatment_status!="ANTI-PD1")]))})
length(which(m4_treatment$p.adj < 0.05 & abs(m4_treatment$dm) > 0.5)) #4

m5_cc_tumors <- MetaData_M5[MetaData_M5$SAMPLE_DESCRIPTION=="TISSUE_TUMOR" & MetaData_M5$HISTOLOGY=="Clear Cell",3] 
m5_tumors <- m5[,which(colnames(m5) %in% m5_cc_tumors)]
m5_treatment_status <- MetaData_M5$Group[match(colnames(m5_tumors),MetaData_M5$SAMPLE_NAME)]
m5_treatment_status <- m5_treatment_status == "UNTREATED_TUMOR"
m5_treatment <- data.frame(p.value = sapply(intersecting_metabolites,function(i){wilcox.test(as.numeric(m5_tumors[i,])~m5_treatment_status)$p.value}))
m5_treatment$p.adj <- p.adjust(m5_treatment$p.value, method = "BH")
m5_treatment$metabolite <- metanno_m5$BIOCHEMICAL[match(intersecting_metabolites,metanno_m5$COMP_IDstr)]
m5_treatment$dm <- sapply(intersecting_metabolites,function(i){mean(as.numeric(m5[i,which(m5_treatment_status==FALSE)]))-mean(as.numeric(m5[i,which(m5_treatment_status==TRUE)]))})
length(which(m5_treatment$p.adj < 0.05 & abs(m5_treatment$dm) > 0.5)) #64


treatment_effect <- data.frame(p.value = sapply(1:602,function(i){combine.test(c(m4_treatment$p.value[i],m5_treatment$p.value[i]), method = "fisher")}),
                               dm = (m5_treatment$dm + m4_treatment$dm)/2)
treatment_effect$p.adj <- p.adjust(treatment_effect$p.value, method = "BH")
treatment_effect$Sig <- ifelse(treatment_effect$p.adj < 0.05 & abs(treatment_effect$dm) > 0.5, "FDR < 0.05", "Not Significant")
treatment_effect$metabolite <- m5_treatment$metabolite 
treatment_effect$pathway <- metanno_m5$SUB_PATHWAY[match(intersecting_metabolites,metanno_m5$COMP_IDstr)]


m4_tumors_pca <- t(m4_tumors)
m4_tumors_pca <- prcomp(m4_tumors_pca,scale. = T,center = T)
m4_tumors_pca <- data.frame(PC1=m4_tumors_pca$x[,1],PC2=m4_tumors_pca$x[,2],Treatment = m4_treatment_status)
ggplot(m4_tumors_pca, aes(x=PC1,y=PC2,col =Treatment)) + geom_point(size=1) + theme_Publication() + ggtitle("M4") + 
  scale_color_manual(values = c("#E41A1C", "#377EB8"))  

m5_tumors_pca <- t(m5_tumors)
m5_tumors_pca <- prcomp(m5_tumors_pca,scale. = T,center = T)
m5_tumors_pca <- data.frame(PC1=m5_tumors_pca$x[,1],PC2=m5_tumors_pca$x[,2],Treatment = m5_treatment_status)
ggplot(m5_tumors_pca, aes(x=PC1,y=PC2,col =Treatment)) + geom_point(size=1) + theme_Publication() + ggtitle("M5") + 
  scale_color_manual(values = c("#377EB8","#E41A1C" ))  




