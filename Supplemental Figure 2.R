load("All Data.RData")

library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(factoextra)
library(survcomp)


m4_scaled <- t(scale(t(m4)))
colnames(m4_scaled) <- MetaData_M4$CLIENT_IDENTIFIER[match(colnames(m4_scaled),MetaData_M4$SAMPLE_NAME)]

raw_counts_m4 <- read.csv("~/Documents/RCC/Data/raw_counts_summed_m4.csv")
raw_counts_m4$Sample <- gsub("-",".",raw_counts_m4$Sample)
raw_counts_m4$ID <- MetaData_M4$CLIENT_IDENTIFIER[match(raw_counts_m4$Sample,MetaData_M4$SAMPLE_NAME)]
raw_counts_m4 <- raw_counts_m4[match(colnames(m4),raw_counts_m4$Sample),]
colnames(m4)==raw_counts_m4$Sample

min(raw_counts_m4$Count) #6744240355
max(raw_counts_m4$Count) #61729011981

col_fun = colorRamp2(c(6744240355, 61729011981), c("white", "red"))

ha4 = HeatmapAnnotation(
  IonCount = raw_counts_m4$Count, 
  col = list(IonCount = col_fun))

m4_hm <- Heatmap(m4_scaled, column_names_gp = gpar(fontsize = 5),show_row_names = FALSE, top_annotation = ha4, show_column_names = TRUE)
draw(m4_hm)
m4_hm_names <- colnames(m4_scaled)[column_order(m4_hm)]
which(m4_hm_names=="611") #173; start of outliers
m4_outliers <- m4_hm_names[173:184]
m4_outliers_id <- MetaData_M4$SAMPLE_NAME[match(m4_outliers, MetaData_M4$CLIENT_IDENTIFIER)]


m4_pca <- prcomp(t(m4), scale. = T, center = T)
m4_pca <- data.frame(PC1=m4_pca$x[,1],PC2=m4_pca$x[,2])
m4_pca$Sample <- colnames(m4_scaled)
m4_pca$Outlier <- m4_pca$Sample %in% m4_outliers
m4_pca$Type <- MetaData_M4$TISSUE_STATUS
m4_pca$Color <- NA
m4_pca$Color <- ifelse(m4_pca$Outlier == TRUE, "Outlier",m4_pca$Color)
m4_pca$Color <- ifelse(m4_pca$Type == "TISSUE_NORMAL", "Normal",m4_pca$Color)
m4_pca$Color[which(is.na(m4_pca$Color))] <- c("Not Outlier")

ha4_2 = HeatmapAnnotation(
  IonCount = raw_counts_m4$Count,
  Type = m4_pca$Color,
  col = list(Type = c("Normal" = "purple", "Not Outlier" = "green", "Outlier" = "orange"),
             IonCount = col_fun))

Heatmap(m4_scaled, column_names_gp = gpar(fontsize = 5),show_row_names = FALSE, top_annotation = ha4_2, show_column_names = FALSE)

ggplot(m4_pca, aes(x=PC1,y=PC2,col =Outlier)) + geom_point(size=0.15) + scale_color_manual(palette = "Pastel1") + 
   theme_Publication(base_size = 7) + ggtitle("M4")#can see NIVOT17_RC is the lone non-outlier that clusters with the outliers
ggsave("~/Documents/RCC Final/Old Figures/October 19, 2021/Figure 1/Outlier Analysis/M4 Outliers PCA Plot.pdf", width = 2, height = 2.5, units = "in")


ggplot(m4_pca, aes(x=PC1,y=PC2,col =Type)) + geom_point(size=0.15) + scale_color_brewer(palette = "Pastel1") + 
  theme_Publication(base_size = 7) + ggtitle("M4")
ggsave("~/Documents/RCC Final/Old Figures/October 19, 2021/Supplemental Figures/M4 tumor vs normal PCA Plot.pdf", width = 2.5, height = 2.5, units = "in")


ggplot(m4_pca, aes(x=PC1,y=PC2,col =Color)) + geom_point(size=0.15) + scale_color_manual(values = c("#CCEBC5","#FBB4AE","#B3CDE3")) + 
  theme_Publication(base_size = 7) 
ggsave("~/Documents/RCC Final/M4 Outliers PCA Plot.pdf", width = 2, height = 2, units = "in")



m5_scaled <- t(scale(t(m5)))
colnames(m5_scaled) <- MetaData_M5$CLIENT_IDENTIFIER[match(colnames(m5_scaled),MetaData_M5$SAMPLE_NAME)]


raw_counts_m5 <- read.csv("~/Documents/RCC/Data/raw_counts_summed_m5.csv")
raw_counts_m5$Sample <- gsub("-",".",raw_counts_m5$Sample)
raw_counts_m5$ID <- MetaData_M5$CLIENT_IDENTIFIER[match(raw_counts_m5$Sample,MetaData_M5$SAMPLE_NAME)]
raw_counts_m5 <- raw_counts_m5[match(colnames(m5),raw_counts_m5$Sample),]
colnames(m5)==raw_counts_m5$Sample

min(raw_counts_m5$Count) #6248401428
max(raw_counts_m5$Count) #58436413245

col_fun2 = colorRamp2(c(6248401428, 58436413245), c("white", "red"))

ha5 = HeatmapAnnotation(
  IonCount = raw_counts_m5$Count, 
  col = list(IonCount = col_fun2))

m5_hm <- Heatmap(m5_scaled, column_names_gp = gpar(fontsize = 5),show_row_names = FALSE, top_annotation = ha5, show_column_names = TRUE)
draw(m5_hm)
m5_hm_names <- colnames(m5_scaled)[column_order(m5_hm)]
which(m5_hm_names=="AH002916") #156; start of outliers
m5_outliers <- m5_hm_names[156:162]
m5_outliers_id <- MetaData_M5$SAMPLE_NAME[match(m5_outliers, MetaData_M5$CLIENT_IDENTIFIER)]
m5_pca <- prcomp(t(m5), scale. = T, center = T)
m5_pca <- data.frame(PC1=m5_pca$x[,1],PC2=m5_pca$x[,2])
m5_pca$Sample <- colnames(m5_scaled)
m5_pca$Outlier <- m5_pca$Sample %in% m5_outliers
m5_pca$Type <- MetaData_M5$SAMPLE_DESCRIPTION

ggplot(m5_pca, aes(x=PC1,y=PC2,col =Outlier)) + geom_point(size=0.15) + scale_color_brewer(palette = "Pastel1") + 
  theme_Publication(base_size = 7) + ggtitle("M5")
ggsave("~/Documents/RCC Final/Old Figures/October 19, 2021/Figure 1/Outlier Analysis/M5 Outliers PCA Plot.pdf", width = 1.75, height = 2, units = "in")

ggplot(m5_pca, aes(x=PC1,y=PC2,col =Type)) + geom_point(size=0.15) + scale_color_brewer(palette = "Pastel1") + 
  theme_Publication(base_size = 7) + ggtitle("M5")
ggsave("~/Documents/RCC Final/Old Figures/October 19, 2021/Supplemental Figures/M5 tumor vs normal PCA Plot.pdf", width = 2.5, height = 2.5, units = "in")


save(m4_outliers,m4_outliers_id,m5_outliers,m5_outliers_id, file = "Outlier Data.Rdata")



#compare raw ion count of outliers to non-outliers
outlier_raw_p4 <- wilcox.test(as.numeric(raw_counts_m4[raw_counts_m4$ID %in% m4_outliers,1]),as.numeric(raw_counts_m4[!raw_counts_m4$ID %in% m4_outliers,1]))
outlier_raw_p5 <- wilcox.test(as.numeric(raw_counts_m5[raw_counts_m5$ID %in% m5_outliers,1]),as.numeric(raw_counts_m5[!raw_counts_m5$ID %in% m5_outliers,1]))
combine.test(c(outlier_raw_p4$p.value,outlier_raw_p5$p.value), method = "fisher")


raw_counts_m4$outlier <- raw_counts_m4$ID %in% m4_outliers
raw_counts_m5$outlier <- raw_counts_m5$ID %in% m5_outliers
raw_counts_bp <- data.frame(raw_count = c(raw_counts_m4$Count,raw_counts_m5$Count),
                            outlier = c(raw_counts_m4$outlier,raw_counts_m5$outlier),
                            batch = c(rep("M4",nrow(raw_counts_m4)),rep("M5",nrow(raw_counts_m5))))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(raw_counts_bp, varname="raw_count", 
                    groupnames=c("outlier", "batch"))
df2$batch=as.factor(df2$batch)

ggplot(df2, aes(x=batch, y=raw_count, fill=outlier)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=raw_count-sd, ymax=raw_count+sd), width=.2,
                position=position_dodge(.9)) + theme_Publication() + scale_fill_brewer(palette = "Pastel1") +
  xlab("Batch") + ylab("Raw Ion Count")
ggsave("~/Documents/RCC Final/Outlier Raw Counts Barplot.pdf", width = 2.5, height = 2.25, units = "in")


#volcano plot of outliers vs non-outliers
outlier_vp <- data.frame(m4p = sapply(intersecting_metabolites,function(i){wilcox.test(as.numeric(m4[i,])~m4_pca$Outlier)$p.value}),
                         m5p = sapply(intersecting_metabolites,function(i){wilcox.test(as.numeric(m5[i,])~m5_pca$Outlier)$p.value}),
                         m4dm = sapply(intersecting_metabolites,function(i){mean(as.numeric(m4[i,which(m4_pca$Outlier)]))-mean(as.numeric(m4[i,which(!m4_pca$Outlier)]))}),
                         m5dm = sapply(intersecting_metabolites,function(i){mean(as.numeric(m5[i,which(m5_pca$Outlier)]))-mean(as.numeric(m5[i,which(!m5_pca$Outlier)]))}),
                         pathway = metanno_m4$SUPER_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                         sub =  metanno_m4$SUB_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)])
outlier_vp$fisherp <- sapply(1:602,function(i){combine.test(outlier_vp[i,1:2], method = "fisher")})
outlier_vp$fisherq <- p.adjust(outlier_vp$fisherp)
outlier_vp$dm <- (outlier_vp$m4dm+outlier_vp$m5dm)/2
outlier_vp$Significant <- ifelse(outlier_vp$fisherq < 0.05 & abs(outlier_vp$dm)> 0.5, outlier_vp$pathway,"Not Significant")
outlier_vp$m4q <- p.adjust(outlier_vp$m4p,method = "BH")
outlier_vp$m5q <- p.adjust(outlier_vp$m5p,method = "BH")
outlier_vp$metabolite <- metanno_m4$BIOCHEMICAL[match(rownames(outlier_vp),metanno_m4$COMP_IDstr)]

write.table(outlier_vp,"Outlier Fold Changes.txt",sep="\t",quote = FALSE,row.names = FALSE)
outlier_vp <- read.delim("Outlier Fold Changes.txt")

ggplot(outlier_vp, aes(x =dm, y = -log10(fisherp)))  +
  geom_point(aes(color = Significant),size = 2) +
  scale_color_manual(values = c("#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072", "#80B1D3", "#D9D9D9", "#FDB462", "#B3DE69","#FCCDE5")) +
  ggtitle("Outliers") + xlab("Outlier - Non-Outlier") +
  theme_Publication(base_size = 7) 
ggsave("~/Documents/RCC Final/Old Figures/October 19, 2021/Figure 1/Outlier Analysis/Outlier Volcano Plo.pdf", width = 2.75, height = 4, units = "in")


outlier_means <- data.frame(outlier = c(sapply(intersecting_metabolites,function(i){mean(as.numeric(m4[i,which(m4_pca$Outlier)]))}),
                                                 sapply(intersecting_metabolites,function(i){mean(as.numeric(m5[i,which(m5_pca$Outlier)]))})),
                                     other = c(sapply(intersecting_metabolites,function(i){mean(as.numeric(m4[i,which(!m4_pca$Outlier)]))}),
                                               sapply(intersecting_metabolites,function(i){mean(as.numeric(m5[i,which(!m5_pca$Outlier)]))})),
                                     batch = c(rep("M4",602),rep("M5",602)),
                            pathway = rep(metanno_m4$SUB_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],2))

outlier_means <- melt(outlier_means)

ggplot(subset(outlier_means,pathway == "Sphingomyelins"), aes(x=variable,y=value,fill=batch)) + geom_boxplot() + theme_Publication() +
  xlab("Outlier Status") + ylab("Spingomyelin Expression") + scale_fill_brewer(palette = "Pastel1")
ggsave("~/Documents/RCC Final/Sphingomyelins Outlier Raw Counts Barplot.pdf", width = 2.5, height = 2.25, units = "in")

ggplot(subset(outlier_means,pathway == "Ceramides"), aes(x=variable,y=value,fill=batch)) + geom_boxplot() + theme_Publication() +
  xlab("Outlier Status") + ylab("Ceramide Expression") + scale_fill_brewer(palette = "Pastel1")
ggsave("~/Documents/RCC Final/Ceramides Outlier Raw Counts Barplot.pdf", width = 2, height = 2.25, units = "in")



#outlier pathway analysis
calculate_da <- function(pathways,data){
  da_scores <- as.data.frame(matrix(NA, nrow = length(pathways),ncol = 5))
  colnames(da_scores) <- c("Metabolite_Count","Up","Down","DA_Score","Pathway")
  da_scores$Metabolite_Count <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i],])[1]})
  da_scores$Up <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i] & data$dm > 0 & data$Significant != "Not Significant",])[1]}) 
  da_scores$Down <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i] & data$dm < 0 & data$Significant != "Not Significant",])[1]}) 
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
outlier_pathways <- names(which(sort(table(subset(outlier_vp,Significant != "Not Significant")[,6]),decreasing = TRUE) > 2)) #only keep pathways with 5+ sig altered metabolites
outlier_pathways_table <- outlier_vp[outlier_vp$sub %in% outlier_pathways,]
outlier_pathways_table$sub_pathway <- as.factor(outlier_pathways_table$sub_pathway)
outlier_da_scores <- calculate_da(outlier_pathways,outlier_pathways_table)
plot_da(outlier_da_scores,"Outliers")
ggsave("~/Documents/RCC Final/Outliers Differential Abundance Pathway Analysis.pdf", width = 4, height = 4, units = "in")

write.table(outlier_da_scores,"Final Tables/Outlier DA Scores.txt",sep="\t",quote = FALSE,row.names = TRUE)





#spingomyelin expression
sphingo <- rownames(outlier_vp[outlier_vp$sub=="Sphingomyelins",])
sapply(sphingo,function(i){as.numeric(m4[i,which(m4_pca$Outlier)])})
sapply(sphingo,function(i){as.numeric(m4[i,which(!m4_pca$Outlier)])})

m4_sphingo <- rbind(melt(sapply(sphingo,function(i){as.numeric(m4[i,which(m4_pca$Outlier)])})),
                    melt(sapply(sphingo,function(i){as.numeric(m4[i,which(!m4_pca$Outlier)])})))
m4_sphingo$outlier <- c(rep("TRUE",312),rep("FALSE",4472))

df3 <- data_summary(m4_sphingo, varname="value", 
                    groupnames=c("outlier", "Var2"))
df3$batch=as.factor(df3$Var2)
df3$Var3 <- rep(sapply(1:26,function(i){df3[df3$Var2 == unique(df3$Var2)[i],3][1]-10}),2) #outlier value -1 so that everything will be normalized to 10
df3$Var4 <- df3$value-df3$Var3
df3$Var2 <- factor(df3$Var2, levels = subset(df3[order(df3$Var4),],outlier=="TRUE")[,2])

ggplot(df3, aes(x=Var2, y=value, fill=outlier)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) + theme_Publication() + scale_fill_brewer(palette = "Pastel1") +
  xlab("Batch") + ylab("Sphingomyelin Expression") + ggtitle("M4") + theme(axis.text.x = element_blank())
ggsave("~/Documents/RCC Final/M4 Outliers Sphingomyelin Barplot.pdf", width = 4, height = 2.25, units = "in")


ggplot(df3, aes(x=Var2, y=Var4, fill=outlier)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Var4-sd, ymax=Var4+sd), width=.2,
                position=position_dodge(.9)) + theme_Publication() + scale_fill_brewer(palette = "Pastel1") +
  xlab("Batch") + ylab("Sphingomyelin Expression") + ggtitle("M4") + theme(axis.text.x = element_blank())
ggsave("~/Documents/RCC Final/M4 Outliers Sphingomyelin Barplot Normalized.pdf", width = 4, height = 2.25, units = "in")


m5_sphingo <- rbind(melt(sapply(sphingo,function(i){as.numeric(m5[i,which(m5_pca$Outlier)])})),
                    melt(sapply(sphingo,function(i){as.numeric(m5[i,which(!m5_pca$Outlier)])})))
m5_sphingo$outlier <- c(rep("TRUE",182),rep("FALSE",4030))

df4 <- data_summary(m5_sphingo, varname="value", 
                    groupnames=c("outlier", "Var2"))
df4$batch=as.factor(df4$Var2)
df4$Var3 <- rep(sapply(1:26,function(i){df4[df4$Var2 == unique(df4$Var2)[i],3][1]-10}),2) #outlier value -1 so that everything will be normalized to 10
df4$Var4 <- df4$value-df4$Var3
df4$Var2 <- factor(df4$Var2, levels = subset(df3[order(df3$Var4),],outlier=="TRUE")[,2])


ggplot(df4, aes(x=Var2, y=value, fill=outlier)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) + theme_Publication() + scale_fill_brewer(palette = "Pastel1") +
  xlab("Batch") + ylab("Sphingomyelin Expression") + ggtitle("M5") + theme(axis.text.x = element_blank())
ggsave("~/Documents/RCC Final/M5 Outliers Sphingomyelin Barplot.pdf", width = 4, height = 2.25, units = "in")

ggplot(df4, aes(x=Var2, y=Var4, fill=outlier)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Var4-sd, ymax=Var4+sd), width=.2,
                position=position_dodge(.9)) + theme_Publication() + scale_fill_brewer(palette = "Pastel1") +
  xlab("Batch") + ylab("Sphingomyelin Expression") + ggtitle("M5") + theme(axis.text.x = element_blank())
ggsave("~/Documents/RCC Final/M5 Outliers Sphingomyelin Barplot Normalized.pdf", width = 4, height = 2.25, units = "in")





