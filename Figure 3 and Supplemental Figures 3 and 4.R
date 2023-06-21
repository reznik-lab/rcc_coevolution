load("All Data.RData")
load("Outlier Data.RData")

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(survcomp)
library(DESeq2)
library(qusage)
library(fgsea)
library(data.table)
library(DescTools)
library(org.Hs.eg.db)
library(annotate)
library(patchwork)

#####RUG PLOT#####
#clean data and set up necessary data; outliers need to be calculated from outlier analysis
m4_mr_tumors <- m4[,MetaData_M4$TISSUE_STATUS[match(colnames(m4),MetaData_M4$SAMPLE_NAME)]=="TISSUE_TUMOR"]
m4_mr_tumors <- m4_mr_tumors[ , !(names(m4_mr_tumors) %in% m4_outliers_id)] #remove outliers
which(MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]=="MR06") #75, this is the last MR sample; rest are single region tumors
m4_mr_tumors <- m4_mr_tumors[,1:75]
which(MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]=="MR04") #66:70 need to remove as it's hlrcc
m4_mr_tumors <- m4_mr_tumors[,-c(66:70)]
m4_mr_tumors <- m4_mr_tumors[intersecting_metabolites,]
m4_mr_tumors_subjectid <- MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]
m4_mr_patients <- unique(m4_mr_tumors_subjectid)
m4_mr_patient_count <- sapply(m4_mr_patients,function(x)sum(m4_mr_tumors_subjectid==x))
m4_mr_patient_sum <- cumsum(m4_mr_patient_count)


m5_mr_tumors <- m5[,MetaData_M5$GROUP[match(colnames(m5),MetaData_M5$SAMPLE_NAME)]=="TUMOR"]
m5_mr_tumors <- m5_mr_tumors[ , !(names(m5_mr_tumors) %in% m5_outliers_id)] #remove outliers
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)]=="SC16") #58, this is the last MR sample; rest are single region tumors
m5_mr_tumors <- m5_mr_tumors[,1:58] #note here that SC01 and SC04 are only 1 region each so they need to be removed 
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)] %in% c("SC01","SC04")) #columns 30 and 34 need to be removed
m5_mr_tumors <- m5_mr_tumors[,-c(30,34)]
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)]=="SC07") #need to remove as sc07 is unclassified
m5_mr_tumors <- m5_mr_tumors[,-c(37:39)]
m5_mr_tumors <- m5_mr_tumors[intersecting_metabolites,]
m5_mr_tumors_subjectid <- MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)]
m5_mr_patients <- unique(m5_mr_tumors_subjectid)
m5_mr_patient_count <- sapply(m5_mr_patients,function(x)sum(m5_mr_tumors_subjectid==x))
m5_mr_patient_sum <- cumsum(m5_mr_patient_count)

median_center <- function(data){
  medians <- apply(data,1,median)
  data <- data - medians
  return(data)
}

m4_mr_tumors_mc <- median_center(m4_mr_tumors)
m5_mr_tumors_mc <- median_center(m5_mr_tumors)

intraheterogeneity<- function(data, npatients, n_p1, sample_sum){ #npatients = number of patients, n_p1 = number of samples in patient 1, #sample sum is vector of cumulative sample count
  intra <- as.data.frame(matrix(NA,ncol=nrow(data),nrow=npatients))
  for(x in 1:nrow(data)){
    distance_matrix <- rdist(as.numeric(data[x,]))
    distance_matrix[lower.tri(distance_matrix, diag = TRUE)] <- NA
    
    pat1 <- median(as.numeric(distance_matrix[1:n_p1,1:n_p1]),na.rm=TRUE)
    rest <- sapply(2:npatients,function(i){median(as.numeric(distance_matrix[as.numeric(sample_sum[i-1]+1):as.numeric(sample_sum[i]),as.numeric(sample_sum[i-1]+1):as.numeric(sample_sum[i])]),na.rm=TRUE)})
    intra[,x] <- c(pat1,rest)}
  return(intra)
}



intra_m4_met <- intraheterogeneity(m4_mr_tumors_mc, 16, 5, m4_mr_patient_sum)
rownames(intra_m4_met) <- as.character(m4_mr_patients)
colnames(intra_m4_met) <- rownames(m4_mr_tumors_mc)

intra_m5_met <- intraheterogeneity(m5_mr_tumors_mc, 14, 5, m5_mr_patient_sum)
rownames(intra_m5_met) <- as.character(m5_mr_patients)
colnames(intra_m5_met) <- rownames(m5_mr_tumors_mc)

m4_heterogeneity <- data.frame(Heterogeneity = c(as.numeric(apply(intra_m4_met,1,median)),as.numeric(apply(inter_m4_met,1,median))),
                               Inter = c(rep("Intratumor",nrow(intra_m4_met)),rep("Intertumor",nrow(inter_m4_met))),
                               Sample = c(rownames(intra_m4_met),rownames(inter_m4_met)))
m5_heterogeneity <- data.frame(Heterogeneity = c(as.numeric(apply(intra_m5_met,1,median)),as.numeric(apply(inter_m5_met,1,median))),
                               Inter = c(rep("Intratumor",nrow(intra_m5_met)),rep("Intertumor",nrow(inter_m5_met))),
                               Sample = c(rownames(intra_m5_met),rownames(inter_m5_met)))


fig3a1 <- qplot(m4_heterogeneity,m5_heterogeneity, data=metabolite_intra_tumor_het, color = col) +theme_Publication(base_size = 7) + xlab("M4 Heterogeneity") + ylab("M5 Heterogeneity") + 
  ylim(0,3) + xlim(0,3) + theme(legend.position = "none") +geom_abline(slope=1,intercept = 0) + scale_color_manual(values=c("#D9D9D9","#FF0000")) + 
  geom_point(data = subset(metabolite_intra_tumor_het,metabolite %in% imp_mets_scatterplot), aes(x=m4_heterogeneity,y=m5_heterogeneity), col = "red") +
  geom_text_repel(data=subset(metabolite_intra_tumor_het, metabolite %in% imp_mets_scatterplot), aes(label=metabolite),size=2,color="black")
ggsave("M4 vs M5 Heterogeneity Scatterplot.pdf", width = 2.5, height = 2.5, units = "in")

fig3a2 <- ggplot(metabolite_intra_tumor_het,aes(x=m4_heterogeneity,y=col2)) + geom_point(shape=124,size=3,color=metabolite_intra_tumor_het$col3) + 
  theme_Publication(base_size = 7) + ylab("") + xlab("") + theme(axis.ticks = element_blank(),axis.text = element_blank()) + xlim(0,3)
ggsave("M4 vs M5 Heterogeneity Rug plot.pdf", width = 2.5, height = 1.5, units = "in")


(fig3a1/fig3a2) + plot_layout(heights = c(3, 1))


#####CLEAN DATA FOR PCA ANALYSIS#####
m4_mr_tumors <- m4[,MetaData_M4$TISSUE_STATUS[match(colnames(m4),MetaData_M4$SAMPLE_NAME)]=="TISSUE_TUMOR"]
m4_mr_tumors <- m4_mr_tumors[ , !(names(m4_mr_tumors) %in% m4_outliers_id)] #remove outliers
which(MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]=="MR06") #75, this is the last MR sample; rest are single region tumors
m4_mr_tumors <- m4_mr_tumors[,1:75]
which(MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]=="MR04") #66:70 need to remove as it's hlrcc
m4_mr_tumors <- m4_mr_tumors[,-c(66:70)]

m5_mr_tumors <- m5[,MetaData_M5$GROUP[match(colnames(m5),MetaData_M5$SAMPLE_NAME)]=="TUMOR"]
m5_mr_tumors <- m5_mr_tumors[ , !(names(m5_mr_tumors) %in% m5_outliers_id)] #remove outliers
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)]=="SC16") #58, this is the last MR sample; rest are single region tumors
m5_mr_tumors <- m5_mr_tumors[,1:58] #note here that SC01 and SC04 are only 1 region each so they need to be removed 
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)] %in% c("SC01","SC04")) #columns 30 and 34 need to be removed
m5_mr_tumors <- m5_mr_tumors[,-c(30,34)]
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)]=="SC07") #need to remove as sc07 is unclassified
m5_mr_tumors <- m5_mr_tumors[,-c(37:39)]

m4_mr_normals <- m4[,MetaData_M4$TISSUE_STATUS[match(colnames(m4),MetaData_M4$SAMPLE_NAME)]=="TISSUE_NORMAL"] #all normals are from MR samples but MR06 only has 1 region so remove it
m4_mr_normals <- m4_mr_normals[,-50]
which(MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_normals),MetaData_M4$SAMPLE_NAME)]=="MR04") #44,45,46
m4_mr_normals <- m4_mr_normals[,-c(44:46)]

m5_mr_normals <- m5[,MetaData_M5$GROUP[match(colnames(m5),MetaData_M5$SAMPLE_NAME)]=="NORMAL"] #all are MR but all the SC patients only have 1 region
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_normals),MetaData_M5$SAMPLE_NAME)]=="SC03") #19 remove 19-28
m5_mr_normals <- m5_mr_normals[,-c(19:28)]

m4_mr_all <- cbind(m4_mr_tumors,m4_mr_normals)
m4_mr_all <- m4[,sort(colnames(m4_mr_all))] #so all samples from the same patient are together
m4_mr_all_ids <- MetaData_M4$SUBJECT_ID[match(colnames(m4[,sort(colnames(m4_mr_all))]),MetaData_M4$SAMPLE_NAME)]
m4_mr_patients <- unique(m4_mr_all_ids)
m4_mr_all_sample_count <- sapply(m4_mr_patients,function(x)sum(m4_mr_all_ids==x))
m4_mr_all_sample_sum <- cumsum(m4_mr_all_sample_count)


m4_patient_centered_all <- as.data.frame(matrix(NA,nrow=nrow(m4_mr_all),ncol = ncol(m4_mr_all)-m4_mr_all_sample_count[1]))
for(x in 1:nrow(m4_mr_all)){
  m4_patient_centered_all[x,] <- unlist(sapply(2:length(m4_mr_patients),function(i){m4_mr_all[x,as.numeric(m4_mr_all_sample_sum[i-1]+1):m4_mr_all_sample_sum[i]]-median(as.numeric(m4_mr_all[x,as.numeric(m4_mr_all_sample_sum[i-1]+1):m4_mr_all_sample_sum[i]]))}))
}
m4_patient_centered_all_p1 <- as.data.frame(matrix(NA,nrow=nrow(m4_mr_all),ncol = m4_mr_all_sample_sum[1]))
for(x in 1:nrow(m4_mr_all)){
  m4_patient_centered_all_p1[x,] <- m4_mr_all[x,1:m4_mr_all_sample_sum[1]]-median(as.numeric(m4_mr_all[x,1:m4_mr_all_sample_sum[1]]))
}
m4_patient_centered_all <- cbind(m4_patient_centered_all_p1,m4_patient_centered_all)
colnames(m4_patient_centered_all) <- colnames(m4_mr_all)
rownames(m4_patient_centered_all) <- rownames(m4_mr_all)
rm(m4_patient_centered_all_p1)



m5_mr_all <- cbind(m5_mr_tumors,m5_mr_normals)
m5_mr_all <- m5[,sort(colnames(m5_mr_all))] #so all samples from the same patient are together
m5_mr_all_ids <- MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_all),MetaData_M5$SAMPLE_NAME)]
m5_mr_patients <- unique(m5_mr_all_ids)
m5_mr_all_sample_count <- sapply(m5_mr_patients,function(x)sum(m5_mr_all_ids==x))
m5_mr_all_sample_sum <- cumsum(m5_mr_all_sample_count)


m5_patient_centered_all <- as.data.frame(matrix(NA,nrow=nrow(m5_mr_all),ncol = ncol(m5_mr_all)-m5_mr_all_sample_count[1]))
for(x in 1:nrow(m5_mr_all)){
  m5_patient_centered_all[x,] <- unlist(sapply(2:length(m5_mr_patients),function(i){m5_mr_all[x,as.numeric(m5_mr_all_sample_sum[i-1]+1):m5_mr_all_sample_sum[i]]-median(as.numeric(m5_mr_all[x,as.numeric(m5_mr_all_sample_sum[i-1]+1):m5_mr_all_sample_sum[i]]))}))
}
m5_patient_centered_all_p1 <- as.data.frame(matrix(NA,nrow=nrow(m5_mr_all),ncol = m5_mr_all_sample_sum[1]))
for(x in 1:nrow(m5_mr_all)){
  m5_patient_centered_all_p1[x,] <- m5_mr_all[x,1:m5_mr_all_sample_sum[1]]-median(as.numeric(m5_mr_all[x,1:m5_mr_all_sample_sum[1]]))
}
m5_patient_centered_all <- cbind(m5_patient_centered_all_p1,m5_patient_centered_all)
colnames(m5_patient_centered_all) <- colnames(m5_mr_all)
rownames(m5_patient_centered_all) <- rownames(m5_mr_all)
rm(m5_patient_centered_all_p1)


#separate just the tumors 
m4_mr_tumors_c <- m4_patient_centered_all[,colnames(m4_mr_tumors)]
m5_mr_tumors_c <- m5_patient_centered_all[,colnames(m5_mr_tumors)]

#separate just the normals 
m4_mr_normals_c <- m4_patient_centered_all[,colnames(m4_mr_normals)]
m5_mr_normals_c <- m5_patient_centered_all[,colnames(m5_mr_normals)]

#####CITRATE HIGH VS LOW FIG 4A#####
#show that citrate high and low have the same features when controlling for patient
mr_all <- cbind(m4_mr_all[intersecting_metabolites,],m5_mr_all[intersecting_metabolites,])
patient_centered_all <- cbind(m4_patient_centered_all[intersecting_metabolites,],m5_patient_centered_all[intersecting_metabolites,])

mr_all_melt <- melt(mr_all["M01564",])
mr_all_melt$patient <- c(m4_mr_all_ids,m5_mr_all_ids)

patient_centered_all_melt <- melt(patient_centered_all["M01564",])
patient_centered_all_melt$patient <- c(m4_mr_all_ids,m5_mr_all_ids)
mr_all_ids <- unique(patient_centered_all_melt$patient)

patient_centered_all_melt$level <- ifelse(patient_centered_all_melt$value > 0,"High","Low") #citrate high vs low
citrate_high <- patient_centered_all_melt[patient_centered_all_melt$level=="High",1]
citrate_low <- patient_centered_all_melt[patient_centered_all_melt$level=="Low",1]

citrate_high_low <- data.frame(ca = as.numeric(c(patient_centered_all["M46173",citrate_high],patient_centered_all["M46173",citrate_low])),
                               gl = as.numeric(c(patient_centered_all["M48152",citrate_high],patient_centered_all["M48152",citrate_low])),
                               fu = as.numeric(c(patient_centered_all["M01643",citrate_high],patient_centered_all["M01643",citrate_low])),
                               lvl = c(rep("high",86),rep("low",101)))
ggplot(citrate_high_low,aes(x=lvl,y=ca)) + geom_boxplot(outlier.size = 0.1) + theme_Publication() + 
  xlab("Citrate Level") + ylab("Expression")
wilcox.test(citrate_high_low$ca~citrate_high_low$lvl) #p-value < 2.2e-16


ggplot(citrate_high_low,aes(x=lvl,y=fu)) + geom_boxplot(outlier.size = 0.1) + theme_Publication() + 
  xlab("Citrate Level") + ylab("Expression")
wilcox.test(citrate_high_low$fu~citrate_high_low$lvl) #p-value = 6.371e-09

#####PCA ANALYSIS#####
#check if the loadings are similar in M4 and M5 too see if you can PCA with both batches
cor.test(prcomp(t(m4_patient_centered_all[intersecting_metabolites,]), center = TRUE, scale. = TRUE)$rotation[,2],
         prcomp(t(m5_patient_centered_all[intersecting_metabolites,]), center = TRUE, scale. = TRUE)$rotation[,2], method = "spearman") #rho = 0.4980312 , p-value < 2.2e-16

pca <- prcomp(t(cbind(m4_patient_centered_all[intersecting_metabolites,],m5_patient_centered_all[intersecting_metabolites,])), scale. = T, center = T)
pca_x <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
pca_x$SampleType <- c(MetaData_M4$TISSUE_STATUS[match(colnames(m4_mr_all),MetaData_M4$SAMPLE_NAME)], 
                      MetaData_M5$SAMPLE_DESCRIPTION[match(colnames(m5_mr_all),MetaData_M5$SAMPLE_NAME)])
pca_x$SubjectID <- c(MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_all),MetaData_M4$SAMPLE_NAME)],
                     MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_all),MetaData_M5$SAMPLE_NAME)])

ggplot(pca_x, aes(x=PC1,y=PC2,col =SampleType)) + geom_point(size=0.2) + theme_Publication(base_size = 7) + scale_color_manual(values=c("#B3CDE3","#FBB4AE")) + geom_hline(yintercept = median(pca_x$PC2), linetype="dashed", color="#737272")
ggsave("PCA Scatterplot.pdf", width = 1.25, height = 1.75, units = "in")

ggplot(subset(pca_x, SubjectID == "NIVOT09"), aes(x=PC1,y=PC2,col =SampleType)) + geom_point(size=1) + theme_Publication(base_size = 7) +  scale_color_manual(values=c("#B3CDE3","#FBB4AE")) + ggtitle("NIVO09") + geom_hline(yintercept = median(pca_x$PC2), linetype="dashed", color="#737272")
ggplot(subset(pca_x, SubjectID == "NIVO24"), aes(x=PC1,y=PC2,col =SampleType)) + geom_point(size=1) + theme_Publication(base_size = 7) + scale_color_brewer(palette="Pastel1") + ggtitle("NIVO24") + geom_hline(yintercept = median(pca_x$PC2), linetype="dashed", color="#737272")


#####ABSOLUTE LOADINGS#####
combined_pca_pc2 <- data.frame(PC2 = pca$rotation[,2],
                               PC2_absolute = abs(pca$rotation[,2]))
combined_pca_pc2$metabolite <- metanno_m4$BIOCHEMICAL[match(rownames(combined_pca_pc2),metanno_m4$COMP_IDstr)]
combined_pca_pc2$super <- metanno_m4$SUPER_PATHWAY[match(rownames(combined_pca_pc2),metanno_m4$COMP_IDstr)]
combined_pca_pc2$sub <- metanno_m4$SUB_PATHWAY[match(rownames(combined_pca_pc2),metanno_m4$COMP_IDstr)]
combined_pca_pc2$rank <- rank(combined_pca_pc2$PC2)
combined_pca_pc2$rank_a <- rank(combined_pca_pc2$PC2_absolute)

#####PCA ON JUST TUMORS#####
pca_x_tumors <- subset(pca_x,SampleType == "TISSUE_TUMOR")
subjects<- unique(pca_x_tumors$SubjectID)
pc2_level <- unlist(sapply(subjects, function(i){pca_x_tumors[pca_x_tumors$SubjectID==i,2]  >= median(pca_x_tumors[pca_x_tumors$SubjectID==i,2])}))
pca_x_tumors$PC2_level <- as.character(ifelse(pc2_level == TRUE, "High","Low"))
pca_x_tumors$ID <- c(MetaData_M4$CLIENT_IDENTIFIER[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)], 
                     MetaData_M5$CLIENT_IDENTIFIER[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)])
pca_x_tumors$purity <- facets_log$Purity[match(pca_x_tumors$ID,facets_log$Sample)]
pca_x_tumors$Batch <- c(rep("M4",ncol(m4_mr_tumors)),rep("M5",ncol(m5_mr_tumors)))
pca_x_tumors$Cysteine <- c(as.numeric(m4_mr_tumors["M01868",]),as.numeric(m5_mr_tumors["M01868",]))
pca_x_tumors$GSSG <- c(as.numeric(m4_mr_tumors["M27727",]),as.numeric(m5_mr_tumors["M27727",]))
pca_x_tumors$GSH <- c(as.numeric(m4_mr_tumors["M02127",]),as.numeric(m5_mr_tumors["M02127",]))
pca_x_tumors$PUFA <- c(as.numeric(m4_mr_tumors["M32417",]),as.numeric(m5_mr_tumors["M32417",]))
pca_x_tumors$Cystine <- c(as.numeric(m4_mr_tumors["M00056",]),as.numeric(m5_mr_tumors["M00056",]))

#####PLOT PC2 HIGH VS LOW METABOLITE LEVELS
cysteine_pc <- data.frame(exp = c(as.numeric(m4_mr_tumors_c["M01868",]),
                                  as.numeric(m5_mr_tumors_c["M01868",])),
                          level = as.factor(pca_x_tumors$PC2_level))
ggplot(cysteine_pc, aes(x=level,y=exp)) + geom_boxplot(outlier.size=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), size = 0.6) +
  theme_Publication() + xlab("PC2 Level") + ylab("Cysteine Abundance") 

cystine_pc <- data.frame(exp = c(as.numeric(m4_mr_tumors_c["M00056",]),
                                 as.numeric(m5_mr_tumors_c["M00056",])),
                         level = as.factor(pca_x_tumors$PC2_level))
ggplot(cystine_pc, aes(x=level,y=exp)) + geom_boxplot(outlier.size = 0.3) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.6) +
  theme_Publication() + xlab("PC2 Level") + ylab("Cystine Abundance")



#####CHECK PURITY EFFECT#####
wilcox.test(as.numeric(pca_x_tumors$purity)~pca_x_tumors$PC2_level) #p.value = 0.8
ggplot(pca_x_tumors,aes(x=PC2_level,y=purity)) + geom_boxplot(outlier.size = 0.3) + theme_Publication() + xlab("PC2 Level") + ylab("Purity") + 
  geom_jitter(position=position_jitter(0.2), size = 0.6)


#####COMPARE PC2 HIGH VS LOW TUMOR REGIONS#####
pc2_pvalues_all <- data.frame(id = intersecting_metabolites,
                              metabolite = metanno_m4$BIOCHEMICAL[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                              pathway = metanno_m4$SUB_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                              m4_p = sapply(1:602,function(i){wilcox.test(as.numeric(m4_mr_tumors_c[intersecting_metabolites[i],])~pca_x_tumors$PC2_level[1:70])$p.value}),
                              m5_p = sapply(1:602,function(i){wilcox.test(as.numeric(m5_mr_tumors_c[intersecting_metabolites[i],])~pca_x_tumors$PC2_level[71:123])$p.value}),
                              m4_dm = sapply(1:602,function(i){mean(as.numeric(m4_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[1:70]=="High")]))-
                                  mean(as.numeric(m4_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[1:70]=="Low")]))}),
                              m5_dm = sapply(1:602,function(i){mean(as.numeric(m5_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[71:123]=="High")]))-
                                  mean(as.numeric(m5_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[71:123]=="Low")]))}))
pc2_pvalues_all$fisher_p <- sapply(1:602,function(i){combine.test(pc2_pvalues_all[i,4:5], method = "fisher")})
pc2_pvalues_all$fisher_q <- p.adjust(pc2_pvalues_all$fisher_p, method = "BH")
pc2_pvalues_all$dm <- (pc2_pvalues_all$m4_dm+pc2_pvalues_all$m5_dm)/2
pc2_pvalues_all$Significant <- ifelse(pc2_pvalues_all$fisher_q < 0.05 & abs(pc2_pvalues_all$dm) > 0.5,"FDR < 0.05","Not Significant")

imp_ids <- c("M31260","M01643","M01303","M48152","M01564","M46173")
pc2_pufas <- c("docosadienoate (22:2n6)","docosapentaenoate (n6 DPA; 22:5n6)","docosatrienoate (22:3n3)","docosapentaenoate (n3 DPA; 22:5n3)","arachidonate (20:4n6)","linoleate (18:2n6)")
imp_mets <- metanno_m4$BIOCHEMICAL[match(imp_ids,metanno_m4$COMP_IDstr)]
imp_mets <- c(imp_mets,"cysteine","alpha-tocopherol","ascorbate (Vitamin C)","glutathione, reduced (GSH)")

ggplot(pc2_pvalues_all, aes(x =dm, y = -log10(fisher_p)))  +
  geom_point(aes(color = Significant),size = 1.5) +
  scale_color_manual(values = c("#FBB4AE","#D9D9D9")) +
  geom_text_repel(data = subset(pc2_pvalues_all, metabolite %in% c(imp_mets, pc2_pufas)), aes(label = metabolite), size=1.5, 
                  box.padding = unit(0.35, "lines"), point.padding = unit(0.35, "lines")) +
  ggtitle("PC2 High vs Low") + xlab("High - Low") +
  theme_Publication(base_size = 7) 


calculate_da <- function(pathways,data){
  da_scores <- as.data.frame(matrix(NA, nrow = length(pathways),ncol = 5))
  colnames(da_scores) <- c("Metabolite_Count","Up","Down","DA_Score","Pathway")
  da_scores$Metabolite_Count <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i],])[1]})
  da_scores$Up <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i] & data$dm > 0 & data$Significant == "FDR < 0.05",])[1]}) 
  da_scores$Down <- sapply(1:length(pathways),function(i){dim(data[data$sub_pathway == pathways[i] & data$dm < 0 & data$Significant == "FDR < 0.05",])[1]}) 
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

pc2_pathways <- names(which(sort(table(subset(pc2_pvalues_all,Significant == "FDR < 0.05")[,3]),decreasing = TRUE) > 2)) #only keep pathways with 3+ sig altered metabolites
pc2_pvalues_pathways <- pc2_pvalues_all[pc2_pvalues_all$pathway %in% pc2_pathways,]
pc2_pvalues_pathways$sub_pathway <- as.factor(pc2_pvalues_pathways$pathway)
pc2_pathways_da <- calculate_da(pc2_pathways,pc2_pvalues_pathways)
plot_da(pc2_pathways_da,"PC2 High vs Low")

#####PC1 HIGH VS LOW FOR ALL SAMPLES INLCUDING NORMALS#####
subjects_n <- unique(pca_x$SubjectID)
pc1_level <- unlist(sapply(subjects_n, function(i){pca_x[pca_x$SubjectID==i,1]  >= median(pca_x[pca_x$SubjectID==i,1])}))
pca_x$PC1_level <- as.character(ifelse(pc1_level == TRUE, "High","Low"))
pc1_pvalues_all <- data.frame(id = intersecting_metabolites,
                              metabolite = metanno_m4$BIOCHEMICAL[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                              pathway = metanno_m4$SUB_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                              m4_p = sapply(1:602,function(i){wilcox.test(as.numeric(m4_mr_all[intersecting_metabolites[i],])~pca_x$PC1_level[1:116])$p.value}),
                              m5_p = sapply(1:602,function(i){wilcox.test(as.numeric(m5_mr_all[intersecting_metabolites[i],])~pca_x$PC1_level[117:187])$p.value}),
                              m4_dm = sapply(1:602,function(i){mean(as.numeric(m4_mr_all[intersecting_metabolites[i],which(pca_x$PC1_level[1:116]=="High")]))-
                                  mean(as.numeric(m4_mr_all[intersecting_metabolites[i],which(pca_x$PC1_level[1:116]=="Low")]))}),
                              m5_dm = sapply(1:602,function(i){mean(as.numeric(m5_mr_all[intersecting_metabolites[i],which(pca_x$PC1_level[117:187]=="High")]))-
                                  mean(as.numeric(m5_mr_all[intersecting_metabolites[i],which(pca_x$PC1_level[117:187]=="Low")]))}))
pc1_pvalues_all$fisher_p <- sapply(1:602,function(i){combine.test(pc1_pvalues_all[i,4:5], method = "fisher")})
pc1_pvalues_all$fisher_q <- p.adjust(pc1_pvalues_all$fisher_p, method = "BH")
pc1_pvalues_all$dm <- (pc1_pvalues_all$m4_dm+pc1_pvalues_all$m5_dm)/2
pc1_pvalues_all$Significant <- ifelse(pc1_pvalues_all$fisher_q < 0.05 & abs(pc1_pvalues_all$dm) > 0.5,"FDR < 0.05","Not Significant")


pc1_pathways <- names(which(sort(table(subset(pc1_pvalues_all,Significant == "FDR < 0.05")[,3]),decreasing = TRUE) > 2)) #only keep pathways with 3+ sig altered metabolites
pc1_pvalues_pathways <- pc1_pvalues_all[pc1_pvalues_all$pathway %in% pc1_pathways,]
pc1_pvalues_pathways$sub_pathway <- as.factor(pc1_pvalues_pathways$pathway)
pc1_pathways_da <- calculate_da(pc1_pathways,pc1_pvalues_pathways)
plot_da(pc1_pathways_da,"PC1 High vs Low")


#####PC2 HIGH VS LOW RNA#####
rna_sample_keep <- intersect(colnames(rna_raw_fixed),rownames(pca_x_tumors)) #samples w rna and metabolomics
rna_raw_tumors <- rna_raw_fixed[,rna_sample_keep]
pca_x_tumor_keep <- pca_x_tumors[rna_sample_keep,]
all(colnames(rna_raw_tumors)==rownames(pca_x_tumor_keep)) #TRUE

dds <- DESeqDataSetFromMatrix(countData = rna_raw_tumors,
                              colData = pca_x_tumor_keep,
                              design= ~ SubjectID + PC2_level)
dds$PC2_level <- factor(dds$PC2_level, levels = c("Low","High")) #change so now it's doing PC2 High - Low
dds <- DESeq(dds)
res <- results(dds)

pca_stat<- res$stat
names(pca_stat) <- rownames(res)

fgseaRes_hm <- fgseaMultilevel(pathways = hallmark_gmt, 
                               stats    = pca_stat,
                               minSize=15,maxSize=500,nPermSimple=2000)
fgseaRes_hm15 <- fgseaRes_hm[order(fgseaRes_hm$NES, decreasing = FALSE)[1:15],]

fgseaRes_hm_padj <- data.frame(pathway= fgseaRes_hm15$pathway,
                               nes = fgseaRes_hm15$NES)
fgseaRes_hm_padj$pathway <- sub("^[^_]*_", "", fgseaRes_hm_padj$pathway)
fgseaRes_hm_padj$pathway <- gsub("_"," ",fgseaRes_hm_padj$pathway)
fgseaRes_hm_padj$pathway <- as.factor(fgseaRes_hm_padj$pathway)
fgseaRes_hm_padj$pathway <- factor(fgseaRes_hm_padj$pathway, levels = rev(fgseaRes_hm_padj$pathway))
fgseaRes_hm_padj$pval <- fgseaRes_hm15$pval

ggplot(fgseaRes_hm_padj, aes(x = nes,y=pathway)) + geom_bar(stat="identity") + xlab("NES") + ylab("") + 
  theme_Publication() + ggtitle("Hallmark Pathways") 

#####PC2 HIGH VS LOW SIGNATURE DIFFERENCES#####
m4_w_immune <- immune_mapping$MetabID[match(rownames(immune4_mr),immune_mapping$ITHID)]
m4_w_immune <- gsub("-",".",m4_w_immune)
m4_met_w_immune <- m4_mr_tumors[,intersect(m4_w_immune,colnames(m4_mr_tumors))]
m4_immune_w_met <- immune4_mr[immune_mapping$ITHID[match(gsub("\\.","-",rownames(pca_x_tumor_keep)[1:66]),immune_mapping$MetabID)],]
m4_immune_w_met <- m4_immune_w_met[,1:94] #remove KEGG pathways
m4_immune_w_met <- t(m4_immune_w_met)


m5_w_immune <- immune_mapping$MetabID[match(rownames(immune5),immune_mapping$ITHID)]
m5_w_immune <- gsub("-",".",m5_w_immune)
m5_met_w_immune <- m5_mr_tumors[,intersect(m5_w_immune, colnames(m5_mr_tumors))]
m5_immune_w_met <- immune5[immune_mapping$ITHID[match(gsub("\\.","-",rownames(pca_x_tumor_keep)[67:112]),immune_mapping$MetabID)],]
m5_immune_w_met <- m5_immune_w_met[,1:94] #remove KEGG pathways
m5_immune_w_met <- t(m5_immune_w_met)

m4_mr_i_ids <- MetaData_M4$SUBJECT_ID[match(colnames(m4[,colnames(m4_met_w_immune)]),MetaData_M4$SAMPLE_NAME)]
m4_mr_patients_i <- unique(m4_mr_i_ids)
m4_mr_i_sample_count <- sapply(m4_mr_patients_i,function(x)sum(m4_mr_i_ids==x))
m4_mr_i_sample_sum <- cumsum(m4_mr_i_sample_count)

m5_mr_i_ids <- MetaData_M5$SUBJECT_ID[match(colnames(m5[,colnames(m5_met_w_immune)]),MetaData_M5$SAMPLE_NAME)]
m5_mr_patients_i <- unique(m5_mr_i_ids)
m5_mr_i_sample_count <- sapply(m5_mr_patients_i,function(x)sum(m5_mr_i_ids==x))
m5_mr_i_sample_sum <- cumsum(m5_mr_i_sample_count)


m4_patient_centered_all_i <- as.data.frame(matrix(NA,nrow=94,ncol = 61))
for(x in 1:nrow(m4_immune_w_met)){
  m4_patient_centered_all_i[x,] <- unlist(sapply(2:length(m4_mr_patients_i),function(i){m4_immune_w_met[x,as.numeric(m4_mr_i_sample_sum[i-1]+1):m4_mr_i_sample_sum[i]]-median(as.numeric(m4_immune_w_met[x,as.numeric(m4_mr_i_sample_sum[i-1]+1):m4_mr_i_sample_sum[i]]))}))
}
m4_patient_centered_all_p1_i <- as.data.frame(matrix(NA,nrow=94,ncol = 5))
for(x in 1:nrow(m4_immune_w_met)){
  m4_patient_centered_all_p1_i[x,] <- m4_immune_w_met[x,1:m4_mr_i_sample_sum[1]]-median(as.numeric(m4_immune_w_met[x,1:m4_mr_i_sample_sum[1]]))
}
m4_patient_centered_all_i <- cbind(m4_patient_centered_all_p1_i,m4_patient_centered_all_i)
colnames(m4_patient_centered_all_i) <- colnames(m4_immune_w_met)
rownames(m4_patient_centered_all_i) <- rownames(m4_immune_w_met)
rm(m4_patient_centered_all_p1_i)


m5_patient_centered_all_i <- as.data.frame(matrix(NA,nrow=94,ncol = 41))
for(x in 1:nrow(m5_immune_w_met)){
  m5_patient_centered_all_i[x,] <- unlist(sapply(2:length(m5_mr_patients_i),function(i){m5_immune_w_met[x,as.numeric(m5_mr_i_sample_sum[i-1]+1):m5_mr_i_sample_sum[i]]-median(as.numeric(m5_immune_w_met[x,as.numeric(m5_mr_i_sample_sum[i-1]+1):m5_mr_i_sample_sum[i]]))}))
}
m5_patient_centered_all_p1_i <- as.data.frame(matrix(NA,nrow=94,ncol = 5))
for(x in 1:nrow(m5_immune_w_met)){
  m5_patient_centered_all_p1_i[x,] <- m5_immune_w_met[x,1:m5_mr_i_sample_sum[1]]-median(as.numeric(m5_immune_w_met[x,1:m5_mr_i_sample_sum[1]]))
}
m5_patient_centered_all_i <- cbind(m5_patient_centered_all_p1_i,m5_patient_centered_all_i)
colnames(m5_patient_centered_all_i) <- colnames(m5_immune_w_met)
rownames(m5_patient_centered_all_i) <- rownames(m5_immune_w_met)
rm(m5_patient_centered_all_p1_i)


m5_pca2_level <- pca_x_tumor_keep$PC2_level[67:112]
pca_x_tumor_keep$PC2_level <- as.factor(pca_x_tumor_keep$PC2_level)

pc2_immune <- data.frame(signature = rownames(m4_patient_centered_all_i),
                         m4_p = sapply(1:94,function(i){wilcox.test(as.numeric(m4_patient_centered_all_i[i,])~pca_x_tumor_keep$PC2_level[1:66])$p.value}),
                         m4_diff = sapply(1:94,function(i){mean(as.numeric(m4_patient_centered_all_i[i,which(pca_x_tumor_keep$PC2_level[1:66]=="High")]))-
                             mean(as.numeric(m4_patient_centered_all_i[i,which(pca_x_tumor_keep$PC2_level[1:66]=="Low")]))}),
                         m5_p = sapply(1:94,function(i){wilcox.test(as.numeric(m5_patient_centered_all_i[i,])~m5_pca2_level)$p.value}),
                         m5_diff = sapply(1:94,function(i){mean(as.numeric(m5_patient_centered_all_i[i,which(m5_pca2_level=="High")]))-
                             mean(as.numeric(m5_patient_centered_all_i[i,which(m5_pca2_level=="Low")]))}))
pc2_immune$wt_p <- sapply(1:94,function(i){tryCatch(combine.test(pc2_immune[i,c(2,4)],method = "fisher"), error=function(e) NA)})
pc2_immune$wt_q <- p.adjust(pc2_immune$wt_p, method = "fdr")
pc2_immune$dm <- (pc2_immune$m4_diff+pc2_immune$m5_diff)/2

pc2_immune_sig <- pc2_immune[which(pc2_immune$wt_q < 0.05),]
pc2_immune_sig$signature <- factor(pc2_immune_sig$signature, levels = pc2_immune_sig[order(pc2_immune_sig$wt_q, decreasing = TRUE),1])
pc2_immune_sig$FoldChange <- ifelse(pc2_immune_sig$dm > 0, "Positive","Negative")

ggplot(pc2_immune_sig, aes(x = signature, y = -log10(wt_q), fill=FoldChange)) + geom_bar(stat="identity") + coord_flip() + theme_Publication() +
  xlab("") + ylab("-Log10 FDR-Corrected P-Value") + scale_fill_manual(values = c("#B3CDE3","#FBB4AE")) 

#plot angiogenesis signature
mcd_angio <- data.frame(exp = c(as.numeric(m4_patient_centered_all_i[59,]),
                                as.numeric(m5_patient_centered_all_i[59,])),
                        level = pca_x_tumor_keep$PC2_level)
ylim <- boxplot.stats(mcd_angio$exp)$stats[c(1, 5)]
c3 <- ggplot(mcd_angio,aes(x=level,y=exp)) + geom_boxplot(outlier.size=0.3) + theme_Publication() + geom_jitter(shape=16, position=position_jitter(0.2), size=0.6) +
  xlab("PC2 Level") + ylab("Angiogenesis Signature Expression") + ylim(c(-250,250))


#####CYSTEINE IS UNIQUELY HETEROGENEOUS#####
amino_acids <- c("arginine","histidine","lysine","aspartate","glutamate","serine","threonine","asparagine","glutamine","cysteine","glycine",      
                 "proline","alanine","valine","isoleucine","leucine","methionine","phenylalanine","tyrosine","tryptophan")
aa_ids <- metanno_m4$COMP_IDstr[match(amino_acids,metanno_m4$BIOCHEMICAL)]

aa_heterogeneity <- rbind(apply(intra_m4_met[,aa_ids],2,median),
                          apply(intra_m5_met[,aa_ids],2,median))
aa_heterogeneity <- as.data.frame(t(aa_heterogeneity))
aa_heterogeneity$metabolite <- amino_acids
aa_heterogeneity$average <- (aa_heterogeneity$V1+aa_heterogeneity$V2)/2
aa_heterogeneity$rank <- rank(aa_heterogeneity$average)

ggplot(aa_heterogeneity, aes(x=rank, y= average)) + geom_point() + 
  theme_Publication() + xlab("") + ylab("Intratumor Heterogeneity") +
  geom_text_repel(data=subset(aa_heterogeneity,metabolite == "cysteine"), aes(label = metabolite), size = 2, color = "black") 

#####CYSTEINE SCATTERPLOTS FIG 4F#####
m4_pca_sp <- data.frame(cysteine = scale(as.numeric(m4_mr_tumors["M01868",])),
                        gssg = scale(as.numeric(m4_mr_tumors["M27727",])),
                        gsh = scale(as.numeric(m4_mr_tumors["M02127",])),
                        glycine = scale(as.numeric(m4_mr_tumors["M00058",])),
                        glutamate = scale(as.numeric(m4_mr_tumors["M00057",])))
m5_pca_sp <- data.frame(cysteine = scale(as.numeric(m5_mr_tumors["M01868",])),
                        gssg = scale(as.numeric(m5_mr_tumors["M27727",])),
                        gsh = scale(as.numeric(m5_mr_tumors["M02127",])),
                        glycine = scale(as.numeric(m5_mr_tumors["M00058",])),
                        glutamate = scale(as.numeric(m5_mr_tumors["M00057",])))

pca_sp <- rbind(m4_pca_sp,m5_pca_sp)
pca_sp$batch <- c(rep("M4",70),rep("M5",53))

ggplot(pca_sp,aes(x=cysteine,y=gssg,color=batch)) + geom_point(size=0.75) + geom_smooth(method=lm, se=FALSE) + 
  scale_color_manual(values=c('#E69F00', '#56B4E9')) + theme_Publication() + xlab("Cysteine") + ylab("Oxidized Glutathione")

ggplot(pca_sp,aes(x=cysteine,y=gsh,color=batch)) + geom_point(size=0.75) + geom_smooth(method=lm, se=FALSE) + 
  scale_color_manual(values=c('#E69F00', '#56B4E9')) + theme_Publication() + xlab("Cysteine") + ylab("Reduced Glutathione")

ggplot(pca_sp,aes(x=cysteine,y=glycine,color=batch)) + geom_point(size=0.75) + geom_smooth(method=lm, se=FALSE) + 
  scale_color_manual(values=c('#E69F00', '#56B4E9')) + theme_Publication() + xlab("Cysteine") + ylab("Glycine")

ggplot(pca_sp,aes(x=cysteine,y=glutamate,color=batch)) + geom_point(size=0.75) + geom_smooth(method=lm, se=FALSE) + 
  scale_color_manual(values=c('#E69F00', '#56B4E9')) + theme_Publication() + xlab("Cysteine") + ylab("Glutamate")

#because these are scaled values can do correlations across cohorts
cor.test(pca_sp$cysteine,pca_sp$gssg,method = "spearman") #rho = -0.02, p = 0.82
cor.test(pca_sp$cysteine,pca_sp$gsh,method = "spearman") #rho = 0.42, p = 1.0x10-6
cor.test(pca_sp$cysteine,pca_sp$glycine,method = "spearman") #rho = 0.66, p < 2.2e-16
cor.test(pca_sp$cysteine,pca_sp$glutamate,method = "spearman") #rho = 0.65, p = 4.2e-16

#####GPX4 STAINING COMPARISONS#####
pca_x_tumors$GPX4 <- GPX4_IHC_New$Score[match(pca_x_tumors$ID,GPX4_IHC_New$Sample)]
pca_x_tumors_ihc <- pca_x_tumors[which(!is.na(pca_x_tumors$GPX4)),]

ihc_samples <- unique(pca_x_tumors_ihc$SubjectID)[1:5]
cysteine_max <- sapply(1:5,function(i){which.max(pca_x_tumors_ihc[pca_x_tumors_ihc$SubjectID==ihc_samples[i],9])})
cysteine_min <- sapply(1:5,function(i){which.min(pca_x_tumors_ihc[pca_x_tumors_ihc$SubjectID==ihc_samples[i],9])})
gpx4_max <- sapply(1:5,function(i){pca_x_tumors_ihc[pca_x_tumors_ihc$SubjectID==ihc_samples[i],14][cysteine_max[i]]})
gpx4_min <- sapply(1:5,function(i){pca_x_tumors_ihc[pca_x_tumors_ihc$SubjectID==ihc_samples[i],14][cysteine_min[i]]})

cysteine_gpx4 <- data.frame(gpx4 = c(gpx4_max,gpx4_min),
                            cysteine = c(rep("max",5),rep("min",5)),
                            sample = rep(ihc_samples,times=2))
ggplot(cysteine_gpx4,aes(x=cysteine,y=gpx4)) + geom_boxplot() + 
  geom_line(aes(x  = cysteine, y = gpx4, group = sample)) +
  theme_Publication() + xlab("Cysteine Level") + ylab("GPX4 H-Score") 

wilcox.test(cysteine_gpx4$gpx4~cysteine_gpx4$cysteine, paired=TRUE) #p = 0.05

######PCA IN NORMAL SAMPLES#####
pca_x_normals <- subset(pca_x,SampleType == "TISSUE_NORMAL")
subjects_n <- unique(pca_x_normals$SubjectID)
pc2_level_n <- unlist(sapply(subjects_n, function(i){pca_x_normals[pca_x_normals$SubjectID==i,2]  >= median(pca_x_normals[pca_x_normals$SubjectID==i,2])}))
pca_x_normals$PC2_level <- as.character(ifelse(pc2_level_n == TRUE, "High","Low"))
pca_x_normals$ID <- c(MetaData_M4$CLIENT_IDENTIFIER[match(colnames(m4_mr_normals),MetaData_M4$SAMPLE_NAME)], 
                      MetaData_M5$CLIENT_IDENTIFIER[match(colnames(m5_mr_normals),MetaData_M5$SAMPLE_NAME)])
pca_x_normals$Batch <- c(rep("M4",ncol(m4_mr_normals)),rep("M5",ncol(m5_mr_normals)))

cysteine_pcn <- data.frame(exp = c(as.numeric(m4_mr_normals_c["M01868",]),
                                   as.numeric(m5_mr_normals_c["M01868",])),
                           level = as.factor(pca_x_normals$PC2_level))
ggplot(cysteine_pcn, aes(x=level,y=exp)) + geom_boxplot(outlier.size = 0.3) + geom_jitter(shape=16, position=position_jitter(0.2),size=0.6) +
  theme_Publication() + xlab("PC2 Level") + ylab("Cysteine Abundance") 


#####PCA VS CD31 IF STAINING IN TUMORS######
if_results_mr <- if_results[1:10,]
if_results_mr$pc2 <- as.factor(c("high","low","high","low","high","low","high","high","low","high"))

wilcox.test(if_results_mr$CD31..um2.~if_results_mr$pc2) #wilcox test p-value = 0.11

c6 <- ggplot(if_results_mr,aes(x=pc2,y=CD31..um2.)) + geom_boxplot(outlier.size = 0.3) + theme_Publication() + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 0.6) +
  xlab("PC2 Level") + ylab("CD31 IF")


######PCA2 HIGH VS LOW LIPIDS ONLY######
#lipid comparison of pc2 high vs low
pc2_pvalues_all <- data.frame(id = intersecting_metabolites,
                              metabolite = metanno_m4$BIOCHEMICAL[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                              super = metanno_m4$SUPER_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                              pathway = metanno_m4$SUB_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)],
                              m4_p = sapply(1:602,function(i){wilcox.test(as.numeric(m4_mr_tumors_c[intersecting_metabolites[i],])~pca_x_tumors$PC2_level[1:70])$p.value}),
                              m5_p = sapply(1:602,function(i){wilcox.test(as.numeric(m5_mr_tumors_c[intersecting_metabolites[i],])~pca_x_tumors$PC2_level[71:123])$p.value}),
                              m4_dm = sapply(1:602,function(i){mean(as.numeric(m4_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[1:70]=="High")]))-
                                  mean(as.numeric(m4_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[1:70]=="Low")]))}),
                              m5_dm = sapply(1:602,function(i){mean(as.numeric(m5_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[71:123]=="High")]))-
                                  mean(as.numeric(m5_mr_tumors_c[intersecting_metabolites[i],which(pca_x_tumors$PC2_level[71:123]=="Low")]))}))
pc2_pvalues_all$fisher_p <- sapply(1:602,function(i){combine.test(pc2_pvalues_all[i,5:6], method = "fisher")})
pc2_pvalues_all <- pc2_pvalues_all[pc2_pvalues_all$super=="Lipid",]
pc2_pvalues_all$fisher_q <- p.adjust(pc2_pvalues_all$fisher_p, method = "BH")
pc2_pvalues_all$dm <- (pc2_pvalues_all$m4_dm+pc2_pvalues_all$m5_dm)/2
pc2_pvalues_all$Significant <- ifelse(pc2_pvalues_all$fisher_q < 0.05 & abs(pc2_pvalues_all$dm) > 0.5, "FDR < 0.05","NS")

pc2_pathways <- names(which(sort(table(pc2_pvalues_all$pathway),decreasing = TRUE)>2)) #only keep pathways with 3+ sig altered metabolites
pc2_pvalues_pathways <- pc2_pvalues_all[pc2_pvalues_all$pathway %in% pc2_pathways,]
pc2_pvalues_pathways$sub_pathway <- as.factor(pc2_pvalues_pathways$pathway)
pc2_pathways_da <- calculate_da(pc2_pathways,pc2_pvalues_pathways)
plot_da(pc2_pathways_da,"PC2 High vs Low")

#####PC2 HIGH VS LOW COMPARING CELL CYCLE AND FERROPTOSIS SIGNATURES#####
new_met_names <- rownames(pca_x_tumor_keep)
new_met_names <- gsub("\\.","-",new_met_names)
new_rna_names <- immune_mapping$RNAID[match(new_met_names,immune_mapping$MetabID)]

ferr_cc <- ferr_cc[new_rna_names,] #only the samples that are mr tumors
ferr_cc$level <- pca_x_tumor_keep$PC2_level

wilcox.test(ferr_cc$WP_FERROPTOSIS~ferr_cc$level) #p.value = 0.9, dm high - low = -0.151
c8 <- ggplot(ferr_cc,aes(x=level,y=WP_FERROPTOSIS)) + geom_boxplot(outlier.size = 0.3) + geom_jitter(shape = 16, size = 0.6) + theme_Publication() + 
  xlab("PC2 Level") + ylab("Ferroptosis Signature Expression")
ggsave("~/Documents/RCC Final/PC2 Ferroptosis Signature Boxplot.pdf", width = 1.5, height = 2, units = "in")

wilcox.test(ferr_cc$Cell.Cycle~ferr_cc$level) #p.value = 0.3, dm high - low = -318.691
c9 <- ggplot(ferr_cc,aes(x=level,y=Cell.Cycle)) + geom_boxplot(outlier.size = 0.3) + geom_jitter(shape = 16, size = 0.6) + theme_Publication() + 
  xlab("PC2 Level") + ylab("Cell Cycle Signature Expression")
ggsave("~/Documents/RCC Final/PC2 Cell Cycle Signature Boxplot.pdf", width = 1.5, height = 2, units = "in")



#####KI67 STAINING#####
ki67$Ki.67..Average <- as.numeric(gsub("%","",ki67$Ki.67..Average))
ki67$Ki.67..Average[10] <- 0
ki67$level <- pca_x_tumors$PC2_level[match(ki67$ID,pca_x_tumors$ID)]

wilcox.test(ki67$Ki.67..Average~ki67$level)
ggplot(ki67,aes(x=level,y=Ki.67..Average)) + geom_boxplot(outlier.size = 0.3) + 
  geom_jitter(shape = 16, size = 0.6) + 
  theme_Publication() + xlab("PC2 Level") + ylab("Average Ki67 Staining")


#####LI ET AL ANALYSIS#####
colnames(met) <- gsub("\\.","-",colnames(met))
met <- met[,-c(1:7)]
colnames(met) <- paste(colnames(met),"-T",sep="")

is <- estimate$ImmuneScore_ESTIMATE

is_cor <- data.frame(metabolite = rownames(met),
                     p.value = sapply(1:183,function(i){cor.test(as.numeric(met[i,]),is,method="spearman")$p.value}),
                     rho = sapply(1:183,function(i){cor.test(as.numeric(met[i,]),is,method="spearman")$estimate}))
is_cor$p.adj <- p.adjust(is_cor$p.value,method = "BH")

rownames(met) <- tolower(rownames(met))
rownames(met)[133] <- c("nicotinamide adenine dinucleotide (NAD+)")
rownames(met)[126] <- c("N6,N6,N6-trimethyllysine")
rownames(met)[33] <- c("alpha-ketoglutarate")
rownames(met)[37] <- c("argininosuccinate")
rownames(met)[78] <- c("fructose 1,6-diphosphate/glucose 1,6-diphosphate/myo-inositol diphosphates")
rownames(met)[87] <- c("glutathione, reduced (GSH)")

int <- intersect(rownames(met),metanno_m4$BIOCHEMICAL[match(intersecting_metabolites,metanno_m4$COMP_IDstr)])

#check which ones overlap with immunescore
is_cor <- data.frame(metabolite = rownames(met),
                     p.value = sapply(1:183,function(i){cor.test(as.numeric(met[i,]),is,method="spearman")$p.value}),
                     rho = sapply(1:183,function(i){cor.test(as.numeric(met[i,]),is,method="spearman")$estimate}))
is_cor$p.adj <- p.adjust(is_cor$p.value,method = "BH")

is_cor <- is_cor[is_cor$metabolite %in% int,]

pb_data = lapply(pb_fs, function(x) read.delim(file=x))

pb_results <- do.call("rbind", pb_data)
pb_results$number <- file_names
pb_results$metID <- javelin_p$metID[as.numeric(pb_results$number)]
pb_results$metabolite <- metanno_m4$BIOCHEMICAL[match(pb_results$metID,metanno_m4$COMP_IDstr)]
pb_results$pathway <- metanno_m4$SUB_PATHWAY[match(pb_results$metID,metanno_m4$COMP_IDstr)]
pb_results$is_adj <- p.adjust(pb_results$is, method = "BH")

pb_results <- pb_results[pb_results$metabolite %in% int,]

rownames(is_cor) <- is_cor$metabolite
rownames(pb_results) <- pb_results$metabolite

is_cor <- is_cor[int,]
pb_results <- pb_results[int,]

is_cor$pb_p <- pb_results$is
is_cor$pb_q <- p.adjust(is_cor$pb_p, method= "BH")

is_cor$p.adj <- p.adjust(is_cor$p.value, method= "BH")

immune_score_cor <- immune_met_correlations(3)
rownames(immune_score_cor) <- immune_score_cor$metabolite
immune_score_cor <- immune_score_cor[int,]

is_cor$mr_r <- immune_score_cor$corR
is_cor$li_sign <- ifelse(is_cor$rho < 0, -1, 1)
is_cor$li_value <- -log10(is_cor$p.adj) * is_cor$li_sign
is_cor$mr_sign <- ifelse(is_cor$mr_r < 0, -1, 1)
is_cor$mr_value <- -log10(is_cor$pb_q) * is_cor$mr_sign

cor.test(is_cor$li_value,is_cor$mr_value, method = "spearman") #rho = 0.5292856 , p = 1.156e-05
ggplot(is_cor, aes(x=li_value,y=mr_value)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + theme_Publication() + 
  xlab("Li et al -Log10 P-Value x Correlation Sign") + ylab("Mixed Effects Bootstrapping -Log10 P-Value x Correlation Sign") +
  ggtitle("Immune Score") + geom_text_repel(data = is_cor, aes(label = metabolite), size = 2)


#####SPATIAL ANALYSES#####
aas = c('alanine','aspartate','asparagine','cysteine','glutamate','glutamine','glycine','arginine','histidine','lysine',
        'serine','threonine','proline','valine','isoleucine','methionine','phenylalanine','tyrosine','tryptophan','leucine')

aa_measured <- intersect(colnames(het),aas)

heterogeneity <- apply(het[,aa_measured],2,median)

heterogeneity <- data.frame(metabolite = aa_measured,
                            het = heterogeneity)
heterogeneity$rank <- rank(heterogeneity$het)

ggplot(heterogeneity, aes(x=rank,y=het)) + geom_point(size = 0.6) + theme_Publication() +
  geom_text_repel(data = heterogeneity, aes(label = metabolite), size = 2) +
  xlab("") + ylab("Intrakidney Heterogeneity")

glutamine_tracing$fisher_p <- sapply(1:3085,function(i){combine.test(glutamine_tracing[i,8:10], method = "fisher")})
glutamine_tracing$corZ <- (FisherZ(glutamine_tracing$cor_0) + FisherZ(glutamine_tracing$cor_1) + FisherZ(glutamine_tracing$cor_2))/3
glutamine_tracing$corR <- FisherZInv(glutamine_tracing$corZ)
glutamine_tracing$fisher_q <- p.adjust(glutamine_tracing$fisher_p, method = "BH")
glutamine_tracing$fisher_p[which(glutamine_tracing$fisher_p == 0)] <- 9.881313e-324
glutamine_tracing$sig <- ifelse(glutamine_tracing$fisher_q < 0.05, "FDR < 0.05", "Not Significant")
glutamine_tracing$col <- ifelse(glutamine_tracing$fisher_q < 0.05, "Other Significant", "Not Significant")
glutamine_tracing$col[199:201] <- c("Cysteine")
glutamine_tracing$col[1920:1929] <- c("Glutathione")
glutamine_tracing$col[1212:1217] <- c("Citrate")

glucose_tracing$fisher_p <- sapply(1:3068,function(i){combine.test(glucose_tracing[i,8:10], method = "fisher")})
glucose_tracing$fisher_p[which(glucose_tracing$fisher_p == 0)] <- 9.881313e-324
glucose_tracing$corZ <- (FisherZ(glucose_tracing$cor_0) + FisherZ(glucose_tracing$cor_1) + FisherZ(glucose_tracing$cor_2))/3
glucose_tracing$corR <- FisherZInv(glucose_tracing$corZ)
glucose_tracing$fisher_q <- p.adjust(glucose_tracing$fisher_p, method = "BH")
glucose_tracing$sig <- ifelse(glucose_tracing$fisher_q < 0.05, "FDR < 0.05", "Not Significant")
glucose_tracing$col <- ifelse(glucose_tracing$fisher_q < 0.05, "Other Significant", "Not Significant")
glucose_tracing$col[196:198] <- c("Cysteine")
glucose_tracing$col[1905:1914] <- c("Glutathione")
glucose_tracing$col[1201:1206] <- c("Citrate")



ggplot(glutamine_tracing, aes(x=corR, y=-log10(fisher_p))) + 
  geom_point(aes(color = col),size = 0.1) +
  geom_point(data = subset(glutamine_tracing, col != 'Other Significant'),
             aes(x = corR, y = -log10(fisher_p), color = col), size = 0.1) +
  scale_color_manual(values = c("#000000","#D9D9D9", "#FBB4AE")) +
  theme_Publication() +
  geom_text_repel(data=subset(glutamine_tracing,col == 'Glutathione'), aes(label = X), size = 2)

ggplot(glucose_tracing, aes(x=corR, y=-log10(fisher_p))) + 
  geom_point(aes(color = col),size = 0.1) +
  geom_point(data = subset(glucose_tracing, col != 'Other Significant'),
             aes(x = corR, y = -log10(fisher_p), color = col), size = 0.1)  +
  scale_color_manual(values = c("#000000","#D9D9D9", "#FBB4AE")) +
  theme_Publication() +
  geom_text_repel(data=subset(glucose_tracing,col == 'Glutathione'), aes(label = X), size = 0.25)






