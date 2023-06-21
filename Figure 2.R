load("All Data.RData")
load("Outlier Data.RData")

library(fields)
library(annotate)
library(org.Hs.eg.db)
library(survcomp)
library(ggplot2)
library(ggrepel)
library(stringr)
library(reshape2)
library(ComplexHeatmap)
library(ape)
library(phangorn)
library(psych)
library(data.table)
library(qusage)
library(fgsea)

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

mr_rna_tumors <- lfiltn_rna_all[,which(rna_mapping$Sample_Type=="Tumor")] #all tumors not just ones with metabolomics
table(rna_mapping$Patient.ID[match(colnames(mr_rna_tumors),rna_mapping$Sample)]) #MR05, SC01 and SC04 only have 1 region, need to remove
which(rna_mapping$Patient.ID[match(colnames(mr_rna_tumors),rna_mapping$Sample)] %in% c("MR05","SC01","SC04","SC07","MR04")) #16, 108, 112, 117, 118, 119
mr_rna_tumors <- mr_rna_tumors[,-c(16, 108,112,117,118,119)]
which(rna_mapping$Sample.Name[match(colnames(mr_rna_tumors),rna_mapping$Sample)] %in% m5_outliers) #no m4 outliers in rna samples; 90,100,125,126
mr_rna_tumors <- mr_rna_tumors[,-c(99,100,125,126)]

rna_mr_tumors_subjectid <- rna_mapping$Patient.ID[match(colnames(mr_rna_tumors),rna_mapping$Sample)]
rna_mr_tumors_patients <- unique(rna_mr_tumors_subjectid)
rna_mr_patient_count <- sapply(rna_mr_tumors_patients,function(x)sum(rna_mr_tumors_subjectid==x))
rna_mr_patient_sum <- cumsum(rna_mr_patient_count)





m4_mr_normals <- m4[,MetaData_M4$TISSUE_STATUS[match(colnames(m4),MetaData_M4$SAMPLE_NAME)]=="TISSUE_NORMAL"] #all normals are from MR samples but MR06 only has 1 region so remove it
m4_mr_normals <- m4_mr_normals[,-50]
which(MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_normals),MetaData_M4$SAMPLE_NAME)]=="MR04") #44,45,46
m4_mr_normals <- m4_mr_normals[,-c(44:46)]
m4_mr_normals <- m4_mr_normals[intersecting_metabolites,]
m4_mr_normals_subjectid <- MetaData_M4$SUBJECT_ID[match(colnames(m4_mr_normals),MetaData_M4$SAMPLE_NAME)]
m4_mr_normals_patients <- unique(m4_mr_normals_subjectid)
m4_mr_normals_count <- sapply(m4_mr_normals_patients,function(x)sum(m4_mr_normals_subjectid==x))
m4_mr_normals_sum <- cumsum(m4_mr_normals_count)

m5_mr_normals <- m5[,MetaData_M5$GROUP[match(colnames(m5),MetaData_M5$SAMPLE_NAME)]=="NORMAL"] #all are MR but all the SC patients only have 1 region
which(MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_normals),MetaData_M5$SAMPLE_NAME)]=="SC03") #19 remove 19-28
m5_mr_normals <- m5_mr_normals[,-c(19:28)]
m5_mr_normals <- m5_mr_normals[intersecting_metabolites,]
m5_mr_normals_subjectid <- MetaData_M5$SUBJECT_ID[match(colnames(m5_mr_normals),MetaData_M5$SAMPLE_NAME)]
m5_mr_normals_patients <- unique(m5_mr_normals_subjectid)
m5_mr_normals_count <- sapply(m5_mr_normals_patients,function(x)sum(m5_mr_normals_subjectid==x))
m5_mr_normals_sum <- cumsum(m5_mr_normals_count)


mr_rna_normals <- lfiltn_rna_all[,which(rna_mapping$Sample_Type=="Normal")] #all tumors not just ones with metabolomics
table(rna_mapping$Patient.ID[match(colnames(mr_rna_normals),rna_mapping$Sample)]) #MR06 SC03 SC04 SC06 SC07 SC08 SC09 SC10 SC12 SC14 only have 1 region, need to remove; no MR04 normals
which(rna_mapping$Patient.ID[match(colnames(mr_rna_normals),rna_mapping$Sample)] %in% c("MR06","SC03","SC04","SC06","SC07","SC08","SC09","SC10","SC12","SC14")) #112 63 64 65 66 67 68 69 70 71
mr_rna_normals <- mr_rna_normals[,-c(12, 63, 64, 65, 66, 67, 68, 69, 70, 71)]
mr_rna_normals_subjectid <- rna_mapping$Patient.ID[match(colnames(mr_rna_normals),rna_mapping$Sample)]
mr_rna_normals_patients <- unique(mr_rna_normals_subjectid)
mr_rna_normals_count <- sapply(mr_rna_normals_patients,function(x)sum(mr_rna_normals_subjectid==x))
mr_rna_normals_sum <- cumsum(mr_rna_normals_count)





#median center the data for heterogeneity calculations
median_center <- function(data){
  medians <- apply(data,1,median)
  data <- data - medians
  return(data)
}

m4_mr_tumors_mc <- median_center(m4_mr_tumors)
m5_mr_tumors_mc <- median_center(m5_mr_tumors)
rna_mr_tumors_mc <- median_center(mr_rna_tumors)

m4_mr_normals_mc <- median_center(m4_mr_normals)
m5_mr_normals_mc <- median_center(m5_mr_normals)
rna_mr_normals_mc <- median_center(mr_rna_normals)




#intrapatient heterogeneity calculation
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

intra_rna <- intraheterogeneity(rna_mr_tumors_mc, 28, 5, rna_mr_patient_sum)
rownames(intra_rna) <- as.character(rna_mr_tumors_patients)
rownames(intra_rna)[1:16] <- gsub("NIVO","NIVOT",rownames(intra_rna)[1:16])
colnames(intra_rna) <- as.character(getSYMBOL(rownames(rna_mr_tumors_mc), data='org.Hs.eg'))



intra_m4_metn <- intraheterogeneity(m4_mr_normals_mc, 16, 3, m4_mr_normals_sum)
rownames(intra_m4_metn) <- as.character(m4_mr_normals_patients)
colnames(intra_m4_metn) <- rownames(m4_mr_normals_mc)

intra_m5_metn <- intraheterogeneity(m5_mr_normals_mc, 6, 3, m5_mr_normals_sum)
rownames(intra_m5_metn) <- as.character(m5_mr_normals_patients)
colnames(intra_m5_metn) <- rownames(m5_mr_normals_mc)

intra_rnan <- intraheterogeneity(rna_mr_normals_mc, 21, 3, mr_rna_normals_sum)
rownames(intra_rnan) <- as.character(mr_rna_normals_patients)
rownames(intra_rnan)[1:16] <- gsub("NIVO","NIVOT",rownames(intra_rnan)[1:16])
colnames(intra_rnan) <- as.character(getSYMBOL(rownames(rna_mr_normals_mc), data='org.Hs.eg'))


#interpatient heterogeneity calculation
interheterogeneity<- function(data, npatients,n_p1,sample_sum){
  inter <- as.data.frame(matrix(NA,ncol=1,nrow=npatients))
  for(x in 1:nrow(data)){
    distance_matrix <- rdist(as.numeric(data[x,]))
    diag(distance_matrix) <- NA
    
    pat1 <- median(as.numeric(distance_matrix[1:n_p1,-c(1:n_p1)]),na.rm=TRUE)
    rest <- sapply(2:npatients,function(i){median(na.omit(as.numeric(unlist(distance_matrix[as.numeric(sample_sum[i-1]+1):sample_sum[i],-c(as.numeric(sample_sum[i-1]+1):sample_sum[i])]))))})
    inter[,x] <- c(pat1,rest)}
  return(inter)
}

inter_m4_met <- interheterogeneity(m4_mr_tumors_mc, 16, 5, m4_mr_patient_sum)
rownames(inter_m4_met) <- as.character(m4_mr_patients)
colnames(inter_m4_met) <- rownames(m4_mr_tumors_mc)

inter_m5_met <- interheterogeneity(m5_mr_tumors_mc, 14, 5, m5_mr_patient_sum)
rownames(inter_m5_met) <- as.character(m5_mr_patients)
colnames(inter_m5_met) <- rownames(m5_mr_tumors_mc)

inter_rna <- interheterogeneity(rna_mr_tumors_mc, 28, 5, rna_mr_patient_sum)
rownames(inter_rna) <- as.character(rna_mr_tumors_patients)
rownames(inter_rna)[1:16] <- gsub("NIVO","NIVOT",rownames(inter_rna)[1:16])
colnames(inter_rna) <- as.character(getSYMBOL(rownames(rna_mr_tumors_mc), data='org.Hs.eg'))



inter_m4_metn <- interheterogeneity(m4_mr_normals_mc, 16, 3, m4_mr_normals_sum)
rownames(inter_m4_metn) <- as.character(m4_mr_normals_patients)
colnames(inter_m4_metn) <- rownames(m4_mr_normals_mc)

inter_m5_metn <- interheterogeneity(m5_mr_normals_mc, 6, 3, m5_mr_normals_sum)
rownames(inter_m5_metn) <- as.character(m5_mr_normals_patients)
colnames(inter_m5_metn) <- rownames(m5_mr_normals_mc)

inter_rnan <- interheterogeneity(rna_mr_normals_mc, 21, 3, mr_rna_normals_sum)
rownames(inter_rnan) <- as.character(mr_rna_normals_patients)
rownames(inter_rnan)[1:16] <- gsub("NIVO","NIVOT",rownames(inter_rnan)[1:16])
colnames(inter_rnan) <- as.character(getSYMBOL(rownames(rna_mr_normals_mc), data='org.Hs.eg'))


#calculate standard deviation, standard error and confidence interval
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



m4_heterogeneity <- data.frame(Heterogeneity = c(as.numeric(apply(intra_m4_met,1,median)),as.numeric(apply(inter_m4_met,1,median))),
                               Inter = c(rep("Intratumor",nrow(intra_m4_met)),rep("Intertumor",nrow(inter_m4_met))),
                               Sample = c(rownames(intra_m4_met),rownames(inter_m4_met)))
m4_heterogeneity_summary <- summarySE(m4_heterogeneity[,1:2],measurevar = "Heterogeneity", groupvars = "Inter")
m4_heterogeneity_summary$median <- c(median(apply(intra_m4_met,1,median)),median(apply(inter_m4_metn,1,median)))


m5_heterogeneity <- data.frame(Heterogeneity = c(as.numeric(apply(intra_m5_met,1,median)),as.numeric(apply(inter_m5_met,1,median))),
                               Inter = c(rep("Intratumor",nrow(intra_m5_met)),rep("Intertumor",nrow(inter_m5_met))),
                               Sample = c(rownames(intra_m5_met),rownames(inter_m5_met)))
m5_heterogeneity_summary <- summarySE(m5_heterogeneity[,1:2],measurevar = "Heterogeneity", groupvars = "Inter")
m5_heterogeneity_summary$median <- c(median(apply(intra_m5_met,1,median)),median(apply(inter_m5_metn,1,median)))


m4_heterogeneityn <- data.frame(Heterogeneity = c(as.numeric(apply(intra_m4_metn,1,median)),as.numeric(apply(inter_m4_metn,1,median))),
                                Inter = c(rep("Intratumor",nrow(intra_m4_metn)),rep("Intertumor",nrow(inter_m4_metn))),
                                Sample = c(rownames(intra_m4_metn),rownames(inter_m4_metn)))
m4_heterogeneity_summaryn <- summarySE(m4_heterogeneityn[,1:2],measurevar = "Heterogeneity", groupvars = "Inter")
m4_heterogeneity_summaryn$median <- c(median(apply(intra_m4_metn,1,median)),median(apply(inter_m4_metn,1,median)))


m5_heterogeneityn <- data.frame(Heterogeneity = c(as.numeric(apply(intra_m5_metn,1,median)),as.numeric(apply(inter_m5_metn,1,median))),
                                Inter = c(rep("Intratumor",nrow(intra_m5_metn)),rep("Intertumor",nrow(inter_m5_metn))),
                                Sample = c(rownames(intra_m5_metn),rownames(inter_m5_metn)))
m5_heterogeneity_summaryn <- summarySE(m5_heterogeneityn[,1:2],measurevar = "Heterogeneity", groupvars = "Inter")
m5_heterogeneity_summaryn$median <- c(median(apply(intra_m5_metn,1,median)),median(apply(inter_m5_metn,1,median)))

#plot barplots comparing heterogeneity
m4_heterogeneity_summary_all <- rbind(m4_heterogeneity_summary,m4_heterogeneity_summaryn)
m4_heterogeneity_summary_all$Tumor <- as.factor(c("Tumor Intertumor","Tumor Intratumor","Normal Intertumor","Normal Intratumor"))
m4_heterogeneity_summary_all$Tumor <- factor(m4_heterogeneity_summary_all$Tumor, levels = c("Tumor Intertumor", "Normal Intertumor", "Tumor Intratumor", "Normal Intratumor"))
ggplot(m4_heterogeneity_summary_all, aes(x=Tumor, y=Heterogeneity)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=Heterogeneity-ci, ymax=Heterogeneity+ci), width=.2, position=position_dodge(.9)) +
  xlab("") + theme_Publication(base_size = 7) + ggtitle("M4 Heterogeneity") + ylim(0,1.4)

m5_heterogeneity_summary_all <- rbind(m5_heterogeneity_summary,m5_heterogeneity_summaryn)
m5_heterogeneity_summary_all$Tumor <- as.factor(c("Tumor Intertumor","Tumor Intratumor","Normal Intertumor","Normal Intratumor"))
m5_heterogeneity_summary_all$Tumor <- factor(m5_heterogeneity_summary_all$Tumor, levels = c("Tumor Intertumor", "Normal Intertumor", "Tumor Intratumor", "Normal Intratumor"))
ggplot(m5_heterogeneity_summary_all, aes(x=Tumor, y=Heterogeneity)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=Heterogeneity-ci, ymax=Heterogeneity+ci), width=.2, position=position_dodge(.9)) +
  xlab("") + theme_Publication(base_size = 7) + ggtitle("M5 Heterogeneity") + ylim(0,1.4)


#linear model to compare heterogeneity and metabolite abundance
m4_ith_level_c <- as.numeric(rep(apply(intra_m4_met,1,median),times = m4_mr_patient_count))
m5_ith_level_c <- as.numeric(rep(apply(intra_m5_met,1,median),times = m5_mr_patient_count))

heterogeneity_table <- cbind(m4_mr_tumors,m5_mr_tumors)
heterogeneity_table <- t(heterogeneity_table)
heterogeneity_table <- as.data.frame(heterogeneity_table)
heterogeneity_table$batch <- c(rep("M4",70),rep("M5",53))
heterogeneity_table$heterogeneity <- c(m4_ith_level_c,m5_ith_level_c)

heterogeneity_met_lm <- function(i){
  m <- lm(as.formula(paste(colnames(heterogeneity_table)[i],"~heterogeneity+batch", collapse = "+", sep = "")), data = heterogeneity_table)
  values <- summary(m)$coefficients[2,3:4]
  return(as.numeric(values))
}

heterogeneity_met_lm_values <- sapply(1:602,heterogeneity_met_lm)
heterogeneity_met_lm_values <- t(heterogeneity_met_lm_values)
heterogeneity_met_lm_values <- as.data.frame(heterogeneity_met_lm_values)
rownames(heterogeneity_met_lm_values) <- intersecting_metabolites
colnames(heterogeneity_met_lm_values) <- c("t.statistic","p.value")
heterogeneity_met_lm_values$metabolite <- metanno_m4$BIOCHEMICAL[match(intersecting_metabolites,metanno_m4$COMP_IDstr)]
heterogeneity_met_lm_values$pathway <- metanno_m4$SUB_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)]
heterogeneity_met_lm_values$p.adj <- p.adjust(heterogeneity_met_lm_values$p.value, method = "fdr")
heterogeneity_met_lm_values$Significant <- ifelse(heterogeneity_met_lm_values$p.adj < 0.05,"FDR < 0.05","Not Significant")
heterogeneity_met_lm_values$super <- metanno_m4$SUPER_PATHWAY[match(intersecting_metabolites,metanno_m4$COMP_IDstr)]
heterogeneity_met_lm_values$Pathway <- ifelse(heterogeneity_met_lm_values$p.adj < 0.05, heterogeneity_met_lm_values$super, "Not Significant")

calculate_da <- function(pathways,data){
  da_scores <- as.data.frame(matrix(NA, nrow = length(pathways),ncol = 5))
  colnames(da_scores) <- c("Metabolite_Count","Up","Down","DA_Score","Pathway")
  da_scores$Metabolite_Count <- sapply(1:length(pathways),function(i){dim(data[data$pathway == pathways[i],])[1]})
  da_scores$Up <- sapply(1:length(pathways),function(i){dim(data[data$pathway == pathways[i] & data$t.statistic > 0 & data$Significant == "FDR < 0.05",])[1]}) 
  da_scores$Down <- sapply(1:length(pathways), function(i){dim(data[data$pathway == pathways[i] & data$t.statistic < 0 & data$Significant == "FDR < 0.05",])[1]})
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

ith_pathways_lm <- names(which(sort(table(subset(heterogeneity_met_lm_values,Significant == "FDR < 0.05")[,4]),decreasing = TRUE) > 2)) #only keep pathways with 3+ sig altered metabolites
ith_pvalues_all_pathways_lm <- heterogeneity_met_lm_values[heterogeneity_met_lm_values$pathway %in% ith_pathways_lm,]
ith_pvalues_all_pathways_lm$pathway <- as.factor(ith_pvalues_all_pathways_lm$pathway)
ith_da_scores_lm <- calculate_da(ith_pathways_lm,ith_pvalues_all_pathways_lm)
plot_da(ith_da_scores_lm,"MetITH Pathway Analysis")


#HLA LOH vs heterogeneity
sample_variable_table <- read.csv("~/Documents/RCC Final/Mahdi RCC Table.csv")
sample_variable_table$Pt[6:17] <- gsub("NIVO","NIVOT",sample_variable_table$Pt[6:17])
sample_variable_table <- sample_variable_table[which(sample_variable_table$Pt %in% c(rownames(intra_m4_met),rownames(intra_m5_met))),]
sample_variable_table$MET.ITH <- NA
sample_variable_table$MET.ITH[1:16] <- sapply(1:16,function(i){m4_ith_level[sample_variable_table$Pt[i]]})
sample_variable_table$MET.ITH[17:28] <- sapply(17:28,function(i){m5_ith_level[sample_variable_table$Pt[i]]})
sapply(6:9,function(i){fisher.test(table(sample_variable_table[,i],sample_variable_table$MET.ITH))$p.value}) #1.00000000 1.00000000 0.01275362 0.12011187

sample_variable_table$MET.ITHV <- NA
sample_variable_table$MET.ITHV[1:16] <- sapply(1:16,function(i){m4_ith_value[sample_variable_table$Pt[i]]})
sample_variable_table$MET.ITHV[17:28] <- sapply(17:28,function(i){m5_ith_value[sample_variable_table$Pt[i]]})

loh_heterogeneity <- data.frame(Heterogeneity = sample_variable_table$MET.ITHV,
                                LOH = as.factor(sample_variable_table$HLALOH))
loh_heterogeneity_summary <- summarySE(loh_heterogeneity[,1:2],measurevar = "Heterogeneity", groupvars = "LOH")
loh_heterogeneity_summary$median <- c(median(loh_heterogeneity[loh_heterogeneity$LOH==0,1]),median(loh_heterogeneity[loh_heterogeneity$LOH==1,1]))
ggplot(loh_heterogeneity_summary, aes(x=LOH, y=Heterogeneity)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=Heterogeneity-ci, ymax=Heterogeneity+ci), width=.2, position=position_dodge(.9)) +
  xlab("") + theme_Publication(base_size = 7) + ggtitle("HLA LOH") + ylim(0,1)


#heatmap of driver genes vs metabolic heterogeneity
drivers <- c("VHL","BAP1","PBRM1","SETD2","KDM5C","TP53","PTEN","PIK3CA","MTOR","TSC2")

multiregion_maf$Patient <- str_extract(multiregion_maf$Sample, "[^_]+")
patients <- unique(multiregion_maf$Patient)

drivers_table <- sapply(patients,function(x){sapply(drivers,function(i){nrow(multiregion_maf[multiregion_maf$Patient==x & multiregion_maf$Hugo_Symbol==i,])>0})})


drivers_table_hm <- t(drivers_table)
drivers_table_hm <- drivers_table_hm[-c(4,23,24,26,28),] #remove samples that aren't in the sample variable table
drivers_table_hm <- ifelse(drivers_table_hm == TRUE, 1,0)
drivers_table_hm <- cbind(drivers_table_hm,sample_variable_table[,c(8,9,16)]) #add CDKN2A/B loss and HLA


sample_hm_annotation <- HeatmapAnnotation(Met.ITH = sample_variable_table$MET.ITHV)

Heatmap(t(drivers_table_hm[,1:12]), top_annotation = sample_hm_annotation)


#RNA linear model vs metabolic heterogeneity (looking at Hallmark pathways via GSEA)
rna_keep <- which(!is.na(gsub("-",".",immune_mapping$MetabID[match(colnames(rna_mr_tumors_mc),immune_mapping$RNAID)])))
rna_mr_tumors_mc_keep <- rna_mr_tumors_mc[,rna_keep]
colnames(rna_mr_tumors_mc_keep) <- gsub("-",".",immune_mapping$MetabID[match(colnames(rna_mr_tumors_mc_keep),immune_mapping$RNAID)])
rna_subjectid_keep <- rna_mr_tumors_subjectid[rna_keep]
rna_subjectid_keep[19:66] <- gsub("NIVO","NIVOT",rna_subjectid_keep[19:66])
m4_het_forrna <- apply(intra_m4_met,1,median)[rna_subjectid_keep[1:66]]
m5_het_forrna <- apply(intra_m5_met,1,median)[rna_subjectid_keep[67:112]]

heterogeneity_table_rna <- t(rna_mr_tumors_mc_keep)
heterogeneity_table_rna <- as.data.frame(heterogeneity_table_rna)
heterogeneity_table_rna$batch <- c(rep("M4",66),rep("M5",46))
heterogeneity_table_rna$heterogeneity <- as.numeric(c(m4_het_forrna,m5_het_forrna))


heterogeneity_rna_lm_rna <- function(i){
  m <- lm(as.formula(paste("`",colnames(heterogeneity_table_rna)[i],"`","~heterogeneity+batch", collapse = "+", sep = "")), data = heterogeneity_table_rna)
  values <- summary(m)$coefficients[2,3:4]
  return(as.numeric(values))
}


heterogeneity_met_lm_values_rna <- sapply(1:19570,heterogeneity_rna_lm_rna)
heterogeneity_met_lm_values_rna <- t(heterogeneity_met_lm_values_rna)
heterogeneity_met_lm_values_rna <- as.data.frame(heterogeneity_met_lm_values_rna)
rownames(heterogeneity_met_lm_values_rna) <- rownames(rna_mr_tumors_mc_keep)
colnames(heterogeneity_met_lm_values_rna) <- c("t.statistic","p.value")
write.table(heterogeneity_met_lm_values_rna, file = "~/Documents/RCC Final/RNA ITH High vs Low Metabolite Differences Linear Model.txt", sep = "\t", quote = FALSE)
heterogeneity_met_lm_values_rna <- read.delim("~/Documents/RCC Final/RNA ITH High vs Low Metabolite Differences Linear Model.txt", row.names = 1)

heterogeneity_met_lm_values_rna$p.adj <- p.adjust(heterogeneity_met_lm_values_rna$p.value, method = "BH")


rna_stat<- heterogeneity_met_lm_values_rna$t.statistic
names(rna_stat) <- rownames(heterogeneity_met_lm_values_rna)

fgseaRes_hallmark_cont <- fgseaMultilevel(pathways = hallmark_gmt, 
                                          stats    = rna_stat,
                                          minSize=15,maxSize=500,nPermSimple = 20000)
fgseaRes_hallmark_cont$sig <- ifelse(fgseaRes_hallmark_cont$padj < 0.05, "sig","not sig")
fgseaRes_hallmark_cont$pathway <- gsub("HALLMARK_","",fgseaRes_hallmark_cont$pathway)
fgseaRes_hallmark_cont$pathway <- gsub("_"," ",fgseaRes_hallmark_cont$pathway)
pathway_order <- fgseaRes_hallmark_cont[order(fgseaRes_hallmark_cont$NES),1]
fgseaRes_hallmark_cont$pathway <- factor(fgseaRes_hallmark_cont$pathway, levels = pathway_order$pathway)

fgseaRes_hallmark_cont_sig <- fgseaRes_hallmark_cont[fgseaRes_hallmark_cont$padj < 0.05,]

ggplot(fgseaRes_hallmark_cont_sig,aes(x=NES,y=pathway)) + geom_bar(stat="identity") + theme_Publication()


#metabolic heterogeneity vs immune signatures
m4_immune_w_met <- immune_all[immune_mapping$ITHID[match(gsub("\\.","-",rna_mr_samples[1:66]),immune_mapping$MetabID)],]
m4_immune_w_met <- m4_immune_w_met[,1:94] #remove KEGG pathways

m5_immune_w_met <- immune_all[immune_mapping$ITHID[match(gsub("\\.","-",rna_mr_samples[67:112]),immune_mapping$MetabID)],]
m5_immune_w_met <- m5_immune_w_met[,1:94] #remove KEGG pathways

immune_w_met <- rbind(m4_immune_w_met,m5_immune_w_met)

rna_patient_names_imm <- rna_mapping$Patient.ID[match(rownames(immune_w_met),rna_mapping$Sample)]
rna_patient_names_imm[1:66] <- gsub("NIVO","NIVOT",rna_patient_names_imm[1:66])
m4_het_forimm <- apply(intra_m4_met,1,median)[rna_patient_names_imm[1:66]]
m5_het_forimm <- apply(intra_m5_met,1,median)[rna_patient_names_imm[67:112]]

heterogeneity_table_imm <- immune_w_met
heterogeneity_table_imm <- as.data.frame(heterogeneity_table_imm)
heterogeneity_table_imm$batch <- c(rep("M4",66),rep("M5",46))
heterogeneity_table_imm$heterogeneity <- as.numeric(c(m4_het_forimm,m5_het_forimm))

test <- lm(CYT~heterogeneity+batch, data = heterogeneity_table_imm) 

heterogeneity_imm_lm <- function(i){
  m <- lm(as.formula(paste(colnames(heterogeneity_table_imm)[i],"~heterogeneity+batch", collapse = "+", sep = "")), data = heterogeneity_table_imm)
  values <- summary(m)$coefficients[2,3:4]
  return(as.numeric(values))
}

heterogeneity_imm_lm_values <- sapply(1:94,heterogeneity_imm_lm)
heterogeneity_imm_lm_values <- t(heterogeneity_imm_lm_values)
heterogeneity_imm_lm_values <- as.data.frame(heterogeneity_imm_lm_values)
rownames(heterogeneity_imm_lm_values) <- colnames(immune_w_met)
colnames(heterogeneity_imm_lm_values) <- c("t.statistic","p.value")
heterogeneity_imm_lm_values$p.adj <- p.adjust(heterogeneity_imm_lm_values$p.value, method = "BH")


heterogeneity_imm_lm_values$Significant <- ifelse(heterogeneity_imm_lm_values$p.adj < 0.05,"FDR < 0.05","Not Significant")

heterogeneity_imm_lm_values <- heterogeneity_imm_lm_values[order(heterogeneity_imm_lm_values$p.adj),]
heterogeneity_imm_lm_values$signature <- rownames(heterogeneity_imm_lm_values)
heterogeneity_imm_lm_values$signature <- factor(heterogeneity_imm_lm_values$signature, levels = rev(heterogeneity_imm_lm_values$signature))

ggplot(heterogeneity_imm_lm_values[1:10,], aes(x=-log10(p.adj), y=signature, fill = t.statistic)) + 
  geom_bar(stat="identity") + theme_Publication(base_size = 7) + ylab("") + xlab("-log10(Adjusted P-Value)") + 
  geom_vline(xintercept = 1.30103, color="red", linetype="dotted") #cutoff is 0.05 



#calculate DNA heterogeneity
m4_tumor_clientID <- MetaData_M4$CLIENT_IDENTIFIER[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]
m5_tumor_clientID <- MetaData_M5$CLIENT_IDENTIFIER[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)]
m4_tumor_clientID_new <- gsub("NIVOT","NIVO",m4_tumor_clientID)
m5_tumor_clientID_new <- gsub("-","_",m5_tumor_clientID)
multiregion_maf <- multiregion_maf[multiregion_maf$Sample %in% c(m4_tumor_clientID_new,m5_tumor_clientID_new),]
multiregion_maf$mutation_cat <- paste(multiregion_maf$Hugo_Symbol,"-",multiregion_maf$ENS_HGVS_c)

maf_patients <- gsub("_.*","",names(table(multiregion_maf$Sample)))
length(unique(maf_patients)) #28 patients
maf_patient_count <- table(maf_patients) #sc08 only has 1 sample so remove it; and SC08_TC
multiregion_maf <- multiregion_maf[-which(multiregion_maf$Sample =="SC08_TC"),]
maf_patients <- gsub("_.*","",names(table(multiregion_maf$Sample)))
length(unique(maf_patients)) #27 patients
maf_patient_count <- table(maf_patients)
maf_patient_sum <- cumsum(maf_patient_count)

length(names(table(multiregion_maf$Sample))) #109 samples in total
mutation_combinations <- combn(names(table(multiregion_maf$Sample)),2)
mutation_combinations <- data.frame(mutation_combinations)

mutation_distance <- function(combination){
  x_muts <- multiregion_maf[multiregion_maf$Sample==mutation_combinations[1,combination],129]
  y_muts <- multiregion_maf[multiregion_maf$Sample==mutation_combinations[2,combination],129]
  shared_muts <- length(intersect(x_muts,y_muts))
  total_muts <- length(unique(c(x_muts,y_muts)))
  sim_scores <- shared_muts/total_muts
  distance <- 1-sim_scores
  return(distance)
}

all_mutation_distance <- sapply(1:5886,function(i){mutation_distance(i)})

mutation_distance_matrix <- matrix(0,nrow=109,ncol=109)
colnames(mutation_distance_matrix) <- names(table(multiregion_maf$Sample))
rownames(mutation_distance_matrix) <- names(table(multiregion_maf$Sample))
mutation_distance_matrix[lower.tri(mutation_distance_matrix, diag=FALSE)] <- all_mutation_distance
mutation_distance_matrix <- t(mutation_distance_matrix)
mutation_distance_matrix[lower.tri(mutation_distance_matrix, diag = TRUE)] <- NA



intraheterogeneity_dna <- function(data, npatients, n_p1, sample_sum){ #npatients = number of patients, n_p1 = number of samples in patient 1, #sample sum is vector of cumulative sample count
  
  pat1 <- median(as.numeric(data[1:n_p1,1:n_p1]),na.rm=TRUE)
  rest <- sapply(2:npatients,function(i){median(as.numeric(data[as.numeric(sample_sum[i-1]+1):as.numeric(sample_sum[i]),as.numeric(sample_sum[i-1]+1):as.numeric(sample_sum[i])]),na.rm=TRUE)})
  intra <- c(pat1,rest)
  return(intra)
}
genetic_intra_het <- intraheterogeneity_dna(mutation_distance_matrix, 27, 5, maf_patient_sum)
genetic_intra_het <- t(data.frame(genetic_intra_het))
colnames(genetic_intra_het) <- unique(maf_patients)


#scatterplots comparing  heterogeneity across different modalities
m4_met_order <- names(sort(apply(intra_m4_met,1,median),decreasing = TRUE))
m4_met_order <- gsub("NIVOT","NIVO",m4_met_order)

m5_met_order <- names(sort(apply(intra_m5_met,1,median),decreasing = TRUE))

m4_annotations <- data.frame(Metabolite = apply(intra_m4_met,1,median)[match(names(sort(apply(intra_m4_met,1,median),decreasing = TRUE)),rownames(intra_m4_met))],
                             RNA = apply(intra_rna,1,median)[match(names(sort(apply(intra_m4_met,1,median),decreasing = TRUE)),rownames(intra_rna))],
                             DNA = genetic_intra_het[1,c(match(m4_met_order,colnames(genetic_intra_het)))])
m5_annotations <- data.frame(Metabolite = apply(intra_m5_met,1,median)[match(names(sort(apply(intra_m5_met,1,median),decreasing = TRUE)),rownames(intra_m5_met))],
                             RNA = as.numeric(apply(intra_rna,1,median)[match(names(sort(apply(intra_m5_met,1,median),decreasing = TRUE)),rownames(intra_rna))]),
                             DNA = as.numeric(genetic_intra_het[match(m5_met_order,colnames(genetic_intra_het))]))

ggplot(m4_annotations,aes(x=Metabolite,y=RNA)) + geom_point() + theme_Publication() + geom_smooth(method=lm, se=FALSE) + ggtitle("M4")
ggplot(m5_annotations,aes(x=Metabolite,y=RNA)) + geom_point() + theme_Publication() + geom_smooth(method=lm, se=FALSE) + ggtitle("M5")

ggplot(m4_annotations,aes(x=Metabolite,y=DNA)) + geom_point() + theme_Publication() + geom_smooth(method=lm, se=FALSE) + ggtitle("M4")
ggplot(m5_annotations,aes(x=Metabolite,y=DNA)) + geom_point() + theme_Publication() + geom_smooth(method=lm, se=FALSE) + ggtitle("M5")

ggplot(m4_annotations,aes(x=RNA,y=DNA)) + geom_point() + theme_Publication() + geom_smooth(method=lm, se=FALSE) + ggtitle("M4")
ggplot(m5_annotations,aes(x=RNA,y=DNA)) + geom_point() + theme_Publication() + geom_smooth(method=lm, se=FALSE) + ggtitle("M5")


cor.test(m4_annotations$Metabolite,m4_annotations$RNA, method = "spearman") #rho = 0.60, p = 0.02
cor.test(m4_annotations$Metabolite,m4_annotations$DNA, method = "spearman") #rho = 0.06, p = 0.83

cor.test(m5_annotations$Metabolite,m5_annotations$RNA, method = "spearman") #rho = 0.87, p = 0.01
cor.test(m5_annotations$Metabolite,m5_annotations$DNA, method = "spearman") #rho = 0.41, p = 0.21

cor.test(m4_annotations$DNA,m4_annotations$RNA, method = "spearman") #rho = 0.16, p = 0.55
cor.test(m5_annotations$RNA,m5_annotations$DNA, method = "spearman") #rho = 0.44, p = 0.18











