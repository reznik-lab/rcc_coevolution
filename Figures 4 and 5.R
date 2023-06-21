#####LOAD DATA#####
library(annotate)
library(org.Hs.eg.db)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(DescTools)
library(survcomp)
library(qusage)
library(fgsea)
library(ComplexHeatmap)
library(textshape)
library(lme4)
library(pbkrtest)
library(performance)
library(patchwork)


load("All Data.RData")
load("Outlier Data.RData")

#####CLEAN DATA#####
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


m4_w_immune <- immune_mapping$MetabID[match(rownames(immune_all),immune_mapping$ITHID)]
m4_w_immune <- gsub("-",".",m4_w_immune)
m4_met_w_immune <- m4_mr_tumors[intersecting_metabolites,intersect(m4_w_immune,colnames(m4_mr_tumors))]
m4_immune_w_met <- immune_all[immune_mapping$ITHID[match(gsub("\\.","-",colnames(m4_met_w_immune)),immune_mapping$MetabID)],]
m4_immune_w_met <- m4_immune_w_met[,1:94] #remove KEGG pathways


m5_w_immune <- immune_mapping$MetabID[match(rownames(immune_all),immune_mapping$ITHID)]
m5_w_immune <- gsub("-",".",m5_w_immune)
m5_met_w_immune <- m5_mr_tumors[intersecting_metabolites,intersect(m5_w_immune, colnames(m5_mr_tumors))]
m5_immune_w_met <- immune_all[immune_mapping$ITHID[match(gsub("\\.","-",colnames(m5_met_w_immune)),immune_mapping$MetabID)],]
m5_immune_w_met <- m5_immune_w_met[,1:94] #remove KEGG pathways


m4_rna_w_met <- lfiltn_rna[,colnames(m4_met_w_immune)]
m5_rna_w_met <- lfiltn_rna[,colnames(m5_met_w_immune)]

#####PREPARE FOR BOOTSTRAPPING ANALYSIS######
is_table <- rbind(t(m4_met_w_immune_mc[intersecting_metabolites,]),t(m5_met_w_immune_mc[intersecting_metabolites,]))
is_table <- cbind(is_table,c(m4_immune_w_met[,3],m5_immune_w_met[,3]),
                  c(m4_immune_w_met[,59],m5_immune_w_met[,59]),
                  c(m4_immune_w_met[,60],m5_immune_w_met[,60]),
                  c(m4_immune_w_met[,53],m5_immune_w_met[,53]),
                  c(m4_immune_w_met[,46],m5_immune_w_met[,46]),
                  c(m4_immune_w_met[,57],m5_immune_w_met[,57]))
colnames(is_table)[603:608] <- c("ImmuneScore","Angiogenesis","Myeloid","JAVELIN","Dendritic","Macrophage")
is_table <- as.data.frame(is_table)
is_table$Batch <- c(rep("M4",66),rep("M5",46))
is_table$Patient <- c(MetaData_M4$SUBJECT_ID[match(colnames(m4_met_w_immune),MetaData_M4$SAMPLE_NAME)],
                      MetaData_M5$SUBJECT_ID[match(colnames(m5_met_w_immune),MetaData_M5$SAMPLE_NAME)])

#####RUN BOOTSTRAPPING --> THIS SHOULD BE DONE ON CLUSTER --> RUN 602 ARGS AS EACH JOB WILL CALCULATE THE P-VALUE FOR A SINGLE METABOLITE#####
args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)

library(lme4)
library(pbkrtest)

is_table <- read.csv("Signature Bootstrapping Table.txt", sep="")

p_calculation_ag <- function(i){
  m <- lmer(as.formula(paste(colnames(is_table)[i],"~Angiogenesis+(1|Patient)", collapse = "+", sep = "")), data = is_table)
  m1 <- update( m, .~. - Angiogenesis)
  pb.b <- PBmodcomp(m, m1, nsim=100000)
  return(pb.b$test$p.value[2]) #returns the fraction of simulated LRT-values that are larger or equal to observed LRT value
}

p_calculation_is <- function(i){
  m <- lmer(as.formula(paste(colnames(is_table)[i],"~ImmuneScore+(1|Patient)", collapse = "+", sep = "")), data = is_table)
  m1 <- update( m, .~. - ImmuneScore)
  pb.b <- PBmodcomp(m, m1, nsim=100000)
  return(pb.b$test$p.value[2]) #returns the fraction of simulated LRT-values that are larger or equal to observed LRT value
}

p_calculation_my <- function(i){
  m <- lmer(as.formula(paste(colnames(is_table)[i],"~Myeloid+(1|Patient)", collapse = "+", sep = "")), data = is_table)
  m1 <- update( m, .~. - Myeloid)
  pb.b <- PBmodcomp(m, m1, nsim=100000)
  return(pb.b$test$p.value[2]) #returns the fraction of simulated LRT-values that are larger or equal to observed LRT value
}

p_calculation_jv <- function(i){
  m <- lmer(as.formula(paste(colnames(is_table)[i],"~JAVELIN+(1|Patient)", collapse = "+", sep = "")), data = is_table)
  m1 <- update( m, .~. - Dendritic)
  pb.b <- PBmodcomp(m, m1, nsim=100000)
  return(pb.b$test$p.value[2]) #returns the fraction of simulated LRT-values that are larger or equal to observed LRT value
}

p_calculation_dc <- function(i){
  m <- lmer(as.formula(paste(colnames(is_table)[i],"~Dendritic+(1|Patient)", collapse = "+", sep = "")), data = is_table)
  m1 <- update( m, .~. - JAVELIN)
  pb.b <- PBmodcomp(m, m1, nsim=100000)
  return(pb.b$test$p.value[2]) #returns the fraction of simulated LRT-values that are larger or equal to observed LRT value
}

p_calculation_mac <- function(i){
  m <- lmer(as.formula(paste(colnames(is_table)[i],"~Macrophage+(1|Patient)", collapse = "+", sep = "")), data = is_table)
  m1 <- update( m, .~. - Macrophage)
  pb.b <- PBmodcomp(m, m1, nsim=100000)
  return(pb.b$test$p.value[2]) #returns the fraction of simulated LRT-values that are larger or equal to observed LRT value
}


null_model_is <- data.frame(is = p_calculation_is(args),
                            ag = p_calculation_ag(args),
                            my = p_calculation_my(args),
                            jv = p_calculation_jv(args),
                            dc = p_calculation_dc(args),
                            mac = p_calculation_mac(args))


write.table(null_model_is,paste(args,"Null Model P-Values.txt"), quote = FALSE, row.names = FALSE, append = FALSE, sep = "\t")

#####ANALYZE THE BOOTSTRAPPING RESULTS#####
#final results of bootstrapping on cluster is saved in a dataframe called pb_results
pb_results$is_adj <- p.adjust(pb_results$is, method = "BH")
pb_results$ag_adj <- p.adjust(pb_results$ag, method = "BH")
pb_results$my_adj <- p.adjust(pb_results$my, method = "BH")
pb_results$jv_adj <- p.adjust(pb_results$jv, method = "BH")
pb_results$dc_adj <- p.adjust(pb_results$dc, method = "BH")
pb_results$mac_adj <- p.adjust(pb_results$mac, method = "BH")


rownames(pb_results) <- pb_results$metID
pb_results <- pb_results[intersecting_metabolites,]

#####FUNCTION TO CALCULATE MARGINAL EXPLAINED VARIANCE#####
signature_lmer <- function(column){
  sig_table <- rbind(t(m4_met_w_immune_mc[intersecting_metabolites,]),t(m5_met_w_immune_mc[intersecting_metabolites,]))
  sig_table <- cbind(sig_table,c(m4_immune_w_met[,column],m5_immune_w_met[,column]))
  colnames(sig_table)[603] <- c("SIG")
  sig_table <- as.data.frame(sig_table)
  sig_table$Batch <- c(rep("M4",66),rep("M5",46))
  sig_table$Patient <- c(MetaData_M4$SUBJECT_ID[match(colnames(m4_met_w_immune),MetaData_M4$SAMPLE_NAME)],
                         MetaData_M5$SUBJECT_ID[match(colnames(m5_met_w_immune),MetaData_M5$SAMPLE_NAME)])
  sig_table$PC1 <- pca_mc$x[,1]
  
  
  r2_calculation_sig <- function(i){
    m <- lmer(as.formula(paste(colnames(sig_table)[i],"~SIG+(1|Patient)", collapse = "+", sep = "")), data = sig_table)
    r2m <- r2(m)
    return(as.numeric(r2m[[2]]))
  }
  
  sig_r2 <- sapply(1:602,r2_calculation_sig)
  sig_r2 <- data.frame(metID = colnames(sig_table)[1:602],
                       metabolite = metanno_m4$BIOCHEMICAL[match(colnames(sig_table)[1:602],metanno_m4$COMP_IDstr)],
                       pathway = metanno_m4$SUB_PATHWAY[match(colnames(sig_table)[1:602],metanno_m4$COMP_IDstr)],
                       r2 = sig_r2)
  sig_r2$rank <- rank(sig_r2$r2)
  
  
  pc1_pred <- lmer(PC1~SIG+(1|Patient), data=sig_table) 
  print(r2(pc1_pred))
  
  return(sig_r2)
  
}



#####IMMUNESCORE ANALYSIS#####
is_r2 <- signature_lmer(3)
is_r2$p.value <- pb_results$is_adj
is_r2$Significant <- c(ifelse(is_r2$p.value < 0.05, "Significant","Not Significant"))

is_labels <- c("glutathione, oxidized (GSSG)","quinolinate","cysteine","nicotinamide","N1-Methyl-2-pyridone-5-carboxamide")


ggplot(is_r2, aes(x=rank, y= r2, color = Significant)) + geom_point() + 
  geom_point(data = subset(is_r2, Significant == 'Significant'),
             aes(x = rank, y = r2, color = Significant)) +
  theme_Publication() + xlab("") + ylab("Marginal Explained Variance") +
  geom_text_repel(data=subset(is_r2, metabolite %in% is_labels), aes(label = metabolite), size = 1, color = "black", max.overlaps = 15) + 
  scale_color_manual(values=c("#000000","#FF0000"))



#pathway analysis
is_sig_pathways <- data.frame(pathway = names(table(pb_results[pb_results$is_adj < 0.05,8])),
                              significant = as.numeric(table(pb_results[pb_results$is_adj < 0.05,8])))
is_sig_pathways$not_sig <- sapply(is_sig_pathways$pathway,function(i){nrow(pb_results[pb_results$pathway==i,])})-is_sig_pathways$significant
is_sig_pathways$sig_other <- 33-is_sig_pathways$significant
is_sig_pathways$not_sig_other <- 569-is_sig_pathways$not_sig
is_sig_pathways$size <- is_sig_pathways$significant+is_sig_pathways$not_sig

is_sig_pathways_3 <- is_sig_pathways[is_sig_pathways$size >=3,]

is_fisher <- data.frame(pathway = is_sig_pathways_3$pathway,
                        fisher_p = sapply(1:19,function(i){fisher.test(matrix(as.numeric(is_sig_pathways_3[i,2:5]),nrow=2))$p.value}))
is_fisher$adj <- p.adjust(is_fisher$fisher_p, method = "BH")
is_fisher <- is_fisher[-c(which(is_fisher$adj==1)),]
is_fisher$pathway <- factor(is_fisher$pathway, levels = is_fisher[order(is_fisher$adj, decreasing = TRUE),1])

ggplot(is_fisher, aes(x = pathway, y = -log10(adj))) + geom_bar(stat = "identity") + coord_flip() + theme_Publication() +
  xlab("") + ylab("-Log 10 FDR-Adjusted P-Value")

#nicotinamide and quinolinate boxplots
#####QUINOLINATE AND NICOTINAMIDE BOXPLOTS FIG 5E#####
#immunescore vs quinolinate and nicotinamide for every patient, scaled
met_is <- data.frame(IS = c(m4_immune_w_met[,3],m5_immune_w_met[,3]),
                     Quinolinate = as.numeric(c(m4_met_w_immune_mc["M01899",],m5_met_w_immune_mc["M01899",])),
                     Nicotinamide = as.numeric(c(m4_met_w_immune_mc["M00594",],m5_met_w_immune_mc["M00594",])),
                     Patient = c(MetaData_M4$SUBJECT_ID[match(colnames(m4_met_w_immune),MetaData_M4$SAMPLE_NAME)],
                                 MetaData_M5$SUBJECT_ID[match(colnames(m5_met_w_immune),MetaData_M5$SAMPLE_NAME)]))
met_is$is_rank <- unlist(sapply(1:28,function(i){rank(met_is[met_is$Patient==unique(met_is$Patient)[i],1])/nrow(met_is[met_is$Patient==unique(met_is$Patient)[i],])}))
met_is$q_rank <- unlist(sapply(1:28,function(i){rank(met_is[met_is$Patient==unique(met_is$Patient)[i],2])/nrow(met_is[met_is$Patient==unique(met_is$Patient)[i],])}))
met_is$n_rank <- unlist(sapply(1:28,function(i){rank(met_is[met_is$Patient==unique(met_is$Patient)[i],3])/nrow(met_is[met_is$Patient==unique(met_is$Patient)[i],])}))


is_max <- sapply(1:28,function(i){which(met_is[met_is$Patient==unique(met_is$Patient)[i],1]==max(met_is[met_is$Patient==unique(met_is$Patient)[i],1]))})
is_min <- sapply(1:28,function(i){which(met_is[met_is$Patient==unique(met_is$Patient)[i],1]==min(met_is[met_is$Patient==unique(met_is$Patient)[i],1]))})

q_max <- sapply(1:28,function(i){which(met_is[met_is$Patient==unique(met_is$Patient)[i],2]==max(met_is[met_is$Patient==unique(met_is$Patient)[i],2]))})
q_min <- sapply(1:28,function(i){which(met_is[met_is$Patient==unique(met_is$Patient)[i],2]==min(met_is[met_is$Patient==unique(met_is$Patient)[i],2]))})

n_max <- sapply(1:28,function(i){which(met_is[met_is$Patient==unique(met_is$Patient)[i],3]==max(met_is[met_is$Patient==unique(met_is$Patient)[i],3]))})
n_min <- sapply(1:28,function(i){which(met_is[met_is$Patient==unique(met_is$Patient)[i],3]==min(met_is[met_is$Patient==unique(met_is$Patient)[i],3]))})


met_is_table <- data.frame(q_is_high = sapply(1:28,function(i){met_is[met_is$Patient==unique(met_is$Patient)[i],2][is_max[[i]]]}),
                           q_is_low = sapply(1:28,function(i){met_is[met_is$Patient==unique(met_is$Patient)[i],2][is_min[[i]]]}),
                           n_is_high = sapply(1:28,function(i){met_is[met_is$Patient==unique(met_is$Patient)[i],3][is_max[[i]]]}),
                           n_is_low = sapply(1:28,function(i){met_is[met_is$Patient==unique(met_is$Patient)[i],3][is_min[[i]]]}),
                           sample = unique(met_is$Patient))
met_is_table$color_q <- c("red")
met_is_table$color_q[which(met_is_table[,1]<met_is_table[,2])] <- c("blue")
met_is_table$color_n <- c("red")
met_is_table$color_n[which(met_is_table[,3]<met_is_table[,4])] <- c("blue")
met_is_table_melt <- melt(met_is_table)
met_is_table_melt$met <- c(rep("q",56),rep("n",56))
met_is_table_melt$variable <- factor(met_is_table_melt$variable, levels = c("q_is_low","n_is_low","q_is_high","n_is_high"))

ggplot(met_is_table_melt[met_is_table_melt$met=="q",],aes(x=variable,y=value)) + 
  geom_boxplot() +
  geom_point(size=2, alpha=0.5,aes(fill=color_q)) +
  geom_line(aes(group=sample), colour="red", linetype="11") + 
  theme_Publication() + xlab("ImmuneScore Level") + ylab("Quinolinate Expression")

ggplot(met_is_table_melt[met_is_table_melt$met=="n",],aes(x=variable,y=value)) + 
  geom_boxplot() +
  geom_point(size=2, alpha=0.5,aes(fill=color_n)) +
  geom_line(aes(group=sample), colour="red", linetype="11") + 
  theme_Publication() + xlab("ImmuneScore Level") + ylab("Nicotinamide Expression")

wilcox.test(as.numeric(met_is_table_melt[met_is_table_melt$met=="q",5])~as.factor(met_is_table_melt[met_is_table_melt$met=="q",4]), paired= TRUE)
wilcox.test(as.numeric(met_is_table_melt[met_is_table_melt$met=="n",5])~as.factor(met_is_table_melt[met_is_table_melt$met=="n",4]), paired= TRUE)


#####COMPARE RNA SIGNATURES TO IF STAINING#####
if_results$Myeloid <- if_results_new$Myeloid_Score[match(if_results$TUMOR_ID,if_results_new$New_ID)]

qplot(if_results$CD31..um2.,if_results$Angio_Score) + theme_Publication() + xlab("CD31 IF") + ylab("Angiogenesis Signature")
qplot(if_results$CD68..cells.mm2.,if_results$Myeloid) + theme_Publication() + xlab("CD68 IF") + ylab("Myeloid Signature")


cor.test(if_results$CD31..um2.,if_results$Angio_Score,method = "spearman") #rho = 0.63, p = 0.000315
cor.test(if_results$CD68..cells.mm2.,if_results$Myeloid,method = "spearman") #rho = 0.58, p = 0.001315


#####SIGNATURES IN FIGURE 6 PLOTS#####
jav_r2_bind <- signature_lmer(53)
jav_r2_bind$p.value <- pb_results$jv_adj
jav_r2_bind$Significant <- c(ifelse(jav_r2_bind$p.value < 0.05, "Significant","Not Significant"))
ggplot(jav_r2_bind, aes(x=rank, y= r2, color = Significant)) + geom_point() + theme_Publication() + xlab("") + ylab("Marginal Explained Variance") +
  geom_text_repel(data=subset(jav_r2_bind,Significant == "Significant"), aes(label = metabolite), size = 2, color = "black") + scale_color_manual(values=c("#000000","#FF0000"))


ang_r2_bind <- signature_lmer(59)
ang_r2_bind$p.value <- pb_results$ag_adj
ang_r2_bind$Significant <- c(ifelse(ang_r2_bind$p.value < 0.05, "Significant","Not Significant"))
ggplot(ang_r2_bind, aes(x=rank, y= r2, color = Significant)) + geom_point() + theme_Publication() + xlab("") + ylab("Marginal Explained Variance") +
  geom_text_repel(data=subset(ang_r2_bind,Significant == "Significant"), aes(label = metabolite), size = 2, color = "black") + scale_color_manual(values=c("#000000","#FF0000"))

mye_r2_bind <- signature_lmer(60) 
mye_r2_bind$p.value <- pb_results$my_adj
mye_r2_bind$Significant <- c(ifelse(mye_r2_bind$p.value < 0.05, "Significant","Not Significant"))

myeloid_mets <- c("N-acetylalanine","N-acetylarginine","N-acetylasparagine","N-acetylaspartate (NAA)","N-acetylglutamine",
                  "N-acetylleucine","N-acetylmethionine","kynurenine","1-methylnicotinamide", "saccharopine","N-palmitoyl-sphingadienine (d18:2/16:0)*")

ggplot(mye_r2_bind, aes(x=rank, y= r2, color = Significant)) + geom_point() + 
  geom_point(data = subset(mye_r2_bind, Significant == 'Significant'),
             aes(x = rank, y = r2, color = Significant)) +
  theme_Publication() + xlab("") + ylab("Marginal Explained Variance") +
  geom_text_repel(data=subset(mye_r2_bind,metabolite %in% myeloid_mets), aes(label = metabolite), size = 2, color = "black") + 
  scale_color_manual(values=c("#000000","#FF0000"))

mac_lmer <- signature_lmer(57) 
mac_lmer$p.value <- pb_results$mac.adj
mac_lmer$Significant <- c(ifelse(mac_lmer$p.value < 0.05, "Significant","Not Significant"))
mac_mets <- c("M52634","M52677","M15497")
ggplot(mac_lmer, aes(x=rank, y= r2, color = Significant)) + geom_point() + 
  geom_point(data = subset(mac_lmer, Significant == 'Significant'),
             aes(x = rank, y = r2, color = Significant)) +
  theme_Publication() + xlab("") + ylab("Marginal Explained Variance") +
  geom_text_repel(data=subset(mac_lmer,metID %in% c(mac_mets,naa) & Significant == "Significant"), aes(label = metabolite), size = 2, color = "black") + 
  scale_color_manual(values=c("#000000","#FF0000"))

dc_lmer <- signature_lmer(46) 
dc_lmer$p.value <- pb_results$dc.adj
dc_lmer$Significant <- c(ifelse(dc_lmer$p.value < 0.05, "Significant","Not Significant"))
dc_mets <- c("M52634","M62100","M27738")
ggplot(dc_lmer, aes(x=rank, y= r2, color = Significant)) + geom_point() + 
  geom_point(data = subset(dc_lmer, Significant == 'Significant'),
             aes(x = rank, y = r2, color = Significant)) +
  theme_Publication() + xlab("") + ylab("Marginal Explained Variance") +
  geom_text_repel(data=subset(dc_lmer,metID %in% c(dc_mets,naa) & Significant == "Significant"), aes(label = metabolite), size = 2, color = "black") + 
  scale_color_manual(values=c("#000000","#FF0000"))


#####CALCULATE CORRELATION RHO BETWEEN N-ACETYL AMINO ACIDS AND NAT GENE EXPRESSION#####
####THIS IS ALL PATIENT CENTERED DATA#####


which(colnames(lfilt_rna)=="MSKC.00349")

rna_m4 <- which(MetaData_M4$TISSUE_STATUS[match(colnames(lfilt_rna)[1:114],MetaData_M4$SAMPLE_NAME)]=="TISSUE_TUMOR")
rna_m4 <- lfilt_rna[,rna_m4]

m4_patient_all_ids <- MetaData_M4$SUBJECT_ID[match(colnames(rna_m4),MetaData_M4$SAMPLE_NAME)]
m4_patients <- unique(m4_patient_all_ids)
m4_all_sample_count <- sapply(m4_patients,function(x)sum(m4_patient_all_ids==x))
m4_all_sample_sum <- cumsum(m4_all_sample_count)



m4_patient_centered_all <- as.data.frame(matrix(NA,nrow=nrow(rna_m4),ncol = ncol(rna_m4)-m4_all_sample_count[1]))
for(x in 1:nrow(rna_m4)){
  m4_patient_centered_all[x,] <- unlist(sapply(2:length(m4_patients),function(i){rna_m4[x,as.numeric(m4_all_sample_sum[i-1]+1):m4_all_sample_sum[i]]-median(as.numeric(rna_m4[x,as.numeric(m4_all_sample_sum[i-1]+1):m4_all_sample_sum[i]]))}))
}
m4_patient_centered_all_p1 <- as.data.frame(matrix(NA,nrow=nrow(rna_m4),ncol = m4_all_sample_sum[1]))
for(x in 1:nrow(rna_m4)){
  m4_patient_centered_all_p1[x,] <- rna_m4[x,1:m4_all_sample_sum[1]]-median(as.numeric(rna_m4[x,1:m4_all_sample_sum[1]]))
}
m4_patient_centered_all <- cbind(m4_patient_centered_all_p1,m4_patient_centered_all)
colnames(m4_patient_centered_all) <- colnames(rna_m4)
rownames(m4_patient_centered_all) <- rownames(rna_m4)
rm(m4_patient_centered_all_p1)


m4_mr_tumors_c <- m4_mr_tumors_c[,colnames(m4_patient_centered_all)]

int <- intersect(colnames(m4_mr_tumors_c),colnames(m4_patient_centered_all))


m4_mr_tumors_c <- m4_mr_tumors_c[,int]
m4_patient_centered_all <- m4_patient_centered_all[,int]

gene_names <- as.character(getSYMBOL(rownames(m4_patient_centered_all), data='org.Hs.eg'))

met_names <- metanno_m4$BIOCHEMICAL[match(intersecting_metabolites,metanno_m4$COMP_IDstr)]
naa <- met_names[grep("N-acetyl",met_names)]
not_naa <- c(1,2,3,4,5,6,12,13,14,21,22,23)
naa <- naa[-not_naa]
naa <- metanno_m4$COMP_IDstr[match(naa,metanno_m4$BIOCHEMICAL)]

nat <- gene_names[grep("NAT",gene_names)]
not_nat <- c(2,3,4,7,8,10,11,13,14,18,19)
nat <- nat[-not_nat]

gene_id <- which(gene_names%in%nat)
gene_id <- rownames(m4_patient_centered_all)[gene_id]

naa_exp4 <- m4_mr_tumors_c[naa,]
nat_exp4 <- m4_patient_centered_all[gene_id,]



naa_nat_corp <- sapply(
  data.frame(t(naa_exp)),
  function(x) Map(function(a,b) cor.test(a,b,method = "spearman")$p.value,
                  list(x),
                  as.data.frame(t(nat_exp))))

naa_nat_corp <- as.data.frame(matrix(as.numeric(naa_nat_corp), nrow = nrow(naa_nat_corp), ncol = ncol(naa_nat_corp), byrow = FALSE))
rownames(naa_nat_corp) <- getSYMBOL(rownames(nat_exp), data='org.Hs.eg')
colnames(naa_nat_corp) <- rownames(naa_exp)


#repeat in m5
rna_m5 <- which(MetaData_M5$SAMPLE_DESCRIPTION[match(colnames(lfilt_rna)[115:190],MetaData_M5$SAMPLE_NAME)]=="TISSUE_TUMOR")
rna_m5_df <- lfilt_rna[,115:190]
rna_m5 <- rna_m5_df[,rna_m5]

m5_patient_all_ids <- MetaData_M5$SUBJECT_ID[match(colnames(rna_m5),MetaData_M5$SAMPLE_NAME)]
m5_patients <- unique(m5_patient_all_ids)
m5_all_sample_count <- sapply(m5_patients,function(x)sum(m5_patient_all_ids==x))
m5_all_sample_sum <- cumsum(m5_all_sample_count)



m5_patient_centered_all <- as.data.frame(matrix(NA,nrow=nrow(rna_m5),ncol = ncol(rna_m5)-m5_all_sample_count[1]))
for(x in 1:nrow(rna_m5)){
  m5_patient_centered_all[x,] <- unlist(sapply(2:length(m5_patients),function(i){rna_m5[x,as.numeric(m5_all_sample_sum[i-1]+1):m5_all_sample_sum[i]]-median(as.numeric(rna_m5[x,as.numeric(m5_all_sample_sum[i-1]+1):m5_all_sample_sum[i]]))}))
}
m5_patient_centered_all_p1 <- as.data.frame(matrix(NA,nrow=nrow(rna_m5),ncol = m5_all_sample_sum[1]))
for(x in 1:nrow(rna_m5)){
  m5_patient_centered_all_p1[x,] <- rna_m5[x,1:m5_all_sample_sum[1]]-median(as.numeric(rna_m5[x,1:m5_all_sample_sum[1]]))
}
m5_patient_centered_all <- cbind(m5_patient_centered_all_p1,m5_patient_centered_all)
colnames(m5_patient_centered_all) <- colnames(rna_m5)
rownames(m5_patient_centered_all) <- rownames(rna_m5)
rm(m5_patient_centered_all_p1)


m5_mr_tumors_c <- m5_mr_tumors_c[,colnames(m5_patient_centered_all)]

int <- intersect(colnames(m5_mr_tumors_c),colnames(m5_patient_centered_all))


m5_mr_tumors_c <- m5_mr_tumors_c[,int]
m5_patient_centered_all <- m5_patient_centered_all[,int]

gene_names5 <- as.character(getSYMBOL(rownames(m5_patient_centered_all), data='org.Hs.eg'))
nat <- gene_names5[grep("NAT",gene_names)]
not_nat <- c(2,3,4,7,8,10,11,13,14,18,19)
nat <- nat[-not_nat]


gene_id <- rownames(m5_patient_centered_all)[gene_id]
gene_id <- which(gene_names%in%nat)
gene_id <- rownames(m5_patient_centered_all)[gene_id]


naa_exp5 <- m5_mr_tumors_c[naa,]
nat_exp5 <- m5_patient_centered_all[gene_id,]



naa_nat_corp5 <- sapply(
  data.frame(t(naa_exp5)),
  function(x) Map(function(a,b) cor.test(a,b,method = "spearman")$p.value,
                  list(x),
                  as.data.frame(t(nat_exp5))))

naa_nat_corp5 <- as.data.frame(matrix(as.numeric(naa_nat_corp5), nrow = nrow(naa_nat_corp5), ncol = ncol(naa_nat_corp5), byrow = FALSE))
rownames(naa_nat_corp5) <- getSYMBOL(rownames(nat_exp5), data='org.Hs.eg')
colnames(naa_nat_corp5) <- rownames(naa_exp5)



naa_nat_cor4 <- sapply(
  data.frame(t(naa_exp4)),
  function(x) Map(function(a,b) cor.test(a,b,method = "spearman")$estimate,
                  list(x),
                  as.data.frame(t(nat_exp4))))

naa_nat_cor4 <- as.data.frame(matrix(as.numeric(naa_nat_cor4), nrow = nrow(naa_nat_cor4), ncol = ncol(naa_nat_cor4), byrow = FALSE))
rownames(naa_nat_cor4) <- getSYMBOL(rownames(nat_exp4), data='org.Hs.eg')
colnames(naa_nat_cor4) <- rownames(naa_exp4)

naa_nat_cor5 <- sapply(
  data.frame(t(naa_exp5)),
  function(x) Map(function(a,b) cor.test(a,b,method = "spearman")$estimate,
                  list(x),
                  as.data.frame(t(nat_exp5))))

naa_nat_cor5 <- as.data.frame(matrix(as.numeric(naa_nat_cor5), nrow = nrow(naa_nat_cor5), ncol = ncol(naa_nat_cor5), byrow = FALSE))
rownames(naa_nat_cor5) <- getSYMBOL(rownames(nat_exp5), data='org.Hs.eg')
colnames(naa_nat_cor5) <- rownames(naa_exp5)

naa_nat_cor4Z <- FisherZ(naa_nat_cor4)
naa_nat_cor5Z <- FisherZ(naa_nat_cor5)

naa_nat_corZ <- (naa_nat_cor4Z+naa_nat_cor5Z)/2
naa_nat_corR <- FisherZInv(naa_nat_corZ)

colnames(naa_nat_corR) <- metanno_m4$BIOCHEMICAL[match(colnames(naa_nat_corR),metanno_m4$COMP_IDstr)]

Heatmap(naa_nat_corR)










