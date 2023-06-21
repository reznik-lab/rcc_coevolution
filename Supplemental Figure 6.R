load("All Data.RData")
load("Outlier Data.RData")

######COMPARE TMB IN BATCH 1 VS BATCH 2######
wes_batch$Sample <- gsub("-","_",wes_batch$Sample)
wes_batch <- wes_batch[-which(is.na(match(wes_batch$Sample,multiregion_maf$Sample))),]

samples <- wes_batch$Sample

n_mut <- sapply(1:142,function(i){nrow(multiregion_maf[multiregion_maf$Sample==samples[i],])})

wes_batch$n_mut <- n_mut

wilcox.test(wes_batch$n_mut~wes_batch$Shipping.Batch) #p-value = 0.9132
ggplot(wes_batch, aes(x=Shipping.Batch,y=n_mut)) + geom_boxplot() + theme_Publication() + 
  xlab("Batch") + ylab("Number of Mutations Per Sample") +
  ggtitle("TMB By Batch")


#PCA to see if there is batch effect in RNA data
rna_batch$Sample.ID <- gsub("-RNA","",rna_batch$Sample.ID)

pca <- prcomp(t(lfiltn_rna_all), center = TRUE, scale. = TRUE)
pca_x <- pca$x
pca_df <- data.frame(sample = rownames(pca_x),
                     PC1 = pca_x[,1],
                     PC2 = pca_x[,2])
pca_df$plate <- rna_batch$Library.Prep.Plate[match(pca_df$sample,rna_batch$Sample.ID)]
pca_df$cell <- rna_batch$Sequencing.Flow.Cell.ID[match(pca_df$sample,rna_batch$Sample.ID)]

ggplot(pca_df, aes(x=PC1,y=PC2,col=plate)) + geom_point(size=0.6) + theme_Publication() + ggtitle("By Library Prep Plate")
ggsave("RNA PCA By Library Prep Plate.pdf",width = 2.5,height = 3,units = "in")
ggplot(pca_df, aes(x=PC1,y=PC2,col=cell)) + geom_point(size=0.6) + theme_Publication() + ggtitle("By Sequencing Flow Cell ID")
ggsave("RNA PCA By Sequencing Flow Cell ID.pdf",width = 2.5,height = 2.933,units = "in")





####DOWNSAMPLING ANALYSIS#####
#downsample to find if sample number is associated with genetic heterogeneity; use all patients with 2+ regions; the trees were with 3+ regions
m4_tumor_clientID <- MetaData_M4$CLIENT_IDENTIFIER[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]
m5_tumor_clientID <- MetaData_M5$CLIENT_IDENTIFIER[match(colnames(m5_mr_tumors),MetaData_M5$SAMPLE_NAME)]

m4_tumor_clientID_new <- gsub("NIVOT","NIVO",m4_tumor_clientID)
m5_tumor_clientID_new <- gsub("-","_",m5_tumor_clientID)

multiregion_maf <- multiregion_maf[multiregion_maf$Sample %in% c(m4_tumor_clientID_new,m5_tumor_clientID_new),]
multiregion_maf$mutation_cat <- paste(multiregion_maf$Hugo_Symbol,"-",multiregion_maf$ENS_HGVS_c)

maf_patients <- gsub("_.*","",names(table(multiregion_maf$Sample)))
length(unique(maf_patients)) #28 patients
maf_patient_count <- table(maf_patients) #sc08 only has 1 sample so remove them; SC08_TC
multiregion_maf <- multiregion_maf[-which(multiregion_maf$Sample %in% c("SC08_TC")),]
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

downsample <- function(patient,n_samples){
  columns <- sample(grep(patient,colnames(mutation_distance_matrix)),n_samples)
  median(as.numeric(mutation_distance_matrix[columns,columns]),na.rm=TRUE)
}

patients <- unique(gsub("_.*","",names(table(multiregion_maf$Sample))))

downsample_all <- function(n_regions){
  patient <- patients[which(table(gsub("_.*","",names(table(multiregion_maf$Sample)))) >= n_regions)]
  as.numeric(sapply(patient,function(i){replicate(10,downsample(i,n_regions))}))
}

downsample_df <- data.frame(n_regions = as.factor(c(rep(2,270),rep(3,230),rep(4,190),rep(5,130))),
                            heterogeneity = unlist(sapply(2:5,function(i){downsample_all(i)})))
ggplot(downsample_df, aes(x=n_regions,y=heterogeneity)) + geom_boxplot() + theme_Publication() + ylab("DNA Heterogeneity") + 
  xlab("Number of Samples") + ggtitle("Spearman p-value = 0.72")
ggsave("DNA Heterogeneity N Samples Boxplot Downsample.pdf", height = 2, width = 2, units = "in")



intraheterogeneity_matrix<- function(data, npatients){
  intra <- as.data.frame(matrix(NA,ncol=nrow(data),nrow=npatients))
  for(x in 1:nrow(data)){
    distance_matrix <- rdist(as.numeric(data[x,]))
    distance_matrix[lower.tri(distance_matrix, diag = TRUE)] <- NA}
  return(distance_matrix)
}
m4_distance_matrix <- intraheterogeneity_matrix(m4_mr_tumors_mc, 16)
colnames(m4_distance_matrix) <- m4_tumor_clientID
m5_distance_matrix <- intraheterogeneity_matrix(m5_mr_tumors_mc, 14)
colnames(m5_distance_matrix) <- m5_tumor_clientID

patients_mets <- str_extract(c(m4_tumor_clientID,m5_tumor_clientID), "[^_]+")
patients_mets <- str_extract(patients_mets, "[^-]+")

m4_patients <- unique(patients_mets)[1:16]
m5_patients <- unique(patients_mets)[17:30]

downsample_met <- function(patient,n_samples){
  data <- NA
  ifelse(patient %in% m4_patients, data <- m4_distance_matrix, data <- m5_distance_matrix)
  columns <- sample(grep(patient,colnames(data)),n_samples)
  median(as.numeric(data[columns,columns]),na.rm=TRUE)
}

downsample_all_met <- function(n_regions){
  patient <- names(which(table(patients_mets) >= n_regions))
  as.numeric(sapply(patient,function(i){replicate(10,downsample_met(i,n_regions))}))
}

downsample_df_met <- data.frame(n_regions = as.factor(c(rep(2,300),rep(3,260),rep(4,220),rep(5,150))),
                                heterogeneity = unlist(sapply(2:5,function(i){downsample_all_met(i)})))
cor.test(as.numeric(downsample_df_met$n_regions),as.numeric(downsample_df_met$heterogeneity), method = "spearman")
ggplot(downsample_df_met, aes(x=n_regions,y=heterogeneity)) + geom_boxplot() + theme_Publication() + ylab("Metabolic Heterogeneity") + 
  xlab("Number of Samples") + ggtitle("Spearman p-value = 0.23")
ggsave("Metabolite Heterogeneity N Samples Boxplot Downsample.pdf", height = 2, width = 2, units = "in")



#do for rna now
rna_distance_matrix <- intraheterogeneity_matrix(rna_mr_tumors_mc, 28)
colnames(rna_distance_matrix) <- rna_mapping$Sample.Name[match(colnames(rna_mr_tumors_mc),rna_mapping$Sample)]
patients_rna <- str_extract(colnames(rna_distance_matrix), "[^-]+")

downsample_rna <- function(patient,n_samples){
  columns <- sample(grep(patient,colnames(rna_distance_matrix)),n_samples)
  median(as.numeric(rna_distance_matrix[columns,columns]),na.rm=TRUE)
}

downsample_all_rna <- function(n_regions){
  patient <- names(which(table(patients_rna) >= n_regions))
  as.numeric(sapply(patient,function(i){replicate(10,downsample_rna(i,n_regions))}))
}

downsample_df_rna <- data.frame(n_regions = as.factor(c(rep(2,280),rep(3,260),rep(4,210),rep(5,150))),
                                heterogeneity = unlist(sapply(2:5,function(i){downsample_all_rna(i)})))
cor.test(as.numeric(downsample_df_rna$n_regions),as.numeric(downsample_df_rna$heterogeneity), method = "spearman")
ggplot(downsample_df_rna, aes(x=n_regions,y=heterogeneity)) + geom_boxplot() + theme_Publication() + ylab("RNA Heterogeneity") + 
  xlab("Number of Samples") + ggtitle("Spearman p-value = 0.15")
ggsave("RNA Heterogeneity N Samples Boxplot Downsample.pdf", height = 2, width = 2, units = "in")






