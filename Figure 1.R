library(stringr)

load("All Data.RData")
load("Outlier Data.RData")

master_mapping$ID_METABOLON_CLIENT[which(master_mapping$RCC_SUBTYPE=="clear cell")][164:392] #all the clear cell tumors in m1 m4 m5
MetaData_M4$SUBJECT_ID[match(master_mapping$ID_METABOLON_CLIENT[which(master_mapping$RCC_SUBTYPE=="clear cell")][164:273],MetaData_M4$CLIENT_IDENTIFIER)] #client id of m4 cc
m4_ids <- MetaData_M4$SUBJECT_ID[match(fig1_hm$ID_METABOLON_CLIENT[1:110],MetaData_M4$CLIENT_IDENTIFIER)]

m5_ids <- master_mapping$ID_METABOLON_CLIENT[which(master_mapping$RCC_SUBTYPE=="clear cell")][274:392]
m5_ids <- gsub("-",".",m5_ids)
m5_ids <- MetaData_M5$SUBJECT_ID[match(m5_ids,MetaData_M5$SAMPLE_NAME)] #client id of m5 clear cell 


m4_m5 <- which(master_mapping$METABOLON_COHORT %in% c("2018_M4","2019_M5"))
cc_tumors <- which(master_mapping$RCC_SUBTYPE=="clear cell")
fig1_hm <- master_mapping[intersect(m4_m5,cc_tumors),]
fig1_hm <- fig1_hm[,c(1,2,9,10,12)]
fig1_hm$Patient <- c(m4_ids,m5_ids)

mr_patients <- c(names(which(table(m4_ids) >1)),names(which(table(m5_ids) >1)))

ifelse(fig1_hm$METABOLON_COHORT == "2018_M4",2,3)
ifelse(fig1_hm$Patient %in% mr_patients, 4,5)
ifelse(!is.na(fig1_hm$ID_METABOLON_CLIENT),1,0)
ifelse(!is.na(fig1_hm$ID_DNASEQ),1,0)
ifelse(!is.na(fig1_hm$ID_RNASEQ),1,0)

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



io_drugs <- c("Ipilimumab/Nivolumab","Nivolumab")
tki_drugs <- c("Cabozantinib")
combo <- c("Ipilimumab/Nivolumab;Lenvatinib/Pembrolizumab","Sunitinib;Ipilimumab/Nivolumab","Sunitinib;Lenvatinib/Pembrolizumab;Nivolumab;Cabozantinib")

treatment_vector <- rep(NA,229)
treatment_vector[which(fig1_hm$THERAPY_EXPOSURE =="Not exposed")] <- 6
treatment_vector[which(fig1_hm$THERAPY_EXPOSURE %in% io_drugs)] <- 7
treatment_vector[which(fig1_hm$THERAPY_EXPOSURE %in% tki_drugs)] <- 8
treatment_vector[which(fig1_hm$THERAPY_EXPOSURE %in% combo)] <- 9



#dataframe to make heatmap
fig1_for_ggplot <- data.frame(Batch = ifelse(fig1_hm$METABOLON_COHORT == "2018_M4",2,3),
                              Sampling = ifelse(fig1_hm$Patient %in% mr_patients, 4,5),
                              Treatment = treatment_vector,
                              Metabolomics = ifelse(!is.na(fig1_hm$ID_METABOLON_CLIENT),1,0),
                              RNA = ifelse(!is.na(fig1_hm$ID_RNASEQ),1,0),
                              DNA = ifelse(!is.na(fig1_hm$ID_DNASEQ),1,0),
                              Patient = fig1_hm$ID_METABOLON_CLIENT)
fig1_for_ggplot <- fig1_for_ggplot[fig1_for_ggplot$Sampling==4,]
fig1_for_ggplot <- fig1_for_ggplot[,-2]
m5_outliers_fixed <- gsub("\\.","-",m5_outliers_id)
fig1_for_ggplot <- fig1_for_ggplot[-which(fig1_for_ggplot$Patient%in%c(m4_outliers,m5_outliers_fixed)),]
fig1_for_ggplot_melt <- melt(fig1_for_ggplot)
fig1_for_ggplot_melt$variable <- factor(fig1_for_ggplot_melt$variable, levels = c("DNA","RNA","Metabolomics","Treatment","Sampling","Batch"))
fig1_for_ggplot_melt$Patient <- factor(fig1_for_ggplot_melt$Patient, levels = fig1_hm$ID_METABOLON_CLIENT)

ggplot(fig1_for_ggplot_melt, aes(x=Patient,y=variable,fill=factor(value))) + geom_tile() +
  scale_fill_manual(values=color_vector) + theme_Publication() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())




m4_cc_client_id <- MetaData_M4$CLIENT_IDENTIFIER[match(colnames(m4_mr_tumors),MetaData_M4$SAMPLE_NAME)]
master_m4_cc <- master_mapping[master_mapping$ID_METABOLON_CLIENT %in% m4_cc_client_id,]
master_m4_cc$Patient <- str_extract(master_m4_cc$ID_METABOLON_CLIENT, "[^_]+")
master_m4_cc$Value <- c("All")
master_m4_cc$Value <- ifelse(is.na(master_m4_cc$ID_RNASEQ),"No RNA",master_m4_cc$Value)
master_mr_cc_melt_m4 <- melt(table(master_m4_cc$Patient,master_m4_cc$Value))


m5_cc_client_id <- gsub("\\.","-",colnames(m5_mr_tumors))
master_m5_cc <- master_mapping[master_mapping$ID_METABOLON_CLIENT %in% m5_cc_client_id,]
master_m5_cc$Patient <- str_extract(master_m5_cc$BIOBANK_MATCH, "[^-]+")
master_m5_cc$Value <- c("All")
m5_no_dna <- which(is.na(master_m5_cc$ID_DNASEQ))
m5_no_rna <- which(is.na(master_m5_cc$ID_RNASEQ))
master_m5_cc$Value[setdiff(m5_no_dna,m5_no_rna)] <- c("Met and RNA")
master_m5_cc$Value[intersect(m5_no_rna,m5_no_dna)] <- c("Met only")
master_mr_cc_melt_m5 <- melt(table(master_m5_cc$Patient,master_m5_cc$Value))




m4_cc_client_id_n <- MetaData_M4$CLIENT_IDENTIFIER[match(colnames(m4_mr_normals),MetaData_M4$SAMPLE_NAME)]
master_m4_cc_n <- master_mapping[master_mapping$ID_METABOLON_CLIENT %in% m4_cc_client_id_n,]
master_m4_cc_n$Patient <- str_extract(master_m4_cc_n$ID_METABOLON_CLIENT, "[^_]+")
master_m4_cc_n$Value <- c("All Normal")
master_m4_cc_n$Value <- ifelse(is.na(master_m4_cc_n$ID_DNASEQ),"No DNA Normal",master_m4_cc_n$Value)
master_mr_cc_melt_m4_n <- melt(table(master_m4_cc_n$Patient,master_m4_cc_n$Value))



m5_cc_client_id_n <- gsub("\\.","-",colnames(m5_mr_normals))
master_m5_cc_n <- master_mapping[master_mapping$ID_METABOLON_CLIENT %in% m5_cc_client_id_n,]
master_m5_cc_n$Patient <- str_extract(master_m5_cc_n$BIOBANK_MATCH, "[^-]+")
master_m5_cc_n$Value <- c("Met only Normal")
master_m5_cc_n$Value <- ifelse(!is.na(master_m5_cc_n$ID_RNASEQ),"No DNA Normal",master_m5_cc_n$Value)
master_mr_cc_melt_m5_n <- melt(table(master_m5_cc_n$Patient,master_m5_cc_n$Value))

master_for_fig1b <- rbind(master_mr_cc_melt_m4,
                          master_mr_cc_melt_m5,
                          master_mr_cc_melt_m4_n,
                          master_mr_cc_melt_m5_n)
master_for_fig1b$SampleType <- c(rep("Tumor",74),rep("Normal",44))
var2_levels <- rev(c("All","Met and RNA","No RNA","Met only",
                     "All Normal","No DNA Normal","Met only Normal"))
master_for_fig1b$Var2 <- factor(master_for_fig1b$Var2, levels = var2_levels)
master_for_fig1b$Var1 <- factor(master_for_fig1b$Var1, levels = rev(mr_patients_cc))


ggplot(master_for_fig1b, aes(x=Var1,y=value,fill=Var2)) + geom_bar(stat = "identity") + coord_flip() + theme_Publication() +
  xlab("") + ylab("Number of Samples") + scale_fill_brewer(palette = "Pastel1")



