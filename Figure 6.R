# Figure 6

### ggplot theme for publication
theme_Publication <- function(base_size=7, family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5, family = family),
            text = element_text(family = family),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(family = family),
            axis.title.y = element_text(angle=90,vjust =2, family = family),
            axis.title.x = element_text(vjust = -0.2, family = family),
            axis.text = element_text(family = family),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            # legend.margin = unit(t=0, unit="cm"),
            legend.title = element_text(face="italic", family = family),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#ffffff",fill="#ffffff"),
            strip.text = element_text()
    ))
  
}

# Panel A left
results_dir = "MIRTH/results_MET_RNA_imputation/tcga"
plot.dir= "MIRTH/results_MET_RNA_imputation/tcga/plots"

########## scatter plot of tumor vs normal differential analysis
results_path =paste0(results_dir,"/wilcox_tn_rc12_tcga_110met_mean_difference.csv")
wilcox_df <- read.csv(results_path)

wilcox_df <- dplyr::rename(wilcox_df,metabolite = X)
wilcox_df <- wilcox_df %>% 
  mutate(sig.status = case_when(
    p_adj>= 0.05 | p_adj_tcga >= 0.05 ~ "ns",
    log2fc > 0 & log2fc_tcga > 0 ~ "significant in both, consistent",
    log2fc < 0 & log2fc_tcga < 0 ~ "significant in both, consistent",
    T ~ "significant in both, inconsistent")) %>%
  mutate(metabolite = case_when(
    metabolite == "glutathione, reduced (GSH)" ~ "GSH",
    metabolite == "nicotinamide adenine dinucleotide (NAD+)" ~ "NAD+",
    metabolite == "glutathione, oxidized (GSSG)" ~ "GSSG",
    metabolite == "gamma-aminobutyrate (GABA)" ~ "GABA",
    metabolite == "fructose 1-phosphate" ~ "F1P",
    metabolite == "fructose-6-phosphate" ~ "F6P",
    metabolite == "glucose 1-phosphate" ~ "G1P",
    T ~ metabolite
  ))

wilcox_df %>% 
  ggplot(aes(x=log2fc,
             y=log2fc_tcga,
             color=sig.status))+
  geom_point()+
  scale_color_manual(values = c("ns" = "grey",
                                "significant in both, consistent" = "#332288",
                                "significant in both, inconsistent" = "red"))+
  scale_x_continuous(breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))+
  scale_y_continuous(breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))+
  coord_equal()+
  theme_Publication()+
  labs(x="RC12 T/N rank difference (measured)",
       y="TCGA T/N rank difference (imputed)",
       title = "TCGA vs RC12")
ggsave(file.path(plot.dir,"rc12_tcga_tn_110met_mean_difference_new.pdf"),width = 6,height = 6)

# Panel A right
####### scatter plot of High vs Low Stage differential analysis
results_path =paste0(results_dir,"/wilcox_stage_rc12_tcga_110met_mean_difference.csv")
wilcox_df <- read.csv(results_path)

wilcox_df <- dplyr::rename(wilcox_df,metabolite = X)
wilcox_df <- wilcox_df %>% 
  mutate(sig.status = case_when(
    p_adj>= 0.05 | p_adj_tcga >= 0.05 ~ "ns",
    log2fc > 0 & log2fc_tcga > 0 ~ "significant in both, consistent",
    log2fc < 0 & log2fc_tcga < 0 ~ "significant in both, consistent",
    T ~ "significant in both, inconsistent")) 
count <- wilcox_df %>% count(sig.status)

wilcox_df %>% 
  ggplot(aes(x=log2fc,
             y=log2fc_tcga,
             color=sig.status))+
  geom_point()+
  scale_color_manual(values = c("ns" = "grey",
                                "significant in both, consistent" = "#332288",
                                "significant in both, inconsistent" = "red"))+
  scale_x_continuous(breaks = c(-0.3,-0.15,0,0.15),limits = c(-0.35,0.17))+
  scale_y_continuous(breaks = c(-0.15,0,0.15),limits = c(-0.17,0.15))+
  theme_Publication()+
  labs(x="RC12 High stage/Low stage rank difference (measured)",
       y="TCGA High stage/Low stage rank difference (imputed)",
       title = "TCGA vs RC12 Stage")
ggsave(file.path(plot.dir,"rc12_tcga_stage_110met_mean_difference_new.pdf"),width = 6,height = 6)

# Panel D
# Meta analysis and Volcano plots Sunitinib Arms
library(metafor)
sub_dir = 'javelin_101'
results_dir = "MIRTH/results_MET_RNA_imputation"
results_path =paste0(results_dir,"/",sub_dir,"/agesex_padj_cox_df_PFS_C_262.csv")
met_cox_df_1 <- read.csv(results_path)

sub_dir = 'IMmotion151'
results_dir = "MIRTH/results_MET_RNA_imputation"
results_path =paste0(results_dir,"/",sub_dir,"/agesex_padj_cox_df_PFS_C_262.csv")
met_cox_df_2 <- read.csv(results_path)

sub_dir = 'Comparz'
results_dir = "MIRTH/results_MET_RNA_imputation"
results_path =paste0(results_dir,"/",sub_dir,"/agesex_padj_cox_df_PFS_C_262.csv")
met_cox_df_3 <- read.csv(results_path)

sub_dir = 'CheckMate214'
results_dir = "MIRTH/results_MET_RNA_imputation"
results_path =paste0(results_dir,"/",sub_dir,"/agesex_padj_cox_df_PFS_C_262.csv")
met_cox_df_4 <- read.csv(results_path)

len_met = 262

meta_df <- data.frame(matrix(ncol = 3, nrow = len_met))
colnames(meta_df) <- c('X','coef', 'p')

for (i in 1:len_met) {
  combined <- rbind(met_cox_df_1[i,],met_cox_df_2[i,],met_cox_df_3[i,],met_cox_df_4[i,])
  res <- rma(combined$coef, sei=combined$se.coef., data=combined)
  meta_df[i,]$X <- combined[1,]$X 
  meta_df[i,]$coef <- res$beta
  meta_df[i,]$p <- res$pval
}
meta_df$p_adj <- p.adjust(meta_df$p,method="BH")
write.csv(meta_df,"MIRTH/results_MET_RNA_imputation/metafor_multi_C_262_PFS_padj.csv", row.names = FALSE)

# Volcano Plot
plot <- EnhancedVolcano(meta_df,
                        lab = meta_df$X,
                        x = 'coef',
                        y = "p",
                        title = 'CPH regression',
                        pCutoff = 0.05/262,
                        FCcutoff = 0,
                        #selectLab = selectlab,
                        legendPosition = 'none',
                        pointSize = 3.0,
                        labSize = 4.0,
                        labCol = 'black',
                        labFace = 'bold',
                        boxedLabels = TRUE,
                        colAlpha = 1,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,
                        colConnectors = 'black',
                        max.overlaps = Inf
)  

plot <- plot + 
  ggplot2::ylab("-log10(p_unadj) ") +
  ggplot2::xlab("log(Hazard Ratio)") +
  ggplot2::scale_x_continuous(breaks = c(-1.5,0,1.5),limits = c(-1.5,1.5))
plot




