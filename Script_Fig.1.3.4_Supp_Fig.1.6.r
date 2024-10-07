################################
#
# September 2024 by MICHEL Elisa
#
# Do Figures in article Founder variants: Fig.3 ; Supp Fig.6 ; Fig.1 ; Fig.4 ; Supp Fig.1
################################
rm(list = ls());
### Packages
library("ggplot2")
library("ggpubr")
library("gridExtra")
library("ggpmisc")
library("stringr")
library("stringi")
library("tidyr")
library("grid")
library("openxlsx")
library("plotly")
library("htmlwidgets")
library("patchwork")

## Set the Directory
setwd("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/")

## File to open with all variants enriched (1304)
df_final <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt", header=T)
## Open file of variants known
variants_known <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/variants_reported_in_SLSJ.txt")

################################################################# Fig.3 ############################################################
## Get a copy of all the info
df_final_founder_RQ_SAG <- df_final
######### Select info about founder variants in SLSJ or UQc or QcP
df_final_founder_RQ_SAG$Founder_SAG <- ifelse(is.na(df_final$Ind_Imput_SAG), 
                                              df_final$founder_WGS_SAG,df_final$founder_Imput_SAG)
df_final_founder_RQ_SAG$Founder_RQ <- ifelse(is.na(df_final$Ind_Imput_RQ), 
                                             df_final$founder_WGS_RQ,df_final$founder_Imput_RQ)
df_final_founder_RQ_SAG$Founder_WQ <- ifelse(is.na(df_final$Ind_Imput_WQ), 
                                             df_final$founder_WGS_WQ,df_final$founder_Imput_WQ)
## Select founder variants SLSJ or UQc or QcP
df_final_founder_RQ_SAG <- df_final_founder_RQ_SAG[df_final_founder_RQ_SAG$Founder_SAG == "Founder" | df_final_founder_RQ_SAG$Founder_RQ == "Founder" | df_final_founder_RQ_SAG$Founder_WQ == "Founder",]
df_final_founder_RQ_SAG <- df_final_founder_RQ_SAG[complete.cases(df_final_founder_RQ_SAG$SNP),]
## Select columns to keep
df_final_founder_RQ_SAG_clean <- df_final_founder_RQ_SAG[,c(1:23)]
df_final_founder_RQ_SAG_clean <- df_final_founder_RQ_SAG_clean[,-c(4,5,7,8,12,13,15,16,21,22)]
## Select carrier rate for SLSJ and UQc (whether Imput or WGS if missing)
df_final_founder_RQ_SAG_clean$CR_SAG <- ifelse(is.na(df_final_founder_RQ_SAG_clean$CR_Imput_SAG), 
                                               df_final_founder_RQ_SAG_clean$CR_WGS_SAG,df_final_founder_RQ_SAG_clean$CR_Imput_SAG)

df_final_founder_RQ_SAG_clean$CR_RQ <- ifelse(is.na(df_final_founder_RQ_SAG_clean$CR_Imput_RQ), 
                                              df_final_founder_RQ_SAG_clean$CR_WGS_RQ,df_final_founder_RQ_SAG_clean$CR_Imput_RQ)
## Select Individuals that correspond in QcP
df_final_founder_RQ_SAG_clean$Ind_WQ <- ifelse(is.na(df_final_founder_RQ_SAG_clean$Ind_Imput_WQ), 
                                               df_final_founder_RQ_SAG_clean$Ind_WGS_WQ,df_final_founder_RQ_SAG_clean$Ind_Imput_WQ)
df_final_founder_RQ_SAG_clean$Source_ind_WQ <- ifelse(is.na(df_final_founder_RQ_SAG_clean$Ind_Imput_WQ), 
                                                     "WGS","Imputed data")
## Say if founder or not 
df_final_founder_RQ_SAG_clean$Founder_SAG <- df_final_founder_RQ_SAG$Founder_SAG
df_final_founder_RQ_SAG_clean$Founder_RQ <- df_final_founder_RQ_SAG$Founder_RQ
## Select rows for graphs
df_for_graphs_RQ_SAG <- df_final_founder_RQ_SAG_clean[,c(1,2,3,14,15,16,17,18,19)]
## Fill all NA values for peak
df_for_graphs_RQ_SAG[is.na(df_for_graphs_RQ_SAG$CR_SAG),4] <- 0
df_for_graphs_RQ_SAG[is.na(df_for_graphs_RQ_SAG$CR_RQ),5] <- 0
## Set missing as not founder
df_for_graphs_RQ_SAG$Founder_RQ = ifelse(is.na(df_for_graphs_RQ_SAG$Founder_RQ), "Not_founder",df_for_graphs_RQ_SAG$Founder_RQ)
df_for_graphs_RQ_SAG$Founder_SAG = ifelse(is.na(df_for_graphs_RQ_SAG$Founder_SAG), "Not_founder",df_for_graphs_RQ_SAG$Founder_SAG)
## Set founder both or none or only one
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Founder'] <- 'Both'
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Not_founder'] <- 'Founder_RQ'
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Not_founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Founder'] <- 'Founder_SAG'
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Not_founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Not_founder'] <- 'None'
## Select With variants known or not
df_for_graphs_RQ_SAG$Status <-ifelse(df_for_graphs_RQ_SAG$ID %in% variants_known[,1], "Known founder variant",
                                     "New founder variant")
## Do the graph for SLSJ
founder_WQ_known <- ggplot(data = df_for_graphs_RQ_SAG) + geom_point(aes(x = CR_SAG, y = CR_RQ,colour =Status,
                                      shape=Source_ind_WQ,text=SNP),position=position_jitter(h=0.00,w=0.15))+
  scale_shape_manual(values=c(19,0))+
  theme(axis.text.x = element_text(face="bold", size = 20),axis.title=element_text(size=22,face="bold"),
        axis.text.y = element_text(face="bold", size = 20),legend.text = element_text(size=14)) + 
  scale_color_manual(values = c("red","blue", "green"))+
  theme(legend.title = element_blank(),legend.position = "bottom")+ guides(color=guide_legend(override.aes = list(size=4)))+
  xlab("Carrier rate in SLSJ") + ylab("Carrier rate in UQc")+
  scale_x_continuous("Carrier rate in SLSJ",breaks=c(0,100,200,300),labels = c(0,"1/100","1/200","1/300")) +
  scale_y_continuous("Carrier rate in UQc",breaks=c(0,2000,4000,6000),
                     labels = c(0,"1/2,000","1/4,000","1/6,000"))

## Save in PDF format
outfile <- paste0("Figure_3_CR_SAG_RQ.pdf")
pdf(outfile,  width = 10, height = 10)
print(founder_WQ_known)
dev.off()
## Save in PNG format
ggsave("Figure_3_CR_SAG_RQ.png", width = 13, height = 13,dpi=300)

################################################################# Supp Fig.6 ############################################################
######### Graph Correlation Imput / WGS 
## Select data for Graph 
df_final_hist <- df_final[,c(1:23)]
## Select variants with WGS and Imput info
df_final_hist_CR_both_complete <- df_final_hist[complete.cases(df_final_hist$Ind_Imput_WQ),]
df_final_hist_CR_both_complete <- df_final_hist_CR_both_complete[complete.cases(df_final_hist_CR_both_complete$Ind_WGS_WQ),]
## Calculate the CR
df_final_hist_CR_both_complete$Prop_WGS_WQ <- df_final_hist_CR_both_complete$Ind_WGS_WQ/1852
df_final_hist_CR_both_complete$Prop_Imput_WQ <- df_final_hist_CR_both_complete$Ind_Imput_WQ/25061
## Do the equation
lm_eqn <- function(df){
  m <- lm(Prop_Imput_WQ~ Prop_WGS_WQ, df);
  eq <- substitute(~~italic(R)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
## Do the figure
figure_method_compare_WGS_Imput <- ggplot(data = df_final_hist_CR_both_complete,aes(x = Prop_WGS_WQ, y = Prop_Imput_WQ)) + geom_point()+
  scale_shape_manual(values=c(1,5))+
  theme(axis.text.x = element_text(face="bold", size = 15),axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(face="bold", size = 15),legend.text = element_text(size=16)) + 
  theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=4)))+
  xlab("Variant frequency in WGS") + ylab("Variant frequency in imputed data")+
  geom_text(x = 0.07, y = 0.26, label = lm_eqn(df_final_hist_CR_both_complete), parse = TRUE,size=8)

## Save in PDF format
outfile <- paste0("Figure_supp_6.pdf")
pdf(outfile,  width = 20, height = 10)
print(figure_method_compare_WGS_Imput)
dev.off()
## Save in PNG format
ggsave("Figure_supp_6.png", width = 8, height = 8,dpi=300)

################################################################# Fig.1 ############################################################
######### Graph Gnomad/(Imput if missing WGS) MAF on enriched Variants (1304) in QcP
## Get Variants with freq of Gnomad for Imput and WGS
NTNFE_WGS <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/enriched_variants/final_table_VTA_variant_enriched_in_WQ_compared_NTNFE_WGS_CaG_all_chr_tresh_0.1.txt", header =T)
NTNFE_WGS$SNP <- paste0("chr",NTNFE_WGS$CHROM,"_",NTNFE_WGS$POS,"_",NTNFE_WGS$REF,"_",NTNFE_WGS$ALT)
NTNFE_WGS <- NTNFE_WGS[,c(5,18,13,15)]
NTNFE_imput <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/enriched_variants/Imput/final_table_VTA_variant_enriched_in_WQ_compared_NTNFE_Imput_CaG_tresh_0.1.txt", header =T)
NTNFE_imput <- NTNFE_imput[,c(1,13,15)]
## Create table
mat = matrix(ncol = 0, nrow = 1304)
df_gnomad=data.frame(mat)
df_gnomad$SNP <- df_final$SNP
df_gnomad <- merge(df_gnomad,NTNFE_imput,by='SNP', all.x=T, all.y=F)
df_gnomad <- merge(df_gnomad,NTNFE_WGS,by='SNP', all.x=T, all.y=F)
df_gnomad <- df_gnomad[,-2]
names(df_gnomad) <- c("SNP","MAF_WQ_Imput","ID","MAF_GNOMAD","MAF_WQ_WGS")
## Get MAF missing
df_gnomad$MAF_WQ <- ifelse(is.na(df_gnomad$MAF_WQ_Imput),df_gnomad$MAF_WQ_WGS,df_gnomad$MAF_WQ_Imput)
## Get the Variants founder in QcP WGS and Imput
variants_WGS_WQ_founder <- read.table('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Graphs_Tables/WGS/Table_founder_variants_in_WQ.txt',header=T)
variants_Imput_WQ_founder <- read.table('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Graphs_Tables/Imput/Table_founder_variants_in_WQ.txt',header=T)
## Get the status of the variant
df_gnomad$founder_Imput_WQ <- df_final$founder_Imput_WQ
df_gnomad$founder_WGS_WQ <- df_final$founder_WGS_WQ
## Get the founder status of either Imput and if missing WGS
df_gnomad$Founder_status<-  ifelse(is.na(df_gnomad$founder_Imput_WQ),df_gnomad$founder_WGS_WQ,df_gnomad$founder_Imput_WQ )
## Replace not_founder by Not founder
df_gnomad$Founder_status<-  gsub("Not_founder", "Non founder with relative difference ≥ 0.1", df_gnomad$Founder_status)
df_gnomad$Founder_status<-  ifelse(is.na(df_gnomad$Founder_status),"Non founder with relative difference ≥ 0.1",df_gnomad$Founder_status)
## Get info if known or new
df_gnomad$Founder_status<-  ifelse((df_gnomad$ID %in% variants_Imput_WQ_founder$ID | (!(df_gnomad$ID %in%variants_Imput_WQ_founder) & df_gnomad$ID %in% variants_WGS_WQ_founder$ID)),
                                    "New founder variant",df_gnomad$Founder_status)
df_gnomad$Founder_status<-  ifelse((df_gnomad$ID %in% variants_known[,1] & df_gnomad$Founder_status == "New founder variant"),
                                   "Known founder variant",df_gnomad$Founder_status)

## Do the figure
gnomad_vs_founder_enriched_only <- ggplot(data = df_gnomad) + geom_point(aes(x = MAF_WQ, y = MAF_GNOMAD))+
  scale_shape_manual(values=c(1,5))+
  theme(axis.text.x = element_text(face="bold", size = 15),axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(face="bold", size = 15),legend.text = element_text(size=16)) + 
  theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=4)))+
  xlab("Frequency in the QcP") + ylab("Frequency in gnomAD nfe") +xlim(0,0.03)+ylim(0,0.03)

######### Graph Gnomad/(Imput if missing WGS) MAF on enriched Variants (1304) in SLSJ 
## Variants with freq of Gnomad
## Get MAF for SLSJ in Imput ad WGS
NTNFE_WGS_SAG <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Data/WGS/SLSJ/mymerged_all_chr_CaG_seq_VTA_SAG_enriched.frq", header =T)
NTNFE_WGS_SAG <- NTNFE_WGS_SAG[,c(2,5)]
NTNFE_imput_SAG <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Data/Imput/SLSJ/mymerged_all_chr_Imput_VTA_SAG_enriched.frq", header =T)
NTNFE_imput_SAG <- NTNFE_imput_SAG[,c(2,5)]
NTNFE_imput_SAG$SNP <- gsub(":", "_", NTNFE_imput_SAG$SNP)
## Get variants freq in Gnomad
Gnomad_freq <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/enriched_variants/final_table_VTA_variant_enriched_in_WQ_compared_NTNFE_WGS_CaG_all_chr_tresh_0.1.txt", header =T)
Gnomad_freq$SNP <- paste0("chr",Gnomad_freq$CHROM,"_",Gnomad_freq$POS,"_",Gnomad_freq$REF,"_",Gnomad_freq$ALT)
Gnomad_freq <- Gnomad_freq[,c(5,13,18)]
## Create table
df_gnomad_SAG=data.frame(mat)
df_gnomad_SAG$SNP <- df_final$SNP
df_gnomad_SAG <- merge(df_gnomad_SAG,NTNFE_imput_SAG,by='SNP', all.x=T, all.y=F)
df_gnomad_SAG <- merge(df_gnomad_SAG,NTNFE_WGS_SAG,by='SNP', all.x=T, all.y=F)
df_gnomad_SAG <- merge(df_gnomad_SAG,Gnomad_freq,by='SNP', all.x=T, all.y=F)
names(df_gnomad_SAG) <- c("SNP","MAF_SAG_Imput","MAF_SAG_WGS","ID","MAF_GNOMAD")
## Select founder variants in SLSJ (Imput and WGS if missing)
## Get Variants founder in SLSJ WGS and Imput
variants_WGS_SAG_founder <- read.table('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Graphs_Tables/WGS/Table_founder_variants_in_SLSJ.txt',header=T)
variants_Imput_SAG_founder <- read.table('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Graphs_Tables/Imput/Table_founder_variants_in_SLSJ.txt',header=T)
## Get the status of the variant
df_gnomad_SAG$founder_Imput_SAG <- df_final$founder_Imput_SAG
df_gnomad_SAG$founder_WGS_SAG <- df_final$founder_WGS_SAG
## Get founder status of either Imput and if missing WGS
df_gnomad_SAG$Founder_status<-  ifelse(is.na(df_gnomad_SAG$founder_Imput_SAG),df_gnomad_SAG$founder_WGS_SAG,df_gnomad_SAG$founder_Imput_SAG)
## Replace not_founder by Not founder
df_gnomad_SAG$Founder_status<-  gsub("Not_founder", "Non founder with relative difference ≥ 0.1", df_gnomad_SAG$Founder_status)
df_gnomad_SAG$Founder_status<-  ifelse(is.na(df_gnomad_SAG$Founder_status),"Non founder with relative difference ≥ 0.1",df_gnomad_SAG$Founder_status)
## Get info if known or new
df_gnomad_SAG$Founder_status<-  ifelse((df_gnomad_SAG$ID %in% variants_Imput_SAG_founder$ID | (!(df_gnomad_SAG$ID %in%variants_Imput_SAG_founder) & df_gnomad_SAG$ID %in% variants_WGS_SAG_founder$ID)),
                                   "New founder variant",df_gnomad_SAG$Founder_status)
df_gnomad_SAG$Founder_status<-  ifelse((df_gnomad_SAG$ID %in% variants_known[,1]& df_gnomad_SAG$Founder_status == "New founder variant"), 
                                       "Known founder variant",df_gnomad_SAG$Founder_status)
## Get MAF (Imput or WGS if missing)
df_gnomad_SAG$MAF_SAG <- ifelse(is.na(df_gnomad_SAG$MAF_SAG_Imput),df_gnomad_SAG$MAF_SAG_WGS,df_gnomad_SAG$MAF_SAG_Imput)

## Do the figure
gnomad_vs_founder_SAG_enriched_only <- ggplot(data = df_gnomad_SAG) + geom_point(aes(x = MAF_SAG, y = MAF_GNOMAD))+
  scale_shape_manual(values=c(1,5))+
  theme(axis.text.x = element_text(face="bold", size = 15),axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(face="bold", size = 15),legend.text = element_text(size=16)) + 
  theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=4)))+
  xlab("Frequency in SLSJ") + ylab("Frequency in gnomAD nfe")+xlim(0,0.03)+ylim(0,0.03)

##### Assemble the figures
figure_1_freq_Gnomad <- ggarrange(gnomad_vs_founder_enriched_only, gnomad_vs_founder_SAG_enriched_only+ rremove("ylab"), ncol=2, common.legend = TRUE, legend="bottom",labels = c("A", "B"))
## Save in PDF format
outfile <- paste0("Figure_1_gnomAD.pdf")
cairo_pdf(outfile,  width = 20, height = 10)
print(figure_1_freq_Gnomad)
dev.off()
## Save in PNG format
ggsave("Figure_1_gnomAD.png", width = 17, height = 8,bg="white",dpi=300)

################################################################# Fig.4 ############################################################
##### Figure with CR > 1/1000 for SLSJ
## Select the columns with info
CR_for_graph <- df_final[,c(1,2,6,9,14,17,22,25,18,19,10,11,26,27)]
## Remove variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
CR_for_graph <- CR_for_graph[!CR_for_graph$SNP %in% variants_to_freq,]
## Select founder in SLSJ and UQc and QcP
CR_for_graph$Founder_SAG <- ifelse(is.na(CR_for_graph$CR_Imput_SAG),
                                   CR_for_graph$founder_WGS_SAG,CR_for_graph$founder_Imput_SAG)
CR_for_graph$Founder_RQ <- ifelse(is.na(CR_for_graph$CR_Imput_RQ),
                                  CR_for_graph$founder_WGS_RQ,CR_for_graph$founder_Imput_RQ)
CR_for_graph$Founder_WQ <- ifelse(is.na(CR_for_graph$CR_Imput_WQ),
                                  CR_for_graph$founder_WGS_WQ,CR_for_graph$founder_Imput_WQ)
CR_for_graph <- CR_for_graph[CR_for_graph$Founder_SAG == "Founder" | CR_for_graph$Founder_RQ == "Founder"| CR_for_graph$Founder_WQ == "Founder",]
## Get CR (either Imput or if missing WGS)
CR_for_graph$CR_SAG <- ifelse(is.na(CR_for_graph$CR_Imput_SAG),CR_for_graph$CR_WGS_SAG,CR_for_graph$CR_Imput_SAG)
CR_for_graph$CR_RQ <- ifelse(is.na(CR_for_graph$CR_Imput_RQ),CR_for_graph$CR_WGS_RQ,CR_for_graph$CR_Imput_RQ)
## Filter CR > 1/1000 in SLSJ
CR_for_graph <- CR_for_graph[CR_for_graph$CR_SAG <= 1000 | is.na(CR_for_graph$CR_SAG),]
CR_for_graph <- CR_for_graph[CR_for_graph$CR_RQ <= 1000| is.na(CR_for_graph$CR_RQ),]
CR_for_graph <- CR_for_graph[complete.cases(CR_for_graph$SNP),]
## Create tables for the graph
CR_for_graph[is.na(CR_for_graph)] <- -1
# Create bins with intervals of 20 for UQc
bins_rq <- seq(0, max(CR_for_graph$CR_RQ) + 100, by = 100)
# Cut the UQc data into bins
CR_for_graph$carrier_rate_intervals <- cut(CR_for_graph$CR_RQ, breaks = bins_rq, right = FALSE)
# Count the frequency of each interval for UQc
frequency_rq <- table(CR_for_graph$carrier_rate_intervals)
# Convert frequency table to dataframe for UQc
frequency_df_rq <- as.data.frame(frequency_rq)
frequency_df_rq$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_rq)))
frequency_df_rq$dataset <- "RQ"
# Calculate the midpoints of the intervals for UQc
frequency_df_rq$midpoints <- (bins_rq[-1] + bins_rq[-length(bins_rq)]) / 2
# Create bins with intervals of 20 for SLSJ
bins_sag <- seq(0, max(CR_for_graph$CR_SAG) + 100, by = 100)
# Cut the SLSJ data into bins
CR_for_graph$carrier_rate_intervals <- cut(CR_for_graph$CR_SAG, breaks = bins_sag, right = FALSE)
# Count the frequency of each interval for SLSJ
frequency_sag <- table(CR_for_graph$carrier_rate_intervals)
# Convert frequency table to dataframe for SLSJ
frequency_df_sag <- as.data.frame(frequency_sag)
frequency_df_sag$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_sag)))
frequency_df_sag$dataset <- "SAG"
# Calculate the midpoints of the intervals for SLSJ
frequency_df_sag$midpoints <- (bins_sag[-1] + bins_sag[-length(bins_sag)]) / 2
# Combine the frequency dataframes for UQc and SLSJ
combined_frequency <- rbind(frequency_df_rq, frequency_df_sag)

########################## Enriched 
## Select the columns with info
CR_for_graph_enriched <- df_final[,c(1,2,6,9,14,17,22,25,18,19,10,11)]
## Remove variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
CR_for_graph_enriched <- CR_for_graph_enriched[!CR_for_graph_enriched$SNP %in% variants_to_freq,]
## Get CR (either Imput or if missing WGS)
CR_for_graph_enriched$CR_SAG <- ifelse(is.na(CR_for_graph_enriched$CR_Imput_SAG),CR_for_graph_enriched$CR_WGS_SAG,CR_for_graph_enriched$CR_Imput_SAG)
CR_for_graph_enriched$CR_RQ <- ifelse(is.na(CR_for_graph_enriched$CR_Imput_RQ),CR_for_graph_enriched$CR_WGS_RQ,CR_for_graph_enriched$CR_Imput_RQ)
## Filter CR > 1/1000 in SLSJ
CR_for_graph_enriched <- CR_for_graph_enriched[CR_for_graph_enriched$CR_SAG <= 1000 | is.na(CR_for_graph_enriched$CR_SAG),]
CR_for_graph_enriched <- CR_for_graph_enriched[CR_for_graph_enriched$CR_RQ <= 1000| is.na(CR_for_graph_enriched$CR_RQ),]
CR_for_graph_enriched <- CR_for_graph_enriched[complete.cases(CR_for_graph_enriched$SNP),]
CR_for_graph_enriched[is.na(CR_for_graph_enriched)] <- -1
# Create bins with intervals of 20 for UQc
bins_rq_enriched <- seq(0, max(CR_for_graph_enriched$CR_RQ) + 100, by = 100)
# Cut the UQc data into bins
CR_for_graph_enriched$carrier_rate_intervals <- cut(CR_for_graph_enriched$CR_RQ, breaks = bins_rq_enriched, right = FALSE)
# Count the frequency of each interval for UQc
frequency_rq_enriched <- table(CR_for_graph_enriched$carrier_rate_intervals)
# Convert frequency table to dataframe for UQc
frequency_df_rq_enriched <- as.data.frame(frequency_rq_enriched)
frequency_df_rq_enriched$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_rq_enriched)))
frequency_df_rq_enriched$dataset <- "RQ_enriched"
# Calculate the midpoints of the intervals for UQc
frequency_df_rq_enriched$midpoints <- (bins_rq_enriched[-1] + bins_rq_enriched[-length(bins_rq_enriched)]) / 2
# Create bins with intervals of 20 for SLSJ
bins_sag_enriched <- seq(0, max(CR_for_graph_enriched$CR_SAG) + 100, by = 100)
# Cut the SLSJ data into bins
CR_for_graph_enriched$carrier_rate_intervals <- cut(CR_for_graph_enriched$CR_SAG, breaks = bins_sag_enriched, right = FALSE)
# Count the frequency of each interval for SLSJ
frequency_sag_enriched <- table(CR_for_graph_enriched$carrier_rate_intervals)
# Convert frequency table to dataframe for SLSJ
frequency_df_sag_enriched <- as.data.frame(frequency_sag_enriched)
frequency_df_sag_enriched$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_sag_enriched)))
frequency_df_sag_enriched$dataset <- "SAG_enriched"
# Calculate the midpoints of the intervals for SLSJ
frequency_df_sag_enriched$midpoints <- (bins_sag_enriched[-1] + bins_sag_enriched[-length(bins_sag_enriched)]) / 2
# Combine the frequency dataframes for UQc and SLSJ
combined_frequency_enriched <- rbind(frequency_df_rq_enriched, frequency_df_sag_enriched)

## Combine the 2 plots
combined_freq_all <- rbind(combined_frequency,combined_frequency_enriched)
combined_freq_all$Dots <- ifelse(stri_detect_fixed(combined_freq_all$dataset,"enriched"),"variants with RFD≥10%","Founder variants")
combined_freq_all$Region <- ifelse(stri_detect_fixed(combined_freq_all$dataset,"RQ"),"UQc","SLSJ")
## Do the figure
plot_final <- ggplot(combined_freq_all, aes(x = midpoints, y = Freq, color = dataset,linetype = Dots)) +
  geom_line() + geom_point(aes(shape=Dots),size = 3) +
  scale_shape_manual(values=c(19,0))+
  scale_x_continuous(breaks = bins_rq_enriched[-length(bins_rq_enriched)],
                     labels = c(0,"1/100","1/200","1/300","1/400","1/500","1/600","1/700","1/800","1/900"))+
  theme(axis.text.x = element_text(face="bold", size = 15),axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(face="bold", size = 15),legend.text = element_text(size=16)) + 
  theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(shape = c(19,0),size=4,linetype = c(1, 3))),shape = "none",linetype ="none")+
  scale_color_manual(labels = c("UQc for founder variants","UQc for variants with RFD≥10%",
                                "SLSJ for founder variants","SLSJ for variants with RFD≥10%"),
                     values = c("green","green","blue","blue"))+
  xlab("Carrier rate") + ylab("Count")

## Save in PDF format
outfile <- paste0("Figure_4.pdf")
cairo_pdf(outfile,  width = 20, height = 10)
print(plot_final)
dev.off()

## Save in PNG format
ggsave("Figure_4.png", width = 17, height = 8,bg="white",dpi=300)

################################### Supp Fig.3 ##############################################
## Copy the table with all info
df_final_CR <- df_final
## Open files with info
References <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/CR_reported_paper.txt", header=T)
References <- References[-c(35,6),]
CR_cumulated <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/CR_cumul/CR_cumulated_all_gene.txt", header=T)
df_final_known_founder <-read.xlsx("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/Table_known_variants.xlsx")
df_final_known_founder$Position.GRCh38 <- gsub(":","_",df_final_known_founder$Position.GRCh38)
## Add info of references
## Change format References
References$Carrier_rate_reported_SLSJ <- as.numeric(str_split_fixed(References$Carrier_rate_reported_SLSJ,"",3)[,3])
df_final_CR_ref <- merge(df_final_CR,References,by="ID",all=F)
## Select known
df_final_CR_ref <- df_final_CR_ref[df_final_CR_ref$SNP %in% df_final_known_founder$Position.GRCh38, ]
## Add info cumulated CR
df_final_CR_ref <- merge(df_final_CR_ref,CR_cumulated,by="ID",all=T)
## Select CR data
df_final_CR_ref$CR_SLSL_data <- ifelse(is.na(df_final_CR_ref$Ind_Imput_SAG), 
                                  df_final_CR_ref$CR_WGS_SAG,df_final_CR_ref$CR_Imput_SAG)
## Select CR cumulated first, if missing CR data
df_final_CR_ref$CR_SLSL_final <- ifelse(is.na(df_final_CR_ref$CR_cumulated_SLSJ), 
                                  df_final_CR_ref$CR_SLSL_data,df_final_CR_ref$CR_cumulated_SLSJ)
df_final_CR_ref <- df_final_CR_ref[complete.cases(df_final_CR_ref$Carrier_rate_reported_SLSJ),]
## Create the equation
lm_eqn <- function(df){
  m <- lm(Carrier_rate_reported_SLSJ~ CR_SLSL_final, df);
  eq <- substitute(~~italic(R)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
## Remove the dot not relevant
df_final_CR_ref <- df_final_CR_ref[df_final_CR_ref$ID != "21439",]
## Do the graph
CR_graph <- ggplot(data = df_final_CR_ref,aes(x = CR_SLSL_final, y = Carrier_rate_reported_SLSJ)) + geom_point()+
  scale_shape_manual(values=c(1,5))+
  theme(axis.text.x = element_text(face="bold", size = 15,angle=45, hjust = 1),axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(face="bold", size = 15),legend.text = element_text(size=16)) + 
  theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=4)))+
  scale_x_continuous("CR found in our analysis in SLSJ",breaks=c(0,25,50,75,100,125,150,175,200),labels = c("0","1/25","1/50","1/75","1/100","1/125","1/150","1/175","1/200")) +
  scale_y_continuous("CR reported in the literature in SLSJ",breaks=c(0,25,50,75,100,125),
                     labels = c("0","1/25","1/50","1/75","1/100","1/125"))+
  geom_text(x = 75, y = 100, label = lm_eqn(df_final_CR_ref), parse = TRUE,size=8)

## Save in PDF format
outfile <- paste0("Figure_supp_1.pdf")
cairo_pdf(outfile,  width = 20, height = 10)
print(CR_graph)
dev.off()

## Save in PNG format
ggsave("Figure_supp_1.png", width = 8, height = 8,bg="white",dpi=300)
