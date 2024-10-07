################################
#
# September 2024 by MICHEL Elisa
#
# Do Figure in article Founder variants: Fig.5
################################
rm(list = ls());
### Packages
library("dplyr")
library("openxlsx")
library("tidyr")
library("stringr")
library("ggpubr")
library("reshape2")
library("ggplot2")

print("Graphs number of carriers (founder variants)")

## Open file with known variants
variants_known <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/variants_reported_in_SLSJ.txt")
## Get known variants positions
variants_enriched <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Remove variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
variants_enriched <- variants_enriched[!variants_enriched$SNP %in% variants_to_freq,]
## Get variants known SNP
variants_known_SNP <- variants_enriched[variants_enriched$ID %in% variants_known[,1],1]

#### Put first row in column names 
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

####################### Graph SLSJ - known variants only and all (known and new) side by side ############################
## Get individuals info in SLSJ
individual_file_SLSJ <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Data/Imput/SLSJ/CaG_Imput_VTA_SAG_enriched_in_WQ_chr_1.tfam'))
## Open file with all info
df_final <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Remove variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
df_final <- df_final[!df_final$SNP %in% variants_to_freq,]
## Keep variants founder in SLSJ, UQc and QcP
df_final$Founder_SAG <- ifelse(is.na(df_final$Ind_Imput_SAG), 
                               df_final$founder_WGS_SAG,df_final$founder_Imput_SAG)
df_final$Founder_RQ <- ifelse(is.na(df_final$Ind_Imput_RQ), 
                              df_final$founder_WGS_RQ,df_final$founder_Imput_RQ)
df_final$Founder_WQ <- ifelse(is.na(df_final$Ind_Imput_WQ), 
                              df_final$founder_WGS_WQ,df_final$founder_Imput_WQ)
## Select founder variants
Founder_SLSJ_RQ <- df_final[df_final$Founder_SAG == "Founder" | df_final$Founder_RQ == "Founder"| df_final$Founder_WQ == "Founder", ]
#### Create table with all variants enriched (compared with gnomAD) - Imput
mat = matrix(ncol = 0, nrow = 5000)
PLPCL_Imput=data.frame(mat)
for (i in 1:22){
  table_ind <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Imput/SLSJ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=F)
  table_ind <- as.data.frame(table_ind)
  if (nrow(table_ind)>5000){
    table_ind <- table_ind[c(-5001:-nrow(table_ind)),]
  } else if (nrow(table_ind)<5000){
    table_ind[nrow(table_ind)+1:(5000-nrow(table_ind)),] <- NA
  }
  PLPCL_Imput <- cbind(PLPCL_Imput,table_ind)
}
## Change : in _ in variant name
PLPCL_Imput[1,] <- gsub(":", "_", PLPCL_Imput[1,])
## Remove NR from imput data
NR_var <- read.table("/lustre03/project/6033529/saguenay_disease/results/correspondance_WGS_imput/snps_to_remove_seuil1_allcriteria_homswitch.txt")
PLPCL_Imput <- PLPCL_Imput[,!PLPCL_Imput[1,] %in% NR_var[,1]]
## Select only variants founders in SLSJ
PLPCL_founder_Imput <- PLPCL_Imput[,PLPCL_Imput[1,] %in% Founder_SLSJ_RQ$SNP]
PLPCL_founder_Imput <- header.true(PLPCL_founder_Imput)
## Found missing variants founder to add with WGS
variants_missing_SAG <- Founder_SLSJ_RQ[!Founder_SLSJ_RQ$SNP %in% colnames(PLPCL_founder_Imput),1]
#### Add info with WGS if imput missing (compared with gnomAD) - WGS for missing variants
PLPCL_WGS=data.frame(mat)
for (i in 1:22){
  table_ind <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/WGS/SLSJ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=F)
  table_ind <- as.data.frame(table_ind)
  if (nrow(table_ind)>5000){
    table_ind <- table_ind[c(-5001:-nrow(table_ind)),]
  } else if (nrow(table_ind)<5000){
    table_ind[nrow(table_ind)+1:(5000-nrow(table_ind)),] <- NA
  }
  PLPCL_WGS <- cbind(PLPCL_WGS,table_ind)
}
## Select only variants founders in SLSJ
PLPCL_founder_WGS <- PLPCL_WGS[,PLPCL_WGS[1,] %in% variants_missing_SAG]
PLPCL_founder_WGS <- header.true(PLPCL_founder_WGS)
## Merge all info between Imput and WGS
PLPCL_founder_all_SAG <- cbind(PLPCL_founder_WGS,PLPCL_founder_Imput)
print(paste0("variants_for_SLSJ : ",ncol(PLPCL_founder_all_SAG)))
## Select info for known variants
PLPCL_founder_all_SAG_known <- PLPCL_founder_all_SAG[,colnames(PLPCL_founder_all_SAG) %in% variants_known_SNP]
######### Do the graph with inds with x variants founder in SLSJ ALL (known and unknown)
#### Search for ind with more than one variant
ind_with_variants=data.frame(mat)
for (i in 1:nrow(individual_file_SLSJ)){
  index=1
  ind <- individual_file_SLSJ[i,1]
  ind_with_variants[index,i] <- ind
  for (l in 1:ncol(PLPCL_founder_all_SAG)){
    if (ind %in% PLPCL_founder_all_SAG[,l]){
      ind_with_variants[index+1,i] <- colnames(PLPCL_founder_all_SAG[l])
      #print(colnames(PLPCL_founder_all_SAG[l]))
      index=index+1
    }
  }
}
#### Select ind with min 2 variants
ind_with_multiple_variants=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants)){
  if (!is.na(ind_with_variants[3,h])){
    index_bis = index_bis +1
    ind_with_multiple_variants[,index_bis] <- ind_with_variants[,h]
  }
}
ind_with_multiple_variants <- header.true(ind_with_multiple_variants)
#### Select ind with min 1 variant
ind_with_one_or_more_variants=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants)){
  if (!is.na(ind_with_variants[2,h])){
    index_bis = index_bis +1
    ind_with_one_or_more_variants[,index_bis] <- ind_with_variants[,h]
  }
}
ind_with_one_or_more_variants <- header.true(ind_with_one_or_more_variants)
#### Select ind with zero variant
ind_with_zero_variants=data.frame(mat)
index_tri=0
for (h in 1:ncol(ind_with_variants)){
  if (is.na(ind_with_variants[2,h])){
    index_tri = index_tri +1
    ind_with_zero_variants[,index_tri] <- ind_with_variants[,h]
  }
}
ind_with_zero_variants <- header.true(ind_with_zero_variants)
#### Do histograms : number of variants by individual
## Calculate how many inds with one variant
ind_with_one_variant <- ind_with_one_or_more_variants[,is.na(ind_with_one_or_more_variants[2,])]
## Calculate how many inds with two variants
ind_with_two_variant <- ind_with_multiple_variants[,is.na(ind_with_multiple_variants[3,])]
## Calculate how many inds with three variants
ind_with_three_or_more_variant <- ind_with_multiple_variants[,!is.na(ind_with_multiple_variants[3,])]
ind_with_three_variant <- ind_with_three_or_more_variant[,is.na(ind_with_three_or_more_variant[4,])]
## Calculate how many inds with four variants
ind_with_four_or_more_variant <- ind_with_three_or_more_variant[,!is.na(ind_with_three_or_more_variant[4,])]
ind_with_four_variant <- ind_with_four_or_more_variant[,is.na(ind_with_four_or_more_variant[5,])]
# Calculate how many inds with five variants
ind_with_five_variant_or_more <- as.data.frame(ind_with_four_or_more_variant[,!is.na(ind_with_four_or_more_variant[5,])])
ind_with_five_variant <- as.data.frame(ind_with_five_variant_or_more[,is.na(ind_with_five_variant_or_more[6,])])
# Calculate how many inds with six variants
ind_with_six_variant_or_more <- as.data.frame(ind_with_five_variant_or_more[,!is.na(ind_with_five_variant_or_more[6,])])
ind_with_six_variant <- as.data.frame(ind_with_six_variant_or_more[,is.na(ind_with_six_variant_or_more[7,])])
# Calculate how many inds with seven variants
ind_with_seven_variant_or_more <- as.data.frame(ind_with_six_variant_or_more[,!is.na(ind_with_six_variant_or_more[7,])])
ind_with_seven_variant <- as.data.frame(ind_with_seven_variant_or_more[,is.na(ind_with_seven_variant_or_more[8,])])

print(paste0("nb var with seven or more SAG :",ncol(ind_with_seven_variant_or_more)))
print(paste0("nb var with seven SAG :",ncol(ind_with_seven_variant)))

######### Do the graph with inds with x variants founder in SLSJ known only
#### Search for ind with more than one variant
ind_with_variants_known=data.frame(mat)
for (i in 1:nrow(individual_file_SLSJ)){
  index=1
  ind <- individual_file_SLSJ[i,1]
  ind_with_variants_known[index,i] <- ind
  for (l in 1:ncol(PLPCL_founder_all_SAG_known)){
    if (ind %in% PLPCL_founder_all_SAG_known[,l]){
      ind_with_variants_known[index+1,i] <- colnames(PLPCL_founder_all_SAG_known[l])
      #print(colnames(PLPCL_founder_all_SAG_known[l]))
      index=index+1
    }
  }
}
#### Select ind with min 2 variants
ind_with_multiple_variants_known=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants_known)){
  if (!is.na(ind_with_variants_known[3,h])){
    index_bis = index_bis +1
    ind_with_multiple_variants_known[,index_bis] <- ind_with_variants_known[,h]
  }
}
ind_with_multiple_variants_known <- header.true(ind_with_multiple_variants_known)
#### Select ind with min 1 variants
ind_with_one_or_more_variants_known=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants_known)){
  if (!is.na(ind_with_variants_known[2,h])){
    index_bis = index_bis +1
    ind_with_one_or_more_variants_known[,index_bis] <- ind_with_variants_known[,h]
  }
}
ind_with_one_or_more_variants_known <- header.true(ind_with_one_or_more_variants_known)
#### Select ind with zero variant
ind_with_zero_variants_known=data.frame(mat)
index_tri=0
for (h in 1:ncol(ind_with_variants_known)){
  if (is.na(ind_with_variants_known[2,h])){
    index_tri = index_tri +1
    ind_with_zero_variants_known[,index_tri] <- ind_with_variants_known[,h]
  }
}
ind_with_zero_variants_known <- header.true(ind_with_zero_variants_known)
#### Do histograms : number of variants by individual
## Calculate how many inds with one variant
ind_with_one_variant_known <- ind_with_one_or_more_variants_known[,is.na(ind_with_one_or_more_variants_known[2,])]
## Calculate how many inds with two variants
ind_with_two_variant_known <- ind_with_multiple_variants_known[,is.na(ind_with_multiple_variants_known[3,])]
## Calculate how many inds with three variants
ind_with_three_or_more_variant_known <- ind_with_multiple_variants_known[,!is.na(ind_with_multiple_variants_known[3,])]
ind_with_three_variant_known <- ind_with_three_or_more_variant_known[,is.na(ind_with_three_or_more_variant_known[4,])]
## Calculate how many inds with four variants
ind_with_four_or_more_variant_known <- as.data.frame(ind_with_three_or_more_variant_known[,!is.na(ind_with_three_or_more_variant_known[4,])])
ind_with_four_variant_known <- as.data.frame(ind_with_four_or_more_variant_known[,is.na(ind_with_four_or_more_variant_known[5,])])
# Calculate how many inds with five variants
ind_with_five_variant_or_more_known <- as.data.frame(ind_with_four_or_more_variant_known[,!is.na(ind_with_four_or_more_variant_known[5,])])
ind_with_five_variant_known <- as.data.frame(ind_with_five_variant_or_more_known[,is.na(ind_with_five_variant_or_more_known[6,])])
# Calculate how many inds with six variants
ind_with_six_variant_or_more_known <- as.data.frame(ind_with_five_variant_or_more_known[,!is.na(ind_with_five_variant_or_more_known[6,])])
ind_with_six_variant_known <- as.data.frame(ind_with_six_variant_or_more_known[,is.na(ind_with_six_variant_or_more_known[7,])])
# Calculate how many inds with seven variants
ind_with_seven_variant_or_more_known <- as.data.frame(ind_with_six_variant_or_more_known[,!is.na(ind_with_six_variant_or_more_known[7,])])
ind_with_seven_variant_known <- as.data.frame(ind_with_seven_variant_or_more_known[,is.na(ind_with_seven_variant_or_more_known[8,])])

print(paste0("nb var with seven or more SAG known :",ncol(ind_with_seven_variant_or_more_known)))
print(paste0("nb var with seven SAG known :",ncol(ind_with_seven_variant_known)))

#### Collect info for Known and all (known and unknown)
mat_2 = matrix(ncol = 3, nrow = 21) 
df_hist_SAG=data.frame(mat_2)
names(df_hist_SAG) <- c("Number_of_variants","Founder variants","Known founder variants")
npop_SAG <- 3589
df_hist_SAG[1,1] <- 0
df_hist_SAG[2,1] <- 1
df_hist_SAG[3,1] <- 2
df_hist_SAG[4,1] <- 3
df_hist_SAG[5,1] <- 4
df_hist_SAG[6,1] <- 5
df_hist_SAG[7,1] <- 6
df_hist_SAG[1,2] <- ncol(ind_with_zero_variants)
df_hist_SAG[2,2] <- ncol(ind_with_one_variant)
df_hist_SAG[3,2] <- ncol(ind_with_two_variant)
df_hist_SAG[4,2] <- ncol(ind_with_three_variant)
df_hist_SAG[5,2] <- ncol(ind_with_four_variant)
df_hist_SAG[6,2] <- ncol(ind_with_five_variant)
df_hist_SAG[7,2] <- ncol(ind_with_six_variant)
df_hist_SAG[8,2] <- ncol(ind_with_seven_variant)
df_hist_SAG[1,3] <- ncol(ind_with_zero_variants_known)
df_hist_SAG[2,3] <- ncol(ind_with_one_variant_known)
df_hist_SAG[3,3] <- ncol(ind_with_two_variant_known)
df_hist_SAG[4,3] <- ncol(ind_with_three_variant_known)
df_hist_SAG[5,3] <- ncol(ind_with_four_variant_known)
df_hist_SAG[6,3] <- ncol(ind_with_five_variant_known)
df_hist_SAG[7,3] <- ncol(ind_with_six_variant_known)
df_hist_SAG <- df_hist_SAG[complete.cases(df_hist_SAG),]
## Transform data for the histograms
flip_results <- melt(df_hist_SAG, id.vars=c('Number_of_variants'))
flip_results$prop <- round(flip_results$value/npop_SAG,digits=3)

#### Do the graph
ID <- 0:6
dodger = position_dodge(width = 0.9)
graph_SAG <- ggplot(flip_results, aes(x=Number_of_variants, y=prop,fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Number of variants in SLSJ", y = "Proportion of individuals") + 
  scale_x_continuous("Number of variants in SLSJ", labels = as.character(ID), breaks = ID)+
  geom_text(aes(label=value,),
            fontface="bold", vjust = -0.5, color="black", size=3,position = dodger) +
  theme(plot.title = element_text(size=18, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=18, face="bold", colour = "black"),    
    axis.title.y = element_text(size=18, face="bold", colour = "black"),    
    axis.text.x = element_text(size=18, face="bold", colour = "black"), 
    axis.text.y = element_text(size=18, face="bold", colour = "black"),
    strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 12, face="bold", colour = "black"),
    legend.title = element_blank(),legend.text = element_text(size=16))+
  scale_fill_manual(values=c("Founder variants"="pink", "Known founder variants"="red")) +
  scale_colour_manual(values=c("Founder variants"="pink", "Known founder variants"="red"), 
                      labels=c("Founder variants", "Known founder variants")) +
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1))
  
print("SAG")
print(flip_results)

####################### Graph UQc - known variants only and all (known and new) side by side ############################
## Get individuals info in UQc
individual_file_RQ <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Data/Imput/RQ/CaG_Imput_VTA_RQ_enriched_in_WQ_chr_1.tfam'))
## Keep founder (SLSJ, UQc and QcP) with specific CR
df_final_CR <- df_final[,c(1,14,17)]
df_final_CR <- df_final_CR[df_final_CR$CR_Imput_RQ <= 5000 | is.na(df_final_CR$CR_Imput_RQ),]
df_final_CR <- df_final_CR[complete.cases(df_final_CR$SNP),]
Founder_SLSJ_RQ <- Founder_SLSJ_RQ[Founder_SLSJ_RQ$SNP %in% df_final_CR$SNP,]
#### Create table with all variants enriched (compared with gnomAD) - Imput
mat = matrix(ncol = 0, nrow = 25000)
PLPCL_Imput_RQ=data.frame(mat)
for (i in 1:22){
  table_ind <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Imput/RQ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=F)
  table_ind <- as.data.frame(table_ind)
  if (nrow(table_ind)>25000){
    table_ind <- table_ind[c(-25001:-nrow(table_ind)),]
  } else if (nrow(table_ind)<25000){
    table_ind[nrow(table_ind)+1:(25000-nrow(table_ind)),] <- NA
  }
  PLPCL_Imput_RQ <- cbind(PLPCL_Imput_RQ,table_ind)
}
## Change : in _ in variant name
PLPCL_Imput_RQ[1,] <- gsub(":", "_", PLPCL_Imput_RQ[1,])
## Remove NR imput
PLPCL_Imput_RQ <- PLPCL_Imput_RQ[,!PLPCL_Imput_RQ[1,] %in% NR_var[,1]]
## Select only variants founders in RQ
PLPCL_founder_Imput_RQ <- PLPCL_Imput_RQ[,PLPCL_Imput_RQ[1,] %in% Founder_SLSJ_RQ$SNP]
PLPCL_founder_Imput_RQ <- header.true(PLPCL_founder_Imput_RQ)
## Found missing variants founder to add with WGS
variants_missing_RQ <- Founder_SLSJ_RQ[!Founder_SLSJ_RQ$SNP %in% colnames(PLPCL_founder_Imput_RQ),1]
#### Add info with WGS if imput missing (compared with gnomAD) - WGS for missing variants
PLPCL_WGS_RQ=data.frame(mat)
for (i in 1:22){
  table_ind <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/WGS/RQ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=F)
  table_ind <- as.data.frame(table_ind)
  if (nrow(table_ind)>25000){
    table_ind <- table_ind[c(-25001:-nrow(table_ind)),]
  } else if (nrow(table_ind)<25000){
    table_ind[nrow(table_ind)+1:(25000-nrow(table_ind)),] <- NA
  }
  PLPCL_WGS_RQ <- cbind(PLPCL_WGS_RQ,table_ind)
}
## Select only variants founders in RQ
PLPCL_founder_WGS_RQ <- PLPCL_WGS_RQ[,PLPCL_WGS_RQ[1,] %in% variants_missing_RQ]
PLPCL_founder_WGS_RQ <- header.true(PLPCL_founder_WGS_RQ)
## Merge all info between Imput and WGS
PLPCL_founder_all_RQ <- cbind(PLPCL_founder_WGS_RQ,PLPCL_founder_Imput_RQ)
print(paste0("variants_for_RQ : ",ncol(PLPCL_founder_all_RQ)))
## Select info for known variants 
PLPCL_founder_all_RQ_known <- PLPCL_founder_all_RQ[,colnames(PLPCL_founder_all_RQ) %in% variants_known_SNP]
######### Do the graph with inds with x variants founder in UQc all (known and unknown)
#### Search for ind with more than one variant
ind_with_variants_RQ=data.frame(mat)
for (i in 1:nrow(individual_file_RQ)){
  index=1
  ind <- individual_file_RQ[i,1]
  ind_with_variants_RQ[index,i] <- ind
  for (l in 1:ncol(PLPCL_founder_all_RQ)){
    if (ind %in% PLPCL_founder_all_RQ[,l]){
      ind_with_variants_RQ[index+1,i] <- colnames(PLPCL_founder_all_RQ[l])
      #print(colnames(PLPCL_founder_all_RQ[l]))
      index=index+1
    }
  }
}
#### Select ind with min 2 variants
ind_with_multiple_variants_RQ=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants_RQ)){
  if (!is.na(ind_with_variants_RQ[3,h])){
    index_bis = index_bis +1
    ind_with_multiple_variants_RQ[,index_bis] <- ind_with_variants_RQ[,h]
  }
}
ind_with_multiple_variants_RQ <- header.true(ind_with_multiple_variants_RQ)
#### Select ind with min 1 variant
ind_with_one_or_more_variants_RQ=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants_RQ)){
  if (!is.na(ind_with_variants_RQ[2,h])){
    index_bis = index_bis +1
    ind_with_one_or_more_variants_RQ[,index_bis] <- ind_with_variants_RQ[,h]
  }
}
ind_with_one_or_more_variants_RQ <- header.true(ind_with_one_or_more_variants_RQ)
#### Select ind with zero variant
ind_with_zero_variants_RQ=data.frame(mat)
index_tri=0
for (h in 1:ncol(ind_with_variants_RQ)){
  if (is.na(ind_with_variants_RQ[2,h])){
    index_tri = index_tri +1
    ind_with_zero_variants_RQ[,index_tri] <- ind_with_variants_RQ[,h]
  }
}
ind_with_zero_variants_RQ <- header.true(ind_with_zero_variants_RQ)
#### Do histograms : number of variants by individual
## Calculate how many inds with one variant
ind_with_one_variant_RQ <- ind_with_one_or_more_variants_RQ[,is.na(ind_with_one_or_more_variants_RQ[2,])]
## Calculate how many inds with two variants
ind_with_two_variant_RQ <- ind_with_multiple_variants_RQ[,is.na(ind_with_multiple_variants_RQ[3,])]
## Calculate how many inds with three variants
ind_with_three_or_more_variant_RQ <- ind_with_multiple_variants_RQ[,!is.na(ind_with_multiple_variants_RQ[3,])]
ind_with_three_variant_RQ <- ind_with_three_or_more_variant_RQ[,is.na(ind_with_three_or_more_variant_RQ[4,])]
## Calculate how many inds with four variants
ind_with_four_or_more_variant_RQ <- ind_with_three_or_more_variant_RQ[,!is.na(ind_with_three_or_more_variant_RQ[4,])]
ind_with_four_variant_RQ <- ind_with_four_or_more_variant_RQ[,is.na(ind_with_four_or_more_variant_RQ[5,])]
# Calculate how many inds with five variants
ind_with_five_variant_or_more_RQ <- as.data.frame(ind_with_four_or_more_variant_RQ[,!is.na(ind_with_four_or_more_variant_RQ[5,])])
ind_with_five_variant_RQ <- as.data.frame(ind_with_five_variant_or_more_RQ[,is.na(ind_with_five_variant_or_more_RQ[6,])])
# Calculate how many inds with six variants
ind_with_six_variant_or_more_RQ <- as.data.frame(ind_with_five_variant_or_more_RQ[,!is.na(ind_with_five_variant_or_more_RQ[6,])])
ind_with_six_variant_RQ <- as.data.frame(ind_with_six_variant_or_more_RQ[,is.na(ind_with_six_variant_or_more_RQ[7,])])

print(paste0("nb var with six or more RQ :",ncol(ind_with_six_variant_or_more_RQ)))
print(paste0("nb var with six RQ :",ncol(ind_with_six_variant_RQ)))

######### Do the graph with inds with x variants founder in UQc KNOWN ONLY
#### Search for ind with more than one variant
ind_with_variants_known_RQ=data.frame(mat)
for (i in 1:nrow(individual_file_RQ)){
  index=1
  ind <- individual_file_RQ[i,1]
  ind_with_variants_known_RQ[index,i] <- ind
  for (l in 1:ncol(PLPCL_founder_all_RQ_known)){
    if (ind %in% PLPCL_founder_all_RQ_known[,l]){
      ind_with_variants_known_RQ[index+1,i] <- colnames(PLPCL_founder_all_RQ_known[l])
      #print(colnames(PLPCL_founder_all_RQ_known[l]))
      index=index+1
    }
  }
}
#### Select ind with min 2 variants
ind_with_multiple_variants_known_RQ=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants_known_RQ)){
  if (!is.na(ind_with_variants_known_RQ[3,h])){
    index_bis = index_bis +1
    ind_with_multiple_variants_known_RQ[,index_bis] <- ind_with_variants_known_RQ[,h]
  }
}
ind_with_multiple_variants_known_RQ <- header.true(ind_with_multiple_variants_known_RQ)
#### Select ind with min 1 variant
ind_with_one_or_more_variants_known_RQ=data.frame(mat)
index_bis=0
for (h in 1:ncol(ind_with_variants_known_RQ)){
  if (!is.na(ind_with_variants_known_RQ[2,h])){
    index_bis = index_bis +1
    ind_with_one_or_more_variants_known_RQ[,index_bis] <- ind_with_variants_known_RQ[,h]
  }
}
ind_with_one_or_more_variants_known_RQ <- header.true(ind_with_one_or_more_variants_known_RQ)
#### Select ind with zero variant
ind_with_zero_variants_known_RQ=data.frame(mat)
index_tri=0
for (h in 1:ncol(ind_with_variants_known_RQ)){
  if (is.na(ind_with_variants_known_RQ[2,h])){
    index_tri = index_tri +1
    ind_with_zero_variants_known_RQ[,index_tri] <- ind_with_variants_known_RQ[,h]
  }
}
ind_with_zero_variants_known_RQ <- header.true(ind_with_zero_variants_known_RQ)
#### Do histograms : number of variants by individual
## Calculate how many inds with one variant
ind_with_one_variant_known_RQ <- ind_with_one_or_more_variants_known_RQ[,is.na(ind_with_one_or_more_variants_known_RQ[2,])]
## Calculate how many inds with two variants
ind_with_two_variant_known_RQ <- ind_with_multiple_variants_known_RQ[,is.na(ind_with_multiple_variants_known_RQ[3,])]
## Calculate how many inds with three variants
ind_with_three_or_more_variant_known_RQ <- ind_with_multiple_variants_known_RQ[,!is.na(ind_with_multiple_variants_known_RQ[3,])]
ind_with_three_variant_known_RQ <- ind_with_three_or_more_variant_known_RQ[,is.na(ind_with_three_or_more_variant_known_RQ[4,])]
## Calculate how many inds with four variants
ind_with_four_or_more_variant_known_RQ <- as.data.frame(ind_with_three_or_more_variant_known_RQ[,!is.na(ind_with_three_or_more_variant_known_RQ[4,])])
ind_with_four_variant_known_RQ <- as.data.frame(ind_with_four_or_more_variant_known_RQ[,is.na(ind_with_four_or_more_variant_known_RQ[5,])])
# Calculate how many inds with five variants
ind_with_five_variant_or_more_known_RQ <- as.data.frame(ind_with_four_or_more_variant_known_RQ[,!is.na(ind_with_four_or_more_variant_known_RQ[5,])])
ind_with_five_variant_known_RQ <- as.data.frame(ind_with_five_variant_or_more_known_RQ[,is.na(ind_with_five_variant_or_more_known_RQ[6,])])
# Calculate how many inds with six variants
ind_with_six_variant_or_more_known_RQ <- as.data.frame(ind_with_five_variant_or_more_known_RQ[,!is.na(ind_with_five_variant_or_more_known_RQ[6,])])
ind_with_six_variant_known_RQ <- as.data.frame(ind_with_six_variant_or_more_known_RQ[,is.na(ind_with_six_variant_or_more_known_RQ[7,])])

print(paste0("nb var with six or more RQ known :",ncol(ind_with_six_variant_or_more_known_RQ)))
print(paste0("nb var with six RQ known :",ncol(ind_with_six_variant_known_RQ)))

#### Collect info for Known and all (known and unknown)
mat_2 = matrix(ncol = 3, nrow = 21) 
df_hist_RQ=data.frame(mat_2)
names(df_hist_RQ) <- c("Number_of_variants","Founder variants","Known founder variants")
npop_RQ <- 21472
df_hist_RQ[1,1] <- 0
df_hist_RQ[2,1] <- 1
df_hist_RQ[3,1] <- 2
df_hist_RQ[4,1] <- 3
df_hist_RQ[5,1] <- 4
df_hist_RQ[1,2] <- ncol(ind_with_zero_variants_RQ)
df_hist_RQ[2,2] <- ncol(ind_with_one_variant_RQ)
df_hist_RQ[3,2] <- ncol(ind_with_two_variant_RQ)
df_hist_RQ[4,2] <- ncol(ind_with_three_variant_RQ)
df_hist_RQ[5,2] <- ncol(ind_with_four_variant_RQ)
df_hist_RQ[6,2] <- ncol(ind_with_five_variant_RQ)
df_hist_RQ[1,3] <- ncol(ind_with_zero_variants_known_RQ)
df_hist_RQ[2,3] <- ncol(ind_with_one_variant_known_RQ)
df_hist_RQ[3,3] <- ncol(ind_with_two_variant_known_RQ)
df_hist_RQ[4,3] <- ncol(ind_with_three_variant_known_RQ)
df_hist_RQ[5,3] <- ncol(ind_with_four_variant_known_RQ)
df_hist_RQ <- df_hist_RQ[complete.cases(df_hist_RQ),]
## Transform data for the histogram
flip_results_RQ <- melt(df_hist_RQ, id.vars=c('Number_of_variants'))
flip_results_RQ$prop <- round(flip_results_RQ$value/npop_RQ,digits=3)
print("RQ")
print(flip_results_RQ)
#### Do the graph
ID <- 0:4
dodger = position_dodge(width = 0.9)
graphs_RQ <- ggplot(flip_results_RQ, aes(x=Number_of_variants, y=prop,fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Number of variants in UQc", y = "Proportion of individuals") + 
  scale_x_continuous("Number of variants in UQc", labels = as.character(ID), breaks = ID)+
  geom_text(aes(label=value,),
            fontface="bold", vjust = -0.5, color="black", size=3,position = dodger) +
  theme(axis.title.x = element_text(size=18, face="bold", colour = "black"),    
        axis.title.y = element_text(size=18, face="bold", colour = "black"),    
        axis.text.x = element_text(size=18, face="bold", colour = "black"), 
        axis.text.y = element_text(size=18, face="bold", colour = "black"),
        strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 12, face="bold", colour = "black"),
        legend.title = element_blank(),legend.text = element_text(size=16))+
  scale_fill_manual(values=c("Founder variants"="pink", "Known founder variants"="red")) +
  scale_colour_manual(values=c("Founder variants"="pink", "Known founder variants"="red"), 
                      labels=c("Founder variants", "Known founder variants")) +
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1))

#### Merge the 2 plots to have figure 5
figure_5 <- ggarrange(graphs_RQ, graph_SAG + rremove("ylab"), ncol=2, common.legend = TRUE, legend="bottom",labels = c("A", "B"))
## Save in JPEG format
ggsave("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/Figure_5.jpeg"
       , width =12, height = 9, dpi=300,bg="white")

print("Done")
