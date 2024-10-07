################################
#
# September 2024 by MICHEL Elisa
#
# Do Figure in article Founder variants: Fig.2
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

print("Graphs number of carriers (enriched variants)")

## Get enriched variants positions
variants_enriched <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Delete variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
variants_enriched <- variants_enriched[!variants_enriched$SNP %in% variants_to_freq,]
#### Put first row in column names 
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

####################### Enriched variants - SLSJ ############################
## Open file with the individuals in SLSJ
individual_file_SLSJ <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Data/Imput/SLSJ/CaG_Imput_VTA_SAG_enriched_in_WQ_chr_1.tfam'))
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
## Remove NR imput
NR_var <- read.table("/lustre03/project/6033529/saguenay_disease/results/correspondance_WGS_imput/snps_to_remove_seuil1_allcriteria_homswitch.txt")
PLPCL_Imput <- PLPCL_Imput[,!PLPCL_Imput[1,] %in% NR_var[,1]]
## Select only variants enriched in SLSJ
PLPCL_Imput <- PLPCL_Imput[,PLPCL_Imput[1,] %in% variants_enriched$SNP]
PLPCL_enriched_Imput <- header.true(PLPCL_Imput)
## Found missing variants enriched to add with WGS
variants_missing_SAG <- variants_enriched[!variants_enriched$SNP %in% colnames(PLPCL_enriched_Imput),1]
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
## Select only variants enriched in SLSJ
PLPCL_enriched_WGS <- PLPCL_WGS[,PLPCL_WGS[1,] %in% variants_missing_SAG]
PLPCL_enriched_WGS <- PLPCL_enriched_WGS[,PLPCL_enriched_WGS[1,] %in% variants_enriched$SNP]
PLPCL_enriched_WGS <- header.true(PLPCL_enriched_WGS)
## Merge all info between Imput and WGS
PLPCL_enriched_all_SAG <- cbind(PLPCL_enriched_WGS,PLPCL_enriched_Imput)
print(paste0("variants_for_SLSJ : ",ncol(PLPCL_enriched_all_SAG)))
######### Do the graph with inds with x variants enriched in SLSJ all enriched
#### Search for ind with more than one variant
ind_with_variants=data.frame(mat)
for (i in 1:nrow(individual_file_SLSJ)){
  index=1
  ind <- individual_file_SLSJ[i,1]
  ind_with_variants[index,i] <- ind
  for (l in 1:ncol(PLPCL_enriched_all_SAG)){
    if (ind %in% PLPCL_enriched_all_SAG[,l]){
      ind_with_variants[index+1,i] <- colnames(PLPCL_enriched_all_SAG[l])
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
# Calculate how many inds with eight variants
ind_with_eight_variant_or_more <- as.data.frame(ind_with_seven_variant_or_more[,!is.na(ind_with_seven_variant_or_more[8,])])
ind_with_eight_variant <- as.data.frame(ind_with_eight_variant_or_more[,is.na(ind_with_eight_variant_or_more[9,])])
# Calculate how many inds with nine variants
ind_with_nine_variant_or_more <- as.data.frame(ind_with_eight_variant_or_more[,!is.na(ind_with_eight_variant_or_more[9,])])
ind_with_nine_variant <- as.data.frame(ind_with_nine_variant_or_more[,is.na(ind_with_nine_variant_or_more[10,])])
# Calculate how many inds with ten variants
ind_with_ten_variant_or_more <- as.data.frame(ind_with_nine_variant_or_more[,!is.na(ind_with_nine_variant_or_more[10,])])
ind_with_ten_variant <- as.data.frame(ind_with_ten_variant_or_more[,is.na(ind_with_ten_variant_or_more[11,])])
# Calculate how many inds with eleven variants
ind_with_eleven_variant_or_more <- as.data.frame(ind_with_ten_variant_or_more[,!is.na(ind_with_ten_variant_or_more[11,])])
ind_with_eleven_variant <- as.data.frame(ind_with_eleven_variant_or_more[,is.na(ind_with_eleven_variant_or_more[12,])])

print(paste0("nb var with eleven or more SAG :",ncol(ind_with_eleven_variant_or_more)))
print(paste0("nb var with eleven SAG :",ncol(ind_with_eleven_variant)))

################################### Enriched variants - UQc ##############################################
## Open file with individuals in UQc
individual_file_RQ <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Data/Imput/RQ/CaG_Imput_VTA_RQ_enriched_in_WQ_chr_1.tfam'))
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
PLPCL_Imput_RQ <- PLPCL_Imput_RQ[,PLPCL_Imput_RQ[1,] %in% variants_enriched$SNP]
PLPCL_enriched_Imput_RQ <- header.true(PLPCL_Imput_RQ)
## Found missing variants enriched to add with WGS
variants_missing_RQ <- variants_enriched[!variants_enriched$SNP %in% colnames(PLPCL_enriched_Imput_RQ),1]
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
PLPCL_enriched_WGS_RQ <- PLPCL_WGS_RQ[,PLPCL_WGS_RQ[1,] %in% variants_missing_RQ]
PLPCL_enriched_WGS_RQ <- PLPCL_enriched_WGS_RQ[,PLPCL_enriched_WGS_RQ[1,] %in% variants_enriched$SNP]
PLPCL_enriched_WGS_RQ <- header.true(PLPCL_enriched_WGS_RQ)
## Merge all info between Imput and WGS
PLPCL_enriched_all_RQ <- cbind(PLPCL_enriched_WGS_RQ,PLPCL_enriched_Imput_RQ)
print(paste0("variants_for_RQ : ",ncol(PLPCL_enriched_all_RQ)))
######### Do the graph with inds with x variants enriched in UQc all enriched
#### Search for ind with more than one variant
ind_with_variants_RQ=data.frame(mat)
for (i in 1:nrow(individual_file_RQ)){
  index=1
  ind <- individual_file_RQ[i,1]
  ind_with_variants_RQ[index,i] <- ind
  for (l in 1:ncol(PLPCL_enriched_all_RQ)){
    if (ind %in% PLPCL_enriched_all_RQ[,l]){
      ind_with_variants_RQ[index+1,i] <- colnames(PLPCL_enriched_all_RQ[l])
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
# Calculate how many inds with seven variants
ind_with_seven_variant_or_more_RQ <- as.data.frame(ind_with_six_variant_or_more_RQ[,!is.na(ind_with_six_variant_or_more_RQ[7,])])
ind_with_seven_variant_RQ <- as.data.frame(ind_with_seven_variant_or_more_RQ[,is.na(ind_with_seven_variant_or_more_RQ[8,])])
# Calculate how many inds with eight variants
ind_with_eight_variant_or_more_RQ <- as.data.frame(ind_with_seven_variant_or_more_RQ[,!is.na(ind_with_seven_variant_or_more_RQ[8,])])
ind_with_eight_variant_RQ <- as.data.frame(ind_with_eight_variant_or_more_RQ[,is.na(ind_with_eight_variant_or_more_RQ[9,])])
# Calculate how many inds with nine variants
ind_with_nine_variant_or_more_RQ <- as.data.frame(ind_with_eight_variant_or_more_RQ[,!is.na(ind_with_eight_variant_or_more_RQ[9,])])
ind_with_nine_variant_RQ <- as.data.frame(ind_with_nine_variant_or_more_RQ[,is.na(ind_with_nine_variant_or_more_RQ[10,])])
# Calculate how many inds with nine variants
ind_with_ten_variant_or_more_RQ <- as.data.frame(ind_with_nine_variant_or_more_RQ[,!is.na(ind_with_nine_variant_or_more_RQ[10,])])
ind_with_ten_variant_RQ <- as.data.frame(ind_with_ten_variant_or_more_RQ[,is.na(ind_with_ten_variant_or_more_RQ[11,])])
# Calculate how many inds with eleven variants
ind_with_eleven_variant_or_more_RQ <- as.data.frame(ind_with_ten_variant_or_more_RQ[,!is.na(ind_with_ten_variant_or_more_RQ[11,])])
ind_with_eleven_variant_RQ <- as.data.frame(ind_with_eleven_variant_or_more_RQ[,is.na(ind_with_eleven_variant_or_more_RQ[12,])])

print(paste0("nb var with eleven or more RQ :",ncol(ind_with_eleven_variant_or_more_RQ)))
print(paste0("nb var with eleven RQ :",ncol(ind_with_eleven_variant_RQ)))

#### Collect info for all (enriched)
mat_2 = matrix(ncol = 3, nrow = 21) 
df_hist_SAG=data.frame(mat_2)
names(df_hist_SAG) <- c("Number_of_variants","Variants with RFD≥10% in SLSJ","Variants with RFD≥10% in UQc")
npop_SAG <- 3589
npop_RQ <- 21472
df_hist_SAG[1,1] <- 0
df_hist_SAG[2,1] <- 1
df_hist_SAG[3,1] <- 2
df_hist_SAG[4,1] <- 3
df_hist_SAG[5,1] <- 4
df_hist_SAG[6,1] <- 5
df_hist_SAG[7,1] <- 6
df_hist_SAG[8,1] <- 7
df_hist_SAG[9,1] <- 8
df_hist_SAG[10,1] <- 9
df_hist_SAG[1,2] <- ncol(ind_with_zero_variants)
df_hist_SAG[2,2] <- ncol(ind_with_one_variant)
df_hist_SAG[3,2] <- ncol(ind_with_two_variant)
df_hist_SAG[4,2] <- ncol(ind_with_three_variant)
df_hist_SAG[5,2] <- ncol(ind_with_four_variant)
df_hist_SAG[6,2] <- ncol(ind_with_five_variant)
df_hist_SAG[7,2] <- ncol(ind_with_six_variant)
df_hist_SAG[8,2] <- ncol(ind_with_seven_variant)
df_hist_SAG[9,2] <- ncol(ind_with_eight_variant)
df_hist_SAG[10,2] <- ncol(ind_with_nine_variant)
df_hist_SAG[11,2] <- ncol(ind_with_ten_variant)
df_hist_SAG[12,2] <- ncol(ind_with_eleven_variant)
df_hist_SAG[1,3] <- ncol(ind_with_zero_variants_RQ)
df_hist_SAG[2,3] <- ncol(ind_with_one_variant_RQ)
df_hist_SAG[3,3] <- ncol(ind_with_two_variant_RQ)
df_hist_SAG[4,3] <- ncol(ind_with_three_variant_RQ)
df_hist_SAG[5,3] <- ncol(ind_with_four_variant_RQ)
df_hist_SAG[6,3] <- ncol(ind_with_five_variant_RQ)
df_hist_SAG[7,3] <- ncol(ind_with_six_variant_RQ)
df_hist_SAG[8,3] <- ncol(ind_with_seven_variant_RQ)
df_hist_SAG[9,3] <- ncol(ind_with_eight_variant_RQ)
df_hist_SAG[10,3] <- ncol(ind_with_nine_variant_RQ)
df_hist_SAG <- df_hist_SAG[complete.cases(df_hist_SAG$Number_of_variants),]
## Transform results for histogram
flip_results <- melt(df_hist_SAG, id.vars=c('Number_of_variants'))
flip_results$prop <- ifelse(flip_results$variable == "Variants with RFD≥10% in SLSJ",
                            round(flip_results$value/npop_SAG,digits=3),round(flip_results$value/npop_RQ,digits=3))
print(flip_results)
#### Do the graph
ID <- 0:9
dodger = position_dodge(width = 0.9)
graph_SAG <- ggplot(flip_results, aes(x=Number_of_variants, y=prop,fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  labs(x = "Number of variants", y = "Proportion of individuals") + 
  scale_x_continuous("Number of variants", labels = as.character(ID), breaks = ID)+
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
  scale_fill_manual(values=c("Variants with RFD≥10% in SLSJ"="blue", "Variants with RFD≥10% in UQc"="green")) +
  scale_colour_manual(values=c("Variants with RFD≥10% in SLSJ"="blue", "Variants with RFD≥10% in UQc"="green"), 
                      labels=c("Variants with RFD≥10% in SLSJ", "Variants with RFD≥10% in UQc"))  +
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.4))
  
## Save in JPEG format
ggsave("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/Figure_2_enriched.jpeg"
       , width =12, height = 9, dpi=300,bg="white")

print("Done")
