################################
#
# GAGNON Laurence August 2023
################################
#
# September 2024 modified by MICHEL Elisa
#
# Do Figure in article Founder variants: Supp Fig.2
################################
rm(list = ls());
### Packages
library('ggplot2')
library("ggpubr")
library("dplyr")
library("stringi")
library("stringr")
library("tidyr")
library("patchwork")

## Get arguments QcP or SLSJ or UQc
cli <- commandArgs(trailingOnly = TRUE)
args <- strsplit(cli, "=", fixed = TRUE)
print(args)
print("Imput")
## Get the directory
work_dir_LPCL <- paste0("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Imput/",args,"/")
setwd(work_dir_LPCL)
## Get info about Variants known
variants_known <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/variants_reported_in_SLSJ.txt")
## Open file with all data
df_final <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Get the known variants ID
SNP_variants_known <- df_final[df_final$ID %in% variants_known[,1],1]
#### Get Table for position GRchr37
df_all_pos_chr37 <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/enriched_variants/Variants_to_keep_Imput_POS_38_and_37.txt",header=F)
df_all_pos_chr37$POS_38 <- str_split_fixed(df_all_pos_chr37$V4, "-",2)[,2]
df_all_pos_chr37 <- df_all_pos_chr37[,c(6,1,2)]
names(df_all_pos_chr37) <- c("POS_38","CHROM","POS_37")
#### Organize the table thanks to the file's name
name_file=list.files()
name_file=as.data.frame(name_file[grep(".sharing.by.pos",name_file)])
name_file <- as.data.frame(str_split_fixed(name_file[,1], "@", 3))
name_file <- cbind(name_file, str_split_fixed(name_file[,3], "_",2)[,1])
name_file <- cbind(name_file, str_split_fixed(name_file[,1], "",12)[,12])
name_file <- cbind(name_file, str_split_fixed(name_file[,2], "",4)[,4])
name_file <- cbind(name_file, str_split_fixed(name_file[,4], "",11)[,11])
name_file <- name_file[,-1:-4]
names(name_file) <- c('var_chr','num_var','num_ind')
name_file_clean <- name_file[!duplicated(name_file), ]
#### Set the colors
color <- as.data.frame(c('blue'))
color_full <- color[rep(seq_len(nrow(color)), each = nrow(name_file_clean)), ]
names(color) <- 'color'
#### List to save the graph
index_known=0
index_unknown=0
plot_list_known_Imput = list()
plot_list_unknown_Imput = list()
#### Number of individual with variant
name_file_clean$num_ind <- as.numeric(name_file_clean$num_ind) 
## Remove ind without IBD - beginning
ind_to_remove <- read.table("/lustre03/project/6033529/cartagene/data/cleanup_2022/Girard432553_axiom.final.rs_geno0.05_hwe10-6_snps_mind0.05.grch38.fam")
ind_to_remove <- ind_to_remove[ind_to_remove$V1 != "11134013" & ind_to_remove$V1 != "11125043",]
ind_to_remove <- ind_to_remove$V1
for (m in 1:22){
  table_ind <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Imput/',args,'/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',m,'.txt'), header=F)
  for (i in 1:length(table_ind)){
    table_ind_snps <- table_ind[,i]
    table_ind_snps <- table_ind_snps[complete.cases(table_ind_snps)]
    table_ind_snps_clean <- table_ind_snps[!table_ind_snps %in% ind_to_remove]
    nb_ind <- length(table_ind_snps_clean) -1
    name_file_clean[name_file_clean$var_chr == m & name_file_clean$num_var == i,4 ] <- nb_ind
  }
}
## Rename the columns
names(name_file_clean) <- c('var_chr','num_var','num_ind','real_num_ind')
## Calculate the pairs
name_file_clean$pairs_vec <- (name_file_clean$real_num_ind * (name_file_clean$real_num_ind - 1) / 2)
## Create table df_info with variants
df_info <- cbind(name_file_clean, color_full)

for (var in 1:nrow(df_info)) {
  ## Add postion in zero with bim file
  bim_file <- read.table("/lustre03/project/6033529/cartagene/data/cleanup_2022/gsa.17k.5300.4224.760.omni.bim")
  ## Clean Chr X and Y
  bim_file <- bim_file[bim_file$V1 <23,]
  bim_file <- as.data.frame(bim_file[,c(1,4)])
  names(bim_file) <- c("chr","Pos_genome")
  print(paste0("info: chr_",df_info[var, 1]," and num_var_",df_info[var, 2]))
  ### Open file with pos to get pos and chr corresponding for redline
  info_pos_VTA <- read.table(paste0('/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/results/Imput/',args,'/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',df_info[var, 1],'.txt'))
  # Get pos and chromosome of the variant
  pos_chr38 <- info_pos_VTA[1,as.numeric(df_info[var,2])]
  pos_chr38 <- as.data.frame(str_split_fixed(pos_chr38, ":", 3))
  chr <- as.numeric(str_split_fixed(pos_chr38[,1], "", 4)[,4])
  my_position_chr38 <- as.numeric(pos_chr38[,2])
  # Get the position in chr37
  pos_chr37 <- as.numeric(df_all_pos_chr37[df_all_pos_chr37$POS_38 == my_position_chr38,3])
  ### Open file IBD
  files = paste0('variant_chr',df_info[var, 1],"@num",df_info[var, 2],"@individual",df_info[var, 3],"_IBD_data_with_indiviudals_in_common_for_patho_variant_chr", c(1:22),"_0.5.sharing.by.pos")
  #data_list = lapply(files, read.table, colClasses="numeric",header=FALSE)
  data_list = lapply(files, function(x) {
    tryCatch(read.table(x, header = F, colClasses="numeric"), error=function(e) NULL)})
  pat_sharing <- do.call(rbind, data_list)
  head(pat_sharing)
  names(pat_sharing) <- c("chr","Pos_genome","pairs")
  ## Add zero positions
  bim_file$SNP <- paste0("chr",bim_file$chr,":",bim_file$Pos_genome)
  pat_sharing$SNP <- paste0("chr",pat_sharing$chr,":",pat_sharing$Pos_genome)
  df <- merge(bim_file,pat_sharing,by='SNP', all.x=T, all.y=F)
  df[is.na(df)] <- 0
  df <- df[,c(2,3,6)]
  names(df) <- c("chr","Pos_genome","pairs")
  df = df[with(df, order(chr, Pos_genome)), ]
  ### Add info
  df$prop <- df$pairs / df_info[var, 5]
  df$pos <- c(1:length(df[,1]))
  ## Get pos on graph on redline
  pos_variant <- as.numeric(which(df$chr == chr & abs(df$Pos_genome-pos_chr37)==min(abs(df[df$chr == chr,]$Pos_genome-pos_chr37))))
  ### Do graphic
  xaxis <- c()
  vlines <- c()
  for (ichr in 1:22){
    mymin <- min(df$pos[df$chr==ichr])
    mymax <- max(df$pos[df$chr==ichr])
    mypos <- ((mymax-mymin) / 2) +  mymin
    xaxis <- c(xaxis, mypos )
    vlines<- c(vlines, mymax)
  }
  SNP_variant <- paste0("chr",chr,"_",my_position_chr38,"_",(str_split_fixed(pos_chr38$V3,':',2)[1]),
                        "_",(str_split_fixed(pos_chr38$V3,':',2)[2]))
  Clinvar_ID <- df_final[df_final$SNP == SNP_variant,2]
  ## Do the graph
  p <- ggplot(df, aes(x=pos, y=prop)) +
    geom_line(linewidth=0.5, col= df_info[var, 6]) +
    theme_grey(base_size = 14) +
    ylim(c(0,1)) +
    scale_x_continuous(name='Chromosome', breaks=xaxis, labels=c(1:22)) +
    geom_vline(xintercept=vlines,linetype="dashed", color = "darkgrey", linewidth=0.25) +
    geom_vline(xintercept=pos_variant, linetype="longdash", color = "red", linewidth=0.25,alpha = 0.70) +
    labs(title=paste0("chr",df_info[var, 1], ":",my_position_chr38,":",str_split_fixed(pos_chr38[,3], ":",2)[1],":",
                      str_split_fixed(pos_chr38[,3], ":",2)[2]," (n = ",df_info[var, 3],") ClinVar ID: ",Clinvar_ID), y = "Proportion of pairs sharing")+
    theme(axis.text.x = element_text(size = 14,face="bold"),
          axis.title.y = element_text(vjust= 1.8, size = 16),
          axis.title.x = element_text(vjust= -0.5, size = 16),
          axis.text.y = element_text(size = 14,face="bold"))
  ## Save the graph
  if (SNP_variant %in% SNP_variants_known){
    index_known=index_known+1
    ### Save file
    plot_list_known_Imput[[index_known]] = p
  } else if (df[pos_variant,4] >= 0.5 & df_info[var,3] >=30) {
    index_unknown=index_unknown+1
    plot_list_unknown_Imput[[index_unknown]] = p
  } else {
    print("Neither known or founder")
  }
}

## Save the list of known and unknown Imput
print("Saving list")
outfile_1=paste0("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/scripts/Scripts_graphs_PEAK/list_plots_Imput_unknown_",args)
outfile_2=paste0("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/scripts/Scripts_graphs_PEAK/list_plots_Imput_known_",args)
saveRDS(plot_list_unknown_Imput,outfile_1)
saveRDS(plot_list_known_Imput,outfile_2)
print("Done")
