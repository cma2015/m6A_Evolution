## current wroking directory
setwd("~/a2z/m6A_13spp/")

## Options Settings
options(stringsAsFactors = F)
options(scipen=200)

## loading libraries
library(ggforce)
library(dplyr)
library(ggplot2)
library(patchwork)

## reload data
load("notebooks/RData/01_m6A_peaks_genes.RData")
load("notebooks/RData/02_detailed_gene_attributes.RData")

comlab_name <- "m6A methylation ratio (%)"

subgenome_divergence <- function(dataMat, col_num, species_list, idname){
  out_tmp <- list()
  dataMat <- dataMat[, col_num]
  dataMat_index <- apply(dataMat, 1, function(x){
    return(x[1]%in%exp_list[[species_list[1]]]&
             x[2]%in%exp_list[[species_list[2]]])
  })
  dataMat_tmp <- dataMat[dataMat_index, ]
  dataMat_m6A_type <- apply(dataMat_tmp, 1, function(x){
    if(x[1]%in%m6A_list[[species_list[1]]]){
      vec_tmp1 <- "1"
    }else{
      vec_tmp1 <- "0"
    }
    if(x[2]%in%m6A_list[[species_list[2]]]){
      vec_tmp2 <- "1"
    }else{
      vec_tmp2 <- "0"
    }
    return(paste0(c(vec_tmp1, vec_tmp2), collapse = ":"))
  })
  dataMat_m6A_type[dataMat_m6A_type=="1:0"|dataMat_m6A_type=="0:1"] <- "DM"
  dataMat_m6A_type[dataMat_m6A_type=="0:0"] <- "NM"
  dataMat_m6A_type[dataMat_m6A_type=="1:1"] <- "IM"
  out_tmp[[1]] <- as.data.frame(table(dataMat_m6A_type))
  colnames(out_tmp[[1]]) <- c("Type", "Count")
  out_tmp[[1]]$Name <- idname
  dataMat_tmp$Type <- dataMat_m6A_type
  dataMat_tmp$Name <- idname
  out_tmp[[2]] <- dataMat_tmp
  return(out_tmp)
}

## evolution origin ####
## read the subgenome information (gene pairs) ####
gar_ghi <- read.table("associated_data/subgenome_clean/gar_ghi_subA_24468.txt", sep = "\t", stringsAsFactors = F, header = T)
ghiAD <- read.table("associated_data/subgenome_clean/ghi_subA_subD_26725.txt", sep = "\t", stringsAsFactors = F, header = T)
pvu_gma <- read.table("associated_data/subgenome_clean/pvu_gma_pairs_20195.txt", sep = "\t", stringsAsFactors = F, header = F)
sbi_zma <- read.table("associated_data/subgenome_clean/zma_sbi_pairs_17898.txt", sep = "\t", stringsAsFactors = F, header = T)
wheat_subA <- read.table("associated_data/subgenome_clean/wheat_subA_17701.txt", sep = "\t", stringsAsFactors = F)
wheat_subB <- read.table("associated_data/subgenome_clean/wheat_subB_17375.txt", sep = "\t", stringsAsFactors = F)
wheat_subD <- read.table("associated_data/subgenome_clean/wheat_subD_18475.txt", sep = "\t", stringsAsFactors = F)
# wheat_ABD <- read.table("associated_data/subgenome_clean/wheat_ABD_triads.txt", sep = "\t", stringsAsFactors = F)

## extract m6A gene pairs ####
Pvu_Gma1_list <- subgenome_divergence(dataMat = pvu_gma, col_num = c(1,2), species_list = c("pvu","gma"), idname="Pvu_Gma1")
Pvu_Gma2_list <- subgenome_divergence(dataMat = pvu_gma, col_num = c(1,3), species_list = c("pvu","gma"), idname="Pvu_Gma2")
# Gma1_Gma2_list <- subgenome_divergence(dataMat = pvu_gma, col_num = c(2,3), species_list = c("gma","gma"), idname="Gma1_Gma2")
Sbi_Zma1_list <- subgenome_divergence(dataMat = sbi_zma, col_num = c(1,2), species_list = c("sbi","zma"), idname="Sbi_Zma1")
Sbi_Zma2_list <- subgenome_divergence(dataMat = sbi_zma, col_num = c(1,3), species_list = c("sbi","zma"), idname="Sbi_Zma2")
# Zma1_Zma2_list <- subgenome_divergence(dataMat = sbi_zma, col_num = c(2,3), species_list = c("zma","zma"), idname="Zma1_Zma2")

Gar_GhiA_list <- subgenome_divergence(dataMat = gar_ghi, col_num = c(1,2), species_list = c("gar","ghi"), idname="Gar_GhiA")
# GhiA_GhiD_list <- subgenome_divergence(dataMat = ghiAD, col_num = c(1,2), species_list = c("ghi","ghi"), idname="GhiA_GhiD")
Ata_TaeD_list <- subgenome_divergence(dataMat = wheat_subD, col_num = c(1,2), species_list = c("ata","tae"), idname="Ata_TaeD")
TdiA_TaeA_list <- subgenome_divergence(dataMat = wheat_subA, col_num = c(1,2), species_list = c("tdi","tae"), idname="TdiA_TaeA")
TdiB_TaeB_list <- subgenome_divergence(dataMat = wheat_subB, col_num = c(1,2), species_list = c("tdi","tae"), idname="TdiB_TaeB")
# TaeA_TaeB_list <- subgenome_divergence(dataMat = wheat_ABD, col_num = c(1,2), species_list = c("tae","tae"), idname="TaeA_TaeB")

rownames(gene_attributes_df) <- gene_attributes_df$Gene
div_m6A_PSs <- div_m6A_ratio <- div_m6A_spp <- vector()

subgenome_df <- rbind(Pvu_Gma1_list[[1]], Pvu_Gma2_list[[1]], Sbi_Zma1_list[[1]], Sbi_Zma2_list[[1]], 
                      Gar_GhiA_list[[1]], Ata_TaeD_list[[1]], TdiB_TaeB_list[[1]], TdiA_TaeA_list[[1]])

tmp_merge <- rep("Identical pattern", nrow(subgenome_df))
tmp_merge[subgenome_df$Type%in%c("DM")] <- "Diverged pattern"
subgenome_df$mergeType <- tmp_merge



## Figure 4B ####
two_spp_seq_var1 <- read.table("associated_data/figure_data/sbi_df.csv", sep = ",", header = T)
two_spp_seq_var2 <- read.table("associated_data/figure_data/pvu_df.csv", sep = ",", header = T)
two_spp_seq_var3 <- read.table("associated_data/figure_data/gar_df.csv", sep = ",", header = T)
two_spp_seq_var4 <- read.table("associated_data/figure_data/ata_df.csv", sep = ",", header = T)
two_spp_seq_var5 <- read.table("associated_data/figure_data/tdi_df.csv", sep = ",", header = T)
colnames(two_spp_seq_var2) <- colnames(two_spp_seq_var3) <- 
  colnames(two_spp_seq_var4) <- colnames(two_spp_seq_var5) <- 
  colnames(two_spp_seq_var1)
two_spp_seq_var <- rbind(two_spp_seq_var1, two_spp_seq_var2, two_spp_seq_var3,
                         two_spp_seq_var4, two_spp_seq_var5)
two_spp_seq_var$allsites <- paste0(two_spp_seq_var$variable, "_", 
                                   two_spp_seq_var$position)
two_spp_seq_var$allsites <- factor(two_spp_seq_var$allsites, levels = unique(two_spp_seq_var$allsites))
two_spp_seq_var$value <- 1-two_spp_seq_var$value

two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="sbi_zma1"] <- "Sbi_Zma1"
two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="sbi_zma2"] <- "Sbi_Zma2"
two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="pvu_gma1"] <- "Pvu_Gma1"
two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="pvu_gma2"] <- "Pvu_Gma2"
two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="garA_ghiA"] <- "Gar_GhiA"
two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="ataD_taeD"] <- "Ata_TaeD"
two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="tdiA_taeA"] <- "TdiA_TaeA"
two_spp_seq_var$subgenome[two_spp_seq_var$subgenome=="tdiB_taeB"] <- "TdiB_TaeB"
two_spp_seq_var$subgenome<-factor(two_spp_seq_var$subgenome, 
                                  levels=c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2", 
                                           "Gar_GhiA", "Ata_TaeD", "TdiA_TaeA", "TdiB_TaeB"))

two_spp_seq_var %>% group_by(allsites, type, subgenome) %>% 
  summarise("medianValues"=median(value), "meanValues"=mean(value)) %>% 
  as.data.frame() ->  two_spp_seq_var_merge
two_spp_seq_var_merge$type[two_spp_seq_var_merge$type==0] <- "NM"
two_spp_seq_var_merge$type[two_spp_seq_var_merge$type==1] <- "DM"
two_spp_seq_var_merge$type[two_spp_seq_var_merge$type==2] <- "IM"

two_spp_seq_var_merge$type <- factor(two_spp_seq_var_merge$type, levels = c("IM", "DM", "NM"))
two_spp_seq_var_merge$meanValues <- two_spp_seq_var_merge$meanValues*100

slideMean<-function(x,windowsize=3,slide=1){
  idx1 <- seq(1,length(x) - windowsize+1, by=slide);
  idx1 + windowsize -> idx2;
  # idx2[idx2>(length(x)+1)] <- length(x)+1;
  c(0, cumsum(x)) -> cx;
  return(c(mean(x[1:2]), (cx[idx2]-cx[idx1])/windowsize, 
           mean(x[(length(x)-1):length(x)])));
}

window_df <- list()
for(ii in unique(two_spp_seq_var_merge$subgenome)){
  for(jj in unique(two_spp_seq_var_merge$type)){
    sub_df <- two_spp_seq_var_merge[two_spp_seq_var_merge$subgenome==ii&two_spp_seq_var_merge$type==jj,]
    sub_df$windowValue <- slideMean(sub_df$meanValues)
    window_df[[paste0(ii, jj)]] <- sub_df
  }
}

window_df_data <- do.call(rbind, window_df) 
# pdf("Figure21.pdf", height = 3, width = 14)
ggplot(data = window_df_data, aes(x=allsites, y=windowValue, color=type, group=type))+
  geom_line(size=0.8) + #geom_point(size=0.8) +
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab("") + 
  facet_grid(.~subgenome, scales = "free_x", space = "free_x") +
  scale_color_brewer(palette = "Dark2") +
  geom_rect(aes(xmin=1, xmax=10.5, ymin=-5, ymax=-1),
            fill="#AAAAAA", color=NA) +
  geom_rect(aes(xmin=10.5, xmax=20.5, ymin=-6, ymax=0),
            fill="#767777", color=NA) +
  geom_rect(aes(xmin=20.5, xmax=30, ymin=-5, ymax=-1),
            fill="#AAAAAA", color=NA)
# dev.off()
# write.table(window_df_data[, c(3,1,2,6)], file = "Figure4E.txt", quote = F, sep = "\t", row.names = F)

## Figure 4C and Supplementary Fig. 5 Ka Ks ####
spe_pairs <- c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2", "Gar_GhiA", "Ata_TaeD", "TdiA_TaeA", "TdiB_TaeB")
two_spp_kaks <- read.table("associated_data/figure_data/Figure22.csv", sep = ",", header = T)
two_spp_kaks <- two_spp_kaks[two_spp_kaks$variable%in%c("ka", "ks"), ]
two_spp_kaks$X3kind <- factor(two_spp_kaks$X3kind, levels = c("IM", "DM", "NM"))
two_spp_kaks$spe[two_spp_kaks$spe=="sbi_zma1"] <- "Sbi_Zma1"
two_spp_kaks$spe[two_spp_kaks$spe=="sbi_zma2"] <- "Sbi_Zma2"
two_spp_kaks$spe[two_spp_kaks$spe=="pvu_gma1"] <- "Pvu_Gma1"
two_spp_kaks$spe[two_spp_kaks$spe=="pvu_gma2"] <- "Pvu_Gma2"
two_spp_kaks$spe[two_spp_kaks$spe=="gar_ghiA"] <- "Gar_GhiA"
two_spp_kaks$spe[two_spp_kaks$spe=="ata_taeD"] <- "Ata_TaeD"
two_spp_kaks$spe[two_spp_kaks$spe=="tdi_taeA"] <- "TdiA_TaeA"
two_spp_kaks$spe[two_spp_kaks$spe=="tdi_taeB"] <- "TdiB_TaeB"
# two_spp_kaks <- two_spp_kaks[two_spp_kaks$spe%in%c(
#   "Gar_GhiA", "Ata_TaeD", "TdiA_TaeA", "TdiB_TaeB"
# ), ]
two_spp_kaks$spe <- factor(two_spp_kaks$spe,
                           levels = intersect(c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2",
                                                "Gar_GhiA", "TdiA_TaeA", "TdiB_TaeB", "Ata_TaeD"),
                                              unique(two_spp_kaks$spe)))

Ka_plot <- ggplot(two_spp_kaks %>% filter(spe%in%spe_pairs, variable=="ka"), aes(x=X3kind, y=value, fill=X3kind)) +
  geom_violin(scale="area") +
  geom_boxplot(width=0.12, fill="#7F7F7F", outlier.colour = NA)+ #
  scale_fill_manual(values = c("#379677", "#C6612D", "#6B6BA2")) +
  theme_bw(base_size = 12)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())+
  xlab("") + ylab("Ka") +
  coord_cartesian(ylim = c(0, 0.27)) + facet_grid(.~spe, scales = "free_x") #facet_grid(.~, scales = "free_x", space = "free_x") 
Ks_plot <- ggplot(two_spp_kaks %>% filter(spe%in%spe_pairs, variable=="ks"), aes(x=X3kind, y=value, fill=X3kind)) +
  geom_violin(scale="area") + ## scale="width" geom_boxplot(outlier.color = "white")+scale_fill_brewer(palette = "Dark2") +
  geom_boxplot(width=0.12, fill="#7F7F7F", outlier.colour = NA)+ #
  scale_fill_manual(values = c("#379677", "#C6612D", "#6B6BA2")) +
  theme_bw(base_size = 12)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())+
  xlab("") + ylab("Ks") +
  coord_cartesian(ylim = c(0, 0.95)) + facet_grid(.~spe, scales = "free_x")

for(ii in unique(two_spp_kaks$spe)){
  sub_tmp_df <- two_spp_kaks[two_spp_kaks$spe==ii, ]
  sub_tmp_df <- sub_tmp_df[sub_tmp_df$variable=="ka", ]
  cat(ii, t.test(sub_tmp_df[sub_tmp_df$X3kind=="IM", "value"], 
                 sub_tmp_df[sub_tmp_df$X3kind=="DM", "value"], alternative = "less")$p.value,
      t.test(sub_tmp_df[sub_tmp_df$X3kind=="IM", "value"], 
             sub_tmp_df[sub_tmp_df$X3kind=="NM", "value"], alternative = "less")$p.value, "\n")
}

# pdf("Figure4C_FigureS5.pdf", height = 7, width = 12)
Ka_plot + Ks_plot + plot_layout(ncol = 1, guides = "collect")
# dev.off()


## Figure 4D and Supplementary Fig. 6 sequence variation ####
two_spp_3utr_var1 <- read.table("associated_data/figure_data/sbi_region_identity.csv", sep = ",", header = T)
two_spp_3utr_var2 <- read.table("associated_data/figure_data/pvu_region_identity.csv", sep = ",", header = T)
two_spp_3utr_var3 <- read.table("associated_data/figure_data/gar_region_identity.csv", sep = ",", header = T)
two_spp_3utr_var4 <- read.table("associated_data/figure_data/ata_region_identity.csv", sep = ",", header = T)
two_spp_3utr_var5 <- read.table("associated_data/figure_data/tdi_region_identity.csv", sep = ",", header = T)
colnames(two_spp_3utr_var2) <- colnames(two_spp_3utr_var3) <- 
  colnames(two_spp_3utr_var4) <- colnames(two_spp_3utr_var5) <- 
  colnames(two_spp_3utr_var1)
two_spp_3utr_var <- rbind(two_spp_3utr_var1, two_spp_3utr_var2, two_spp_3utr_var3,
                          two_spp_3utr_var4, two_spp_3utr_var5)

two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="sbi_zma1"] <- "Sbi_Zma1"
two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="sbi_zma2"] <- "Sbi_Zma2"
two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="pvu_gma1"] <- "Pvu_Gma1"
two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="pvu_gma2"] <- "Pvu_Gma2"
two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="garA_ghiA"] <- "Gar_GhiA"
two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="ataD_taeD"] <- "Ata_TaeD"
two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="tdiA_taeA"] <- "TdiA_TaeA"
two_spp_3utr_var$subgenome[two_spp_3utr_var$subgenome=="tdiB_taeB"] <- "TdiB_TaeB"
# two_spp_3utr_var <- two_spp_3utr_var[two_spp_3utr_var$subgenome%in%c(
#   "Gar_GhiA", "Ata_TaeD", "TdiA_TaeA", "TdiB_TaeB"
# ), ]
two_spp_3utr_var$subgenome <- factor(two_spp_3utr_var$subgenome, 
                                     levels = intersect(c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2",   
                                                "Gar_GhiA", "TdiA_TaeA", "TdiB_TaeB", "Ata_TaeD"),
                                                unique(two_spp_3utr_var$subgenome)))
two_spp_3utr_var$type[two_spp_3utr_var$type==0] <- "NM"
two_spp_3utr_var$type[two_spp_3utr_var$type==1] <- "DM"
two_spp_3utr_var$type[two_spp_3utr_var$type==2] <- "IM"
two_spp_3utr_var <- two_spp_3utr_var[two_spp_3utr_var$type%in%c("IM", "DM"), ]
two_spp_3utr_var$type <- factor(two_spp_3utr_var$type, 
                                levels = intersect(c("IM", "DM", "NM"), unique(two_spp_3utr_var$type)))
two_spp_3utr_var$three_prime_UTR <- two_spp_3utr_var$three_prime_UTR*100

for(ii in unique(two_spp_3utr_var$subgenome)){
  sub_tmp_df <- two_spp_3utr_var[two_spp_3utr_var$subgenome==ii, ]
  cat(ii, t.test(sub_tmp_df[sub_tmp_df$type=="IM", "three_prime_UTR"], 
                 sub_tmp_df[sub_tmp_df$type=="DM", "three_prime_UTR"], alternative = "greater")$p.value,
      t.test(sub_tmp_df[sub_tmp_df$type=="IM", "three_prime_UTR"], 
             sub_tmp_df[sub_tmp_df$type=="NM", "three_prime_UTR"], alternative = "greater")$p.value, "\n")
}

# pdf("Figure4D_FigureS6.pdf", height = 4, width = 10)
ggplot(data = two_spp_3utr_var, aes(x=type, y=three_prime_UTR, fill=type))+
  geom_violin(scale="width") +
  geom_boxplot(width=0.12, fill="#7F7F7F", outlier.colour = NA)+
  scale_fill_manual(values = c("#379677", "#C6612D")) + #"#6B6BA2"
  theme_bw(base_size = 12)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())+
  xlab("") + ylab("3' UTR sequence Identity (%)") + facet_grid(.~subgenome)
# dev.off()


### Figure 4F autopolyploidization and allopolyploidization merge IM, DM, NM ####
poly_type <- rep("autopolyploidization", nrow(two_spp_kaks))
poly_type[two_spp_kaks$spe%in%c("Gar_GhiA", "Ata_TaeD","TdiA_TaeA", "TdiB_TaeB")] <- "allopolyploidization"
two_spp_kaks$poly_type <- poly_type
two_spp_kaks$poly_type <- factor(two_spp_kaks$poly_type, 
                                 levels = c("autopolyploidization", "allopolyploidization"))

two_spp_kaks %>% filter(spe%in%spe_pairs, variable=="ka") -> two_spp_ka
ka_merge_plot <- ggplot(two_spp_ka, aes(x=poly_type, y=value, fill=poly_type)) +
geom_violin() + geom_boxplot(width=0.05, fill="#7F7F7F", outlier.colour = NA)+ #
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 12)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "none")+
  xlab("") + ylab("Ka") +
  coord_cartesian(ylim = c(0, 0.2)) #facet_grid(.~, scales = "free_x", space = "free_x") 

two_spp_kaks %>% filter(spe%in%spe_pairs, variable=="ks") -> two_spp_ks
ks_merge_plot <- ggplot(two_spp_ks, aes(x=poly_type, y=value, fill=poly_type)) +
  geom_violin() + geom_boxplot(width=0.05, fill="#7F7F7F", outlier.colour = NA)+ #
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 12)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "none")+
  xlab("") + ylab("Ks") +
  coord_cartesian(ylim = c(0, 0.8))

wilcox.test(two_spp_ka[two_spp_ka$poly_type%in%c("autopolyploidization"), 
                       "value"],
            two_spp_ka[two_spp_ka$poly_type%in%c("allopolyploidization"),
                       "value"])

wilcox.test(two_spp_ks[two_spp_ks$poly_type%in%c("autopolyploidization"), 
                       "value"],
            two_spp_ks[two_spp_ks$poly_type%in%c("allopolyploidization"),
                       "value"])

# pdf("Figure4F.pdf", height = 3, width = 4)
ka_merge_plot + ks_merge_plot + plot_layout(ncol = 1, guides = "collect")
# dev.off()


## Figure 4E ####
pie_stats <- function(df, x0, y0, r0, r1, amount, explode, label_perc) {
  x0 <- enquo(x0)
  y0 <- enquo(y0)
  r0 <- enquo(r0)
  r1 <- enquo(r1)
  amount <- enquo(amount)
  explode <- enquo(explode)
  df %>%
    mutate(
      `x0` = !! x0,
      `y0` = !! y0,
      `r0` = !! r0,
      `r1` = !! r1,
      `explode` = !! explode
    ) %>%
    group_by(x0, y0) %>%
    mutate(end = cumsum(!! amount) / sum(!! amount) * 2 * pi + pi-pi) %>% #pi/2
    mutate(start = lag(end, default = pi-pi)) %>% #pi/2
    ungroup() %>%
    mutate(
      x_lab = (!! x0) +
        ((!! r0) + label_perc * ((!! r1) - (!! r0)) + (!! explode)) *
        sin((end + start) / 2),
      y_lab = (!! y0) +
        ((!! r0) + label_perc * ((!! r1) - (!! r0)) + (!! explode)) *
        cos((end + start) / 2)
    )
}

subgenome_plots <- list()
for(ii in unique(subgenome_df$Name)){
  subgenome_df_sub <- subgenome_df %>% filter(Name==ii) %>% 
    group_by(Name, mergeType)%>% summarise("Value"=sum(Count)) %>% 
    as.data.frame()
  # subgenome_df_sub$Type <- as.character(subgenome_df_sub$Type)
  # subgenome_df_sub <- subgenome_df_sub[c(2,3,1), ]
  # subgenome_df_sub$Type <- factor(as.character(subgenome_df_sub$Type), levels = c("NM", "DM","IM"))
  subgenome_plots[[ii]] <- ggplot(
    data = subgenome_df_sub %>% pie_stats(0, 0, 0, 1, Value, .1 * (mergeType == "Diverged pattern"), .5)
    ) +
    ggforce::geom_arc_bar(aes(x0 = x0, y0 = y0, r0 = r0, r = r1,
                              start = start, end = end, explode = explode,
        fill = factor(mergeType, levels = c("Diverged pattern", "Identical pattern")))) +
    geom_text(aes(x = x_lab, y = y_lab, label = scales::percent(Value/sum(Value), accuracy = 0.1)), size = 3) +
    coord_equal() + theme_void() + 
    labs(fill = expression("m"^"6"*"A patterns"), title = ii) + 
    scale_fill_manual(values = c("#F0861B", "#9FA0A0"))  #c("#FF8900", "#086FA1")
}

# pdf("Figure4E.pdf", height = 5, width = 10)
wrap_plots(subgenome_plots, ncol = 2, guides = "collect")
# dev.off()


## 3'UTR length ####
region_length_plot <- function(gene_pairs, species_name, first_region, second_region){
  pvu_gma_df_tmp <- apply(gene_pairs,1,function(x){
    tmp_df1 <- first_region[first_region$V1%in%x[1], 3:4]
    tmp_df1$Species <- species_name[1]
    tmp_df2 <- second_region[second_region$V1%in%x[2], 3:4]
    tmp_df2$Species <- species_name[2]
    tmp_index <- tmp_df1[,2]>100&tmp_df2[,2]>100
    rbind(tmp_df1[tmp_index, ], tmp_df2[tmp_index, ])
  })
  ##
  pvu_gma_df_tmp <- do.call(rbind, pvu_gma_df_tmp)
  colnames(pvu_gma_df_tmp) <- c("Region", "Length", "Species")
  pvu_gma_df_tmp$Length <- log2(pvu_gma_df_tmp$Length+1)
  pvu_gma_df_tmp$Species <- factor(pvu_gma_df_tmp$Species, levels = species_name)
  pvu_gma_df_tmp$Region <- factor(pvu_gma_df_tmp$Region, levels = c("5UTR", "CDS", "3UTR"))
  ##
  p <- ggplot(pvu_gma_df_tmp, aes(x=Region, y=Length, fill=Species)) +
    geom_boxplot(outlier.color = "white")+
    scale_fill_manual(values = c("#0072B2", "#E69F00"))+
    geom_point(position=position_jitterdodge(jitter.width = 0.05)) +
    theme_bw(base_size = 14)+
    theme(axis.text.x= element_text(angle=15,hjust = 1,vjust = 1, colour="black",size=10),
          axis.line = element_line(colour = "black",size=0.5),
          legend.text=element_text(colour="black", size=16),
          legend.title=element_text(colour="black", size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    stat_compare_means(method = 'wilcox.test',label.y = max(pvu_gma_df_tmp$Length)+1, size=2) +
    xlab("") + ylab("Percentage")
  p
}

region_difference  <- function(gene_pairs, species_name, first_region, second_region, type_name){
  pvu_gma_df_tmp <- apply(gene_pairs,1,function(x){
    tmp_df1 <- first_region[first_region$V1%in%x[1], ]
    tmp_df2 <- second_region[second_region$V1%in%x[2], ]
    tmp_index <- tmp_df1[,4]>0&tmp_df2[,4]>0
    tmp_df1 <- tmp_df1[tmp_index, ]
    tmp_df2 <- tmp_df2[tmp_index, ]
    if("test"=="notest"){
      tmp <- "tmp"
      ### test
      # if(x[1]%in%m6A_list[[species_name[1]]]&x[2]%in%m6A_list[[species_name[2]]]){
      #   tmp_df1[,4] <- abs(tmp_df1[,4]-tmp_df2[,4])
      # } else if (x[1]%in%m6A_list[[species_name[1]]]){
      #   tmp_df1[,4] <- tmp_df1[,4]-tmp_df2[,4]
      # } else if (x[2]%in%m6A_list[[species_name[2]]]){
      #   tmp_df1[,4] <- tmp_df2[,4] - tmp_df1[,4]
      # } else {
      #   tmp_df1[,4] <- abs(tmp_df1[,4]-tmp_df2[,4])
      # }
    }else{
      tmp_df1[,4] <- tmp_df1[,4]-tmp_df2[,4]
    }
    tmp_df1[,2] <- tmp_df2[,1]
    tmp_df1$Type <- type_name
    tmp_df1
  })
  ##
  pvu_gma_df_tmp <- do.call(rbind, pvu_gma_df_tmp)
  colnames(pvu_gma_df_tmp) <- c("Gene1", "Gene2", "Region", "Difference", "Type")
  pvu_gma_df_tmp
}

region_length_change <- function(gene_pair_df, species_name, 
                                 first_region, second_region, type_name){
  yes_yes <- gene_pair_df[gene_pair_df$Type=="IM",]
  yes_no <- gene_pair_df[gene_pair_df$Type=="DM",]
  no_no <- gene_pair_df[gene_pair_df$Type=="NM",]
  pp_out <- list()
  pp_out[[1]] <- region_difference(gene_pairs = yes_yes, species_name, first_region, second_region, type_name = "IM")
  pp_out[[2]] <- region_difference(yes_no, species_name, first_region, second_region, type_name = "DM")
  pp_out[[3]] <- region_difference(no_no, species_name, first_region, second_region, type_name = "NM")
  pp_out <- do.call(rbind, pp_out)
  pp_out$Species <- type_name
  pp_out
}

## 
gar_region_df <- read.table("associated_data/pep_cds/region_length/gar.txt", sep = "\t", stringsAsFactors = F)
ghi_region_df <- read.table("associated_data/pep_cds/region_length/ghi.txt", sep = "\t", stringsAsFactors = F)

ata_region_df <- read.table("associated_data/pep_cds/region_length/ata.txt", sep = "\t", stringsAsFactors = F)
tdi_region_df <- read.table("associated_data/pep_cds/region_length/tdi.txt", sep = "\t", stringsAsFactors = F)
tae_region_df <- read.table("associated_data/pep_cds/region_length/tae.txt", sep = "\t", stringsAsFactors = F)

pvu_region_df <- read.table("associated_data/pep_cds/region_length/pvu.txt", sep = "\t", stringsAsFactors = F)
gma_region_df <- read.table("associated_data/pep_cds/region_length/gma.txt", sep = "\t", stringsAsFactors = F)

sbi_region_df <- read.table("associated_data/pep_cds/region_length/sbi.txt", sep = "\t", stringsAsFactors = F)
zma_region_df <- read.table("associated_data/pep_cds/region_length/zma.txt", sep = "\t", stringsAsFactors = F)

pvu_gma1_rlc <- region_length_change(gene_pair_df = Pvu_Gma1_list[[2]], species_name = c("pvu", "gma"),
                                     first_region = pvu_region_df, second_region = gma_region_df, 
                                     type_name = "Pvu_Gma1")
pvu_gma2_rlc <- region_length_change(gene_pair_df = Pvu_Gma2_list[[2]], species_name = c("pvu", "gma"),
                                     first_region = pvu_region_df, second_region = gma_region_df, 
                                     type_name = "Pvu_Gma2")
# gma1_gma2_rlc <- region_length_change(gene_pair_df = Gma1_Gma2_list[[2]], species_name = c("gma", "gma"),
#                                       first_region = gma_region_df, second_region = gma_region_df,
#                                       type_name = "Gma1_Gma2")
sbi_zma1_rlc <- region_length_change(gene_pair_df = Sbi_Zma1_list[[2]], species_name = c("sbi", "zma"),
                                     first_region = sbi_region_df, second_region = zma_region_df,
                                     type_name = "Sbi_Zma1")
sbi_zma2_rlc <- region_length_change(gene_pair_df = Sbi_Zma2_list[[2]], species_name = c("sbi", "zma"),
                                     first_region = sbi_region_df, second_region = zma_region_df,
                                     type_name = "Sbi_Zma2")
# zma1_zma2_rlc <- region_length_change(gene_pair_df = Zma1_Zma2_list[[2]], species_name = c("zma", "zma"),
#                                      first_region = zma_region_df, second_region = zma_region_df,
#                                      type_name = "Zma1_Zma2")
gar_ghiA_rlc <- region_length_change(gene_pair_df = Gar_GhiA_list[[2]], species_name = c("gar", "ghi"),
                                    first_region = gar_region_df, second_region = ghi_region_df, 
                                    type_name = "Gar_GhiA")
# ghiA_ghiD_rlc <- region_length_change(gene_pair_df = GhiA_GhiD_list[[2]], species_name = c("ghi", "ghi"),
#                                     first_region = ghi_region_df, second_region = ghi_region_df, 
#                                     type_name = "GhiA_GhiD")
ataD_taeD_rlc <- region_length_change(gene_pair_df = Ata_TaeD_list[[2]], species_name = c("ata", "tae"),
                                      first_region = ata_region_df, second_region = tae_region_df,
                                      type_name = "Ata_TaeD")
tdiA_taeA_rlc <- region_length_change(gene_pair_df = TdiA_TaeA_list[[2]], species_name = c("tdi", "tae"),
                                      first_region = tdi_region_df, second_region = tae_region_df, 
                                      type_name = "TdiA_TaeA")
tdiB_taeB_rlc <- region_length_change(gene_pair_df = TdiB_TaeB_list[[2]], species_name = c("tdi", "tae"),
                                      first_region = tdi_region_df,second_region = tae_region_df,
                                      type_name = "TdiB_TaeB")
# taeA_taeB_rlc <- region_length_change(gene_pair_df = TaeA_TaeB_list[[2]], species_name = c("tae", "tae"),
#                                       first_region = tae_region_df,second_region = tae_region_df,
#                                       type_name = "TaeA_TaeB")

region_plot <- rbind(pvu_gma1_rlc, pvu_gma2_rlc,
                     sbi_zma1_rlc, sbi_zma2_rlc,
                     gar_ghiA_rlc, ataD_taeD_rlc,
                     tdiA_taeA_rlc, tdiB_taeB_rlc)

region_plot$Type <- factor(region_plot$Type, levels = c("IM","DM", "NM"))
region_plot$Species <- factor(region_plot$Species, levels = unique(region_plot$Species) )
# write.table(region_plot, "region_difference_in_two_species.txt", quote = F, sep = "\t")
region_plot2  <- region_plot[region_plot$Region%in%"3UTR"&region_plot$Type%in%c("IM","DM"),]
region_plot2 <- region_plot2[region_plot2$Difference>= -500& region_plot2$Difference<= 500, ]


# pdf("Figure24.pdf", width = 25, height = 3.5)
#, "Gar_GhiA", "Ata_TaeD", "TdiA_TaeA", "TdiB_TaeB"
ggplot(data = region_plot2[region_plot2$Species%in%c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2"),],
       mapping = aes(Difference, fill=Type)) +
  geom_density(size=0.7) + facet_grid(.~Species, scales = "free")  +coord_cartesian(xlim = c(-500,500)) +
  theme_bw(base_size = 14)+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"))+
  xlab("The difference of 3'UTR")+ylab("Density") +
  scale_fill_manual(values = c("#547B77", "#AB8472"))
  #scale_color_manual(values = rev(c("grey", "#EA2427","#3A53A3"))) +
  # scale_linetype_manual(values = rev(c('solid', 'solid', 'twodash')))
# dev.off()
# write.table(region_plot2[, c(6,1,2,3,5,4)], file = "Figure4D.txt", quote = F, sep = "\t", row.names = F)


for(ii in unique(region_plot2$Species)){
  sub_df <- region_plot2[region_plot2$Species==ii,]
  type1_value <- sub_df[sub_df$Type=="IM", "Difference"]
  type2_value <- sub_df[sub_df$Type=="DM", "Difference"]
  type3_value <- sub_df[sub_df$Type=="NM", "Difference"]
  cat(ii, wilcox.test(type1_value, type2_value)$p.value, 
      wilcox.test(type2_value, type3_value)$p.value, "\n")
}