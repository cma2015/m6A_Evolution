# gene information
# Time: 2021-06-03

## current wroking directory
setwd("E:/a2z_projects/03_m6A_13spp") #~/a2z/m6A_13spp/
options(stringsAsFactors = F)

## loading libraries
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## reload data
load("notebooks/RData/01_m6A_peaks_genes.RData")
load("notebooks/RData/02_detailed_gene_attributes.RData")

## Phylogenetic relationships among 34 species ####
SPECIES34 <- c("ath", "Bstricta", "Thassleriana", "gar", "ghi", "Egrandis", "pvu", "gma", 
               "Mtruncatula", "Ptrichocarpa", "Vvinifera", "Stuberosum", "sly", 
               "Sitalica", "Sviridis", "Pvirgatum", "sbi", "zma", "ata", "tdi", "tae",
               "Hvulgare", "Bdistachyon", "osa", "Eguineensis", "Spolyrhiza", "Atrichopoda", 
               "Paab","Gibi","Smoellendorffii", "ppa", "Mpolymorpha", "Creinhardtii", 
               "MpusillaCCMP1545")

OGs_34spp <- read.table("associated_data/Results_Jan16_34spp/Orthogroups/Orthogroups.tsv", 
                        sep = "\t", header = T, row.names = 1, stringsAsFactors = F)[, SPECIES34]
OGs_34spp_count <- t(apply(OGs_34spp, 1, function(x){
  return(sapply(strsplit(x, split = ", "), length))
}))

Brassicaceae <- c("ath", "Bstricta")
Brassicales <- list(Brassicaceae, "Thassleriana")
Malvaceae <- c("gar", "ghi")
BM <- list(unlist(Brassicales), unlist(Malvaceae))
Malvids <- list(unlist(BM), "Egrandis")
Phaseoleae <- c("pvu", "gma")
Fabaceae <- list(Phaseoleae, "Mtruncatula")
Fabids <- list("Ptrichocarpa", unlist(Fabaceae))
Eurosids <- list(unlist(Malvids), unlist(Fabids))
Rosids <- list(unlist(Eurosids), "Vvinifera")
Solanaceae <- c("Stuberosum", "sly")
Eudicots <- list(unlist(Rosids), Solanaceae)
Setaria <- c("Sitalica", "Sviridis")
Paniceae <- list(Setaria, "Pvirgatum")
Andropogoneae <- c("sbi", "zma")
PACMAD <- list(unlist(Paniceae), Andropogoneae)
Triticum <- c("tdi", "tae")
AT <- list(Triticum, "ata")
# Hordeeae (formerly Triticeae)
Hordeeae <- list(unlist(AT), "Hvulgare")
Pooideae <- list(unlist(Hordeeae), "Bdistachyon")
BOP <- list(unlist(Pooideae), "osa")
Poaceae <- list(unlist(PACMAD), unlist(BOP))
Commelinids <- list(unlist(Poaceae), "Eguineensis")
Monocots <- list(unlist(Commelinids), "Spolyrhiza")
Mesangiospermae <- list(unlist(Eudicots), unlist(Monocots))
Angiosperm <- list(unlist(Mesangiospermae), "Atrichopoda")
Acrogymnospermae <- c("Paab","Gibi")
Spermatophyta <- list(unlist(Angiosperm), Acrogymnospermae)
Tracheophyta <- list(unlist(Spermatophyta), "Smoellendorffii")
Stomatophyta <- list(unlist(Tracheophyta),"ppa")
Embryophyta <- list(unlist(Stomatophyta),"Mpolymorpha")
Chlorophyta <- c("Creinhardtii", "MpusillaCCMP1545")
Chloroplastida <- list(Chlorophyta, unlist(Embryophyta))

PS_details_names <- c("Chlorophyta", "Embryophyta", "Stomatophyta", "Tracheophyta", 
                      "Spermatophyta", "Acrogymnospermae", "Angiosperm", "Mesangiospermae", 
                      "Monocots", "Commelinids", "Poaceae", "BOP", "Pooideae", "Hordeeae", 
                      "AT", "Triticum", "PACMAD", "Andropogoneae", "Paniceae", "Setaria", 
                      "Eudicots", "Rosids", "Eurosids","Malvids", "BM", "Brassicales", 
                      "Brassicaceae", "Malvaceae", "Fabids","Fabaceae","Phaseoleae","Solanaceae")

all_species <- unlist(Chloroplastida)

PS_details <- rep("Chloroplastida", nrow(OGs_34spp_count))
for(PS_name in PS_details_names){
  cat(PS_name, "\n")
  tmp_PS <- eval(parse(text = PS_name))
  tmp_PS_others <- setdiff(all_species, unlist(tmp_PS))
  tmp_PS_value <- rep(0, nrow(OGs_34spp_count))
  index2 <- apply(OGs_34spp_count[,tmp_PS_others], 1, function(x){sum(x>0)})==0
  if(is.character(tmp_PS)){
    if(length(tmp_PS)>=2){
      index1 <- apply(OGs_34spp_count[,tmp_PS], 1, function(x){sum(x>0)})>=2
      PS_details[index1&index2] <- PS_name
    }
  } else{
    for(PS_each in 1:length(tmp_PS)){
      if(length(tmp_PS[[PS_each]])==1){
        PS_each_value <- OGs_34spp_count[,tmp_PS[[PS_each]]]>0
      } else {
        PS_each_value <- apply(OGs_34spp_count[,tmp_PS[[PS_each]]], 1, function(x){sum(x>0)})>0
      }
      tmp_PS_value <- tmp_PS_value + PS_each_value
    }
    PS_details[tmp_PS_value>=2&index2] <- PS_name
  }
}

## Get the PS level for each OGs
table(PS_details)

## species-specific OGs
for(tmp_ss in all_species){
  tmp_PS_others <- setdiff(all_species, tmp_ss)
  index2 <- apply(OGs_34spp_count[,tmp_PS_others], 1, function(x){sum(x>0)})==0
  index1 <- OGs_34spp_count[,tmp_ss]!=0
  PS_details[index1&index2] <- "LS_OGs" 
}

table(PS_details)[intersect(c(PS_details_names, "LS_OGs"), unique(PS_details))]

## add the PS level to data frame
OGs_34spp <- cbind(OGs_34spp, "PSnames" = PS_details)
table(OGs_34spp$PSnames)
# Chloroplastida: 6059 Embryophyta: 2864 Stomatophyta: 739 Tracheophyta: 628 Spermatophyta:4006
OGs_34spp_count <- cbind(OGs_34spp_count, "PSnames" = PS_details)
# write.table(OGs_34spp_count, file = "OG_gene_count_and_PS_levels.txt", quote = F, sep = "\t")

for(ii in SPECIES34){
  tmp_genes <- OGs_34spp[OGs_34spp[,ii]!="-", ii]
  tmp_genes <- unlist(strsplit(tmp_genes, ", "))
  cat(ii, length(tmp_genes)==length(unique(tmp_genes)), "\n")
}

## 13 species ####
OGs_13spp <- OGs_34spp[, c(SPECIES13, "PSnames")]
tmp_rm_index <- apply(OGs_13spp[,1:13], 1, function(x){sum(x=="")})<13
table(tmp_rm_index)
OGs_13spp <- OGs_13spp[tmp_rm_index, ]
OGs_13spp_count <- OGs_34spp_count[tmp_rm_index, SPECIES13]

OGs_cross_spp <- OGs_13spp[apply(OGs_13spp_count, 1, function(x){sum(x>0)>1}), ]
for(ii in 1:13){
  tmp_genes <- unlist(strsplit(OGs_cross_spp[,ii], ", "))
  tmp_genes <- tmp_genes[!tmp_genes%in%""]
  cat(colnames(OGs_cross_spp)[ii], length(tmp_genes)==length(unique(tmp_genes)), "\n")
}


## gene attributes data frame ####
gene_attributes_list <- list()

for(species in SPECIES13){
  ## gene names
  tmp_df <- read.table(paste0("associated_data/pep_cds/gene_names/", 
                              species, "_gene_cDNA_pep.txt"), sep = "\t")[,1:2]
  rownames(tmp_df) <- tmp_df[,1]
  colnames(tmp_df) <- c("Gene", "Transcript")
  tmp_exp_val <- round(exp_list$exp[[species]]$MeanExp, 2)
  names(tmp_exp_val) <- exp_list$exp[[species]]$Gene
  tmp_df$Exp <- tmp_exp_val[rownames(tmp_df)]
  tmp_m6A_val <- round(m6A_list$m6A[[species]]$Ratio, 2)
  names(tmp_m6A_val) <- m6A_list$m6A[[species]]$Gene
  tmp_df$m6A <- tmp_m6A_val[rownames(tmp_df)]
  cat(species, nrow(tmp_df), sum(tmp_df[,1]%in%exp_list[[species]]), sum(tmp_df[,1]%in%m6A_list[[species]]), "\n")
  ## Add OGs
  OG_length <- sapply(strsplit(OGs_13spp[, species], ", "), length)
  gene_names <- unlist(strsplit(OGs_13spp[OGs_13spp[, species]!="", species], ", "))
  OG_names <- rep(rownames(OGs_13spp), OG_length)
  PS_names <- rep(OGs_13spp$PSnames, OG_length)
  names(OG_names) <- names(PS_names) <- gene_names
  tmp_df$OGs <- OG_names[rownames(tmp_df)]
  tmp_df$PSs <- PS_names[rownames(tmp_df)]
  LS_one_index <- is.na(tmp_df$PSs)
  tmp_df$PSs[LS_one_index] <- "LS_OC"
  tmp_df$OGs[LS_one_index] <- paste0("LS", 1:sum(LS_one_index))
  ## Add duplication
  tmp_df$Duplication <- "singleton"
  for(dup in c("wgd", "tandem", "proximal", "transposed", "dispersed")){
    dup_tmp_df <- read.table(paste0("associated_data/Duplication_Genes/", 
                                    species, "/pep_DupGen_results/", species, ".", 
                                    dup, ".genes-unique"), sep = "\t", 
                             stringsAsFactors = F, header = T)
    tmp_df[dup_tmp_df[,1], "Duplication"] <- dup 
  }
  ## Add Transcription and Kinase
  tmp_tf <- as.matrix(read.table(paste0("associated_data/TF/", species, 
                                        ".fa_output/tf_classification.txt"), 
                                 stringsAsFactors = F, sep = "\t", row.names = 1))
  tmp_names <- tmp_tf[,1]
  names(tmp_names) <- rownames(tmp_tf)
  tmp_df$Species <- species
  tmp_df$TFs <- tmp_names[rownames(tmp_df)]
  tmp_kinase <- read.table(paste0("associated_data/TF/", species, 
                                  ".fa_output/shiu_classification.txt"), 
                           stringsAsFactors = F, sep = "\t", row.names = 1)
  tmp_names <- tmp_kinase[,1]
  names(tmp_names) <- rownames(tmp_kinase)
  tmp_df$Kinases <- tmp_names[rownames(tmp_df)]
  tmp_df[is.na(tmp_df)] <- ""
  ## final table
  gene_attributes_list[[species]] <- tmp_df
}
gene_attributes_df <- do.call(rbind, gene_attributes_list)
# write.table(gene_attributes_df, "gene_attributes_df.txt", quote = F, sep = "\t", row.names = F)

PS_all_names <-  c("Chloroplastida", PS_details_names, "LS_OGs", "LS_OC")


## The OGs and SGs in 13 species ####
gene_group_13spp <- read.table("associated_data/Results_Oct22_13spp/Orthogroups/Orthogroups.tsv",
                               sep = "\t", header = T, row.names = 1, stringsAsFactors = F)[, SPECIES13]
gene_group_index <- apply(gene_group_13spp[,1:13], 1, function(x){sum(x!="")})

OGs_list_13spp <- apply(gene_group_13spp[gene_group_index >1, ], 2, function(x){
  unlist(strsplit(x, ", "))
})

## save data
save(OGs_13spp, OGs_list_13spp, OGs_34spp, OGs_34spp_count,
     PS_all_names, gene_attributes_df, file = "notebooks/RData/02_detailed_gene_attributes.RData")


## plots ####
comlab_name <- "m6A methylation ratio (%)"

## Figure 2A ####
### The proportion of m6A modification in monocots and dicots
DM_data <- read.table("associated_data/figure_data/Figure6_8_9_10.csv", sep = ",", header = T)
DM_data$ratio1 <- DM_data$ratio1*100
DM_data$ratio2 <- DM_data$ratio2*100
DM_data$species2 <- apply(DM_data, 1, function(x){x[1:2][!x[1:2]%in%x[3]]})

monodi_type <- c("dicots_dicots", "dicots_monocots", "dicots_Bryophyta",
                 "monocots_monocots","monocots_dicots",  "monocots_Bryophyta")
monodi_abbr <- c("DDs", "DMs", "DPs", "MMs", "MDs", "MPs")
names(monodi_abbr) <- monodi_type
DM_data_filter <- DM_data[DM_data$type%in%monodi_type, ]
DM_data_filter$type <- monodi_abbr[DM_data_filter$type]
DM_data_filter$type <- factor(DM_data_filter$type, levels = monodi_abbr)

# pdf("Figure2A.pdf", width = 5, height = 5)
ggplot(DM_data_filter, aes(x=type, y=ratio1, fill=type)) +
  geom_violin(trim=FALSE) + # geom_jitter(width = 0.2) +
  geom_boxplot(width=0.2, fill="#7F7F7F")+
  scale_fill_manual(values = c("#237E7F", "#309388", "#47A082",
                               "#826F74", "#A78889", "#D0B9B9")) +
  theme_bw(base_size = 14)+ 
  theme(axis.text= element_text(colour="black", size=14),
        axis.line = element_line(colour = "black", size=0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # coord_cartesian(ylim = c(20,85)) +
  stat_compare_means(method = 'wilcox.test', method.args = list(alternative = "two.sided"),
                     label.y =85, size=4,
                     comparisons = list( c("DDs", "DMs"), c("DDs", "DPs"), c("DMs", "DPs"),
                                         c("MMs", "MDs"), c("MMs", "MPs"), c("MDs", "MPs"))) +
  xlab("") + ylab(comlab_name)
# dev.off()


## Figure 2B ####
gene_count <- table(gene_attributes_df$Species)
PS_details_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(PS_all_names))
names(PS_details_colors) <- PS_all_names

gene_attributes_df %>% filter(Species%in%SPECIES13) %>% 
  group_by(PSs, Species) %>% 
  summarise("Count"=length(PSs), 
            "GeneRatio"=length(PSs)/gene_count[unique(Species)]*100, 
            "m6ARatio"=sum(m6A!="")/sum(Exp!="")*100, ## 20210729 -- sum(m6A!="")/length(Gene)*100
            "ExpRatio"=sum(Exp!="")/length(Gene)) %>% 
  as.data.frame() -> gene_num_species

PS_plot <- intersect(PS_all_names, gene_attributes_df$PSs)
PS_plot <- PS_plot[!PS_plot%in%c("LS_OGs", "LS_OC")]
PS_trend <- list(PS_plot[c(1:7, 18:24)], # "Brassicaceae"
                 PS_plot[c(1:7, 18:22, 25)], # "Malvaceae"
                 PS_plot[c(1:7, 18:20, 26:28)], # "Phaseoleae"
                 PS_plot[c(1:7, 18, 29)], # "Solanaceae"
                 PS_plot[c(1:7, 8:10, 16:17)], # "Andropogoneae"
                 PS_plot[c(1:7, 8:15)] # "Triticum"
)

count_ratio_plot1 <- count_exp_plot2 <-  list()
for(ii in 1:length(PS_trend)){
  count_ratio_df <- gene_num_species[gene_num_species$PSs%in%PS_trend[[ii]], ]
  tmp_cr_df <- count_ratio_df[1:2,]
  for(tmp_pss in intersect(PS_all_names, count_ratio_df$PSs)){
    tmp_cr_df <- rbind(tmp_cr_df, count_ratio_df[count_ratio_df$PSs==tmp_pss, ])
  }
  tmp_cr_df <- tmp_cr_df[-c(1,2),]
  # write.table(tmp_cr_df[, c(1,2,5)], paste0("Figure2B_", ii, ".txt"), 
  #             quote = F, sep = "\t", col.names = T, row.names = F)
  count_ratio_df$PSs <- factor(count_ratio_df$PSs, levels = intersect(PS_all_names, count_ratio_df$PSs))
  # count_ratio_df$Species <- factor(count_ratio_df$Species, levels = intersect(SPECIES13, unique(count_ratio_df$Species)))
  # tmp_dm_type <- rep("Monocots", nrow(count_ratio_df))
  # tmp_dm_type[count_ratio_df$Species%in%c("ath", "gar", "ghi", "pvu", "gma", "sly")] <- "Dicots"
  # count_ratio_df$DMtype <- factor(tmp_dm_type, levels = c("Dicots", "Monocots"))
  count_ratio_plot1[[ii]] <- count_ratio_df %>% 
    ggplot(aes(x = PSs, y = m6ARatio))+
    geom_boxplot(outlier.color = NA) + 
    # geom_smooth(method="loess",se=FALSE, color="#2F2F2F") +
    theme_bw(base_size = 12) + 
    theme(axis.line = element_line(colour = "black", size = 0.5),
          axis.text.x = element_text(colour = "black", 
                                     angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("") + ylab(comlab_name) +
    scale_fill_manual(values = PS_details_colors)
  ## 
  # exp_genes <- gene_attributes_df[gene_attributes_df$Exp!="", ]
  # exp_genes <- exp_genes[exp_genes$PSs%in%PS_trend[[ii]], ]
  # exp_genes$Exp <- as.numeric(exp_genes$Exp)
  # exp_genes$logExp <- log2(as.numeric(exp_genes$Exp))
  # tmp_exp_genes_PSs <- intersect(PS_all_names, exp_genes$PSs)
  # for(m6A_gene_index in 1:(length(tmp_exp_genes_PSs)-1)){
  #   cat(tmp_exp_genes_PSs[m6A_gene_index], 
  #       t.test(exp_genes[exp_genes$PSs==tmp_exp_genes_PSs[m6A_gene_index], "logExp"],
  #              exp_genes[exp_genes$PSs==tmp_exp_genes_PSs[m6A_gene_index+1], "logExp"], 
  #              alternative = "greater")$p.value,
  #       "\n")
  # }
  # exp_classify <- rep("non-m6A", nrow(exp_genes))
  # exp_classify[exp_genes$m6A!=""] <- "m6A"
  # exp_genes$Type <- exp_classify
  # exp_genes$PSs <- factor(exp_genes$PSs, levels = tmp_exp_genes_PSs)
  # count_exp_plot2[[ii]] <- ggplot(exp_genes, aes(x = PSs, y = logExp))+ #as.numeric(Exp)
  #   geom_boxplot(outlier.size = 0.6) +
  #   theme_bw(base_size = 12) + 
  #   theme(axis.line = element_line(colour = "black", size = 0.5),
  #         axis.text.x = element_text(colour = "black", 
  #                                    angle = 90, vjust = 0.5, hjust = 1), #angle = 30, vjust = 1, hjust = 1),
  #         axis.text.y = element_text(colour = "black"),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   xlab("") + ylab(expression("Log"[2]*"(FPKM)")) +
  #   scale_fill_manual(values = c("#FF9500", "#CCCCCC"))
}

# pdf("Figure15.pdf", width = 20, height = 4)
count_ratio_plot1[[1]] + count_ratio_plot1[[2]] + count_ratio_plot1[[3]] + count_ratio_plot1[[4]] + count_ratio_plot1[[5]] +
  count_ratio_plot1[[6]] + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position='none')
# dev.off()
# count_exp_plot2[[1]] + count_exp_plot2[[2]] + count_exp_plot2[[3]] + count_exp_plot2[[4]] + count_exp_plot2[[5]] + count_exp_plot2[[6]]


## Supplementary Fig. 2 ####
#"associated_data/figure_data/Figure8.csv"
raw_data <- read.table("associated_data/figure_exp_data/SupFig2.csv", 
                       sep = ",", header = T, stringsAsFactors = F)
colnames(raw_data)[6] <- "species"
raw_data <- raw_data[raw_data$species!="ppa", ]
raw_data <- raw_data[raw_data$type!="P.patens", ] # Embryophyta
raw_data$ratio <- raw_data$m6ARatio*100 # ratio
raw_data_level <- raw_data %>% group_by(species, type) %>% 
  summarise("mean"=mean(ratio), 
            "sd"=sd(ratio), 
            "se"=sd(ratio)/sqrt(length(ratio))) %>% 
  as.data.frame()

raw_data_level$species <- firstup(raw_data_level$species)
raw_data_level$species <- factor(raw_data_level$species, levels = firstup(SPECIES13[-13]))
colnames(raw_data_level)[2] <- "Type"
# pdf("Figure8.pdf", width = 9, height = 3)
ggplot(data = raw_data_level, aes(x=species, y=mean, fill=Type))+
  geom_bar(stat='identity', position="dodge") +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, position = position_dodge(0.9)) +
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab(expression("m"^"6"*"A levels (%)")) + 
  scale_fill_manual(values=c("#3F756D", "#B78D67"))
# dev.off()

for(i in unique(raw_data$species)){
  cat(i, 
      t.test(raw_data[raw_data$species%in%i&raw_data$type%in%"dicot", "m6ARatio"], 
             raw_data[raw_data$species%in%i&raw_data$type%in%"monocot", "m6ARatio"])$p.value, "\n")
}

## Supplementary Fig. 3 ####
gene_PSs_count <- dcast(gene_attributes_df[,c("Species", "PSs", "Gene")], Species~PSs)
rownames(gene_PSs_count) <- gene_PSs_count[,1]
gene_PSs_count <- gene_PSs_count[,-1]
gene_PSs_count_out <- t(gene_PSs_count[SPECIES13, intersect(PS_all_names, colnames(gene_PSs_count))])
write.table(gene_PSs_count_out, "gene_PSs_count.txt", quote = F, sep = "\t")

# tmp_gene_attributes_df_tmp <- gene_attributes_df
# tmp_gene_attributes_df_tmp[tmp_gene_attributes_df_tmp$PSs%in%c("LS_OGs", "LS_OC"), "PSs"] <- "SGs"
# tmp_gene_attributes_df_tmp %>% filter(Species%in%SPECIES13) %>% 
#   group_by(PSs, Species) %>% 
#   summarise("Count"=length(PSs), "GeneRatio"=length(PSs)/gene_count[unique(Species)]*100, "m6ARatio"=sum(m6A!="")/length(Gene)*100, "ExpRatio"=sum(Exp!="")/length(Exp), 
#             "m6A_genes"=sum(m6A!=""), "Gene2"=length(Exp)) %>% 
#   as.data.frame() -> tmp_gene_num_species
# write.table(tmp_gene_num_species, "tmp_gene_num_species.txt", quote = F, sep = "\t", row.names = F, col.names = T)
