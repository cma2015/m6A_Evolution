## current wroking directory
setwd("~/a2z/m6A_13spp/")

## Options Settings
options(stringsAsFactors = F)
options(scipen=200)

## loading libraries
library(ggplot2)
library(ggpubr)
library(patchwork)
library(progress)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

## reload data
load("notebooks/RData/01_m6A_peaks_genes.RData")
load("notebooks/RData/02_detailed_gene_attributes.RData")

comlab_name <- "m6A methylation ratio (%)"
system("echo 6 > proc")

## Figure3A ####
### Species divergenece time and m6A divergence
data3a <- read.table("associated_data/figure_data/Figure14_1.csv", 
                     sep = ",", header = T)
data3a$Spearman <- 1-data3a$Spearman
data3a[data3a$Spearman>1, "Spearman"] <- 1
data3a$type <- factor(data3a$type, 
                      levels = c("monocots_monocots", "dicots_dicots", 
                                 "dicots_monocots", "bryophyta_monocots", 
                                 "bryophyta_dicots"))

# pdf("Figure.pdf", width = 7.2, height = 5)
ggplot(data3a, mapping = aes(x=time, y=Spearman, color=type)) +
  geom_point(size=2.8) + 
  scale_color_manual(values = c("#25926F", "#D05E15", "#AC7FB6",
                                "#C73279", "#275BAA")) +
  theme_bw(base_size = 14)+ 
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black",size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Species divergenece time (MYA)") + ylab(comlab_name) + 
  ylim(c(0.2, 1.05)) +
  geom_smooth(method="lm", se=F, color="#A3A2A2") #method="lm", se=T, color="#DA2416"
# dev.off()

## Supplementary Fig. 2 ####
raw_data <- read.table("associated_data/figure_data/Figure8.csv", sep = ",", header = T, stringsAsFactors = F)
raw_data <- raw_data[raw_data$type!="Embryophyta", ]
raw_data$ratio <- raw_data$ratio*100
raw_data_level <- raw_data %>% group_by(species, type) %>% summarise("mean"=mean(ratio), "sd"=sd(ratio), "se"=sd(ratio)/sqrt(length(ratio))) %>% as.data.frame()
# raw_data_level <- melt(raw_data_level)
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
  # scale_fill_manual(values=c("#F74275", "#EDDA66"))
  # dev.off()
  
  for(i in unique(raw_data$species)){
    cat(i, t.test(raw_data[raw_data$species%in%i&raw_data$type%in%"Dicots", "ratio"], 
                  raw_data[raw_data$species%in%i&raw_data$type%in%"Monocots", "ratio"])$p.value, "\n")
  }


## Figure3C plot single-copy OGs in 7 diploid ###
summary_table_list1 <- summary_table_list2 <- list()
output_list1 <- output_list2 <- list()

kaks_list <- list(c("ath", "sly"), c("ath", "sbi"), c("ath", "osa"),
                  c("ath", "gar"), c("ath", "pvu"), c("ath", "ata"), 
                  c("sly", "sbi"), c("sly", "osa"), c("sbi", "osa"),
                  c("sly", "gar"), c("sly", "pvu"), c("sly", "ata"),
                  c("sbi", "gar"), c("sbi", "pvu"), c("sbi", "ata"),
                  c("osa", "gar"), c("osa", "pvu"), c("osa", "ata"),
                  c("gar", "pvu"), c("gar", "ata"), c("pvu", "ata"))

PS_index_names <- c("Chloroplastida", "Embryophyta", "Spermatophyta", "Angiosperm") #c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae")
limit_PSs <- c("Chloroplastida", "Embryophyta", "Spermatophyta", "Angiosperm")

evolutionary_rate_cor <- m6A_divergence_cor <- vector()

species_type <- vector()
for(tmp_spe in kaks_list){
  if(tmp_spe[1]%in% c("ath", "gar", "pvu", "sly")){
    species_type1 <- "Dicots"
  }else{
    species_type1 <- "Monocots"
  }
  if(tmp_spe[2] %in% c("ath", "gar", "pvu", "sly")){
    species_type2 <- "Dicots"
  }else{
    species_type2 <- "Monocots"
  }
  tmp_type <- sort(c(species_type1, species_type2))
  species_type <- c(species_type, paste0(tmp_type[1], " vs ", tmp_type[2]))
}

pairs_number <- data.frame("species_pair" = rep(sapply(kaks_list, function(x){paste0(x, collapse = "-")}),
                                                each=14),
                           "PS_names" = rep(rep(PS_index_names, each=2), 21), 
                           "Type" = rep(c("Count", "Exp"), 147), "count_value" = 0)

## main
pairs_number_count <- vector()
pb <- progress_bar$new(total = length(kaks_list), clear = FALSE)
identity_score_list <- list()
for(two_spp in kaks_list){
  pb$tick()
  two_spp_name <- paste0(two_spp, collapse = "-")
  div_val_list <- kaks_val_list <- list()
  m6a_levels <- exp_levels <- list()
  two_spp_orthologues <- read.table(paste0("associated_data/Results_Jan16_34spp/Orthologues/Orthologues_",
                                           two_spp[1], "/", two_spp[1], "__v__", two_spp[2],".tsv"), 
                                    header = T, sep = "\t")
  # sp1_cdna <- read.fasta(paste0("associated_data/pep_cds/cDNA/", two_spp[1], ".fa"),
  #                        seqtype = "DNA", as.string = T, forceDNAtolower =F)
  # sp2_cdna <- read.fasta(paste0("associated_data/pep_cds/cDNA/", two_spp[2], ".fa"), 
  #                        seqtype = "DNA", as.string = T, forceDNAtolower =F)
  one_divergence <- one_kaks <- list()
  identity_score_list[[two_spp_name]] <- list()
  for(PS_index_name in PS_index_names) {
    tmp_ele_name <- paste0(two_spp_name, "_", PS_index_name)
    PS_two_spp <- OGs_13spp[OGs_13spp$PSnames==PS_index_name, two_spp]
    speciesone <- unlist(strsplit(PS_two_spp[,1], split = ", "))
    speciestwo <- unlist(strsplit(PS_two_spp[,2], split = ", "))
    speciesot_index <- apply(two_spp_orthologues, 1, function(x){
      tmp_index1 <- all(unlist(strsplit(x[2], split = ", "))%in%speciesone)
      tmp_index2 <- all(unlist(strsplit(x[3], split = ", "))%in%speciestwo)
      return(tmp_index1&tmp_index2)
    })
    two_spp_ortho_filter <- two_spp_orthologues[speciesot_index, 2:3]
    PS_two_spp_type <- apply(two_spp_ortho_filter, 1, function(x){
      tmp1 <- unlist(strsplit(x[1], split = ", "))
      tmp2 <- unlist(strsplit(x[2], split = ", "))
      if(length(tmp1)>1){
        a <- "m"
      }else{
        a <- "1"
      }
      if(length(tmp2)>1){
        b <- "m"
      }else{
        b <- "1"
      }
      return(paste0(c(a,b), collapse = ":"))
    })
    # table(PS_two_spp_type)
    two_spp_ortho_filter <- two_spp_ortho_filter[PS_two_spp_type=="1:1", ]
    pairs_number_count <- c(pairs_number_count, nrow(two_spp_ortho_filter))
    PS_two_spp_exp <- apply(two_spp_ortho_filter, 1, function(x){
      if(x[1]%in%exp_list[[two_spp[1]]]){
        a <- "1"
      }else{
        a <- "0"
      }
      if(x[2]%in%exp_list[[two_spp[2]]]){
        b <- "1"
      }else{
        b <- "0"
      }
      return(paste0(c(a,b), collapse = ":"))
    })
    # table(PS_two_spp_exp)
    two_spp_ortho_filter <- two_spp_ortho_filter[PS_two_spp_exp=="1:1", ]
    pairs_number_count <- c(pairs_number_count, nrow(two_spp_ortho_filter))
    PS_two_spp_m6A <- apply(two_spp_ortho_filter, 1, function(x){
      if(x[1]%in%m6A_list[[two_spp[1]]]){
        a <- "1"
      }else{
        a <- "0"
      }
      if(x[2]%in%m6A_list[[two_spp[2]]]){
        b <- "1"
      }else{
        b <- "0"
      }
      return(paste0(c(a,b), collapse = ":"))
    })
    # summary_table_list1[[tmp_ele_name]] <- table(PS_two_spp_type)
    # summary_table_list2[[tmp_ele_name]] <- table(PS_two_spp_m6A)
    # tmp_PS_two_spp_m6A <- PS_two_spp_m6A
    # tmp_PS_two_spp_m6A[tmp_PS_two_spp_m6A=="0:1"] <- "1:0"
    # uniden_part <- table(tmp_PS_two_spp_m6A)[c("1:0")]/sum(table(tmp_PS_two_spp_m6A))*100
    # m6A_divergence_cor <- c(m6A_divergence_cor, uniden_part)
    # div_val_list[[PS_index_name]] <- c(as.numeric(table(tmp_PS_two_spp_m6A)[c("1:0")]), 
    #                                    sum(table(tmp_PS_two_spp_m6A))-as.numeric(table(tmp_PS_two_spp_m6A)[c("1:0")]))
    # iden_part <- 100-uniden_part
    # one_divergence[[tmp_ele_name]] <- data.frame("Species"=two_spp_name, "PSs"=PS_index_name,
    #                                              "Type"=c("Identical pattern", "Diverged pattern"),
    #                                              "value"=c(iden_part, uniden_part))
    # m6A_divergence_index <- PS_two_spp_m6A%in%c("1:0", "0:1")
    # if(sum(m6A_divergence_index) >0 ){
    #   m6A_divergence_two_spp <- two_spp_ortho_filter[m6A_divergence_index, ]
    #   identity_score <- apply(m6A_divergence_two_spp, 1, function(x){
    #     write.fasta(sequences = sp1_cdna[[x[1]]], names = x[1], file.out = paste0("tmp_seq/", x[1],".fa"))
    #     write.fasta(sequences = sp2_cdna[[x[2]]], names = x[2], file.out = paste0("tmp_seq/", x[2],".fa"))
    #     system(command = paste0("~/miniconda3/envs/metaPlants/bin/_needle -gapopen=10 -gapextend=0.5 -asequence ", 
    #                             paste0("tmp_seq/", x[1],".fa"), " -bsequence ", 
    #                             paste0("tmp_seq/", x[2],".fa"), " ", 
    #                             paste0("tmp_seq/ali_", x[1],"_", x[2],".fa")))
    #     system(command = paste0("grep 'Identity:' ", paste0("tmp_seq/ali_", x[1],"_", x[2],".fa | grep -oP '[0-9]+\\.[0-9]+%'")),intern = T)
    #   })
    #   identity_score_list[[two_spp_name]][[PS_index_name]] <- gsub("%", "", identity_score)
    #   ## ka ks
    #   if("run" == "nonrun"){
    #     system(paste0("mkdir -p associated_data/KaKs_m6A/",tmp_ele_name))
    #     write.table(m6A_divergence_two_spp, paste0("associated_data/KaKs_m6A/", tmp_ele_name,
    #                                                "/gene_pairs.txt"), quote = F, sep = "\t",
    #                 col.names = F, row.names = F)
    #     write(paste0("cat associated_data/pep_cds/cds/", two_spp[1],
    #                  ".fa associated_data/pep_cds/cds/", two_spp[2],
    #                  ".fa >associated_data/KaKs_m6A/", tmp_ele_name, "/spp_merge_cds.fa"),
    #           file = "KaKs_m6A_divergence.sh", append = T)
    # 
    #     write(paste0("cat associated_data/pep_cds/pep/", two_spp[1],
    #                  ".fa associated_data/pep_cds/pep/", two_spp[2],
    #                  ".fa >associated_data/KaKs_m6A/", tmp_ele_name, "/spp_merge_pep.fa"),
    #           file = "KaKs_m6A_divergence.sh", append = T)
    # 
    #     write(paste0("tools/ParaAT2.0/ParaAT.pl -h associated_data/KaKs_m6A/", tmp_ele_name,
    #                  "/gene_pairs.txt -a associated_data/KaKs_m6A/",
    #                  tmp_ele_name, "/spp_merge_pep.fa ",
    #                  " -n associated_data/KaKs_m6A/", tmp_ele_name,
    #                  "/spp_merge_cds.fa -p proc -o associated_data/KaKs_m6A/",
    #                  tmp_ele_name, "/output -f axt"),
    #           file = "KaKs_m6A_divergence.sh", append = T)
    # 
    #     write(paste0("ls associated_data/KaKs_m6A/", tmp_ele_name,
    #                  "/output | grep '.axt' | parallel -I% tools/KaKs_Calculator ",
    #                  "-i associated_data/KaKs_m6A/", tmp_ele_name,
    #                  "/output/%  -o associated_data/KaKs_m6A/", tmp_ele_name,
    #                  "/output/%.kaks -m MA"), file = "KaKs_m6A_divergence.sh", append = T)
    # 
    #     write(paste0("cat associated_data/KaKs_m6A/", tmp_ele_name,
    #                  "/output/*.kaks | awk '{if(!(NR>1&&$1~/Sequence/)){print $0}}'",
    #                  " >associated_data/KaKs_m6A/", tmp_ele_name,
    #                  "/kaks_", tmp_ele_name, "_out.txt"),
    #           file = "KaKs_m6A_divergence.sh", append = T)
    #   } else {
    #     tmp_kaks <- read.table(paste0("associated_data/KaKs_m6A/", tmp_ele_name, "/kaks_", tmp_ele_name, "_out.txt"),
    #                            sep = "\t", header = T, stringsAsFactors = F)
    #     one_kaks[[tmp_ele_name]] <-  data.frame("Species"=two_spp_name, "PSs"=PS_index_name,
    #                                             "value"=tmp_kaks$Ka.Ks)
    #     kaks_val_list[[PS_index_name]] <- tmp_kaks$Ka.Ks
    #     evolutionary_rate_cor <- c(evolutionary_rate_cor, median(one_kaks[[tmp_ele_name]]$value))
    #   }
    }
  }
  ## statistics
  test_chi <- function(inputdata, var1, var2){
    chisq.test(matrix(c(inputdata[[var1]][1], inputdata[[var1]][2], inputdata[[var2]][1], inputdata[[var2]][2]),
                      nrow = 2, byrow = T))
  }
  test_wilcox <- function(inputdata, var1, var2){
    wilcox.test(inputdata[[var1]], inputdata[[var1]], alternative = "greater")
  }
  
  # test_chi(inputdata = div_val_list, var1 = limit_PSs[1], var2 = limit_PSs[2])
  # test_chi(inputdata = div_val_list, var1 = limit_PSs[2], var2 = limit_PSs[3])
  # test_chi(inputdata = div_val_list, var1 = limit_PSs[3], var2 = limit_PSs[4])
  # 
  # test_wilcox(inputdata = kaks_val_list, var1 = limit_PSs[1], var2 = limit_PSs[2])
  # test_wilcox(inputdata = kaks_val_list, var1 = limit_PSs[2], var2 = limit_PSs[3])
  # test_wilcox(inputdata = kaks_val_list, var1 = limit_PSs[3], var2 = limit_PSs[4])
  ## percentage barplot
  one_divergence_df <- do.call(rbind, one_divergence)
  one_divergence_df <- one_divergence_df[one_divergence_df$Type=="Diverged pattern", ]
  one_divergence_df <- one_divergence_df[one_divergence_df$PSs%in%limit_PSs, ]
  one_divergence_df$PSs <- factor(one_divergence_df$PSs, levels = limit_PSs)
  div_val <- max(one_divergence_df[, "value"])
  output_list1[[two_spp_name]] <- ggplot(data = one_divergence_df, aes(x=PSs, y=value, fill=PSs))+
    geom_bar(stat = "identity", width=0.6) +
    theme_bw(base_size = 14) +
    theme(axis.line = element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    xlab("") + ylab("Percentage of m6A divergence") + 
    coord_cartesian(ylim = c(0, div_val+5)) +
    scale_fill_manual(values = c("#036F6F", "#B06A6A", 
                                  "#9B70A7", "#C88F17"))
  ### Ka/Ks boxplot
  one_kaks_df <- do.call(rbind, one_kaks)
  limval <- boxplot(one_kaks_df[one_kaks_df$PSs=="Angiosperm", ]$value, plot=F)[[1]][5,1] + 0.05
  one_kaks_df <- one_kaks_df[one_kaks_df$PSs%in%limit_PSs, ]
  one_kaks_df$PSs <- factor(one_kaks_df$PSs, levels = limit_PSs)
  output_list2[[two_spp_name]] <- ggplot(data = one_kaks_df, aes(x=PSs, y=value, fill=PSs))+
    geom_boxplot(outlier.color = NA) +
    theme_bw(base_size = 14) +
    theme(axis.ticks.x = element_blank(), 
      axis.line = element_line(colour = "black", size = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    xlab("") + ylab("Ka/Ks") +
    coord_cartesian(ylim = c(0, limval+0.05)) +
    scale_fill_manual(values = c("#036F6F", "#B06A6A", 
                                 "#9B70A7", "#C88F17"))
}

## Figure 3C ath gar pvu sly sbi ata osa ####
# pdf("Figure3C.pdf", width = 20, height = 15)
(plot_spacer()             | output_list1[["ath-gar"]] | output_list1[["ath-pvu"]] | output_list1[["ath-sly"]] | output_list1[["ath-sbi"]] | output_list1[["ath-ata"]] | output_list1[["ath-osa"]]) /
  (output_list2[["ath-gar"]] | plot_spacer()             | output_list1[["gar-pvu"]] | output_list1[["sly-gar"]] | output_list1[["sbi-gar"]] | output_list1[["gar-ata"]] | output_list1[["osa-gar"]]) /
  (output_list2[["ath-pvu"]] | output_list2[["gar-pvu"]] | plot_spacer()             | output_list1[["sly-pvu"]] | output_list1[["sbi-pvu"]] | output_list1[["pvu-ata"]] | output_list1[["osa-pvu"]]) /
  (output_list2[["ath-sly"]] | output_list2[["sly-gar"]] | output_list2[["sly-pvu"]] | plot_spacer()             | output_list1[["sly-sbi"]] | output_list1[["sly-ata"]] | output_list1[["sly-osa"]]) /
  (output_list2[["ath-sbi"]] | output_list2[["sbi-gar"]] | output_list2[["sbi-pvu"]] | output_list2[["sly-sbi"]] | plot_spacer()             | output_list1[["sbi-ata"]] | output_list1[["sbi-osa"]]) /
  (output_list2[["ath-ata"]] | output_list2[["gar-ata"]] | output_list2[["pvu-ata"]] | output_list2[["sly-ata"]] | output_list2[["sbi-ata"]] | plot_spacer()             | output_list1[["osa-ata"]]) /
  (output_list2[["ath-osa"]] | output_list2[["osa-gar"]] | output_list2[["osa-pvu"]] | output_list2[["sly-osa"]] | output_list2[["sbi-osa"]] | output_list2[["osa-ata"]] | plot_spacer()) +
  plot_layout(guides = 'collect')
# dev.off()

## Figure 3D ####
evo_div_df <- data.frame(
  "m6A_divergence" = m6A_divergence_cor,
  "evolutionary_rate" = evolutionary_rate_cor,
  "type" = species_type,
  "Type" = rep(PS_index_names, 21)
)

evo_div_df <- evo_div_df[evo_div_df$Type%in%limit_PSs,]
cor.test(evo_div_df$m6A_divergence, evo_div_df$evolutionary_rate)

fit_stat2 <- summary(lm(as.numeric(evo_div_df[["m6A_divergence"]])~as.numeric(evo_div_df[["evolutionary_rate"]])))
cat(fit_stat2$adj.r.squared, fit_stat2$coefficients[2,4], "\n")

evo_div_df$type <- factor(evo_div_df$type, levels = c("Monocots vs Monocots", "Dicots vs Dicots", "Dicots vs Monocots"))
evo_div_df$Type <- factor(evo_div_df$Type, levels = limit_PSs)

# pdf("Figure3D.pdf", width = 6, height = 4, useDingbats = F)
ggplot(data = evo_div_df, mapping = aes(x=as.numeric(m6A_divergence), 
                                        y=as.numeric(evolutionary_rate),
                                        color=Type)) +
  geom_point(size=2) +
  scale_color_manual(values = c("#036F6F", "#B06A6A", "#9B70A7", "#C88F17")) +
  labs(x="m6A_divergence", y=expression("evolutionary_rate")) +
  theme_bw(base_size = 12)+ 
  theme(axis.line = element_line(colour = "black",size=0.5),
        axis.text = element_text(colour="black"),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  # expand_limits(x=0, y=0) +
  labs(x="Percentage of m6A divergence (%)", y = "Ka/Ks") +
  geom_smooth(method="lm", se=F, color="#A3A2A2", size=1)  #+
# dev.off()


## subgenome m6A levels ####
two_spp <- c("gar", "ghi")
m6a_levels <- exp_levels <- list()
two_spp_name <- paste0(two_spp, collapse = "-")
two_spp_orthologues <- gar_ghi[gar_ghi[,2]!="-",1:2]
for(PS_index_name in PS_index_names) {
  tmp_ele_name <- paste0(two_spp_name, "_", PS_index_name)
  PS_two_spp <- OGs_13spp[OGs_13spp$PSnames==PS_index_name, two_spp]
  speciesone <- unlist(strsplit(PS_two_spp[,1], split = ", "))
  speciestwo <- unlist(strsplit(PS_two_spp[,2], split = ", "))
  speciesot_index <- apply(two_spp_orthologues, 1, function(x){
    tmp_index1 <- all(unlist(strsplit(x[1], split = ", "))%in%speciesone)
    tmp_index2 <- all(unlist(strsplit(x[2], split = ", "))%in%speciestwo)
    return(tmp_index1&tmp_index2)
  })
  two_spp_ortho_filter <- two_spp_orthologues[speciesot_index, 1:2]
  PS_two_spp_type <- apply(two_spp_ortho_filter, 1, function(x){
    tmp1 <- unlist(strsplit(x[1], split = ", "))
    tmp2 <- unlist(strsplit(x[2], split = ", "))
    if(length(tmp1)>1){
      a <- "m"
    }else{
      a <- "1"
    }
    if(length(tmp2)>1){
      b <- "m"
    }else{
      b <- "1"
    }
    return(paste0(c(a,b), collapse = ":"))
  })
  # table(PS_two_spp_type)
  two_spp_ortho_filter <- two_spp_ortho_filter[PS_two_spp_type=="1:1", ]
  pairs_number_count <- c(pairs_number_count, nrow(two_spp_ortho_filter))
  PS_two_spp_exp <- apply(two_spp_ortho_filter, 1, function(x){
    if(x[1]%in%exp_list[[two_spp[1]]]){
      a <- "1"
    }else{
      a <- "0"
    }
    if(x[2]%in%exp_list[[two_spp[2]]]){
      b <- "1"
    }else{
      b <- "0"
    }
    return(paste0(c(a,b), collapse = ":"))
  })
  # table(PS_two_spp_exp)
  two_spp_ortho_filter <- two_spp_ortho_filter[PS_two_spp_exp=="1:1", ]
  pairs_number_count <- c(pairs_number_count, nrow(two_spp_ortho_filter))
  PS_two_spp_m6A <- apply(two_spp_ortho_filter, 1, function(x){
    if(x[1]%in%m6A_list[[two_spp[1]]]){
      a <- "1"
    }else{
      a <- "0"
    }
    if(x[2]%in%m6A_list[[two_spp[2]]]){
      b <- "1"
    }else{
      b <- "0"
    }
    return(paste0(c(a,b), collapse = ":"))
  })
  ## cor
  IM_df <- two_spp_ortho_filter[PS_two_spp_m6A=="1:1", ]
  m6a_levels <- c(m6a_levels, cor(m6A_list$m6A[[two_spp[1]]][IM_df[,1], "Ratio"], 
                                  m6A_list$m6A[[two_spp[2]]][IM_df[,2], "Ratio"], 
                                  method = "spearman"))
  exp_levels <- c(exp_levels, cor(exp_list$exp[[two_spp[1]]][IM_df[,1], "MeanExp"], 
                                  exp_list$exp[[two_spp[2]]][IM_df[,2], "MeanExp"], 
                                  method = "spearman"))
}
