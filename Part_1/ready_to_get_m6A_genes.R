# peak and expression values
# Time: 2021-06-03
# Raw path: Cal-Server

setwd("~/a2z/m6A_13spp/Part_1/")
options(stringsAsFactors = F)


# library for manipulating data frame ####
library(dplyr)
library(parallel)


# data frames for expression and m6A information ####
species12.names <- c("ath", "gar", "ghi", "pvu", "gma", "sbi", "zma", "ata", 
                     "tdi", "tae", "osa", "ppa")
species13.names <- c(species12.names[1:5], "sly", species12.names[6:12])
data.path <- "workflow/data/"
exp.list <- m6A.list <- list()
exp.list$exp <- m6A.list$m6A <- m6A.list$peaks <- list()


# importing m6A peaks and genes from 12 species ####
for (sp.tmp in species12.names){ # sp - species
  # expression matrix    
  exp.df <- read.table(paste0(data.path, sp.tmp, "/04Gene/", sp.tmp, "_gene_feature.txt"), 
                       sep = "\t", stringsAsFactors = F)
  colnames(exp.df) <- c("Chr", "Start", "End", "Gene", "Type", "Strand", "MeanExp")
  exp.df <- exp.df %>% filter(MeanExp >= 1)
  rownames(exp.df) <- exp.df$Gene
  # m6A modification matrix
  mod.df <- read.table(paste0(data.path, sp.tmp, "/05Ratio/", sp.tmp, "_gene_ratio.txt"), 
                       sep = "\t", stringsAsFactors = F)
  colnames(mod.df) <- c("Chr", "Start", "End", "Gene", "Type", "Strand", 
                        "IP1", "IP2", "IP3", "Input1", "Input2", "Input3")
  # peak information
  peak.df <- mod.df %>% filter(Gene%in%exp.df$Gene) %>%
    mutate(IPval = IP1+IP2+IP3+1, Inputval = Input1+Input2+Input3+1) %>%
    mutate(Ratio = IPval/Inputval) %>% 
    filter(Ratio >= 2 ) # IP/Input > 2
  # print(summary(peak.df$End-peak.df$Start+1))
  # merge peak to gene
  mod.df <- peak.df %>% 
    mutate(IPval = IP1+IP2+IP3, Inputval = Input1+Input2+Input3) %>%
    group_by(Gene) %>%
    summarise(IPsum = sum(IPval)+1,
              Inputsum = sum(Inputval)+1) %>%
    mutate(Ratio = IPsum/Inputsum) %>%
    as.data.frame()
  rownames(mod.df) <- mod.df$Gene
  # save data
  out.peak.df <- peak.df[, 1:6]
  out.peak.df[,5] <- "."
  out.peak.df[,2] <- out.peak.df[,2]-1
  write.table(out.peak.df, paste0(data.path, sp.tmp, "/05Ratio/R_output_peaks.bed"), 
              quote = F, sep = "\t",
              row.names = F, col.names = F)
  # distinguish with public data
  if(sp.tmp%in%c("ath", "zma", "osa")){sp.tmp <- paste0(sp.tmp, "_own")}
  # output data
  exp.list[[sp.tmp]] <- unique(exp.df$Gene)
  m6A.list[[sp.tmp]] <- unique(peak.df$Gene)
  print(length(unique(peak.df$Gene))==nrow(mod.df))
  exp.list[['exp']][[sp.tmp]] <- exp.df
  m6A.list[['m6A']][[sp.tmp]] <- mod.df
  m6A.list[['peaks']][[sp.tmp]] <- peak.df
}
rm(exp.df, mod.df, peak.df, out.peak.df)

# public data with three replication ####
species.others <- c("maizePpMc", "slyGbQgz1", "slyGbQgz2")
for (sp.tmp in species.others){
  # expression matrix
  exp.df <- read.table(paste0(data.path, sp.tmp, "/04Gene/", sp.tmp, "_gene_feature.txt"), 
                       sep = "\t", stringsAsFactors = F)
  colnames(exp.df) <- c("Chr", "Start", "End", "Gene", "Type", "Strand", "MeanExp")
  exp.df <- exp.df %>% 
    filter(Type == "PCG") %>% 
    filter(!Chr %in% c("Mt", "Pt")) %>% 
    filter(MeanExp >= 1)
  rownames(exp.df) <- exp.df$Gene
  # m6A modification matrix
  mod.df <- read.table(paste0(data.path, sp.tmp, "/05Ratio/", sp.tmp, "_gene_ratio.txt"), 
                       sep = "\t", stringsAsFactors = F)
  colnames(mod.df) <- c("Chr", "Start", "End", "Gene", "Type", "Strand", 
                        "IP1", "IP2", "IP3", "Input1", "Input2", "Input3")
  # peak information
  peak.df <- mod.df %>% filter(Gene%in%exp.df$Gene) %>%
    filter(Type == "PCG") %>% 
    filter(!Chr %in% c("Mt", "Pt")) %>% 
    mutate(IPval = IP1+IP2+IP3+1, Inputval = Input1+Input2+Input3+1) %>%
    mutate(Ratio = IPval/Inputval) %>% 
    filter(Ratio >= 2)
  print(summary(peak.df$End-peak.df$Start+1))
  # merge peak to gene
  mod.df <- peak.df %>% 
    mutate(IPval = IP1+IP2+IP3, Inputval = Input1+Input2+Input3) %>%
    group_by(Gene) %>%
    summarise(IPsum = sum(IPval)+1,
              Inputsum = sum(Inputval)+1) %>%
    mutate(Ratio = IPsum/Inputsum) %>%
    as.data.frame()
  rownames(mod.df) <- mod.df$Gene
  # save data
  out.peak.df <- peak.df[, 1:6]
  out.peak.df[,5] <- "."
  out.peak.df[,2] <- out.peak.df[,2]-1
  write.table(out.peak.df, paste0(data.path, sp.tmp, "/05Ratio/R_output_peaks.bed"), 
              quote = F, sep = "\t",
              row.names = F, col.names = F)
  # output data
  exp.list[[sp.tmp]] <- unique(exp.df$Gene)
  m6A.list[[sp.tmp]] <- unique(peak.df$Gene)
  exp.list[['exp']][[sp.tmp]] <- exp.df
  m6A.list[['m6A']][[sp.tmp]] <- mod.df
  m6A.list[['peaks']][[sp.tmp]] <- peak.df
  # cat(sp.tmp, nrow(mod.df), nrow(unique(mod.df[,1:3])), nrow(peak.df), nrow(unique(peak.df[,1:3])), length(unique(peak.df$Gene)), "\n")
}
rm(exp.df, mod.df, peak.df, out.peak.df)


# public data with two replication ####
species.others <- c("maizePgLyj1", "maizePgLyj2", "maizePgLyj3")
for (sp.tmp in species.others){
  # expression matrix
  exp.df <- read.table(paste0(data.path, sp.tmp, "/04Gene/", sp.tmp, "_gene_feature.txt"), 
                       sep = "\t", stringsAsFactors = F)
  colnames(exp.df) <- c("Chr", "Start", "End", "Gene", 
                        "Type", "Strand", "MeanExp")
  exp.df <- exp.df %>% filter(MeanExp >= 1)
  rownames(exp.df) <- exp.df$Gene
  # m6A modification matrix
  mod.df <- read.table(paste0(data.path, sp.tmp, "/05Ratio/", sp.tmp, "_gene_ratio.txt"), sep = "\t", stringsAsFactors = F)
  colnames(mod.df) <- c("Chr", "Start", "End", "Gene", "Type", "Strand", 
                        "IP1", "IP2", "Input1", "Input2")
  # peak information
  peak.df <- mod.df %>% filter(Gene%in%exp.df$Gene) %>%
    mutate(IPval = IP1+IP2+1, Inputval = Input1+Input2+1) %>%
    mutate(Ratio = IPval/Inputval) %>% 
    filter(Ratio >= 2)
  print(summary(peak.df$End-peak.df$Start+1))
  # merge peak to gene
  mod.df <- peak.df %>% 
    mutate(IPval = IP1+IP2, Inputval = Input1+Input2) %>%
    group_by(Gene) %>%
    summarise(IPsum = sum(IPval)+1,
              Inputsum = sum(Inputval)+1) %>%
    mutate(Ratio = IPsum/Inputsum) %>%
    as.data.frame()
  rownames(mod.df) <- mod.df$Gene
  # save data
  out.peak.df <- peak.df[, 1:6]
  out.peak.df[,5] <- "."
  out.peak.df[,2] <- out.peak.df[,2]-1
  write.table(out.peak.df, paste0(data.path, sp.tmp, "/05Ratio/R_output_peaks.bed"), 
              quote = F, sep = "\t",
              row.names = F, col.names = F)
  # output data
  exp.list[[sp.tmp]] <- unique(exp.df$Gene)
  m6A.list[[sp.tmp]] <- unique(peak.df$Gene)
  exp.list[['exp']][[sp.tmp]] <- exp.df
  m6A.list[['m6A']][[sp.tmp]] <- mod.df
  m6A.list[['peaks']][[sp.tmp]] <- peak.df
  # cat(sp.tmp, nrow(mod.df), nrow(unique(mod.df[,1:3])), nrow(peak.df), nrow(unique(peak.df[,1:3])), length(unique(peak.df$Gene)),"\n")
}
rm(exp.df, mod.df, peak.df, out.peak.df)


species.others <- c("maizePpHy", "riceNcbJgf1", "riceNcbJgf2", "athCrBgd1", "athCrBgd2")
for (sp.tmp in species.others){
  # expression matrix
  exp.df <- read.table(paste0(data.path, sp.tmp, "/04Gene/", sp.tmp, "_gene_feature.txt"), sep = "\t", stringsAsFactors = F)
  colnames(exp.df) <- c("Chr", "Start", "End", "Gene", "Type", "Strand", "MeanExp")
  exp.df <- exp.df %>% 
    filter(Type == "PCG") %>% 
    filter(!Chr %in% c("Mt", "Pt")) %>% 
    filter(MeanExp >= 1)
  rownames(exp.df) <- exp.df$Gene
  # m6A modification matrix
  mod.df <- read.table(paste0(data.path, 
                              sp.tmp, "/05Ratio/", sp.tmp, "_gene_ratio.txt"), 
                       sep = "\t", stringsAsFactors = F)
  colnames(mod.df) <- c("Chr", "Start", "End", "Gene", "Type", "Strand", 
                        "IP1", "IP2", "Input1", "Input2")
  # peak information
  peak.df <- mod.df %>% filter(Gene%in%exp.df$Gene) %>%
    filter(Type == "PCG") %>% 
    filter(!Chr %in% c("Mt", "Pt")) %>% 
    mutate(IPval = IP1+IP2+1, Inputval = Input1+Input2+1) %>%
    mutate(Ratio = IPval/Inputval) %>% 
    filter(Ratio >= 2)
  print(summary(peak.df$End-peak.df$Start+1))
  # merge peak to gene
  mod.df <- peak.df %>% 
    mutate(IPval = IP1+IP2, Inputval = Input1+Input2) %>%
    group_by(Gene) %>%
    summarise(IPsum = sum(IPval)+1,
              Inputsum = sum(Inputval)+1) %>%
    mutate(Ratio = IPsum/Inputsum) %>%
    as.data.frame()
  rownames(mod.df) <- mod.df$Gene
  # save data
  out.peak.df <- peak.df[, 1:6]
  out.peak.df[,5] <- "."
  out.peak.df[,2] <- out.peak.df[,2]-1
  write.table(out.peak.df, paste0(data.path, sp.tmp, "/05Ratio/R_output_peaks.bed"), 
              quote = F, sep = "\t",
              row.names = F, col.names = F)
  # output data
  exp.list[[sp.tmp]] <- unique(exp.df$Gene)
  m6A.list[[sp.tmp]] <- unique(peak.df$Gene)
  exp.list[['exp']][[sp.tmp]] <- exp.df
  m6A.list[['m6A']][[sp.tmp]] <- mod.df
  m6A.list[['peaks']][[sp.tmp]] <- peak.df
  # cat(sp.tmp, nrow(mod.df), nrow(unique(mod.df[,1:3])), nrow(peak.df), nrow(unique(peak.df[,1:3])), length(unique(peak.df$Gene)),"\n")
}
rm(exp.df, mod.df, peak.df, out.peak.df)


# merge maize rice tomato ####
exp.list$zma <- sort(unique(c(exp.list[["zma_own"]], exp.list[['maizePgLyj1']], 
                              exp.list[["maizePgLyj2"]], exp.list[["maizePgLyj3"]], 
                              exp.list[["maizePpHy"]], exp.list[["maizePpMc"]])))
m6A.list$zma <- sort(unique(c(m6A.list[["zma_own"]], m6A.list[["maizePgLyj1"]], 
                              m6A.list[["maizePgLyj2"]], m6A.list[["maizePgLyj3"]], 
                              m6A.list[["maizePpHy"]], m6A.list[["maizePpMc"]])))

exp.list$osa <- sort(unique(c(exp.list[["osa_own"]], exp.list[["riceNcbJgf1"]], 
                              exp.list[["riceNcbJgf2"]])))
m6A.list$osa <- sort(unique(c(m6A.list[["osa_own"]], m6A.list[["riceNcbJgf1"]], 
                              m6A.list[["riceNcbJgf2"]])))

exp.list$sly <- sort(unique(c(exp.list[["slyGbQgz1"]], exp.list[["slyGbQgz2"]])))
m6A.list$sly <- sort(unique(c(m6A.list[["slyGbQgz1"]], m6A.list[["slyGbQgz2"]])))

exp.list$ath <- sort(unique(c(exp.list[["ath_own"]], exp.list[["athCrBgd1"]], 
                              exp.list[["athCrBgd2"]])))
m6A.list$ath <- sort(unique(c(m6A.list[["ath_own"]], m6A.list[["athCrBgd1"]], 
                              m6A.list[["athCrBgd2"]])))

tmp.df <- rbind(exp.list$exp$zma_own, exp.list$exp$maizePgLyj1, 
                exp.list$exp$maizePgLyj2, exp.list$exp$maizePgLyj3, 
                exp.list$exp$maizePpHy, exp.list$exp$maizePpMc)
exp.list$exp$zma <- tmp.df %>% group_by(Gene) %>% 
  filter(MeanExp == max(MeanExp)) %>% as.data.frame()
rownames(exp.list$exp$zma) <- exp.list$exp$zma$Gene
tmp.df <- rbind(m6A.list$m6A$zma_own, m6A.list$m6A$maizePgLyj1, 
                m6A.list$m6A$maizePgLyj2, m6A.list$m6A$maizePgLyj3, 
                m6A.list$m6A$maizePpHy, m6A.list$m6A$maizePpMc)
m6A.list$m6A$zma <- tmp.df %>% group_by(Gene) %>% 
  filter(Ratio == max(Ratio)) %>% as.data.frame()
rownames(m6A.list$m6A$zma) <- m6A.list$m6A$zma$Gene

tmp.df <- rbind(exp.list$exp$osa_own, exp.list$exp$riceNcbJgf1, exp.list$exp$riceNcbJgf2)
exp.list$exp$osa <- tmp.df %>% group_by(Gene) %>% 
  filter(MeanExp == max(MeanExp)) %>% as.data.frame()
rownames(exp.list$exp$osa) <- exp.list$exp$osa$Gene
tmp.df <- rbind(m6A.list$m6A$osa_own, m6A.list$m6A$riceNcbJgf1, m6A.list$m6A$riceNcbJgf2)
m6A.list$m6A$osa <- tmp.df %>% group_by(Gene) %>% 
  filter(Ratio == max(Ratio)) %>% as.data.frame()
rownames(m6A.list$m6A$osa) <- m6A.list$m6A$osa$Gene

tmp.df <- rbind(exp.list$exp$slyGbQgz1, exp.list$exp$slyGbQgz2)
exp.list$exp$sly <- tmp.df %>% group_by(Gene) %>% 
  filter(MeanExp == max(MeanExp)) %>% as.data.frame()
rownames(exp.list$exp$sly) <- exp.list$exp$sly$Gene
tmp.df <- rbind(m6A.list$m6A$slyGbQgz1, m6A.list$m6A$slyGbQgz2)
m6A.list$m6A$sly <- tmp.df %>% group_by(Gene) %>% 
  filter(Ratio == max(Ratio)) %>% as.data.frame()
rownames(m6A.list$m6A$sly) <- m6A.list$m6A$sly$Gene

tmp.df <- rbind(exp.list$exp$ath_own, exp.list$exp$athCrBgd1, exp.list$exp$athCrBgd2)
exp.list$exp$ath <- tmp.df %>% group_by(Gene) %>% 
  filter(MeanExp == max(MeanExp)) %>% as.data.frame()
rownames(exp.list$exp$ath) <- exp.list$exp$ath$Gene
tmp.df <- rbind(m6A.list$m6A$ath_own, m6A.list$m6A$athCrBgd1, m6A.list$m6A$athCrBgd2)
m6A.list$m6A$ath <- tmp.df %>% group_by(Gene) %>% 
  filter(Ratio == max(Ratio)) %>% as.data.frame()
rownames(m6A.list$m6A$ath) <- m6A.list$m6A$ath$Gene
rm(tmp.df)


for (sp.tmp in species13.names){
  write.table(m6A.list$m6A[[sp.tmp]][, c("Gene", "Ratio")], 
              file = paste0("workflow/gene_list/", sp.tmp, "_m6A_abundances.txt"),
              quote = F, row.names = F, col.names = F, sep = "\t")
  # write.table(exp.list[[sp.tmp]], file = paste0("workflow/gene_list/", sp.tmp, "_exp_gene.txt"),
  #             quote = F, row.names = F, col.names = F, sep = "\t")
  # write.table(m6A.list[[sp.tmp]], file = paste0("workflow/gene_list/", sp.tmp, "_m6A_gene.txt"),
  #             quote = F, row.names = F, col.names = F, sep = "\t")
  # cat(sp.tmp, length(m6A.list[[sp.tmp]]), length(exp.list[[sp.tmp]]),
  #     sum(!m6A.list[[sp.tmp]]%in%exp.list[[sp.tmp]]), "\n")
}


# save variates
save(exp.list, m6A.list, species12.names, species13.names, 
     file = "part1_m6A_peaks_genes.RData")
