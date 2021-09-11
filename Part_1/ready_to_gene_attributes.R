# gene information
# Time: 2021-06-03

# 
options(stringsAsFactors = F)

# reload data
load("data1_m6A_peaks_genes.RData")
link.data.path <- "../associated_data/"

# Phylogenetic relationships among 34 species ####
species34.names <- c("ath", "Bstricta", "Thassleriana", "gar", "ghi", "Egrandis", 
               "pvu", "gma", "Mtruncatula", "Ptrichocarpa", "Vvinifera", 
               "Stuberosum", "sly", "Sitalica", "Sviridis", "Pvirgatum", "sbi", 
               "zma", "ata", "tdi", "tae", "Hvulgare", "Bdistachyon", "osa", 
               "Eguineensis", "Spolyrhiza", "Atrichopoda", "Paab","Gibi",
               "Smoellendorffii", "ppa", "Mpolymorpha", "Creinhardtii", 
               "MpusillaCCMP1545")

OGs.34spp <- read.table(paste0(link.data.path, "Results_Jan16_34spp/Orthogroups/Orthogroups.tsv"), 
                        sep = "\t", header = T, row.names = 1, 
                        stringsAsFactors = F)[, species34.names]
OGs.34spp.count <- t(apply(OGs.34spp, 1, function(x){
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
Hordeeae <- list(unlist(AT), "Hvulgare") # Hordeeae (formerly Triticeae)
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

PS.level.names <- c("Chlorophyta", "Embryophyta", "Stomatophyta", "Tracheophyta", 
                      "Spermatophyta", "Acrogymnospermae", "Angiosperm", 
                      "Mesangiospermae", "Monocots", "Commelinids", "Poaceae", 
                      "BOP", "Pooideae", "Hordeeae", "AT", "Triticum", "PACMAD", 
                      "Andropogoneae", "Paniceae", "Setaria", "Eudicots", 
                      "Rosids", "Eurosids","Malvids", "BM", "Brassicales", "Brassicaceae", 
                      "Malvaceae", "Fabids","Fabaceae","Phaseoleae","Solanaceae")
all.species.names <- unlist(Chloroplastida)

PS.details.names <- rep("Chloroplastida", nrow(OGs.34spp.count))
for (tmp.PS in PS.level.names){
  tmp.PS.list <- eval(parse(text = tmp.PS))
  tmp.PS.others <- setdiff(all.species.names, unlist(tmp.PS.list))
  tmp.PS.values <- rep(0, nrow(OGs.34spp.count))
  index2 <- apply(OGs.34spp.count[, tmp.PS.others], 1, function(x){sum(x>0)}) == 0
  if (is.character(tmp.PS.list)){
    if (length(tmp.PS.list) >= 2){
      index1 <- apply(OGs.34spp.count[,tmp.PS.list], 1, function(x){sum(x>0)}) >= 2
      PS.details.names[index1 & index2] <- tmp.PS
    }
  } else {
    for (PS_each in 1:length(tmp.PS.list)){
      if (length(tmp.PS.list[[PS_each]]) == 1){
        each.PS.value <- OGs.34spp.count[, tmp.PS.list[[PS_each]]] > 0
      } else {
        each.PS.value <- apply(OGs.34spp.count[, tmp.PS.list[[PS_each]]], 1, 
                               function(x){sum(x > 0)}) > 0
      }
      tmp.PS.values <- tmp.PS.values + each.PS.value
    }
    PS.details.names[tmp.PS.values >= 2 & index2] <- tmp.PS
  }
  cat(tmp.PS, "\n")
}

# species-specific OGs
for(tmp.sp in all.species.names){
  tmp.PS.others <- setdiff(all.species.names, tmp.sp)
  index2 <- apply(OGs.34spp.count[,tmp.PS.others], 1, function(x){sum(x>0)}) == 0
  index1 <- OGs.34spp.count[, tmp.sp] != 0
  PS.details.names[index1 & index2] <- "LS_OGs" 
}

# Get the PS level for each OGs
table(PS.details.names)[intersect(c("Chloroplastida", PS.level.names, "LS_OGs"), 
                                  unique(PS.details.names))]


# add the PS level to data frame
OGs.34spp <- cbind(OGs.34spp, "PSnames" = PS.details.names)
# Chloroplastida: 6059 Embryophyta: 2864 
# Stomatophyta: 739 Tracheophyta: 628 Spermatophyta:4006

OGs.34spp.count <- cbind(OGs.34spp.count, "PSnames" = PS.details.names)
# write.table(OGs.34spp.count, file = "OG_gene_count_and_PS_levels.txt", quote = F, sep = "\t")

for(ii in species34.names){
  tmp.genes <- OGs.34spp[OGs.34spp[, ii]!="-", ii]
  tmp.genes <- unlist(strsplit(tmp.genes, ", "))
  cat(ii, length(tmp.genes)==length(unique(tmp.genes)), "\n")
}


# 13 species ####
OGs.13spp <- OGs.34spp[, c(species13.names, "PSnames")]
tmp.retain.index <- apply(OGs.13spp[, 1:13], 1, function(x){sum(x == "")}) < 13
table(tmp.retain.index)

OGs.13spp <- OGs.13spp[tmp.retain.index, ]
OGs.13spp.count <- OGs.34spp.count[tmp.retain.index, species13.names]


# gene attributes data frame ####
gene.attributes.list <- list()

for (tmp.sp in species13.names){
  tmp.gene.df <- read.table(paste0(link.data.path, "pep_cds/gene_names/", 
                                   tmp.sp, "_gene_cDNA_pep.txt"), 
                       sep = "\t")[,1:2]
  # gene names
  rownames(tmp.gene.df) <- tmp.gene.df[,1]
  colnames(tmp.gene.df) <- c("Gene", "Transcript")
  tmp.exp.val <- round(exp.list$exp[[tmp.sp]]$MeanExp, 2)
  names(tmp.exp.val) <- exp.list$exp[[tmp.sp]]$Gene
  tmp.gene.df$Exp <- tmp.exp.val[rownames(tmp.gene.df)]
  tmp.m6A.val <- round(m6A.list$m6A[[tmp.sp]]$Ratio, 2)
  names(tmp.m6A.val) <- m6A.list$m6A[[tmp.sp]]$Gene
  tmp.gene.df$m6A <- tmp.m6A.val[rownames(tmp.gene.df)]
  cat(tmp.sp, nrow(tmp.gene.df), sum(tmp.gene.df[,1]%in%exp.list[[tmp.sp]]), 
      sum(tmp.gene.df[,1]%in%m6A.list[[tmp.sp]]), "\n")
  # Add OGs
  OGs.length <- sapply(strsplit(OGs.13spp[, tmp.sp], ", "), length)
  tmp.gene.names <- unlist(strsplit(OGs.13spp[OGs.13spp[, tmp.sp] != "", tmp.sp], ", "))
  OG.rep.names <- rep(rownames(OGs.13spp), OGs.length)
  PS.tmp.names <- rep(OGs.13spp$PSnames, OGs.length)
  names(OG.rep.names) <- names(PS.tmp.names) <- tmp.gene.names
  tmp.gene.df$OGs <- OG.rep.names[rownames(tmp.gene.df)]
  tmp.gene.df$PSs <- PS.tmp.names[rownames(tmp.gene.df)]
  # one copy genes
  LS.one.index <- is.na(tmp.gene.df$PSs)
  tmp.gene.df$PSs[LS.one.index] <- "LS_OC"
  tmp.gene.df$OGs[LS.one.index] <- paste0("LS", 1:sum(LS.one.index))
  # Add duplication
  tmp.gene.df$Duplication <- "singleton"
  for(dup in c("wgd", "tandem", "proximal", "transposed", "dispersed")){
    tmp.dup.df <- read.table(paste0(link.data.path, "Duplication_Genes/", 
                                    tmp.sp, "/pep_DupGen_results/", tmp.sp, ".", 
                                    dup, ".genes-unique"), sep = "\t", 
                             stringsAsFactors = F, header = T)
    tmp.gene.df[tmp.dup.df[,1], "Duplication"] <- dup 
  }
  # Add Transcription and Kinase
  tmp.tf.df <- as.matrix(read.table(paste0(link.data.path, "TF/", tmp.sp, 
                                        ".fa_output/tf_classification.txt"), 
                                 stringsAsFactors = F, sep = "\t", row.names = 1))
  tmp.tf.names <- tmp.tf.df[, 1]
  names(tmp.tf.names) <- rownames(tmp.tf.df)
  tmp.gene.df$Species <- tmp.sp
  tmp.gene.df$TFs <- tmp.tf.names[rownames(tmp.gene.df)]
  tmp.kinase.df <- read.table(paste0(link.data.path, "TF/", tmp.sp, 
                                  ".fa_output/shiu_classification.txt"), 
                           stringsAsFactors = F, sep = "\t", row.names = 1)
  tmp.k.names <- tmp.kinase.df[, 1]
  names(tmp.k.names) <- rownames(tmp.kinase.df)
  tmp.gene.df$Kinases <- tmp.k.names[rownames(tmp.gene.df)]
  tmp.gene.df[is.na(tmp.gene.df)] <- ""
  # final table
  gene.attributes.list[[tmp.sp]] <- tmp.gene.df
}

gene.attributes.df <- do.call(rbind, gene.attributes.list)
# write.table(gene.attributes.df, "gene.attributes.df.txt", quote = F, sep = "\t", row.names = F)
PS.all.names <-  c("Chloroplastida", PS.level.names, "LS_OGs", "LS_OC")


# The OGs and SGs in 13 species ####
gene.13spp.group <- read.table(paste0(link.data.path, 
                                      "Results_Oct22_13spp/Orthogroups/Orthogroups.tsv"),
                               sep = "\t", header = T, row.names = 1, 
                               stringsAsFactors = F)[, species13.names]
gene.group.index <- apply(gene.13spp.group[,1:13], 1, function(x){sum(x != "")})
OGs.only.13spp.list <- apply(gene.13spp.group[gene.group.index >1, ], 2, function(x){
  unlist(strsplit(x, ", "))
})


# save data
save(OGs.13spp, OGs.13spp.count, OGs.34spp, OGs.34spp.count, PS.all.names, 
     OGs.only.13spp.list, gene.attributes.df, 
     file = "data2_gene_attributes.RData")
