# Divergence of m6A modification between homologous genes in the context of polyploidization in plants
#

# current wroking directory
rm(list=ls())
setwd("~/a2z/m6A_13spp/Part_4/")
options(stringsAsFactors = F)

## loading libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggforce)

# reload data
load("../Part_1/data1_m6A_peaks_genes.RData")
load("../Part_1/data2_gene_attributes.RData")

link.data.path <- "../associated_data/"
comlab.name <- "m6A methylation ratio (%)"
system(paste0("mkdir -p figures"))
rownames(gene.attributes.df) <- gene.attributes.df$Gene

# define functions
FirstUp <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

BaseTheme <- function(){
  theme(axis.line = element_line(colour = "black",size=0.5),
        axis.text = element_text(colour="black"),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        axis.ticks = element_line(colour = "black",size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}

SubgenomeDivergence <- function(dataMat, colNum, 
                                speciesList, idName){
  out.tmp <- list()
  dataMat <- dataMat[, colNum]
  dataMat.index <- apply(dataMat, 1, function(x){
    return(x[1]%in%exp.list[[speciesList[1]]] & x[2]%in%exp.list[[speciesList[2]]])
  })
  dataMat.tmp <- dataMat[dataMat.index, ]
  dataMat.m6A.type <- apply(dataMat.tmp, 1, function(x){
    if(x[1]%in%m6A.list[[speciesList[1]]]){vec.tmp1 <- "1"}else{vec.tmp1 <- "0"}
    if(x[2]%in%m6A.list[[speciesList[2]]]){vec.tmp2 <- "1"}else{vec.tmp2 <- "0"}
    return(paste0(c(vec.tmp1, vec.tmp2), collapse = ":"))
  })
  dataMat.m6A.type[dataMat.m6A.type == "1:0"|dataMat.m6A.type == "0:1"] <- "DM"
  dataMat.m6A.type[dataMat.m6A.type == "0:0"] <- "NM"
  dataMat.m6A.type[dataMat.m6A.type == "1:1"] <- "IM"
  out.tmp[[1]] <- as.data.frame(table(dataMat.m6A.type))
  colnames(out.tmp[[1]]) <- c("Type", "Count")
  out.tmp[[1]]$Name <- idName
  dataMat.tmp$Type <- dataMat.m6A.type
  dataMat.tmp$Name <- idName
  out.tmp[[2]] <- dataMat.tmp
  return(out.tmp)
}

SlideMean <- function(x, windowsize=3, slide=1){
  idx1 <- seq(1, length(x) - windowsize+1, by=slide);
  idx1 + windowsize -> idx2;
  # idx2[idx2>(length(x)+1)] <- length(x)+1;
  c(0, cumsum(x)) -> cx;
  return(c(mean(x[1:2]), (cx[idx2]-cx[idx1])/windowsize, 
           mean(x[(length(x)-1):length(x)])));
}


PieStats <- function(df, x0, y0, r0, r1, amount, explode, label_perc) {
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


## 3'UTR length ####
RegionLengthPlot <- function(gene_pairs, species_name, first_region, second_region){
  tmp.spp.df <- apply(gene_pairs,1,function(x){
    tmp.df1 <- first_region[first_region$V1%in%x[1], 3:4]
    tmp.df1$Species <- species_name[1]
    tmp.df2 <- second_region[second_region$V1%in%x[2], 3:4]
    tmp.df2$Species <- species_name[2]
    tmp.index <- tmp.df1[,2]>100&tmp.df2[,2]>100
    rbind(tmp.df1[tmp.index, ], tmp.df2[tmp.index, ])
  })
  ##
  tmp.spp.df <- do.call(rbind, tmp.spp.df)
  colnames(tmp.spp.df) <- c("Region", "Length", "Species")
  tmp.spp.df$Length <- log2(tmp.spp.df$Length+1)
  tmp.spp.df$Species <- factor(tmp.spp.df$Species, levels = species_name)
  tmp.spp.df$Region <- factor(tmp.spp.df$Region, levels = c("5UTR", "CDS", "3UTR"))
  ##
  p <- ggplot(tmp.spp.df, aes(x=Region, y=Length, fill=Species)) +
    geom_boxplot(outlier.color = "white")+
    scale_fill_manual(values = c("#0072B2", "#E69F00"))+
    geom_point(position=position_jitterdodge(jitter.width = 0.05)) +
    theme_bw(base_size = 14)+
    theme(axis.text.x= element_text(angle=15,hjust = 1,vjust = 1, colour="black",size=10))+
    BaseTheme() +
    stat_compare_means(method = 'wilcox.test',label.y = max(tmp.spp.df$Length)+1, size=2) +
    xlab("") + ylab("Percentage")
  p
}

RegionDifference  <- function(gene_pairs, species_name, 
                              first_region, second_region, type_name){
  tmp.spp.df <- apply(gene_pairs,1,function(x){
    tmp.df1 <- first_region[first_region$V1%in%x[1], ]
    tmp.df2 <- second_region[second_region$V1%in%x[2], ]
    tmp.index <- tmp.df1[,4]>0 & tmp.df2[,4]>0
    tmp.df1 <- tmp.df1[tmp.index, ]
    tmp.df2 <- tmp.df2[tmp.index, ]
    if("test"=="notest"){
      # if(x[1]%in%m6A.list[[species_name[1]]]&x[2]%in%m6A.list[[species_name[2]]]){
      #   tmp.df1[,4] <- abs(tmp.df1[,4]-tmp.df2[,4])
      # } else if (x[1]%in%m6A.list[[species_name[1]]]){
      #   tmp.df1[,4] <- tmp.df1[,4]-tmp.df2[,4]
      # } else if (x[2]%in%m6A.list[[species_name[2]]]){
      #   tmp.df1[,4] <- tmp.df2[,4] - tmp.df1[,4]
      # } else {
      #   tmp.df1[,4] <- abs(tmp.df1[,4]-tmp.df2[,4])
      # }
    }else{
      tmp.df1[,4] <- tmp.df1[,4]-tmp.df2[,4]
    }
    tmp.df1[,2] <- tmp.df2[,1]
    tmp.df1$Type <- type_name
    tmp.df1
  })
  ##
  tmp.spp.df <- do.call(rbind, tmp.spp.df)
  colnames(tmp.spp.df) <- c("Gene1", "Gene2", "Region", "Difference", "Type")
  tmp.spp.df
}

RegionLenDiff <- function(gene_pair_df, species_name, 
                                 first_region, second_region, type_name){
  yes_yes <- gene_pair_df[gene_pair_df$Type=="IM",]
  yes_no <- gene_pair_df[gene_pair_df$Type=="DM",]
  no_no <- gene_pair_df[gene_pair_df$Type=="NM",]
  pp_out <- list()
  pp_out[[1]] <- RegionDifference(yes_yes, species_name, first_region, 
                                  second_region, type_name = "IM")
  pp_out[[2]] <- RegionDifference(yes_no, species_name, first_region, 
                                  second_region, type_name = "DM")
  pp_out[[3]] <- RegionDifference(no_no, species_name, first_region, 
                                  second_region, type_name = "NM")
  pp_out <- do.call(rbind, pp_out)
  pp_out$Species <- type_name
  pp_out
}


# read the subgenome information (gene pairs) ####
Pvu.vs.Gma <- read.table(paste0(link.data.path, "subgenome_clean/pvu_gma_pairs_20195.txt"), 
                      sep = "\t", header = F)
Sbi.vs.Zma <- read.table(paste0(link.data.path, "subgenome_clean/zma_sbi_pairs_17898.txt"), 
                      sep = "\t", header = T)
Gar.vs.GhiA <- read.table(paste0(link.data.path, "subgenome_clean/gar_ghi_subA_24468.txt"), 
                      sep = "\t", header = T)
TdiA.vs.TaeA <- read.table(paste0(link.data.path, "subgenome_clean/wheat_subA_17701.txt"), 
                         sep = "\t")
TdiB.vs.TaeB <- read.table(paste0(link.data.path, "subgenome_clean/wheat_subB_17375.txt"), 
                         sep = "\t")
Ata.vs.TaeD <- read.table(paste0(link.data.path, "subgenome_clean/wheat_subD_18475.txt"), 
                         sep = "\t")


# extract m6A gene pairs ####
Pvu.vs.Gma1.list <- SubgenomeDivergence(
  dataMat = Pvu.vs.Gma, colNum = c(1,2),
  speciesList = c("pvu","gma"), idName="Pvu_Gma1")
Pvu.vs.Gma2.list <- SubgenomeDivergence(
  dataMat = Pvu.vs.Gma, colNum = c(1,3),
  speciesList = c("pvu","gma"), idName="Pvu_Gma2")
Sbi.vs.Zma1.list <- SubgenomeDivergence(
  dataMat = Sbi.vs.Zma, colNum = c(1,2),
  speciesList = c("sbi","zma"), idName="Sbi_Zma1")
Sbi.vs.Zma2.list <- SubgenomeDivergence(
  dataMat = Sbi.vs.Zma, colNum = c(1,3),
  speciesList = c("sbi","zma"), idName="Sbi_Zma2")
Gar.vs.GhiA.list <- SubgenomeDivergence(
  dataMat = Gar.vs.GhiA, colNum = c(1,2),
  speciesList = c("gar","ghi"), idName="Gar_GhiA")
TdiA.vs.TaeA.list <- SubgenomeDivergence(
  dataMat = TdiA.vs.TaeA, colNum = c(1,2),
  speciesList = c("tdi","tae"), idName="TdiA_TaeA")
TdiB.vs.TaeB.list <- SubgenomeDivergence(
  dataMat = TdiB.vs.TaeB, colNum = c(1,2),
  speciesList = c("tdi","tae"), idName="TdiB_TaeB")
Ata.vs.TaeD.list <- SubgenomeDivergence(
  dataMat = Ata.vs.TaeD, colNum = c(1,2),
  speciesList = c("ata","tae"), idName="Ata_TaeD")

subgenome.df <- rbind(Pvu.vs.Gma1.list[[1]], Pvu.vs.Gma2.list[[1]], 
                      Sbi.vs.Zma1.list[[1]], Sbi.vs.Zma2.list[[1]], 
                      Gar.vs.GhiA.list[[1]], Ata.vs.TaeD.list[[1]], 
                      TdiB.vs.TaeB.list[[1]], TdiA.vs.TaeA.list[[1]])
tmp.merge <- rep("Identical pattern", nrow(subgenome.df))
tmp.merge[subgenome.df$Type%in%c("DM")] <- "Diverged pattern"
subgenome.df$mergeType <- tmp.merge


# Figure 4A ####
# Figure 4B ####
seq.var1 <- read.table(paste0(link.data.path, "figure_data/Fig4B/pvu_gma.csv"), 
                       sep = ",", header = T)
seq.var2 <- read.table(paste0(link.data.path, "figure_data/Fig4B/sbi_zma.csv"), 
                       sep = ",", header = T)
seq.var3 <- read.table(paste0(link.data.path, "figure_data/Fig4B/gar_ghi.csv"),
                       sep = ",", header = T)
seq.var4 <- read.table(paste0(link.data.path, "figure_data/Fig4B/tdi_tae.csv"), 
                       sep = ",", header = T)
seq.var5 <- read.table(paste0(link.data.path, "figure_data/Fig4B/ata_tae.csv"), 
                       sep = ",", header = T)
colnames(seq.var1) <- colnames(seq.var2) <- colnames(seq.var3) <- 
  colnames(seq.var4) <- colnames(seq.var5)
seq.var <- rbind(seq.var1, seq.var2, seq.var3, seq.var4, seq.var5)
seq.var$all.sites <- paste0(seq.var$region, "_", seq.var$position)
seq.var$all.sites <- factor(seq.var$all.sites, levels = unique(seq.var$all.sites))
seq.var$subgenome <- sapply(strsplit(seq.var$subgenome, "_"), function(x){
  paste0(FirstUp(x[1]), "_", FirstUp(x[2]))
})

seq.var.merge <- seq.var %>% 
  group_by(all.sites, type, subgenome) %>% 
  summarise("medianValues"=median(mutationRatio), "meanValues"=mean(mutationRatio)) %>% 
  as.data.frame()
table(seq.var.merge$type)
table(seq.var.merge$subgenome)
seq.var.merge$meanValues <- seq.var.merge$meanValues*100

window.df <- list()
for (ii in unique(seq.var.merge$subgenome)){
  for (jj in unique(seq.var.merge$type)){
    tmp.sub.df <- seq.var.merge[seq.var.merge$subgenome==ii & seq.var.merge$type==jj,]
    tmp.sub.df$windowValue <- SlideMean(tmp.sub.df$meanValues)
    window.df[[paste0(ii, jj)]] <- tmp.sub.df
  }
}

window.df.data <- do.call(rbind, window.df) 
window.df.data$type <- factor(window.df.data$type, levels = c("IM", "DM", "NM"))
window.df.data$subgenome <- factor(window.df.data$subgenome,
                                   levels=c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2",
                                            "GarA_GhiA", "TdiA_TaeA", "TdiB_TaeB", "AtaD_TaeD"))

pdf("figures/Figure4B.pdf", height = 3, width = 16)
ggplot(data = window.df.data, aes(x=all.sites, y=windowValue, color=type, group=type))+
  geom_line(size=0.8) + #geom_point(size=0.8) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom", 
        legend.direction = "horizontal") +
  BaseTheme() +
  xlab("") + ylab("") + 
  facet_grid(.~subgenome, scales = "free_x", space = "free_x") +
  scale_color_brewer(palette = "Dark2") +
  geom_rect(aes(xmin=1, xmax=10.5, ymin=-5, ymax=-1),
            fill="#AAAAAA", color=NA) +
  geom_rect(aes(xmin=10.5, xmax=20.5, ymin=-6, ymax=0),
            fill="#767777", color=NA) +
  geom_rect(aes(xmin=20.5, xmax=30, ymin=-5, ymax=-1),
            fill="#AAAAAA", color=NA)
dev.off()


## Figure 4C and Supplementary Fig. 12 Ka Ks ####
kaks.df <- read.table(paste0(link.data.path, "figure_data/Fig4C.csv"),
                           sep = ",", header = T)
kaks.df <- kaks.df[kaks.df$variable%in%c("ka", "ks"), ]
# colnames(kaks.df)[6] <- "divType"
kaks.df$divType <- factor(kaks.df$divType, levels = c("IM", "DM", "NM"))
kaks.df$spe <- sapply(strsplit(kaks.df$spe, "_"), function(x){
  paste0(FirstUp(x[1]), "_", FirstUp(x[2]))
})

subtype1 <- c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2")
subtype2 <- c("Gar_GhiA", "TdiA_TaeA", "TdiB_TaeB", "AtaD_TaeD")
for (tmp.type in c("subtype1", "subtype2")){
  tmp.type.value <- eval(parse(text = tmp.type))
  tmp.kaks.df <- kaks.df[kaks.df$spe%in%tmp.type.value, ]
  tmp.kaks.df$spe <- factor(tmp.kaks.df$spe, levels = tmp.type.value)
  if(tmp.type=="subtype1"){ssvalue <- "area";box.wid <- 0.12} else {ssvalue <- "width";box.wid <- 0.08}
  # Ka plot
  Ka_plot <- ggplot(tmp.kaks.df %>% 
                      filter(variable == "ka"), 
                    aes(x = divType, y = value, fill = divType)) +
    geom_violin(scale = ssvalue) +
    geom_boxplot(width = box.wid, fill = "#7F7F7F", outlier.colour = NA)+ #
    scale_fill_manual(values = c("#379677", "#C6612D", "#6B6BA2")) +
    theme_bw(base_size = 12)+
    theme(legend.title = element_blank())+
    BaseTheme() +
    xlab("") + ylab("Ka") +
    facet_grid(.~spe, scales = "free_x") +
    coord_cartesian(ylim = c(0, 0.27))
  # Ks plot
  Ks_plot <- ggplot(tmp.kaks.df %>% 
                      filter(variable=="ks"), 
                    aes(x = divType, y = value, fill = divType)) +
    geom_violin(scale = ssvalue) +
    geom_boxplot(width = box.wid, fill = "#7F7F7F", outlier.colour = NA)+ #
    scale_fill_manual(values = c("#379677", "#C6612D", "#6B6BA2")) +
    theme_bw(base_size = 12)+
    theme(legend.title = element_blank())+
    BaseTheme() +
    xlab("") + ylab("Ks") +
    coord_cartesian(ylim = c(0, 0.95)) + 
    facet_grid(.~spe, scales = "free_x")
  for(ii in unique(tmp.kaks.df$spe)){
    tmp.tmp.df <- tmp.kaks.df[tmp.kaks.df$spe==ii, ]
    sub.tmp.df <- tmp.tmp.df[tmp.tmp.df$variable=="ka", ]
    cat(ii, "ka", t.test(sub.tmp.df[sub.tmp.df$divType=="IM", "value"], 
                   sub.tmp.df[sub.tmp.df$divType=="DM", "value"], alternative = "less")$p.value,
        t.test(sub.tmp.df[sub.tmp.df$divType=="IM", "value"], 
               sub.tmp.df[sub.tmp.df$divType=="NM", "value"], alternative = "less")$p.value,
        t.test(sub.tmp.df[sub.tmp.df$divType=="DM", "value"], 
               sub.tmp.df[sub.tmp.df$divType=="NM", "value"], alternative = "less")$p.value,
        "\n")
    sub.tmp.df <- tmp.tmp.df[tmp.tmp.df$variable=="ks", ]
    cat(ii, "ks", t.test(sub.tmp.df[sub.tmp.df$divType=="IM", "value"], 
                   sub.tmp.df[sub.tmp.df$divType=="DM", "value"], alternative = "less")$p.value,
        t.test(sub.tmp.df[sub.tmp.df$divType=="IM", "value"], 
               sub.tmp.df[sub.tmp.df$divType=="NM", "value"], alternative = "less")$p.value,
        t.test(sub.tmp.df[sub.tmp.df$divType=="DM", "value"], 
               sub.tmp.df[sub.tmp.df$divType=="NM", "value"], alternative = "less")$p.value,
        "\n")
    
  }
  if(tmp.type=="subtype1"){
     outfile <- "figures/Figure4C.pdf"
    wid.value <- 7.5
  } else {
      outfile <- "figures/FigureS12.pdf"
      wid.value <- 12
  }
  # pdf(outfile, height = 7, width = wid.value )
  # print(Ka_plot + Ks_plot + plot_layout(ncol = 1, guides = "collect"))
  # dev.off()
}


## Figure 4D and Supplementary Fig. 13 sequence variation ####
utr3.var1 <- read.table(paste0(link.data.path, "figure_data/Fig4D/pvu_gma.csv"), 
                        sep = ",", header = T)
utr3.var2 <- read.table(paste0(link.data.path, "figure_data/Fig4D/sbi_zma.csv"), 
                        sep = ",", header = T)
utr3.var3 <- read.table(paste0(link.data.path, "figure_data/Fig4D/gar_ghi.csv"),
                        sep = ",", header = T)
utr3.var4 <- read.table(paste0(link.data.path, "figure_data/Fig4D/tdi_tae.csv"), 
                        sep = ",", header = T)
utr3.var5 <- read.table(paste0(link.data.path, "figure_data/Fig4D/ata_tae.csv"),
                        sep = ",", header = T)
colnames(utr3.var2) <- colnames(utr3.var3) <- colnames(utr3.var4) <- 
  colnames(utr3.var5) <- colnames(utr3.var1)
utr3.var <- rbind(utr3.var1, utr3.var2, utr3.var3, utr3.var4, utr3.var5)
utr3.var$subgenome <- sapply(strsplit(utr3.var$subgenome, "_"), function(x){
  paste0(FirstUp(x[1]), "_", FirstUp(x[2]))
})
utr3.var <- utr3.var[utr3.var$type%in%c("IM", "DM"), ]
table(utr3.var$subgenome)

subtype1 <- c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2")
subtype2 <- c("GarA_GhiA", "TdiA_TaeA", "TdiB_TaeB", "AtaD_TaeD")
for (tmp.type in c("subtype1", "subtype2")){
  tmp.type.value <- eval(parse(text = tmp.type))
  tmp.utr3.var <- utr3.var[utr3.var$subgenome%in%tmp.type.value, ]
  tmp.utr3.var$subgenome <- factor(tmp.utr3.var$subgenome, levels = tmp.type.value)
  tmp.utr3.var$type <- factor(tmp.utr3.var$type, levels = c("IM", "DM"))
  tmp.utr3.var$three_prime_UTR <- tmp.utr3.var$three_prime_UTR*100
  if(tmp.type=="subtype1"){
    outfile <- "figures/Figure4D.pdf"
    wid.value <- 5
    box.width <- 0.12
  } else {
    outfile <- "figures/FigureS13.pdf"
    wid.value <- 8
    box.width <- 0.08
  }
  pdf(outfile, height = 4, width = wid.value)
  print(ggplot(data = tmp.utr3.var, aes(x = type, y = three_prime_UTR, fill = type))+
          geom_violin(scale = "width") +
          geom_boxplot(width = box.width, fill = "#7F7F7F", outlier.colour = NA)+
          scale_fill_manual(values = c("#379677", "#C6612D")) +
          theme_bw(base_size = 12) +
          theme(legend.title = element_blank()) +
          BaseTheme() +
          xlab("") + ylab("3' UTR sequence Identity (%)") + 
          facet_grid(.~subgenome))
  dev.off()
  for(ii in unique(tmp.utr3.var$subgenome)){
    sub.tmp.df <- utr3.var[utr3.var$subgenome == ii, ]
    cat(ii, t.test(sub.tmp.df[sub.tmp.df$type == "IM", "three_prime_UTR"], 
                   sub.tmp.df[sub.tmp.df$type == "DM", "three_prime_UTR"], 
                   alternative = "greater")$p.value, "\n")
  }
}


## Figure 4E ####
mergeType <- rep("Identical pattern", nrow(subgenome.df))
mergeType[subgenome.df$Type=="DM"] <- "Diverged pattern"
subgenome.df$mergeType <- mergeType

subgenome.plots <- list()
for (ii in unique(subgenome.df$Name)){
  subgenome.ID.df <- subgenome.df %>% filter(Name == ii) %>% 
    group_by(Name, mergeType)%>% summarise("Value" = sum(Count)) %>% 
    as.data.frame()
  subgenome.plots[[ii]] <- ggplot(data = subgenome.ID.df %>% 
                                    PieStats(0, 0, 0, 1, Value, .1 * (mergeType == "Diverged pattern"), .5)) +
    ggforce::geom_arc_bar(aes(x0 = x0, y0 = y0, r0 = r0, r = r1,
                              start = start, end = end, explode = explode,
                              fill = factor(mergeType, levels = c("Diverged pattern", "Identical pattern")))) +
    geom_text(aes(x = x_lab, y = y_lab, label = scales::percent(Value/sum(Value), accuracy = 0.1)), size = 3) +
    coord_equal() + theme_void() + 
    labs(fill = expression("m"^"6"*"A patterns"), title = ii) + 
    scale_fill_manual(values = c("#F0861B", "#9FA0A0"))
}

pdf("figures/Figure4E.pdf", height = 5, width = 10)
wrap_plots(subgenome.plots, ncol = 2, guides = "collect")
dev.off()
wilcox.test(c(15.2, 15.9, 20.4, 20.0), c(10.3, 9.8, 11.1, 9.2))


# Figure 4F autopolyploidization and allopolyploidization merge IM, DM, NM ####
subtype2 <- c("Gar_GhiA", "TdiA_TaeA", "TdiB_TaeB", "AtaD_TaeD")
poly.type <- rep("autopolyploidization", nrow(kaks.df))
poly.type[kaks.df$spe%in%subtype2] <- "allopolyploidization"
kaks.df$poly.type <- poly.type
kaks.df$poly.type <- factor(kaks.df$poly.type, levels = c("autopolyploidization", "allopolyploidization"))

merge.ka.plot <- ggplot(kaks.df %>% filter(variable=="ka"), 
                        aes(x=poly.type, y=value, fill=poly.type)) +
  geom_violin() + geom_boxplot(width=0.05, fill="#7F7F7F", outlier.colour = NA)+ #
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 12)+
  theme(legend.title = element_blank(), legend.position = "none") +
  BaseTheme() +
  xlab("") + ylab("Ka") +
  coord_cartesian(ylim = c(0, 0.2))

merge.ks.plot <- ggplot(kaks.df %>%  filter(variable=="ks"),
                        aes(x=poly.type, y=value, fill=poly.type)) +
  geom_violin() + geom_boxplot(width=0.05, fill="#7F7F7F", outlier.colour = NA)+ #
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 12)+
  theme(legend.title = element_blank(), legend.position = "none")+
  BaseTheme() +
  xlab("") + ylab("Ks") +
  coord_cartesian(ylim = c(0, 0.8))

pdf("figures/Figure4F.pdf", height = 3, width = 7)
merge.ka.plot + merge.ks.plot + plot_layout(ncol = 2, guides = "collect")
dev.off()


# Supplementary Figure 14: Differnece of 3'UTR lengths ####
gar.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/gar.txt"), 
                            sep = "\t")
ghi.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/ghi.txt"), 
                            sep = "\t")
ata.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/ata.txt"), 
                            sep = "\t")
tdi.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/tdi.txt"), 
                            sep = "\t")
tae.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/tae.txt"), 
                            sep = "\t")
pvu.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/pvu.txt"), 
                            sep = "\t")
gma.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/gma.txt"), 
                            sep = "\t")
sbi.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/sbi.txt"), 
                            sep = "\t")
zma.region.df <- read.table(paste0(link.data.path, "pep_cds/region_length/zma.txt"), 
                            sep = "\t")

Pva.vs.Gma1.rlc <- RegionLenDiff(gene_pair_df = Pvu.vs.Gma1.list[[2]], species_name = c("pvu", "gma"),
                                 first_region = pvu.region.df, second_region = gma.region.df, 
                                 type_name = "Pvu_Gma1")
Pva.vs.Gma2.rlc <- RegionLenDiff(gene_pair_df = Pvu.vs.Gma2.list[[2]], species_name = c("pvu", "gma"),
                                 first_region = pvu.region.df, second_region = gma.region.df, 
                                 type_name = "Pvu_Gma2")
Sbi.vs.Zma1.rlc <- RegionLenDiff(gene_pair_df = Sbi.vs.Zma1.list[[2]], species_name = c("sbi", "zma"),
                                 first_region = sbi.region.df, second_region = zma.region.df,
                                 type_name = "Sbi_Zma1")
Sbi.vs.Zma2.rlc <- RegionLenDiff(gene_pair_df = Sbi.vs.Zma2.list[[2]], species_name = c("sbi", "zma"),
                                 first_region = sbi.region.df, second_region = zma.region.df,
                                 type_name = "Sbi_Zma2")
Gar.vs.GhiA.rlc <- RegionLenDiff(gene_pair_df = Gar.vs.GhiA.list[[2]], species_name = c("gar", "ghi"),
                                 first_region = gar.region.df, second_region = ghi.region.df, 
                                 type_name = "Gar_GhiA")
TdiA.vs.TaeA.rlc <- RegionLenDiff(gene_pair_df = TdiA.vs.TaeA.list[[2]], species_name = c("tdi", "tae"),
                                  first_region = tdi.region.df, second_region = tae.region.df, 
                                  type_name = "TdiA_TaeA")
TdiB.vs.TaeB.rlc <- RegionLenDiff(gene_pair_df = TdiB.vs.TaeB.list[[2]], species_name = c("tdi", "tae"),
                                  first_region = tdi.region.df,second_region = tae.region.df,
                                  type_name = "TdiB_TaeB")
Ata.vs.TaeD.rlc <- RegionLenDiff(gene_pair_df = Ata.vs.TaeD.list[[2]], species_name = c("ata", "tae"),
                                 first_region = ata.region.df, second_region = tae.region.df,
                                 type_name = "Ata_TaeD")
region.plot <- rbind(Pva.vs.Gma1.rlc, Pva.vs.Gma2.rlc, Sbi.vs.Zma1.rlc, Sbi.vs.Zma2.rlc,
                     Gar.vs.GhiA.rlc, Ata.vs.TaeD.rlc, TdiA.vs.TaeA.rlc, TdiB.vs.TaeB.rlc)

type8.names <- c("Pvu_Gma1", "Pvu_Gma2", "Sbi_Zma1", "Sbi_Zma2",
                "Gar_GhiA", "TdiA_TaeA", "TdiB_TaeB", "Ata_TaeD")
region.plot$Type <- factor(region.plot$Type, levels = c("IM","DM", "NM"))
region.plot$Species <- factor(region.plot$Species, 
                              levels = type8.names)
region.plot2  <- region.plot[region.plot$Region%in%"3UTR" & region.plot$Type%in%c("IM", "DM"),]
region.plot2 <- region.plot2[region.plot2$Difference >= -500 & region.plot2$Difference <= 500, ]


pdf("figures/Figure14A.pdf", width = 15, height = 7)
ggplot(data = region.plot2,
       mapping = aes(Difference, fill=Type)) +
  geom_density(size=0.7) + 
  facet_wrap(.~Species, scales = "free_x", ncol = 4) +
  coord_cartesian(xlim = c(-500,500)) +
  theme_bw(base_size = 14)+
  BaseTheme()+
  xlab("The difference of 3'UTR")+ylab("Density") +
  scale_fill_manual(values = c("#547B77", "#AB8472"))
dev.off()


for(ii in type8.names){
  tmp.sub.df <- region.plot2[region.plot2$Species == ii,]
  aIM <- tmp.sub.df[tmp.sub.df$Type=="IM", "Difference"]
  bDM <- tmp.sub.df[tmp.sub.df$Type=="DM", "Difference"]
  df <- merge(
    as.data.frame(density(aIM, from = min(c(aIM,bDM)), to = max(c(aIM,bDM)))[c("x", "y")]),
    as.data.frame(density(bDM, from = min(c(aIM,bDM)), to = max(c(aIM,bDM)))[c("x", "y")]),
    # as.data.frame(density(aIM, from = -500, to = 500)[c("x", "y")]),
    # as.data.frame(density(bDM, from = -500, to = 500)[c("x", "y")]),
    by = "x", suffixes = c(".aIM", ".bDM")
  )
  df$comp <- as.numeric(df$y.aIM > df$y.bDM)
  df$cross <- c(NA, diff(df$comp))
  plot(density(aIM), xlim=c(-500,500), col="green", main =ii)
  lines(density(bDM), col="red")
  points(df[which(df$cross != 0), c("x", "y.aIM")])
  point_values <- df[which(df$cross != 0), 1]
  point_right <- which(point_values>0)[1]
  if(point_right==1){
    break()
  }
  ## chisq.test
  var1 <- sum(aIM >= point_values[point_right-1] & aIM <= point_values[point_right])
  var2 <- length(aIM)
  var3 <- sum(bDM >= point_values[point_right-1] & bDM <= point_values[point_right])
  var4 <- length(bDM)
  cat(ii, point_values[point_right-1], point_values[point_right], 
      chisq.test(matrix(c(var1, var2-var1, var3, var4-var3), ncol=2, 
                        byrow = T))$p.value, "\n")
  # cat(ii, mean(df[df_index, 2]), mean(df[df_index, 3]), tmp_tval, "\n")
}


# Supplementary Figure 15 ####
soybean.duplicates <- Pvu.vs.Gma %>% filter(V2!="-", V3!="-")
soybean.duplicates <- c(soybean.duplicates[, 2], soybean.duplicates[, 3])
soybean.singletons <- c(Pvu.vs.Gma[, 2], Pvu.vs.Gma[, 3])
soybean.singletons <- setdiff(setdiff(soybean.singletons, "-"), soybean.duplicates)
maize.duplicates <- Sbi.vs.Zma %>% filter(B73maize1_Gene!="-", B73maize2_Gene!="-")
maize.duplicates <- c(maize.duplicates[, 2], maize.duplicates[, 3])
maize.singletons <- c(Sbi.vs.Zma[, 2], Sbi.vs.Zma[, 3])
maize.singletons <- setdiff(setdiff(maize.singletons, "-"), maize.duplicates)
sd.m6A <- sum(soybean.duplicates%in%m6A.list[["gma"]])/length(soybean.duplicates)*100
ss.m6A <- sum(soybean.singletons%in%m6A.list[["gma"]])/length(soybean.singletons)*100
md.m6A <- sum(maize.duplicates%in%m6A.list[["zma"]])/length(maize.duplicates)*100
ms.m6A <- sum(maize.singletons%in%m6A.list[["zma"]])/length(maize.singletons)*100
sm.df <- data.frame("type1" = rep(c("gma.d", "gma.s", "zma.d", "zma.s"), each=2),
                    "type2" = rep(c("m6A", "non-m6A"), 4),
                    "value" = c(sd.m6A, 100-sd.m6A, ss.m6A, 100-ss.m6A,
                                md.m6A, 100-md.m6A, ms.m6A, 100-ms.m6A)
                    )
chisq.test(matrix(c(8430, 12822-8430, 960, 1603-960), nrow = 2))

sm.df$type1 <- factor(sm.df$type1, levels = unique(sm.df$type1))
sm.df$type2 <- factor(sm.df$type2, levels = c("non-m6A", "m6A"))
pdf("figures/FigureS15.pdf", height = 4, width = 5)
ggplot(data = sm.df, aes(x=type1, y=value, fill=type2))+
  geom_bar(stat='identity', width=0.7) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle=30,hjust = 1,vjust = 1, colour="black",size=10)) +
  BaseTheme() +
  xlab("") + ylab(comlab.name) + coord_cartesian(ylim = c(25, 75)) +
  scale_fill_manual(values = c( "#CCCCCC", "#547975"))
dev.off()
