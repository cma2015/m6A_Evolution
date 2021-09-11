# Transcriptome-wide mapping of m6A modifications for 13 plant species
#

# current wroking directory
rm(list=ls())
setwd("~/a2z/m6A_13spp/Part_1/")
options(stringsAsFactors = F)

# loading libraries ####
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(DiffBind)
library(rtracklayer)

# reload data
load("data1_m6A_peaks_genes.RData")
load("data2_gene_attributes.RData")

link.data.path <- "../associated_data/"
comlab.name <- "m6A methylation ratio (%)"
system(paste0("mkdir -p figures"))

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

# plot for overlap between replication
CorrelationScatterPlot <- function(x, y, pointCol = rgb(0,0,0,0.7),
                                   legendPos = "topleft", legendCex = 1,
                                   Xlab = "",Ylab = "", ... ){
  ifval <- is.infinite(x)|is.infinite(y)
  x <- x[!ifval]
  y <- y[!ifval]
  valna <- is.na(x)|is.na(y)
  r2 <- cor(x[!valna], y[!valna], method = "pearson")
  pVal <- cor.test(x[!valna], y[!valna], method = "pearson")$p.value
  rp <- vector('expression',1)
  rp[1] <- substitute(expression(PCC == valueA),
                      list(valueA = format(r2,dig = 3)))[2]
  plot(x, y, pch=19, cex=0.3, las=1, xlab = Xlab,ylab = Ylab,col = pointCol,...)
  legend(legendPos, inset=-0.01,legend = rp, bty = 'n', 
         cex=legendCex)
}

CCvaluePlot <- function(df_input, var1, var2, text_pos, xlab, ylab){
  df_input[, var1] <- as.numeric(df_input[, var1])
  df_input[, var2] <- as.numeric(df_input[, var2])
  ##
  print(cor.test(df_input[, var1], df_input[, var2]))
  fit_stat <- summary(lm(df_input[, var1] ~ df_input[, var2]))
  R2_adj <- fit_stat$adj.r.squared
  p_value <- fit_stat$coefficients[2,4]
  out_plot <- ggplot(df_input, aes_string(x=var1, y=var2)) +
    geom_point(size = 3) +
    labs(x= xlab, y= ylab ) +
    theme_bw(base_size = 12)+ 
    theme(axis.line = element_line(colour = "black",size=0.5),
          axis.text = element_text(colour="black"),
          legend.text=element_text(colour="black"),
          legend.title=element_text(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_smooth(method="lm", se=F, color="#A3A2A2") +
    annotate("text", label = sprintf('Pearson\'s r = %.2f', 
                                     cor(df_input[, var1], df_input[, var2])), 
             x = text_pos[1], y = text_pos[2], size = 4) +
    annotate("text", label = sprintf('P-value = %.2g', p_value), 
             x = text_pos[3], y = text_pos[4], size = 4)
  out_plot
}


# Supplementary Data 3 peaks and genes ####
for (tmp.sp in c("maizePgLyj1", "maizePgLyj2", "maizePgLyj3")){ #species12.names
  for (rep in 1:2){ #1:3
    peak.each.rep <- read.table(paste0("workflow/data/", tmp.sp, "/", tmp.sp, "_R", rep, "_cut.bed"),
                                sep = "\t", stringsAsFactors = F)
    cat(tmp.sp, paste0("rep", rep), nrow(peak.each.rep), "\n")
  }
  # if(tmp.sp%in%c("ath", "zma", "osa")){tmp.sp.add <- paste0(tmp.sp, "_own")}else{tmp.sp.add <- tmp.sp}
  # cat(tmp.sp, nrow(m6A.list$peaks[[tmp.sp.add]]),
  #     length(unique(m6A.list$peaks[[tmp.sp.add]]$Gene)), length(m6A.list[[tmp.sp.add]]), "\n")
}


# Supplementary Data 4 summary table ####
summary.table <- as.data.frame(matrix(nrow = 13, ncol = 17))
rownames(summary.table) <- species13.names
colnames(summary.table) <- c("PCGs", "Exps", "Exp.PCGs","m6A.PCGs", "m6A.ratio", 
                             "OGs", "OGs.PCGs","Exp.OGs", "ExpOGs.OGs", "m6A.OGs", 
                             "m6AOGs.ratio","SGs",	"SGs.PCGs", "Exp.SGs", 
                             "ExpSGs.SGs", "m6A.SGs", "m6ASGs.ratio")
# c(1:2,4:6,8,10:12,14,16:17)
for (sp.tmp in species13.names){
  gene_names <- gene.attributes.df[gene.attributes.df$Species==sp.tmp, "Gene"]
  tmp.exp.genes <- exp.list[[sp.tmp]]
  tmp.m6A.genes <- m6A.list[[sp.tmp]]
  OLGs <- OGs.only.13spp.list[[sp.tmp]]
  LSGs <- setdiff(gene_names, OLGs)
  summary.table[sp.tmp, ] <- c(length(gene_names), 
                               length(tmp.exp.genes),
                               round(length(tmp.exp.genes)/length(gene_names)*100, 2),
                               length(tmp.m6A.genes),
                               round(length(tmp.m6A.genes)/length(gene_names)*100, 2),
                               length(OLGs),
                               round(length(OLGs)/length(gene_names)*100, 2),
                               length(intersect(OLGs, tmp.exp.genes)),
                               round(length(intersect(OLGs, tmp.exp.genes))/length(OLGs)*100, 2),
                               length(intersect(OLGs, tmp.m6A.genes)),
                               round(length(intersect(OLGs, tmp.m6A.genes))/length(OLGs)*100, 2),
                               length(LSGs),
                               round(length(LSGs)/length(gene_names)*100, 2),
                               length(intersect(LSGs, tmp.exp.genes)),
                               round(length(intersect(LSGs, tmp.exp.genes))/length(LSGs)*100, 2),
                               length(intersect(LSGs, tmp.m6A.genes)),
                               round(length(intersect(LSGs, tmp.m6A.genes))/length(LSGs)*100, 2)
                               )
}
summary.table$Type <- c(rep("Dicots", 6), rep("Monocots",6), "ppa")
output.summary.table <- summary.table[, c(1,4:6,10:12,16:17)]


# Figure 1 ####
# Percentage of expressed genes and m6A genes
summary.df <- output.summary.table
summary.df$species <- factor(rownames(summary.df), levels = rev(species13.names))
m6A.ratio.bar <- ggplot(summary.df, aes(x=species, y=m6A.ratio, fill="A")) + 
  scale_fill_manual(values = "#8A6060") +
  geom_bar(stat='identity', position=position_dodge(width=0.5)) +
  theme_bw(base_size = 12)+ 
  theme(axis.text = element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("") + ylab("m6A methylation ratio (%)") + coord_flip ()

pdf("figures/Figure1A.pdf", width = 3, height = 7.5)
m6A.ratio.bar
dev.off()


# m6A genes and background ####
gg.df <- read.table(paste0(link.data.path, "figure_data_old/Figure4.csv"), 
                       sep = ",", header = T, stringsAsFactors = F, row.names = 1)
summary.df$GenomeSize <- gg.df[as.character(summary.df$species), "GenomeSize"]
summary.df$GeneNum <- gg.df[as.character(summary.df$species), "GeneNum"]

genome.plot <- CCvaluePlot(df_input = summary.df, var1 ="GenomeSize", var2 = "m6A.ratio", 
           text_pos=c(5000,30,5000,27), xlab = "Genome size (MB)", ylab = comlab.name )

gene.plot <- CCvaluePlot(df_input = summary.df, var1 ="GeneNum", var2 = "m6A.ratio", 
           text_pos=c(50000,30,50000,27), xlab = "Gene number", ylab = comlab.name )

pdf("figures/Figure1B_C.pdf", width = 5, height = 8, useDingbats = F)
genome.plot + gene.plot + plot_layout(ncol = 1, guides = "collect") 
dev.off()


# Supplementary Fig. 6: m6A levels and YTH-domain proteins ####
tmp.species <- c("ath_own", "gar", "ghi", "pvu", "gma", "slyGbQgz1", "sbi", 
                 "zma_own", "ata", "tdi", "tae", "osa_own", "ppa")

level.domain.df <- data.frame(
  "spp"=tmp.species,
  "scale.levels"=sapply(m6A.list$m6A, function(x){
    sum(log2(x$Ratio))*sum(log2(x$Inputsum))
    })[tmp.species],
  "total.levels"=sapply(m6A.list$m6A, function(x){sum(x$Ratio)})[tmp.species],
  "YTHdomain"=c(13, 13, 26, 11, 19, 9, 12, 24, 13, 26, 39, 12, 4),
  "Writer.number"=c(6, 7, 14, 6, 15, 6, 6, 9, 8, 14, 22, 9, 10)
  )

Writer.plot <- CCvaluePlot(df_input = level.domain.df, var1 ="total.levels", 
                           var2 = "Writer.number", 
                           text_pos=c(85000, 18, 85000, 16), 
                           xlab = "Global m6A level", ylab = "Number of writer genes" )
YTH.plot <- CCvaluePlot(df_input = level.domain.df, var1 ="total.levels", 
                        var2 = "YTHdomain", 
                        text_pos=c(75000, 35, 75000, 32), 
                        xlab = "Global m6A level", ylab = "Number of reader genes")

pdf("figures/FigureS6.protein_m6A_levels.pdf", width = 8, height = 4, useDingbats = F)
Writer.plot + YTH.plot 
dev.off()


# Figure 1D ####
summary.table.df <- data.frame(
  "Type"=rep(c("OGs", "SGs"), each=13),
  "Values"=c(output.summary.table$m6AOGs.ratio,
             output.summary.table$m6ASGs.ratio))

pdf("figures/Figure1D.pdf", width = 3.3, height = 3)
ggplot(summary.table.df, aes(x=Type, y=Values, fill=Type)) +
  geom_violin() +
  geom_boxplot(width=0.2, fill="#7F7F7F", outlier.color = NA)+
  scale_fill_manual(values = c("#3D7970", "#BB8F68"))+
  theme_bw(base_size = 12)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none")+
  stat_compare_means(method = 'wilcox.test',label.y = 83, size=2) +
  xlab("") + ylab(comlab.name)
dev.off()


# Supplementary Figure 3 ####
for (sample.name in species12.names){
  cat(paste0("====", sample.name, "===="))
  tamoxifen <- dba(sampleSheet=paste0("workflow/data/", sample.name, "/DiffBind_corrplot.csv"))
  peak.granges <- lapply(tamoxifen$samples$Peaks, import)
  peak.granges.list <- GRangesList(peak.granges)
  png(paste0("workflow/venn_corrplot/", sample.name, "_corplot.png"), 
      res = 300, width = 2200, height = 1550)
  layout(matrix(c(1:6), 2, 3, byrow = FALSE))
  for (cycVec in list(c(1,2), c(1,3), c(2,3))){
    peak.coverage.res <- coverage(peak.granges.list[cycVec])
    covered.ranges.res <- IRanges::slice(peak.coverage.res, lower=1, rangesOnly=T)
    covered.granges.res <- GRanges(covered.ranges.res)
    cat(sample.name, "_", cycVec[1], "_", cycVec[2], "DBA_count_start!!\n")
    tamoxifen.count <- dba.count(tamoxifen, score = "DBA_SCORE_RPKM_FOLD", 
                                 peaks = covered.granges.res,
                                 filter = NULL, bUseSummarizeOverlaps = T)
    CorrelationScatterPlot(log2(tamoxifen.count$peaks[[cycVec[1]]][,5]+1), 
                           log2(tamoxifen.count$peaks[[cycVec[2]]][,5]+1), 
                           Xlab="", Ylab="")
    mtext(side=1, text=paste0(sample.name, "_IP_Rep", cycVec[1]), line=2, cex = 0.6)
    mtext(side=2, text=paste0(sample.name, "_IP_Rep", cycVec[2]), line=2, cex = 0.6)
    CorrelationScatterPlot(log2(tamoxifen.count$peaks[[cycVec[1]]][,7]+1), 
                           log2(tamoxifen.count$peaks[[cycVec[2]]][,7]+1), 
                           Xlab="", Ylab="")
    mtext(side=1, text=paste0(sample.name, "_Input_Rep", cycVec[1]), line=2, cex = 0.6)
    mtext(side=2, text=paste0(sample.name, "_Input_Rep", cycVec[2]), line=2, cex = 0.6)
  }
  dev.off()
}


# Supplementary Fig. 4: enrichment plot ####
input.df <- read.table(paste0(link.data.path, "/figure_data/FigureS4.csv"),
                       sep = ",", stringsAsFactors = F, header = T)
input.df$region <- factor(input.df$region, 
                            levels = c("5UTR", "star_codon", 
                                       "CDS", "stop_codon", "3UTR"))
input.df$sp <- factor(input.df$sp, levels=species13.names)

enrich.bar.plot <- ggplot(input.df, aes(x=region, y=relativeEnrichment, fill=region)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c('#BDD1C5','#4D6372','#E7D86E','#558ED5',"#A984A8"))+
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid =element_blank(),
        legend.position = "bottom") +
  facet_wrap(sp~., ncol = 3) + 
  labs(x="", y="Relative enrichment") +
  scale_y_continuous(breaks = seq(0, 10, 2))

pdf("figures/FigureS4.relative enrichment of m6A peaks.pdf", width = 8.15, height = 9.5)
enrich.bar.plot 
dev.off()


# Supplementary Data 5 ####
# two.spp.list <- list(c("ath", "gar"), c("ath", "pvu"), c("ath", "sly"),
#                      c("ath", "sbi"), c("ath", "ata"), c("ath", "osa"), 
#                      c("gar", "pvu"), c("gar", "sly"), c("gar", "sbi"),
#                      c("gar", "ata"), c("gar", "osa"), c("pvu", "sly"),
#                      c("pvu", "sbi"), c("pvu", "ata"), c("pvu", "osa"),
#                      c("sly", "sbi"), c("sly", "ata"), c("sly", "osa"),
#                      c("sbi", "ata"), c("sbi", "osa"), c("ata", "osa"))
two.spp.mat <- t(combn(species13.names, 2))
share.m6A.df <- as.data.frame(matrix(0, nrow = nrow(two.spp.mat), ncol = 8), 
                              stringsAsFactors = F)
rownames(share.m6A.df) <- apply(two.spp.mat, 1, function(x){paste0(x, collapse = "-")})

for (tmp.num in 1:nrow(two.spp.mat)){
  cat(tmp.two.spp, "\n")
  tmp.two.spp <- two.spp.mat[tmp.num, ]
  two_spp_name <- paste0(tmp.two.spp, collapse = "-")
  two.spp.ortho <- read.table(paste0(link.data.path, "Results_Jan16_34spp/Orthologues/Orthologues_",
                                     tmp.two.spp[1], "/", tmp.two.spp[1], "__v__", tmp.two.spp[2],".tsv"), 
                              header = T, sep = "\t")
  spp.ortho.filter <- two.spp.ortho[, 2:3]
  PS.two.spp.type <- apply(spp.ortho.filter, 1, function(x){
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
  # table(PS.two.spp.type)
  spp.ortho.filter <- spp.ortho.filter[PS.two.spp.type=="1:1", ]
  PS.two.spp.exp <- apply(spp.ortho.filter, 1, function(x){
    if(x[1]%in%exp.list[[tmp.two.spp[1]]]){
      a <- "1"
    }else{
      a <- "0"
    }
    if(x[2]%in%exp.list[[tmp.two.spp[2]]]){
      b <- "1"
    }else{
      b <- "0"
    }
    return(paste0(c(a,b), collapse = ":"))
  })
  # table(PS.two.spp.exp)
  spp.ortho.filter <- spp.ortho.filter[PS.two.spp.exp=="1:1", ]
  spp1.m6a <- sum(spp.ortho.filter[,1]%in%m6A.list[[tmp.two.spp[1]]])
  spp2.m6a <- sum(spp.ortho.filter[,2]%in%m6A.list[[tmp.two.spp[2]]])
  oo.count <- nrow(spp.ortho.filter)
  expect.value <- (spp1.m6a*spp2.m6a)/oo.count
  PS.two.spp.m6A <- apply(spp.ortho.filter, 1, function(x){
    if(x[1]%in%m6A.list[[tmp.two.spp[1]]]){
      a <- "1"
    }else{
      a <- "0"
    }
    if(x[2]%in%m6A.list[[tmp.two.spp[2]]]){
      b <- "1"
    }else{
      b <- "0"
    }
    return(paste0(c(a,b), collapse = ":"))
  })
  m6A.share <- sum(PS.two.spp.m6A=="1:1")
  pvalue <- phyper(m6A.share-1, spp1.m6a, oo.count-spp1.m6a, spp2.m6a, lower.tail=F)
  share.m6A.df[two_spp_name, ] <- c(FirstUp(tmp.two.spp[1]), FirstUp(tmp.two.spp[2]), 
                                    oo.count, spp1.m6a, spp2.m6a, 
                                    round(expect.value,digits = 2), 
                                    m6A.share, pvalue)
}
share.m6A.df$V8 <- format(as.numeric(share.m6A.df$V8), scientific = T, digits = 3)

abbr.names <- c("Ath", "Gar", "Ghi", "Pvu", "Gma", "Sly", "Sbi", "Zma", "Ata",
                "Tdi", "Tae", "Osa", "Ppa")
bi.names <- c("Arabidopsis thaliana", "Gossypium arboreum","Gossypium hirsutum",
              "Phaseolus vulgaris", "Glycine max", "Solanum lycopersicum",
              "Sorghum bicolor", "Zea mays", "Aegilops tauschii", "Triticum dicoccoides",
              "Triticum aestivum", "Oryza sativa", "Physcomitrella patens")
name.change.df <- data.frame("abbr" = abbr.names, 
                             "biname" = bi.names)
rownames(name.change.df) <- name.change.df$abbr
share.m6A.df$V1 <- name.change.df[share.m6A.df$V1, 2]
share.m6A.df$V2 <- name.change.df[share.m6A.df$V2, 2]
write.table(share.m6A.df, "figures/share_m6A_OGs.txt" ,sep = "\t", quote = F, row.names = F)

# 
length(unique(intersect(m6A.list$ath_own, m6A.list$athCrBgd1)))/length(m6A.list$athCrBgd1)*100
length(unique(intersect(m6A.list$ath_own, m6A.list$athCrBgd2)))/length(m6A.list$athCrBgd2)*100
length(unique(intersect(m6A.list$zma_own, m6A.list$maizePgLyj1)))/length(m6A.list$maizePgLyj1)*100
length(unique(intersect(m6A.list$zma_own, m6A.list$maizePgLyj2)))/length(m6A.list$maizePgLyj2)*100
length(unique(intersect(m6A.list$zma_own, m6A.list$maizePgLyj3)))/length(m6A.list$maizePgLyj3)*100
length(unique(intersect(m6A.list$zma_own, m6A.list$maizePpMc)))/length(m6A.list$maizePpMc)*100
length(unique(intersect(m6A.list$zma_own, m6A.list$maizePpHy)))/length(m6A.list$maizePpHy)*100
length(unique(intersect(m6A.list$osa_own, m6A.list$riceNcbJgf1)))/length(m6A.list$riceNcbJgf1)*100
length(unique(intersect(m6A.list$osa_own, m6A.list$riceNcbJgf2)))/length(m6A.list$riceNcbJgf2)*100
