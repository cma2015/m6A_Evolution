# Divergence of m6A modification in orthologous genes with different evolutionary origins
#

# current wroking directory
rm(list=ls())
setwd("~/a2z/m6A_13spp/Part_2/")
options(stringsAsFactors = F)


# loading libraries
library(dplyr)
library(ggpubr)
library(ggplot2)
library(patchwork)


# reload data
load("../Part_1/data1_m6A_peaks_genes.RData")
load("../Part_1/data2_gene_attributes.RData")

link.data.path <- "../associated_data/"
comlab.name <- "m6A methylation ratio (%)"
system(paste0("mkdir -p figures"))


# define functions
FirstUp <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

BaseTheme <- function(){
  theme(panel.border = element_rect(colour = "black"),
        axis.line = element_line(colour = "black",size=0.5),
        axis.text = element_text(colour="black"),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        axis.ticks = element_line(colour = "black",size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}


# Figure 2A ####
# The proportion of m6A modification in monocots and dicots
DM.df <- read.table(paste0(link.data.path, "figure_data/Figure2A.csv"),
                    sep = ",", header = T)
DM.df$m6ARatio <- DM.df$m6ARatio*100
mono.di.abbr <- c("DDs", "DMs", "DPs", "MMs", "MDs", "MPs")
all(DM.df$kind%in%mono.di.abbr)
DM.df$kind <- factor(DM.df$kind, levels = mono.di.abbr)

DM.df %>% group_by(kind) %>% summarise("mean"=mean(m6ARatio), "sd"=sd(m6ARatio))

pdf("figures/Figure2A.pdf", width = 5, height = 6.5) # 42
ggplot(DM.df, aes(x=kind, y=m6ARatio, fill=kind)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.2, fill="#7F7F7F")+
  scale_fill_manual(values = c("#237E7F", "#309388", "#47A082",
                               "#826F74", "#A78889", "#D0B9B9")) +
  theme_bw(base_size = 14)+ 
  theme(legend.position = "none") + BaseTheme() +
  # coord_cartesian(ylim = c(20,85)) +
  stat_compare_means(method = 't.test',  # wilcox.test
                     method.args = list(alternative = "less"),
                     label.y =85, size=4,
                     comparisons = list( c("DDs", "DMs"), 
                                         c("DDs", "DPs"), 
                                         c("DMs", "DPs"),
                                         c("MMs", "MDs"), 
                                         c("MMs", "MPs"), 
                                         c("MDs", "MPs"))) +
  xlab("") + ylab(comlab.name)
dev.off()


# Figure 2B ####
gene_count <- table(gene.attributes.df$Species)
gene.attributes.df %>% filter(Species%in%species13.names) %>% 
  group_by(PSs, Species) %>% 
  summarise("Count"=length(PSs), 
            "GeneRatio"=length(PSs)/gene_count[unique(Species)]*100, 
            "m6ARatio"=sum(m6A!="")/length(Gene)*100, # m6ARatio
            "ExpRatio"=sum(Exp!="")/length(Gene)) %>% 
  as.data.frame() -> ratios.PS.spp

limit.PS.level <- intersect(PS.all.names, gene.attributes.df$PSs)
limit.PS.level <- limit.PS.level[!limit.PS.level%in%c("LS_OGs", "LS_OC")]
PS.trend <- list(limit.PS.level[c(1:7, 18:24)], # "Brassicaceae"
                 limit.PS.level[c(1:7, 18:22, 25)], # "Malvaceae"
                 limit.PS.level[c(1:7, 18:20, 26:28)], # "Phaseoleae"
                 limit.PS.level[c(1:7, 18, 29)], # "Solanaceae"
                 limit.PS.level[c(1:7, 8:10, 16:17)], # "Andropogoneae"
                 limit.PS.level[c(1:7, 8:15)] # "Triticum"
                 )

ratio.PSs.box <-  list()
color.list <- c("#E83828", "#CD6315", "#C5950F", "#3AA568", "#2A68B2", "#9C76B2")

for (ii in 1:length(PS.trend)){
  ratios.PS.df <- ratios.PS.spp[ratios.PS.spp$PSs%in%PS.trend[[ii]], ]
  ratios.PS.df$PSs <- factor(ratios.PS.df$PSs, 
                             levels = intersect(PS.all.names, ratios.PS.df$PSs))
  ratio.PSs.box[[ii]] <- ratios.PS.df %>% ggplot(aes(x = PSs, y = m6ARatio))+
    geom_boxplot(outlier.color = NA, fill=color.list[ii]) + 
    theme_bw(base_size = 12) + 
    theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust = 1)) +
    BaseTheme() + coord_cartesian(ylim = c(0,80)) +
    xlab("") + ylab(comlab.name)
}

pdf("figures/Figure2B.pdf", width = 20, height = 4)
ratio.PSs.box[[1]] + ratio.PSs.box[[2]] + ratio.PSs.box[[3]] + 
  ratio.PSs.box[[4]] + ratio.PSs.box[[5]] +
  ratio.PSs.box[[6]] + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position='none')
dev.off()


# Supplementary Fig. 7 ####
ratio.each.sp.df <- read.table(paste0(link.data.path, "figure_data/FigureS7.csv"), 
                               sep = ",", header = T, stringsAsFactors = F)
colnames(ratio.each.sp.df)[6] <- "species"
ratio.each.sp.df <- ratio.each.sp.df[ratio.each.sp.df$species!="ppa", ]
ratio.each.sp.df <- ratio.each.sp.df[ratio.each.sp.df$type!="P.patens", ] # Embryophyta
ratio.each.sp.df$ratio <- ratio.each.sp.df$m6ARatio*100 # ratio
ratio.each.sp.level <- ratio.each.sp.df %>% group_by(species, type) %>% 
  summarise("mean"=mean(ratio), 
            "sd"=sd(ratio), 
            "se"=sd(ratio)/sqrt(length(ratio))) %>% 
  as.data.frame()

ratio.each.sp.level$species <- FirstUp(ratio.each.sp.level$species)
ratio.each.sp.level$species <- factor(ratio.each.sp.level$species, 
                                      levels = FirstUp(species13.names[-13]))
colnames(ratio.each.sp.level)[2] <- "Type"

pdf("figures/FigureS7.pdf", width = 9, height = 3)
ggplot(data = ratio.each.sp.level, aes(x=species, y=mean, fill=Type))+
  geom_bar(stat='identity', position="dodge") +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, position = position_dodge(0.9)) +
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab(comlab.name) + 
  scale_fill_manual(values=c("#3F756D", "#B78D67"))
dev.off()

for(i in unique(ratio.each.sp.df$species)){
  cat(i, 
      t.test(ratio.each.sp.df[ratio.each.sp.df$species%in%i&ratio.each.sp.df$type%in%"dicot", "m6ARatio"], 
             ratio.each.sp.df[ratio.each.sp.df$species%in%i&ratio.each.sp.df$type%in%"monocot", "m6ARatio"])$p.value, "\n")
}


# Supplementary Fig. 9 ####
gene.PSs.count <- dcast(gene.attributes.df[, c("Species", "PSs", "Gene")], Species~PSs)
rownames(gene.PSs.count) <- gene.PSs.count[, 1]
gene.PSs.count <- gene.PSs.count[, -1]

gene.PSs.out <- t(gene.PSs.count[species13.names, 
                                       intersect(PS.all.names, 
                                                 colnames(gene.PSs.count))])
write.table(gene.PSs.out, "figures/FigureS9.gene_PSs_count.txt", quote = F, sep = "\t")
