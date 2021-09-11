# Effects of m6A modification on translational efficiency in genomic duplications
#

# current wroking directory
rm(list=ls())
setwd("~/a2z/m6A_13spp/Part_6/")
options(stringsAsFactors = F)

## Options Settings
options(stringsAsFactors = F)

## loading libraries
library(reshape2)
library(dplyr)
library(igraph)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggplotify)
library(rgl)

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
  theme(axis.line = element_line(colour = "black",size=0.5),
        axis.text = element_text(colour="black"),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        axis.ticks = element_line(colour = "black",size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}

BarFig <- function(tmp.input.df){
  tmp.input.df %>% ggplot()+
    geom_bar(mapping = aes(x=Type, y=value, fill=Type), 
             stat='identity', position="dodge") +
    theme_bw(base_size = 12) +
    theme(axis.line = element_line(colour = "black", size = 0.5),
          axis.text = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("") + ylab("") + scale_fill_manual(values = c("#557571", "#A37F70"))
}

BarErrorFig <- function(tmp.input.df){
  tmp.input.df %>% ggplot()+
    geom_bar(mapping = aes(x=Type, y=value, fill=Type), 
             stat='identity', position="dodge") +
    geom_errorbar(aes(x=Type, y=value, ymin = value - se, ymax = value + se), 
                  width = 0.2, position = position_dodge(0.9)) +
    theme_bw(base_size = 12) +
    theme(axis.line = element_line(colour = "black", size = 0.5),
          axis.text = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("") + ylab("") + 
    scale_fill_manual(values = c("#557571", "#A37F70"))
}


## load PP Data 
m6A.translation.df <- read.table(paste0(link.data.path, 
                                        "m6A_translation/m6A_translation.txt"), 
                                 header = T, sep = "\t")
m6A.trans.df <- m6A.translation.df %>% group_by(Gene.ID) %>% 
  summarise("m6A"=max(m6A.level)) %>% 
  as.data.frame()
rownames(m6A.trans.df) <- m6A.trans.df$Gene.ID
rm(m6A.translation.df)

### GSM3536193	B_IP_1   GSM3536194	B_IP_2
### GSM3536197	B_CK_1   GSM3536198	B_CK_2
### GSM3536201	B_PP_1   GSM3536202	B_PP_2
RNA.df1 <- read.table(paste0(link.data.path, "m6A_translation/GSM3536197_B_CK_1.genes_FPKM.txt"), 
                      sep = "\t")
cat(sum(RNA.df1$V6), "\n")
RNA.df1$TPM <- RNA.df1$V6/sum(RNA.df1$V6)*10^6
RNA.df1$Type <- "RNA_rep1"
RNA.df2 <- read.table(paste0(link.data.path, "m6A_translation/GSM3536198_B_CK_2.genes_FPKM.txt"), 
                      sep = "\t")
cat(sum(RNA.df2$V6), "\n")
RNA.df2$TPM <- RNA.df2$V6/sum(RNA.df2$V6)*10^6
RNA.df2$Type <- "RNA_rep2"
TT.df1 <- read.table(paste0(link.data.path, "m6A_translation/GSM3536201_B_PP_1.genes_FPKM.txt"), 
                     sep = "\t")
cat(sum(TT.df1$V6), "\n")
TT.df1$TPM <- TT.df1$V6/sum(TT.df1$V6)*10^6
TT.df1$Type <- "Trans_rep1"
TT.df2 <- read.table(paste0(link.data.path, "m6A_translation/GSM3536202_B_PP_2.genes_FPKM.txt"), 
                     sep = "\t")
TT.df2$TPM <- TT.df2$V6/sum(TT.df2$V6)*10^6
TT.df2$Type <- "Trans_rep2"

RNA.df <- rbind(RNA.df1[, c(4, 7, 8)], RNA.df2[, c(4, 7, 8)], 
                TT.df1[, c(4, 7, 8)], TT.df2[, c(4, 7, 8)])
colnames(RNA.df) <- c("Gene", "exp", "Type")
exp.trans.df <- dcast(RNA.df, Gene ~ Type, value.var = "exp")
rownames(exp.trans.df) <- exp.trans.df$Gene

## Gene number: 45578 --> 39005
zma.gene.names <- gene.attributes.df[gene.attributes.df$Species == "zma", "Gene"]
exp.trans.df <- exp.trans.df[exp.trans.df$Gene%in%zma.gene.names, ]
exp.trans.df$m6A <- m6A.trans.df[exp.trans.df$Gene, "m6A"]

exp.trans.df$ts.value <- apply(exp.trans.df, 1, function(x){
  if ( max(as.numeric(x[2:5])) == 0 ){
    return(0)
  } else if (mean(as.numeric(x[2:3]))==0){
    return((mean(as.numeric(x[4:5]))+1)/(mean(as.numeric(x[2:3]))+1))
  } else {
    # return((mean(as.numeric(x[4:5]))+1)/(mean(as.numeric(x[2:3]))+1))
    return((mean(as.numeric(x[4:5])))/(mean(as.numeric(x[2:3]))))
  }
})

# Add attributes
# m6A non-m6A
m6A.type <- rep("non-m6A", nrow(exp.trans.df))
m6A.type[!is.na(exp.trans.df$m6A)] <- "m6A"
exp.trans.df$met.type <- m6A.type

### duplicates singletons
sbi.zma.ortho <- read.table(paste0(link.data.path, 
                                   "subgenome_clean/zma_sbi_pairs_17898.txt"),
                      sep = "\t", stringsAsFactors = F, header = T)
for (ii in 1:nrow(sbi.zma.ortho)){
  if(!sbi.zma.ortho[ii, 2]%in%exp.trans.df$Gene){
    sbi.zma.ortho[ii, 2] <- "-"
  }
  if(!sbi.zma.ortho[ii,3]%in%exp.trans.df$Gene){
    sbi.zma.ortho[ii, 3] <- "-"
  }
}
zma_duplicates <- sbi.zma.ortho[sbi.zma.ortho$B73maize1_Gene != "-" & sbi.zma.ortho$B73maize2_Gene != "-", ]
zma1_dup <- zma_duplicates$B73maize1_Gene
zma2_dup <- zma_duplicates$B73maize2_Gene
zma1_sing <- sbi.zma.ortho[sbi.zma.ortho$B73maize2_Gene=="-", "B73maize1_Gene"]
zma2_sing <- sbi.zma.ortho[sbi.zma.ortho$B73maize1_Gene=="-", "B73maize2_Gene"]

ds.type <- rep("Others", nrow(exp.trans.df))
ds.type[exp.trans.df$Gene%in%zma1_dup] <- "zma1_Duplicates"
ds.type[exp.trans.df$Gene%in%zma2_dup] <- "zma2_Duplicates"
ds.type[exp.trans.df$Gene%in%zma1_sing] <- "zma1_Singletons"
ds.type[exp.trans.df$Gene%in%zma2_sing] <- "zma2_Singletons"
table(ds.type)
exp.trans.df$DS.type <- ds.type

# duplicate types
remove.wgd.type <- c("Others", "zma1_Singletons", "zma2_Singletons")
specific.genes <- exp.trans.df[exp.trans.df$DS.type%in%remove.wgd.type, "Gene"]

local.dup.gene.list <- list()
for(dup in c("wgd", "tandem", "proximal", "dispersed", "transposed")){
  tmp.dup.df <- read.table(paste0(link.data.path, "Duplication_Genes/", 
                                  "zma", "/pep_DupGen_results/", 
                                  "zma", ".", dup, ".pairs-unique"), 
                           sep = "\t", stringsAsFactors = F, header = T)
  tmp.dup.df$type <- dup
  cat(nrow(tmp.dup.df), length(unique(c(tmp.dup.df[,1], tmp.dup.df[,3]))), "\n")
  colnames(tmp.dup.df) <- c("Duplicate.1", "Location", "Duplicate.2", 
                            "Location.1", "E.value", "type")
  local.dup.gene.list[[dup]] <- tmp.dup.df
}
sapply(local.dup.gene.list, nrow)
local.dup.df <- do.call(rbind, local.dup.gene.list[c("tandem", "proximal", "dispersed")]) # "transposed"
dim(local.dup.df)
local.dup.index <- which(local.dup.df[, 1]%in%specific.genes&local.dup.df[, 3]%in%specific.genes)
local.dup.df <- local.dup.df[local.dup.index, ]
dim(local.dup.df)

# duplicate types
merge.dup <- rep("Others", nrow(exp.trans.df))
merge.dup[exp.trans.df$DS.type%in%c("zma1_Duplicates", "zma2_Duplicates")] <- "wgd"
merge.dup[exp.trans.df$Gene%in%c(local.dup.df[,1], local.dup.df[,3])] <- "Local"
exp.trans.df$merge.dup <- merge.dup

exp.trans.df$DS.mtype <- gsub("zma._", "", exp.trans.df$DS.type)
exp.trans.df$expression <- (exp.trans.df$RNA_rep1 + exp.trans.df$RNA_rep2)/2

# Filtered by expression values ####
exp.trans.TPM1.df <- exp.trans.df[exp.trans.df$expression >= 1, ]

## Figure 6A and 6B ####
m6A.genes <- exp.trans.df[!is.na(exp.trans.df$m6A), "Gene"]
cat(length(exp.trans.df$Gene), length(m6A.genes), "\n")

g <- graph.data.frame(local.dup.df[,c(1,3)], directed = F)
gene.graph <- cluster_infomap(g)

loc.dup.list <- list()
length(unique(gene.graph$membership)) #2566
for(ii in unique(gene.graph$membership)){
  name_index <- gene.graph$membership%in%ii
  loc.dup.list[[ii]] <- names(g[1])[name_index]
}

# log file
family.attr <- matrix(0, nrow = length(loc.dup.list), ncol = 7)
class(family.attr) <- "numeric"
family.attr <- as.data.frame(family.attr)
cluster_raw <- list()

for(gene_index in 1:length(loc.dup.list)){
  tmp_genes <- loc.dup.list[[gene_index]]
  m6A_level <- sum(tmp_genes%in%m6A.genes)/length(tmp_genes)*100 # 2021-08-01
  ts.value <- exp.trans.df[tmp_genes, "ts.value"]
  ts.value[is.na(ts.value)] <- 0
  ts_mean <- mean(ts.value)
  exp_value <- exp.trans.df[tmp_genes, "expression"]
  exp_value[is.na(exp_value)] <- 0
  exp_mean <- mean(exp_value)
  cluster_raw[[paste0("c_", gene_index)]] <- data.frame(
    "trans_type"=c("gene expression", "Translational status", "m6A level(%)"),
    "abundance"=c(exp_mean, ts_mean, m6A_level),
    "Type"=paste0("M_", length(tmp_genes)))
  family.attr[gene_index, ] <- c(length(tmp_genes), sum(tmp_genes%in%m6A.genes),
                                 sum(tmp_genes%in%m6A.genes)/length(tmp_genes)*100,  # 2021-08-01
                                 exp_mean, ts_mean, sum(exp_value), sum(ts.value))
}

final.clusters.df <- do.call(rbind, cluster_raw)
# final.clusters.df <- final.clusters.df[!final.clusters.df$Type%in%c("M_0", "M_1"), ]

output.plots <- list()
for (tmp.type in 1:3){
  tmp.types <- c("m6A level(%)", "Translational status", "gene expression")
  # pdf(file = paste0("Figure_34_", tmp.type, ".pdf"), width = 5, height = 3, useDingbats = FALSE)
  vtype.index <- final.clusters.df$trans_type%in%tmp.types[tmp.type]
  bardata.df <- final.clusters.df[vtype.index, ]
  bardata.df <- bardata.df[which(bardata.df$abundance > 0), ]
  # errorbar plot
  ebardata.df <- bardata.df %>% group_by(Type) %>% 
    summarise("mean"=mean(abundance), 
              "se"=sd(abundance)/sqrt(length(abundance))) %>% 
    as.data.frame()
  ebardata.df <- ebardata.df[!is.na(ebardata.df$se), ]
  ebardata.df$count <- as.numeric(sapply(strsplit(as.character(ebardata.df$Type), "_"), function(x){x[2]}))
  ebardata.df <- ebardata.df[order(ebardata.df$count, decreasing = T), ]
  ebardata.df$Type <- factor(ebardata.df$Type, levels = ebardata.df$Type)
  if (tmp.type == 1){m6A.df <- ebardata.df} else if (tmp.type == 2){te.df <- ebardata.df}
  output.plots[[tmp.type]] <- ebardata.df %>% ggplot()+
    geom_bar(mapping = aes(x=Type, y=mean), stat='identity', 
             position="dodge", fill="#A07D6F") + # A07D6F 55736F
    geom_errorbar(data = ebardata.df, 
                  aes(x=Type, y=mean, ymin = mean - se, ymax = mean + se), 
                  width = 0.2, position = position_dodge(0.9)) +
    theme_bw(base_size = 12) +
    theme(axis.ticks.x = element_blank()) +
    BaseTheme() +
    xlab("") + ylab(tmp.types[tmp.type])
}

pdf("figures/Figure6AB_S22.pdf", width = 5, height = 10)
patchwork::wrap_plots(output.plots, ncol = 1)
dev.off()


## Figure 6C and 6D ####
dup_sin_df <- exp.trans.TPM1.df[exp.trans.TPM1.df$DS.mtype!="Others", ]

pdf("figures/Figure6C.pdf", width = 4, height = 4)
ggplot(dup_sin_df, aes(x=DS.mtype, y=log2(ts.value+1),fill=DS.mtype)) +
  geom_violin() + geom_boxplot(width=0.3, fill="#7F7F7F", outlier.color = NA)+
  scale_fill_manual(values = c("#4B7972", "#B1886D"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 't.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational efficiency + 1)"))
dev.off()

pdf("figures/Figure6D.pdf", width = 4, height = 4)
ggplot(dup_sin_df, 
       aes(x=met.type, y=log2(ts.value+1), fill=met.type)) +
  geom_violin() + geom_boxplot(width=0.3, fill="#7F7F7F", outlier.color = NA)+
  scale_fill_manual(values = c("#9B70A7", "#C88F17"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 't.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational efficiency + 1)")) + facet_grid(.~DS.mtype)
dev.off()

# Supplementary figure 21 ####
ts.plot <- ggplot(exp.trans.TPM1.df, aes(x=met.type, y=ts.value, fill=met.type)) +
  geom_boxplot(outlier.color = NA) +
  scale_fill_manual(values = c("#557571", "#A37F70")) +
  theme_bw(base_size = 14) +
  BaseTheme() +
  stat_compare_means(method = 't.test', size=3, label.y = 5) +
  xlab("") + ylab("Translational efficiency") + coord_cartesian(ylim = c(0, 5))

expression.plot <- ggplot(exp.trans.TPM1.df, aes(x=met.type, y=log2(expression+1), fill=met.type)) + 
  geom_boxplot(outlier.color = NA) +
  scale_fill_manual(values = c("#557571", "#A37F70")) +
  theme_bw(base_size = 14) +
  BaseTheme() +
  stat_compare_means(method = 't.test', size=3, label.y = 7) +
  xlab("") + ylab("Expression level") + coord_cartesian(ylim = c(1, 6.7))

pdf( "figures/FigureS21.pdf", width = 6, height = 5  )
ts.plot + expression.plot + plot_layout(guides = "collect")
dev.off()


# Supplementary figure 22 ####
pdf( "figures/FigureS22.pdf", width = 5, height = 4)
ggplot(dup_sin_df, aes(x=met.type, y=log2(expression+1), fill=met.type)) +  # ts.value
  geom_violin(trim = T) + geom_boxplot(width=0.3, fill="#7F7F7F", outlier.color = NA)+
  scale_fill_manual(values = c("#9B70A7", "#C88F17")) +
  theme_bw(base_size = 14) +
  BaseTheme() +
  stat_compare_means(method = 't.test', size=3, label.y = 8) +
  xlab("") + ylab("log2(expression+1)") + 
  facet_grid(.~DS.mtype) +
  coord_cartesian(ylim = c(1, 8))
dev.off()

# Supplementary figure 23 example ####
small.df <- exp.trans.df[loc.dup.list[[564]], ]
small.df[order(small.df$Gene), ]
small.m6A.exp <- sum(small.df[small.df$met.type=="m6A", "expression"])/sum(small.df$met.type=="m6A")
small.nonm6A.exp <- sum(small.df[small.df$met.type=="non-m6A", "expression"])/sum(small.df$met.type=="non-m6A")

small.m6A.ts <- sum(small.df[small.df$met.type=="m6A", "ts.value"])/sum(small.df$met.type=="m6A")
small.nonm6A.ts <- sum(small.df[small.df$met.type=="non-m6A", "ts.value"])/sum(small.df$met.type=="non-m6A")

large.df <- exp.trans.df[loc.dup.list[[111]], ] #124
large.df[order(large.df$Gene), ]
large.m6A.exp <- sum(large.df[large.df$met.type=="m6A", "expression"])/sum(large.df$met.type=="m6A")
large.nonm6A.exp <- sum(large.df[large.df$met.type=="non-m6A", "expression"])/sum(large.df$met.type=="non-m6A")

large.m6A.ts <- sum(large.df[large.df$met.type=="m6A", "ts.value"])/sum(large.df$met.type=="m6A")
large.nonm6A.ts <- sum(large.df[large.df$met.type=="non-m6A", "ts.value"])/sum(large.df$met.type=="non-m6A")

plot1 <- BarFig(tmp.input.df = data.frame("value" = c(sum(small.df$expression), 
                                                      sum(large.df$expression)),
                                          "Type" = c("few", "many")))
plot2 <- BarFig(tmp.input.df = data.frame("value" = c(sum(small.df$ts.value), 
                                                      sum(large.df$ts.value)),
                                          "Type" = c("few", "many")))
plot3 <- BarFig(tmp.input.df = data.frame("value" = c(3/3*100, 3/7*100), 
                                          "Type" = c("few", "many")))
plot4 <- BarErrorFig(tmp.input.df = data.frame("value" = c(mean(small.df$expression), 
                                                           mean(large.df$expression)),
                                               "se" = c(sd(small.df$expression)/sqrt(3), 
                                                        sd(large.df$expression)/sqrt(7)), 
                                               "Type" = c("few", "many")))
plot5 <- BarErrorFig(tmp.input.df = data.frame("value" = c(mean(small.df$ts.value), 
                                                           mean(large.df$ts.value)),
                                               "se" = c(sd(small.df$ts.value)/sqrt(3), 
                                                        sd(large.df$ts.value)/sqrt(7)),
                                               "Type" = c("few", "many")))

pdf("figures/FigureS23 example.pdf", width = 8, height = 2.5)
plot1  + plot2 + plot3  + plot4 + plot5 + plot_layout(guides = "collect", nrow = 1)
dev.off()

# tmp_compare_df <- data.frame("value" = c(small.m6A.exp, small.nonm6A.exp, large.m6A.exp, large.nonm6A.exp),
#                              "Type" = c("few_m6A", "few_nonm6A", "many_m6A", "many_nonm6A"))
# tmp_compare_df <- data.frame("value" = c(small.m6A.ts, small.nonm6A.ts, large.m6A.ts, large.nonm6A.ts), #mean(small.df$ts.value), mean(large.df$ts.value)),
#                              "Type" = c("few_m6A", "few_nonm6A", "many_m6A", "many_nonm6A"))


# Supplementary figure 24 rgl 3D PCA #### 
merged.data1 <- full_join(m6A.df, te.df,by = "count")
merged.data1 <- merged.data1[!is.na(merged.data1$Type.x), ]
merged.data2 <- merged.data1
merged.data2[, "count"] <- 0
merged.data3 <- merged.data2
lm(merged.data2$mean.y ~ merged.data2$mean.x)
merged.data3$mean.y <- 0.10798 + merged.data3$mean.x * 0.02373
cor.test(merged.data2$mean.x, merged.data2$mean.y)

plot3d(merged.data1[,c(2,6,4)], 
       xlab = "", ylab = "", zlab = "", axes=F,
       size =1.5, col="#B1886D",
       type = "s", lit=T)
par3d(lty=3)
plot3d(merged.data1[,c(2,6,4)], size =9, col="#B1886D", type = "h",add = T)
plot3d(merged.data2[,c(2,6,4)], size =1, col="grey", type = "s", lit=F, add = T)
plot3d(merged.data3[,c(2,6,4)],
       lwd =2, col="grey", type = "l", lty=3, add = T)
axes3d(edges = c("x--", "y+-", "z--"),
       ntick = 6,
       cex = .75, expand =6)
mtext3d("m6A methylation ratios", edge = "x--", line = 2)
mtext3d("Translational efficiency", edge = "y+-", line = 3)
mtext3d("Number of family member", edge = "z--", line = 3)
box3d()
rgl.snapshot( "figures/FigureS24.png", fmt = "png")
