## current wroking directory
setwd("~/a2z/m6A_13spp/")

## Options Settings
options(stringsAsFactors = F)
options(scipen=200)

## loading libraries
library(reshape2)
library(dplyr)
library(igraph)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggplotify)
library(rgl)

## reload data
load("notebooks/RData/01_m6A_peaks_genes.RData")
load("notebooks/RData/02_detailed_gene_attributes.RData")

## load PP Data 
m6A_translation_df <- read.table("associated_data/m6A_translation/m6A_translation.txt", 
                                 header = T, sep = "\t")
m6A_trans_df <- m6A_translation_df %>% group_by(Gene.ID) %>% 
  summarise("m6A"=max(m6A.level), "ts"=unique(Translational.status)) %>% 
  as.data.frame()
rownames(m6A_trans_df) <- m6A_trans_df$Gene.ID
rm(m6A_translation_df)
### GSM3536193	B_IP_1   GSM3536194	B_IP_2
### GSM3536197	B_CK_1   GSM3536198	B_CK_2
### GSM3536201	B_PP_1   GSM3536202	B_PP_2
RNA_df1 <- read.table("associated_data/m6A_translation/GSM3536197_B_CK_1.genes_FPKM.txt", sep = "\t")
RNA_df1$Type <- "RNA_rep1"
RNA_df2 <- read.table("associated_data/m6A_translation/GSM3536198_B_CK_2.genes_FPKM.txt", sep = "\t")
RNA_df2$Type <- "RNA_rep2"
TT_df1 <- read.table("associated_data/m6A_translation/GSM3536201_B_PP_1.genes_FPKM.txt", sep = "\t")
TT_df1$Type <- "Trans_rep1"
TT_df2 <- read.table("associated_data/m6A_translation/GSM3536202_B_PP_2.genes_FPKM.txt", sep = "\t")
TT_df2$Type <- "Trans_rep2"
RNA_df <- rbind(RNA_df1[,c(4,6,7)], RNA_df2[,c(4,6,7)], 
                TT_df1[,c(4,6,7)], TT_df2[,c(4,6,7)])
colnames(RNA_df) <- c("Gene", "exp", "Type")
gene_exp_trans <- dcast(RNA_df, Gene~Type, value.var = "exp")
rownames(gene_exp_trans) <- gene_exp_trans$Gene

## Gene number: 45578 --> 39005
zma_gene_names <- gene_attributes_df[gene_attributes_df$Species=="zma", "Gene"]
gene_exp_trans <- gene_exp_trans[gene_exp_trans$Gene%in%zma_gene_names, ]
gene_exp_trans$m6A <- m6A_trans_df[gene_exp_trans$Gene, "m6A"]
gene_exp_trans$ts <- m6A_trans_df[gene_exp_trans$Gene, "ts"]

gene_exp_trans$ts_value <- apply(gene_exp_trans, 1, function(x){
  if(max(as.numeric(x[2:5]))==0){
    return(0)
  }else if(mean(as.numeric(x[2:3]))==0){
    return((mean(as.numeric(x[4:5]))+1)/(mean(as.numeric(x[2:3]))+1))
  }else{
    return(mean(as.numeric(x[4:5]))/mean(as.numeric(x[2:3])))
  }
})

## Add attributes
### m6A non-m6A
m6A_type <- rep("Non-m6A", nrow(gene_exp_trans))
m6A_type[!is.na(gene_exp_trans$m6A)] <- "m6A"
gene_exp_trans$met_type <- m6A_type

# m6A_genes <- gene_exp_trans[gene_exp_trans$met_type=="m6A", "Gene"]
# exp_genes <- gene_exp_trans[(gene_exp_trans$RNA_rep1 + gene_exp_trans$RNA_rep2)/2>1, "Gene"]
# exp_values <- (gene_exp_trans$RNA_rep1 + gene_exp_trans$RNA_rep2)/2
# names(exp_values) <- rownames(gene_exp_trans)
# names_vec <- c("Zm00001d018386", "Zm00001d010190", "Zm00001d016198", "Zm00001d042685",
#       "Zm00001d043382")
# names_vec <- c("Zm00001d034851", "Zm00001d050133", "Zm00001d019988")
# names_vec <- c("Zm00001d013223")
# names_vec <- c("Zm00001d050133", "Zm00001d019988", "Zm00001d023278", "Zm00001d050134",
#                "Zm00001d005535")
# names_vec <- c("Zm00001d018555")
# names_vec <- c("Zm00001d017430", "Zm00001d037768", "Zm00001d043780", "Zm00001d029812",
#                "Zm00001d038904") 
# cat(sum(names_vec%in%rownames(gene_exp_trans)), sum(names_vec%in%exp_genes), sum(names_vec%in%m6A_genes)) #rownames(gene_exp_trans)
# cat(sum(exp_values[rownames(gene_exp_trans)%in%names_vec]), 
#     sum(exp_values[rownames(gene_exp_trans)%in%intersect(names_vec, exp_genes)]),
#     sum(exp_values[rownames(gene_exp_trans)%in%intersect(names_vec, m6A_genes)]))
# 
# df <- data.frame("auth_num"=c(exp_values[rownames(gene_exp_trans)%in%intersect(names_vec, exp_genes)],
#                               exp_values[rownames(gene_exp_trans)%in%intersect(names_vec, m6A_genes)]),
#                  "type"=rep(c("exp","m6A"), c(5,4)))
# ggplot(df, aes(x=auth_num, color=type)) + 
#   stat_ecdf(geom="smooth", se=F, size=1.2) + 
#   theme_bw() 
# 
# names_vec <- c("Zm00001d009927", "Zm00001d038146", "Zm00001d045417", "Zm00001d027731", 
#                "Zm00001d048424", "Zm00001d034005", "Zm00001d007890", "Zm00001d018896", 
#                "Zm00001d015194", "Zm00001d054001", "Zm00001d013962", "Zm00001d032600", 
#                "Zm00001d017566", "Zm00001d027671", "Zm00001d013829")
# which(names_vec%in%exp_genes)
# which(names_vec%in%m6A_genes)
# exp_values[intersect(names_vec, exp_genes)]
# exp_values[intersect(names_vec, m6A_genes)]


### duplicates singletons
sbi_zma <- read.table("associated_data/subgenome_clean/zma_sbi_pairs_17898.txt", 
                      sep = "\t", stringsAsFactors = F, header = T)
for(ii in 1:nrow(sbi_zma)){
  if(!sbi_zma[ii,2]%in%gene_exp_trans$Gene){
    sbi_zma[ii,2] <- "-"
  }
  if(!sbi_zma[ii,3]%in%gene_exp_trans$Gene){
    sbi_zma[ii,3] <- "-"
  }
}
zma_duplicates <- sbi_zma[sbi_zma$B73maize1_Gene!="-"&sbi_zma$B73maize2_Gene!="-", ]
zma1_dup <- zma_duplicates$B73maize1_Gene
zma2_dup <- zma_duplicates$B73maize2_Gene
zma1_sing <- sbi_zma[sbi_zma$B73maize2_Gene=="-", "B73maize1_Gene"]
zma2_sing <- sbi_zma[sbi_zma$B73maize1_Gene=="-", "B73maize2_Gene"]

ds_type <- rep("Others", nrow(gene_exp_trans))
# ds_type[gene_exp_trans$Gene%in%zma_duplicate_genes] <- "Duplicates"
# ds_type[gene_exp_trans$Gene%in%zma_singleton_genes] <- "Singleton"
ds_type[gene_exp_trans$Gene%in%zma1_dup] <- "zma1_Duplicates"
ds_type[gene_exp_trans$Gene%in%zma2_dup] <- "zma2_Duplicates"
ds_type[gene_exp_trans$Gene%in%zma1_sing] <- "zma1_Singletons"
ds_type[gene_exp_trans$Gene%in%zma2_sing] <- "zma2_Singletons"
table(ds_type)
gene_exp_trans$DS_type <- ds_type

## duplicate types
specific_gene <- gene_exp_trans[gene_exp_trans$DS_type%in%c("Others", "zma1_Singletons", "zma2_Singletons"), "Gene"]
gene_local_dup_list <- list()
for(dup in c("wgd", "tandem", "proximal", "dispersed", "transposed")){
  dup_tmp_df <- read.table(paste0("./associated_data/Duplication_Genes/", 
                                  "zma", "/pep_DupGen_results/", 
                                  "zma", ".", dup, ".pairs-unique"), 
                           sep = "\t", stringsAsFactors = F, header = T)
  dup_tmp_df$type <- dup
  cat(nrow(dup_tmp_df), length(unique(c(dup_tmp_df[,1], dup_tmp_df[,3]))), "\n")
  colnames(dup_tmp_df) <- c("Duplicate.1", "Location", "Duplicate.2", 
                            "Location.1", "E.value", "type")
  gene_local_dup_list[[dup]] <- dup_tmp_df
}
sapply(gene_local_dup_list, nrow)
gene_local_dup <- do.call(rbind, gene_local_dup_list[c("tandem", "proximal", "dispersed", "transposed")]) #
dim(gene_local_dup)
local_dup_index <- which(gene_local_dup[,1]%in%specific_gene&gene_local_dup[,3]%in%specific_gene)
gene_local_dup <- gene_local_dup[local_dup_index, ]
dim(gene_local_dup)

## duplicate types
merge_dup <- rep("Others", nrow(gene_exp_trans))
merge_dup[gene_exp_trans$DS_type%in%c("zma1_Duplicates", "zma2_Duplicates")] <- "wgd"
merge_dup[gene_exp_trans$Gene%in%c(gene_local_dup[,1], gene_local_dup[,3])] <- "Local"
gene_exp_trans$merge_dup <- merge_dup

gene_exp_trans$DS_mtype <- gsub("zma._", "", gene_exp_trans$DS_type)

## Filtered by expression values ####
gene_exp_trans_exp1 <- gene_exp_trans[(gene_exp_trans$RNA_rep1 + gene_exp_trans$RNA_rep2)/2>1, ]


## Figure 6A and 6B ####
all_genes <- gene_exp_trans$Gene
exp_genes <- gene_exp_trans_exp1$Gene
m6A_genes <- gene_exp_trans[!is.na(gene_exp_trans$m6A), "Gene"]
cat(length(all_genes), length(exp_genes), length(m6A_genes), "\n")

g=graph.data.frame(gene_local_dup[,c(1,3)], directed = F)
gene_graph <- cluster_infomap(g)

genes_list <- list()
length(unique(gene_graph$membership)) #2566
for(ii in unique(gene_graph$membership)){
  name_index <- gene_graph$membership%in%ii
  genes_list[[ii]] <- names(g[1])[name_index]
}

cluster_raw <- list()
# pk_cluster_raw <- list()
# tf_cluster_raw <- list()
# tf_all_site <- pk_all_site <- vector()
# ppi_cluster_raw <- list()

for(gene_index in 1:length(genes_list)){
  tmp_genes <- genes_list[[gene_index]]
  tmp_genes <- tmp_genes[tmp_genes%in%exp_genes]
  tmp_m6A_level <- sum(tmp_genes%in%m6A_genes)/sum(tmp_genes%in%exp_genes)*100 ## 2021-08-02 length(tmp_genes)
  m6A_value <- gene_exp_trans[tmp_genes, "m6A"]
  # m6A_value[is.na(m6A_value)] <- 0
  # m6A_mean <- mean(m6A_value)
  m6A_mean <- mean(m6A_value[!is.na(m6A_value)])
  ts_value <- gene_exp_trans[tmp_genes, "ts_value"]
  ts_value[is.na(ts_value)] <- 0
  ts_mean <- mean(ts_value) #, na.rm = T
  # ts_mean <- mean(ts_value[!is.na(ts_value)])
  cluster_raw[[paste0("c_", gene_index)]] <- data.frame(
    "trans_type"=c("m6A abundance", "Translational status", "m6A level(%)"),
    "abundance"=c(m6A_mean, ts_mean, tmp_m6A_level),
    "Type"=paste0("M_", sum(tmp_genes%in%exp_genes))) ## length(tmp_genes) 
  ## cluster in ppi
  # if(any(tmp_genes %in%names(zma_edge_num))){cat(gene_index, "\n")}
  # if(sum(tmp_genes%in%names(zma_edge_num))>0){
  #   zma_ppi_genes <- tmp_genes[tmp_genes%in%names(zma_edge_num)]
  #   ppi_m6A_level <- sum(zma_ppi_genes%in%m6A_genes)/length(zma_ppi_genes)*100
  #   edge_num <- zma_edge_num[zma_ppi_genes]
  #   edge_mean <- mean(edge_num)
  #   ppi_cluster_raw[[paste0("c_", gene_index)]] <- data.frame(
  #     "trans_type"=c("m6A level(%)", "edge number"),
  #     "abundance"=c(ppi_m6A_level, edge_mean),
  #     "Type"=paste0("M_", length(zma_ppi_genes)))
  # }
  # if(sum(gene_exp_trans[tmp_genes,"PK_index"]=="PK")>(length(tmp_genes)/2)){
  #   pk_all_site <- c(pk_all_site, gene_index)
  #   pk_cluster_raw[[paste0("c_",gene_index)]] <- data.frame(
  #     "trans_type"=c("m6A abundance", "Translational status", "m6A level(%)"),
  #     "abundance"=c(m6A_mean, ts_mean, tmp_m6A_level),
  #     "Type"=paste0("M_", length(tmp_genes)))
  # }
  # if(sum(gene_exp_trans[tmp_genes,"TF_index"]=="TF")>(length(tmp_genes)/2)){
  #   tf_all_site <- c(tf_all_site, gene_index)
  #   tf_cluster_raw[[paste0("c_",gene_index)]] <- data.frame(
  #     "trans_type"=c("m6A abundance", "Translational status", "m6A level(%)"),
  #     "abundance"=c(m6A_mean, ts_mean, tmp_m6A_level),
  #     "Type"=paste0("M_", length(tmp_genes)))
  # }
}

###
# cluster_range_raw <- list()
# range_vec <- c(100,40,30,20,15,10,5,0)
# for(range_index in 2:length(range_vec)){ #sum(tmp_cluster_num>4)
#   range_one <- range_vec[range_index]
#   range_two <- range_vec[range_index-1]
#   tmp_genes_list <- genes_list[sapply(genes_list, function(x)
#     {length(x)>range_one&length(x)<=range_two})]
#   tmp_genes <- unlist(tmp_genes_list)
#   tmp_m6A_level <- sum(tmp_genes%in%m6A_genes)/length(tmp_genes)*100
#   m6A_value <- gene_exp_trans[tmp_genes, "m6A"]
#   m6A_value[is.na(m6A_value)] <- 0
#   m6A_mean <- mean(m6A_value)
#   ts_value <- gene_exp_trans[tmp_genes, "ts_value"]
#   ts_value[is.na(ts_value)] <- 0
#   ts_mean <- mean(ts_value) #, na.rm = T
#   if(any(tmp_genes %in%names(zma_edge_num))){cat(range_one, "\n")}
#   cluster_range_raw[[paste0("c_", range_one)]] <- data.frame(
#     "trans_type"=c("Translational status", "m6A level(%)"),
#     "abundance"=c(-ts_mean, tmp_m6A_level),
#     "Type"=paste0("M_", range_one))
# }

final_cluster_raw <- cluster_raw # pk_cluster_raw # tf_cluster_raw ppi_cluster_raw cluster_raw
final_cluster_raw <- do.call(rbind, final_cluster_raw)
final_cluster_raw <- final_cluster_raw[!final_cluster_raw$Type%in%c("M_0", "M_1"), ]
final_cluster_raw$Type <- factor(final_cluster_raw$Type, 
                                 levels = unique(final_cluster_raw$Type))

tmp_type <- 1
tmp_types <- c("m6A level(%)", "Translational status", 
               "m6A abundance", "edge number")

# pdf(file = paste0("Figure_34_", tmp_type, ".pdf"), width = 5, height = 3, useDingbats = FALSE)
value_type_index <- final_cluster_raw$trans_type%in%tmp_types[tmp_type]
m6A_abundance_df <- final_cluster_raw[value_type_index, ]
m6A_abundance_df <- m6A_abundance_df[which(m6A_abundance_df$abundance>0), ]

bar_m6A_abundance_df <- m6A_abundance_df %>% group_by(Type) %>% 
  summarise("mean"=mean(abundance), "sd"=sd(abundance), 
            "se"=sd(abundance)/sqrt(length(abundance))) %>% 
  as.data.frame()
bar_m6A_abundance_df <- bar_m6A_abundance_df[!is.na(bar_m6A_abundance_df$sd), ]
bar_m6A_abundance_df$count <- as.numeric(sapply(strsplit(as.character(bar_m6A_abundance_df$Type), "_"), function(x){x[2]}))
bar_m6A_abundance_df <- bar_m6A_abundance_df[order(bar_m6A_abundance_df$count, decreasing = T), ]
bar_m6A_abundance_df$Type <- factor(bar_m6A_abundance_df$Type, levels = bar_m6A_abundance_df$Type)

bar_m6A_abundance_df %>% ggplot()+
  geom_bar(mapping = aes(x=Type, y=mean), stat='identity', position="dodge", fill="#A07D6F") + # A07D6F 55736F
  geom_errorbar(data =bar_m6A_abundance_df, 
                aes(x=Type, y=mean, ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(0.9)) +
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("") + ylab(tmp_types[tmp_type])
# geom_line(mapping = aes(x = Type, y = count/15, group=1), size = 0.7) + 
# geom_point(mapping = aes(x = Type, y = count/15), size = 2, shape = 21, fill = "#C18C1D") +
# scale_y_continuous(name = tmp_types[tmp_type], sec.axis = sec_axis(~ . * 1*15, name = "Number of family members"))
# dev.off()

## 
merged_data1 <- full_join(m6A_df, te_df,by = "count")
merged_data2 <- merged_data1
merged_data2[, "count"] <- 0
merged_data3 <- merged_data2
lm(merged_data2$mean.y ~ merged_data2$mean.x)
merged_data3$mean.y <- 0.08344+merged_data3$mean.x*0.01838
# cor.test(merged_data2$mean.x, merged_data2$mean.y)

plot3d(merged_data1[,c(2,7,5)], 
       xlab = "", ylab = "", zlab = "", axes=F,
       size =1.5, col="#B1886D",
       type = "s", lit=T)
par3d(lty=3)
plot3d(merged_data1[,c(2,7,5)], size =9, col="#B1886D", type = "h",add = T)
plot3d(merged_data2[,c(2,7,5)], 
       size =1, col="grey",type = "s", lit=F, add = T)
plot3d(merged_data3[,c(2,7,5)], 
       lwd =2, col="grey", type = "l", lty=3, add = T)
axes3d(edges = c("x--", "y+-", "z--"),
       ntick = 6,
       cex = .75, expand =6)
mtext3d("m6A methylation ratios", edge = "x--", line = 2)
mtext3d("Translational efficiency", edge = "y+-", line = 3)
mtext3d("Number of family member", edge = "z--", line = 3)
box3d()
rgl.postscript("00plot3d.ps", fmt = "ps")
rgl.snapshot( "00plot3d.png", fmt = "png")



## Figure 6C and 6D ####
dup_sin_df <- gene_exp_trans_exp1[gene_exp_trans_exp1$DS_mtype!="Others", ]

# pdf("Figure6A.pdf", width = 4, height = 5)
# pdf("Figure6B_6D.pdf", width = 6, height = 4)
ggplot(dup_sin_df, aes(x=DS_mtype, y=log2(ts_value+1),fill=DS_mtype)) +
  geom_violin() + geom_boxplot(width=0.3, fill="#7F7F7F", outlier.color = NA)+
  scale_fill_manual(values = c("#4B7972", "#B1886D"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 'wilcox.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational efficiency + 1)"))
# dev.off()

ggplot(dup_sin_df, 
       aes(x=met_type, y=log2(ts_value+1), fill=met_type)) +
  geom_violin() + geom_boxplot(width=0.3, fill="#7F7F7F", outlier.color = NA)+
  scale_fill_manual(values = c("#9B70A7", "#C88F17"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 'wilcox.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational efficiency + 1)")) + facet_grid(.~DS_mtype)
# dev.off()


## Others ####
tmp_gene_exp_trans <- gene_exp_trans_exp1 # gene_exp_trans
singleton_duplicates <- tmp_gene_exp_trans[tmp_gene_exp_trans$DS_mtype%in%c("Duplicates", "Singleton"), ]
ggplot(singleton_duplicates, #[!is.na(gene_exp_trans_exp1$m6A), ], 
       aes(x=DS_mtype, y=log2(ts_value+1), fill=PPI_index)) +
  geom_boxplot(outlier.color = "white")+ facet_grid(.~m6A_index) +
  # scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 'wilcox.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational status)")) + 
  coord_cartesian(ylim = c(0, 3))

ggplot(singleton_duplicates, #[!is.na(gene_exp_trans_exp1$m6A), ], 
       aes(x=DS_mtype, y=log2(ts_value+1), fill=PPI_index)) +
  geom_boxplot(outlier.color = "white")+ facet_grid(.~m6A_index) +
  # scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 'wilcox.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational status)")) + 
  coord_cartesian(ylim = c(0, 3))


ggplot(singleton_duplicates, #[!is.na(gene_exp_trans_exp1$m6A), ], 
       aes(x=PPI_index, y=log2(ts_value+1), fill=m6A_index)) +
  geom_boxplot(outlier.color = "white")+ facet_grid(.~DS_type) +
  # scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 'wilcox.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational status)")) + 
  coord_cartesian(ylim = c(0, 3))


ggplot(gene_exp_trans[!is.na(gene_exp_trans$m6A), ], 
       aes(x=PPI_index, y=m6A, fill=PPI_index)) +
  geom_boxplot(outlier.color = "black")+
  # scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 'wilcox.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("m6A level"))

cor.test(gene_exp_trans_exp1$m6A, log2(gene_exp_trans_exp1$ts_value+1))

inside_df <- gene_exp_trans_exp1[gene_exp_trans_exp1$PPI_index=="Inside", ]
sum(inside_df$m6A_index=="m6A")/nrow(inside_df)
inside_df_m6A <- inside_df[inside_df$m6A_index=="m6A", ]
cor.test(inside_df_m6A$m6A, log2(inside_df_m6A$ts_value+1))
plot(inside_df_m6A$m6A, log2(inside_df_m6A$ts_value+1))

outside_df <- gene_exp_trans_exp1[gene_exp_trans_exp1$PPI_index=="Outside", ]
sum(outside_df$m6A_index=="m6A")/nrow(outside_df)
outside_df_m6A <- outside_df[outside_df$m6A_index=="m6A", ]
cor.test(outside_df_m6A$m6A, log2(outside_df_m6A$ts_value+1))

chisq.test(matrix(c(3213, 5251-3213, 5052, 12535-5052), nrow = 2))

side_level <- data.frame(
  "Type"=rep(c("Inside", "Outside"), 2), 
  "value"=c(0.6118835, 0.4030315, 
            1-0.6118835, 1-0.4030315),
  "mtf"=rep(c("m6A", "non-m6A"), each=2))
side_level$Type <- factor(side_level$Type, levels = unique(side_level$Type))
side_level$mtf <- factor(side_level$mtf, levels = c("non-m6A", "m6A"))
# dev.off()

pdf("Figure36.pdf", height = 4, width = 3.5)
ggplot(data = side_level, aes(x=Type, y=value, fill=mtf))+
  geom_bar(stat='identity', position = "fill", width=0.9) +
  theme_bw(base_size = 12) + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(angle=30,hjust = 1,vjust = 1, 
                                   colour="black",size=10),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(suffix = "")) +
  xlab("") + ylab(expression("m"^"6"*"A level (%)")) + #coord_cartesian(ylim = c(50, 70)) +
  scale_fill_manual(values = c( "#CCCCCC", "#FF9500"))
dev.off()


sum(gene_exp_trans$PPI_index=="Inside"&gene_exp_trans$m6A_index=="m6A")/sum(gene_exp_trans$PPI_index=="Inside")
sum(gene_exp_trans$PPI_index=="Outside"&gene_exp_trans$m6A_index=="m6A")/sum(gene_exp_trans$PPI_index=="Outside")


ts_value_df1 <- sort(sample(1:nrow(gene_exp_trans_exp1), 5000))
ts_value_df2 <- sort(sample(1:nrow(gene_exp_trans_exp1), 5000))
length(intersect(ts_value_df1, ts_value_df2))



p1 <- ggplot(gene_exp_trans_exp1, aes(x=met_type, y=log2(ts_value+1))) +
  geom_boxplot(outlier.color = "white")+
  # scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme_bw(base_size = 14)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(method = 'wilcox.test', size=3, label.y = 3) +
  xlab("") + ylab(expression("Log"["2"]*"(Translational status)")) + 
  coord_cartesian(ylim = c(0, 3))


plot(density(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="Non-m6A", "ts_value"]))
lines(density(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="m6A", "ts_value"]), 
      col="red")

boxplot(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="m6A", "ts_value"], 
        gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="Non-m6A", "ts_value"])

## m6A 0.7755 Non-m6A 0.6795
var(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="Non-m6A", "ts_value"])
var(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="m6A", "ts_value"])
plot(density(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="Non-m6A", "ts_value"]), 
     lwd=3, 
     main = "",col="#D8C573", 
     xlab = "Translational efficiency", las=1)
# axis(1, pos=0)
# axis(2, pos=-0.5, las=1)
lines(density(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="m6A", "ts_value"]), 
      lwd=3, col="#D54472")
legend("topright",legend = c("m6A", "Non-m6A"),
       lty=1, col = c("#D54472", "#D8C573"), box.col = "white", lwd = 3)
text(6,0.4, labels = "m6A:0.7755;Non-m6A:0.6795")


sd(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="Non-m6A", "ts_value"])
sd(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="m6A", "ts_value"])

gene_exp_trans_exp1 %>% group_by(met_type) %>% summarise("sd"=sd(ts_value))

summary(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="Non-m6A", "ts_value"])
summary(gene_exp_trans_exp1[gene_exp_trans_exp1$met_type=="m6A", "ts_value"])