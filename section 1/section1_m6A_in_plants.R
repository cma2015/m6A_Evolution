# Transcriptome-wide mapping of m6A modifications for 13 plant species

## current wroking directory
setwd("E:/a2z_projects/03_m6A_13spp") #~/a2z/m6A_13spp/

## Options Settings
options(stringsAsFactors = F)
options(scipen=200)

## loading libraries
library(ggplot2)
library(patchwork)
library(ggpubr)

## reload data
load("notebooks/RData/01_m6A_peaks_genes.RData")
load("notebooks/RData/02_detailed_gene_attributes.RData")

comlab_name <- "m6A methylation ratio (%)"

## Table S4 summary table ####
summary_table <- as.data.frame(matrix(nrow = 13, ncol = 17))
rownames(summary_table) <- SPECIES13
colnames(summary_table) <- c("PCGs", "Exps",	"Exps_PCGs",
                             "m6APCGs", "m6A_ratio", "OGs", "OGs_PCGs",
                             "ExpOGs", "ExpOGs_OGs", "m6AOGs", "m6AOGs_ExpOGs",
                             "SGs",	"SGs_PCGs", "ExpSGs", "ExpSGs_SGs",
                             "m6ASGs", "m6ASGs_ExpSGs")

for(oness in SPECIES13){
  gene_names <- gene_attributes_df[gene_attributes_df$Species==oness, "Gene"]
  exp_genes <- exp_list[[oness]]
  m6A_genes <- m6A_list[[oness]]
  OLGs <- OGs_list_13spp[[oness]]
  LSGs <- setdiff(gene_names, OLGs)
  summary_table[oness, ] <- c(length(gene_names), length(exp_genes),
                              paste0(round(length(exp_genes)/length(gene_names)*100, 2), "%"),
                              length(m6A_genes),
                              paste0(round(length(m6A_genes)/length(exp_genes)*100, 2), "%"),
                              length(OLGs),
                              paste0(round(length(OLGs)/length(gene_names)*100, 2), "%"),
                              length(intersect(OLGs, exp_genes)),
                              paste0(round(length(intersect(OLGs, exp_genes))/length(OLGs)*100, 2), "%"),
                              length(intersect(OLGs, m6A_genes)),
                              paste0(round(length(intersect(OLGs, m6A_genes))/length(intersect(OLGs, exp_genes))*100, 2), "%"),
                              length(LSGs),
                              paste0(round(length(LSGs)/length(gene_names)*100, 2), "%"),
                              length(intersect(LSGs, exp_genes)),
                              paste0(round(length(intersect(LSGs, exp_genes))/length(LSGs)*100, 2), "%"),
                              length(intersect(LSGs, m6A_genes)),
                              paste0(round(length(intersect(LSGs, m6A_genes))/length(intersect(LSGs, exp_genes))*100, 2), "%"))
}
summary_table$Type <- c(rep("Dicots", 6), rep("Monocots",6), "ppa")

output_summary_table <- summary_table[, c(1:2,4:6,8,10:12,14,16:17)]

## Figure 1 ####
## Percentage of expressed genes and m6A genes
tmp_summary_df <- summary_table
tmp_summary_df$species <- factor(rownames(tmp_summary_df), levels = rev(SPECIES13))
tmp_summary_df$m6A_ratio <- as.numeric(gsub("%", "", tmp_summary_df$m6A_ratio))
plot1a <- ggplot(tmp_summary_df, aes(x=species, y=m6A_ratio, fill="A")) + 
  scale_fill_manual(values = "#8A6060") +
  geom_bar(stat='identity', position=position_dodge(width=0.5)) +
  theme_bw(base_size = 12)+ 
  theme(axis.text = element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("") + ylab("m6A methylation ratio (%)") + coord_flip ()

# pdf("Figure3.pdf", width = 3, height = 7.5)
plot1a
# dev.off()



#### m6A genes and background
gg_table <- read.table("associated_data/figure_data/Figure4.csv", sep = ",", 
                         header = T, stringsAsFactors = F, row.names = 1)

tmp_summary_df$GenomeSize <- gg_table[as.character(tmp_summary_df$species), "GenomeSize"]
tmp_summary_df$GeneNum <- gg_table[as.character(tmp_summary_df$species), "GeneNum"]

xycor_plot <- function(df_input, var1, var2, text_pos, xlab, ylab){
  df_input[, var1] <- as.numeric(df_input[, var1])
  df_input[, var2] <- as.numeric(df_input[, var2])
  ##
  print(cor(df_input[, var1], df_input[, var2]))
  
  fit_stat <- summary(lm(df_input[, var1] ~ df_input[, var2]))
  R2_adj <- fit_stat$adj.r.squared
  p_value <- fit_stat$coefficients[2,4]

  out_plot <- ggplot(df_input, aes_string(x=var1, y=var2)) +
    geom_point() +
    labs( x= xlab, y= ylab ) +
    theme_bw(base_size = 12)+ 
    theme(axis.line = element_line(colour = "black",size=0.5),
          axis.text = element_text(colour="black"),
          legend.text=element_text(colour="black"),
          legend.title=element_text(colour="black")) + 
    # geom_smooth(method="lm", se=T, color="#DA2416") +
    geom_smooth(method="lm", se=F, color="#A3A2A2") +
    annotate("text", label = sprintf('Pearson\'s r = %.2f', cor(df_input[, var1], df_input[, var2])), 
             x = text_pos[1], y = text_pos[2], size = 4) +
    annotate("text", label = sprintf('P-value = %.2g', p_value), 
             x = text_pos[3], y = text_pos[4], size = 4)
  out_plot
}

xycor_plot(df_input = tmp_summary_df, var1 ="GenomeSize", var2 = "m6APCGs_Exps", 
           text_pos=c(5000,30,5000,27), xlab = "Genome size (MB)", ylab = comlab_name )

xycor_plot(df_input = tmp_summary_df, var1 ="GeneNum", var2 = "m6APCGs_Exps", 
           text_pos=c(50000,30,50000,27), xlab = "Gene number", ylab = comlab_name )

xycor_plot(df_input = tmp_summary_df, var1 ="Exps", var2 = "m6APCGs_Exps", 
           text_pos=c(30000,30,30000,27), xlab = "Gene number", ylab = comlab_name )


# pdf("Figure4.pdf", width = 5, height = 8)
plot1 + plot2 + plot_layout(ncol = 1, guides = "collect") 
# dev.off()


### m6A levels and YTH-domain proteins
level_domain_df <- data.frame("spp"=rownames(summary_table),
                              "level"=sapply(m6A_list$m6A, function(x){sum(x$Ratio)})[rownames(summary_table)],
                              "ythpro"=c(13, 13, 26, 11, 19, 9, 12, 24, 13, 26, 39, 12, 4),
                              "writerpro"=c(6, 7, 14, 6, 15, 6, 6, 9, 7, 13, 20, 9, 7))
xycor_plot(df_input = level_domain_df, var1 ="level", var2 = "ythpro", 
           text_pos=c(60000, 10, 60000, 8), xlab = "m6A level", ylab = "YTH_protein" )

xycor_plot(df_input = level_domain_df, var1 ="level", var2 = "writerpro", 
           text_pos=c(60000, 10, 60000, 8), xlab = "m6A level", ylab = "YTH_protein" )


### Figure 5
summary_table_df <- data.frame("Type"=rep(c("OGs", "SGs"), each=13), 
                               "Values"=as.numeric(gsub("%", "", c(tmp_summary_df$m6AOGs_ExpOGs, 
                                          summary_table$m6ASGs_ExpSGs))))

ggplot(summary_table_df, aes(x=Type, y=Values, fill=Type)) +
  geom_violin() +
  geom_boxplot(width=0.2, fill="#7F7F7F")+
  scale_fill_manual(values = c("#3D7970", "#BB8F68"))+
  # geom_point(position=position_jitterdodge(jitter.width = 0.1)) +
  theme_bw(base_size = 12)+
  theme(axis.text= element_text(colour="black"),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black"),
        legend.title=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none")+
  stat_compare_means(method = 'wilcox.test',label.y = 83, size=2) +
  xlab("") + ylab(comlab_name)

# pdf("Figure5.pdf", width = 3, height = 4)
# ggplot(summary_table_df, aes(x=Type, y=Values, fill=Type)) +
#   geom_boxplot(outlier.color = "white")+
#   scale_fill_manual(values = c("#0072B2", "#E69F00"))+
#   geom_point(position=position_jitterdodge(jitter.width = 0.1)) +
#   theme_bw(base_size = 12)+
#   theme(axis.text= element_text(colour="black"),
#         axis.line = element_line(colour = "black",size=0.5),
#         legend.text=element_text(colour="black"),
#         legend.title=element_text(colour="black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), legend.position = "none")+
#   stat_compare_means(method = 'wilcox.test',label.y = 83, size=2) +
#   xlab("") + ylab(expression("m"^"6"*"A levels (%)"))
# dev.off()


## Supplementary Fig. 1 Pie and enrichment plot #### 
data_input <- read.table("associated_data/figure_data/Figure2.csv", sep = ",", 
                         stringsAsFactors = F, header = T)
data_input$position <- factor(data_input$position, 
                              levels = c("5UTR", "star_codon", "CDS", "stop_codon", "3UTR"))
data_input$spe <- factor(data_input$spe, levels = SPECIES13)
plot1 <- ggplot(data = data_input, mapping = aes(x = "Count", y = peak_num, fill = position)) +
  scale_fill_manual(values = c('#BDD1C5','#4D6372','#E7D86E','#558ED5',"#A984A8"), guide = FALSE) + #scale_fill_brewer(palette = "Set2")+
  geom_bar(stat = 'identity', position = 'fill', width = 0.5) +
  labs(x = '', y = '', title = '') + 
  coord_polar(theta = 'y', direction = -1) + theme_void(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) + 
  facet_grid(spe~., )

plot2 <- ggplot(data_input, aes(x = position, y = enrichment, fill=position)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c('#BDD1C5','#4D6372','#E7D86E','#558ED5',"#A984A8"))+
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid =element_blank()) +
  facet_grid(spe~., ) + labs(x="", y="Relative enrichment")

pdf("Figure2 relative enrichment of m6A peaks.pdf", width = 6, height = 15)
plot1 + plot2 + plot_layout(ncol = 2, guides = "collect") 
dev.off()


library(UpSetR)

merge_genes <- unique(c(m6A_list$ath_own, m6A_list$athCrBgd1, m6A_list$athCrBgd2))
ath_mat <- as.data.frame(matrix(0, nrow = length(merge_genes), ncol = 3))
rownames(ath_mat) <- merge_genes
colnames(ath_mat) <- c("This study", "CrBgd1", "CrBgd2")
ath_mat[m6A_list$ath_own, 1] <- 1
ath_mat[m6A_list$athCrBgd1, 2] <- 1
ath_mat[m6A_list$athCrBgd2, 3] <- 1
pdf("Fig. Ath_overlap_upset.pdf", height = 4, width = 5)
upset(ath_mat, sets = c("CrBgd2", "CrBgd1", "This study"), 
      order.by = "freq", empty.intersections = "on", keep.order=T)
dev.off()
