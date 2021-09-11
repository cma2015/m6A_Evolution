# GO enrichment analysis
#

# libraries
library(clusterProfiler)
library(GO.db)
GO.name.list <- as.list(GOTERM)

# reload data
load("data1_m6A_peaks_genes.RData")
load("data2_gene_attributes.RData")

link.data.path <- "../associated_data/"

# m6A OGs
empty.count.index <-  apply(OGs.13spp, 1, function(x){sum(x=="")})
all.orthos <- OGs.13spp[empty.count.index==0, 1:13]
all.m6A.orthos <- t(apply(all.orthos, 1, function(x){
  for(ii in 1:13){
    tmp.genes <- strsplit(x[ii], split=", ")[[1]]
    tmp.genes <- tmp.genes[tmp.genes%in%m6A.list[[colnames(all.orthos)[ii]]]]
    x[ii] <- paste0(tmp.genes, collapse = ", ")
  }
  x
}))
m6A.empty.count <-  apply(all.m6A.orthos, 1, function(x){sum(x=="")})
m6A.orthos.highc <- all.m6A.orthos[m6A.empty.count==0, ]
m6A.orthos.lowc <- all.m6A.orthos[m6A.empty.count==12, ]

# GO results
go.results <- list()
go.results[["hc"]] <- go.results[["lc"]] <- list()

GObarPlot <- function(ttdf, ttsp){
  colnames(ttdf)[2] <- "Term"
  colnames(ttdf)[7] <- "elimFisher"
  if(nrow(ttdf)>20){
    ttdf <- ttdf[1:20,]
  }
  ttdf$Term <- factor(ttdf$Term, levels = rev(unique(ttdf$Term)))
  ggplot(ttdf, aes(x=Term, y=-log10(as.numeric(elimFisher)))) +
    stat_summary(geom = "bar", position = "dodge") +
    xlab("") +
    ylab("-log10(FDR)") +
    labs(title = ttsp) +
    theme_bw(base_size=12) +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
      axis.title=element_text(size=12, face="bold"),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=12),  #Text size
      title=element_text(size=12)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
}


go.output.plots <- list()
for (ii in c("highc", "lowc")){
  for ( spp in c("ath", "ghi", "zma", "ata", "tdi", "tae", "ppa")){
    cat(ii, spp, "\n")
    out_hc_genes <- unique(unlist(strsplit(eval(parse(text = paste0("m6A.orthos.", ii)))[, spp], ", ")))
    if (length(out_hc_genes)>0){
      GO.BG <- read.table(paste0(link.data.path, "GO_terms/", spp, ".txt"), sep = "\t", stringsAsFactors = F)
      GO.BG <- GO.BG[, 2:1]
      colnames(GO.BG) <- c("Term", "Gene")
      GO.name.BG <- GO.BG
      sum(GO.name.BG[,1]%in%names(GO.name.list))/length(GO.name.BG[,1])
      GO.name.vec <- sapply(GO.name.list[GO.name.BG$Term], function(x){
        if (length(x)==1){return(x@Term)}else{return(NA)}
      })
      GO.name.BG[,2] <- GO.name.vec
      if(length(out_hc_genes) <= 10){
        next();
      }
      go.output <- enricher(gene = out_hc_genes, TERM2GENE = GO.BG, TERM2NAME = GO.name.BG)
      go.output.df <- go.output@result
      go.output.df <- go.output.df[go.output.df$p.adjust < 0.05, ]
      go.output.df$p.adjust <- format(go.output.df$p.adjust, scientific = T, digits = 3)
      go.output.df$pvalue <- format(go.output.df$pvalue, scientific = T, digits = 3)
      GO_type_vec <- sapply(GO.name.list[go.output.df$ID], function(x){
        if(length(x)==1){return(x@Ontology)}else{return(NA)}
      })
      go.output.df$Type <- GO_type_vec
      go.output.df <- go.output.df[!is.na(go.output.df$Type), ]
      go.output.df <- go.output.df[, c(1,2,10,3,4,5,6,8)]
      go.output.df <- go.output.df[order(go.output.df$Type), ]
      # write.table(go.output.df, file = paste0("GO_output/", spp, "_", ii, ".txt"), quote = F, row.names = F, sep = "\t")
      if(sum(go.output.df$Type=="BP")>0){
        tmp.df <- go.output.df[go.output.df$Type=="BP", ]
        go.output.plots[[paste0(spp, "_", ii)]] <- GObarPlot(tmp.df, paste0(spp, "_", ii))
      }
    }
  }
}

# pdf("GO_high_low_conserved_m6A.pdf", width = 17, height = 25)
# patchwork::wrap_plots(go.output.plots, ncol = 2, heights = c(1,1,1,1,1,1,1,0.5,0.2))
# dev.off()

GO.plant.df <- rbind(go.output.plots$ath_highc$data, go.output.plots$ghi_highc$data, 
                  go.output.plots$zma_highc$data,go.output.plots$ata_highc$data, 
                  go.output.plots$tdi_highc$data, go.output.plots$tae_highc$data,
                  go.output.plots$ppa_highc$data)
GO.plant.df$elimFisher <- as.numeric(GO.plant.df$elimFisher)
GO.plant.result <- GO.plant.df %>% filter(Type == "BP")  %>% 
  group_by(Term) %>% 
  filter(elimFisher == min(elimFisher))
GO.plant.result <- GO.plant.result[order(GO.plant.result$elimFisher), ]
GO.plant.result$elimFisher <- format(GO.plant.result$elimFisher, 
                                  scientific = T, digits = 3)
cat(nrow(GO.plant.result), unique(length(GO.plant.result$Term)), "\n")
hc.bar.output <- GObarPlot(GO.plant.result, paste0("Orthologous genes with high-conservation of m6A methylation"))

GO.plant.df <- rbind(go.output.plots$zma_lowc$data, go.output.plots$ppa_lowc$data)
GO.plant.df$elimFisher <- as.numeric(GO.plant.df$elimFisher)
GO.plant.result <- GO.plant.df %>% filter(Type == "BP")  %>% 
  group_by(Term) %>% 
  filter(elimFisher == min(elimFisher))
GO.plant.result <- GO.plant.result[order(GO.plant.result$elimFisher), ]
GO.plant.result$elimFisher <- format(GO.plant.result$elimFisher, 
                                  scientific = T, digits = 3)
cat(nrow(GO.plant.result), unique(length(GO.plant.result$Term)), "\n")
lc.bar.output <- GObarPlot(GO.plant.result, paste0("Orthologous genes with species-specific m6A methylation"))

pdf("figures/FigureS5.GO_high_low_conserved_m6A.pdf", width = 12, height = 15)
patchwork::wrap_plots(list(hc.bar.output, lc.bar.output), ncol = 1, heights = c(1,0.6))
dev.off()
