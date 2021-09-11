## Significantly expanded gene families #####

# current wroking directory
rm(list=ls())
setwd("~/a2z/m6A_13spp/Part_5/")
options(stringsAsFactors = F)

## loading libraries
library(ggplot2)
library(corrplot)
library(raster)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(ggplotify)
library(ggtree)


# reload data
load("../Part_1/data1_m6A_peaks_genes.RData")
load("../Part_1/data2_gene_attributes.RData")

link.data.path <- "../associated_data/"
comlab.name <- "m6A methylation ratio (%)"
system(paste0("mkdir -p figures"))
species7.names <- c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa")

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
    BaseTheme() +
    geom_smooth(method="lm", se=F, color="#A3A2A2") +
    annotate("text", label = sprintf('Pearson\'s r = %.2f', 
                                     cor(df_input[, var1], df_input[, var2])), 
             x = text_pos[1], y = text_pos[2], size = 4) +
    annotate("text", label = sprintf('P-value = %.2g', p_value), 
             x = text_pos[3], y = text_pos[4], size = 4)
  out_plot
}


# Computational Analysis of gene Family Evolution ####
species8.names <- c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa", "ppa")
ec.PS.names <- c("Chlorophyta", "Chloroplastida", "Embryophyta", "Stomatophyta", 
                 "Tracheophyta", "Spermatophyta", "Acrogymnospermae",
                 "Angiosperm", "Mesangiospermae", "Monocots", "Commelinids", 
                 "Poaceae", "BOP", "Pooideae", "Hordeeae", "PACMAD", 
                 "Paniceae", "Setaria","Eudicots", "Solanaceae", "Rosids", 
                 "Eurosids", "Malvids", "subMalvids", "Brassicales", 
                 "Brassicaceae", "Fabids", "Fabaceae")
species34.names <- colnames(OGs.34spp)[1:34]
ec.PS.count <- matrix(nrow = length(c(ec.PS.names, species34.names)) , 
                      ncol = length(ec.PS.names))
colnames(ec.PS.count) <- ec.PS.names
rownames(ec.PS.count) <- c(ec.PS.names, species34.names)

for(tmp.PS in ec.PS.names){
  tmp.exp.con <- read.table(paste0(link.data.path, "cafe_result/",
                                   tmp.PS, "/", tmp.PS, "_all_fams.txt"), sep = "\t", 
                            row.names = 1, stringsAsFactors = F, comment.char = "")
  cat("\n### ", tmp.PS, " ###\n", tmp.exp.con[1,1], "\n")
}

name.change.list <- list(
  "Chlorophyta" = data.frame("Nodes"=c("MpusillaCCMP1545<0>:", "Creinhardtii<2>:"), "Names"=c("MpusillaCCMP1545", "Creinhardtii")),
  "Chloroplastida" = data.frame(
    "Nodes"=c("MpusillaCCMP1545<0>:", "<1>:", "Creinhardtii<2>:", "Mpolymorpha<4>:", "<5>:", "ppa<6>:", "<7>:", "Smoellendorffii<8>:", 
              "<9>:","osa<10>:", "<11>:", "Hvulgare<12>:", "<13>:", "ata<14>:", "<15>:", "Bdistachyon<16>:", "<17>:","Sitalica<18>:",
              "<19>:", "Sviridis<20>:", "<21>:", "sbi<22>:", "<23>:", "Eguineensis<24>:", "<25>:", "Spolyrhiza<26>:", "<27>:",
              "Vvinifera<28>:", "<29>:", "Egrandis<30>:", "<31>:", "gar<32>:", "<33>:", "ath<34>:", "<35>:", "Bstricta<36>:", 
              "<37>:", "pvu<38>:", "<39>:", "Mtruncatula<40>:", "<41>:", "Ptrichocarpa<42>:", "<43>:", "Stuberosum<44>:", 
              "<45>:", "sly<46>:", "<47>:", "Atrichopoda<48>:", "<49>:", "Gibi<50>:", "<51>:","Paab<52>:"),
    "Names"=c("MpusillaCCMP1545", "Chlorophyta", "Creinhardtii", "Mpolymorpha", "Embryophyta", "ppa", "Stomatophyta", 
              "Smoellendorffii", "Tracheophyta", "osa", "BOP", 
              "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae","Sitalica", "Setaria", "Sviridis", "PACMAD", 
              "sbi", "Commelinids", "Eguineensis", "Monocots", "Spolyrhiza", "Mesangiospermae","Vvinifera", "Rosids", "Egrandis", 
              "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", "Fabaceae", "Mtruncatula", 
              "Fabids", "Ptrichocarpa", "Eudicots", "Stuberosum", "Solanaceae", "sly", "Angiosperm", 
              "Atrichopoda", "Spermatophyta", "Gibi", "Acrogymnospermae","Paab")),
  "Embryophyta" = data.frame(
    "Nodes"=c("Mpolymorpha<0>:", "ppa<2>:", "<3>:", "Smoellendorffii<4>:", "<5>:", "osa<6>:", "<7>:", "Hvulgare<8>:", "<9>:", "ata<10>:",
              "<11>:", "Bdistachyon<12>:", "<13>:","Sitalica<14>:", "<15>:", "Sviridis<16>:", "<17>:", "sbi<18>:", "<19>:", 
              "Eguineensis<20>:", "<21>:", "Spolyrhiza<22>:", "<23>:","Vvinifera<24>:", "<25>:", "Egrandis<26>:", "<27>:", 
              "gar<28>:", "<29>:", "ath<30>:", "<31>:", "Bstricta<32>:", "<33>:", "pvu<34>:", "<35>:", "Mtruncatula<36>:", 
              "<37>:", "Ptrichocarpa<38>:", "<39>:", "Stuberosum<40>:", "<41>:", "sly<42>:", "<43>:", "Atrichopoda<44>:", 
              "<45>:", "Gibi<46>:", "<47>:","Paab<48>:"),
    "Names"=c("Mpolymorpha", "ppa", "Stomatophyta", "Smoellendorffii", "Tracheophyta", "osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Commelinids", "Eguineensis", "Monocots", "Spolyrhiza", "Mesangiospermae",
              "Vvinifera", "Rosids", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa", "Eudicots", "Stuberosum", "Solanaceae", "sly", "Angiosperm", 
              "Atrichopoda", "Spermatophyta", "Gibi", "Acrogymnospermae","Paab")),
  "Stomatophyta" = data.frame(
    "Nodes"=c("ppa<0>:", "Smoellendorffii<2>:", "<3>:", "osa<4>:", "<5>:", "Hvulgare<6>:", "<7>:", "ata<8>:", "<9>:", "Bdistachyon<10>:",
              "<11>:","Sitalica<12>:", "<13>:", "Sviridis<14>:", "<15>:", "sbi<16>:", "<17>:", "Eguineensis<18>:", "<19>:", 
              "Spolyrhiza<20>:", "<21>:","Vvinifera<22>:", "<23>:", "Egrandis<24>:", "<25>:", "gar<26>:", "<27>:", "ath<28>:", "<29>:", 
              "Bstricta<30>:", "<31>:", "pvu<32>:", "<33>:", "Mtruncatula<34>:", "<35>:", "Ptrichocarpa<36>:", "<37>:", 
              "Stuberosum<38>:", "<39>:", "sly<40>:", "<41>:", "Atrichopoda<42>:", "<43>:", "Gibi<44>:", "<45>:","Paab<46>:"),
    "Names"=c("ppa", "Smoellendorffii", "Tracheophyta", "osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
             "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Commelinids", "Eguineensis", "Monocots", "Spolyrhiza", "Mesangiospermae",
             "Vvinifera", "Rosids", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
             "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa", "Eudicots", "Stuberosum", "Solanaceae", "sly", "Angiosperm", 
             "Atrichopoda", "Spermatophyta", "Gibi", "Acrogymnospermae","Paab")),
  "Tracheophyta" = data.frame(
    "Nodes"=c("Smoellendorffii<0>:", "osa<2>:", "<3>:", "Hvulgare<4>:", "<5>:", "ata<6>:", "<7>:", "Bdistachyon<8>:",
              "<9>:","Sitalica<10>:", "<11>:", "Sviridis<12>:", "<13>:", "sbi<14>:", "<15>:", "Eguineensis<16>:", "<17>:", 
              "Spolyrhiza<18>:", "<19>:","Vvinifera<20>:", "<21>:", "Egrandis<22>:", "<23>:", "gar<24>:", "<25>:", "ath<26>:", "<27>:", 
              "Bstricta<28>:", "<29>:", "pvu<30>:", "<31>:", "Mtruncatula<32>:", "<33>:", "Ptrichocarpa<34>:", "<35>:", 
              "Stuberosum<36>:", "<37>:", "sly<38>:", "<39>:", "Atrichopoda<40>:", "<41>:", "Gibi<42>:", "<43>:","Paab<44>:"),
    "Names"=c("Smoellendorffii", "osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Commelinids", "Eguineensis", "Monocots", "Spolyrhiza", "Mesangiospermae",
              "Vvinifera", "Rosids", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa", "Eudicots", "Stuberosum", "Solanaceae", "sly", "Angiosperm", 
              "Atrichopoda", "Spermatophyta", "Gibi", "Acrogymnospermae","Paab")),
  "Spermatophyta" = data.frame(
    "Nodes"=c("osa<0>:", "<1>:", "Hvulgare<2>:", "<3>:", "ata<4>:", "<5>:", "Bdistachyon<6>:", "<7>:","Sitalica<8>:", "<9>:", 
              "Sviridis<10>:", "<11>:", "sbi<12>:", "<13>:", "Eguineensis<14>:", "<15>:", 
              "Spolyrhiza<16>:", "<17>:","Vvinifera<18>:", "<19>:", "Egrandis<20>:", "<21>:", "gar<22>:", "<23>:", "ath<24>:", "<25>:", 
              "Bstricta<26>:", "<27>:", "pvu<28>:", "<29>:", "Mtruncatula<30>:", "<31>:", "Ptrichocarpa<32>:", "<33>:", 
              "Stuberosum<34>:", "<35>:", "sly<36>:", "<37>:", "Atrichopoda<38>:", "Gibi<40>:", "<41>:","Paab<42>:"),
    "Names"=c("osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Commelinids", "Eguineensis", "Monocots", "Spolyrhiza", "Mesangiospermae",
              "Vvinifera", "Rosids", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa", "Eudicots", "Stuberosum", "Solanaceae", "sly", "Angiosperm", 
              "Atrichopoda", "Gibi", "Acrogymnospermae","Paab")),
  "Acrogymnospermae" =  data.frame("Nodes"=c("Gibi<0>:", "Paab<2>:"), "Names"=c("Gibi", "Paab")),
  "Angiosperm" = data.frame(
    "Nodes"=c("osa<0>:", "<1>:", "Hvulgare<2>:", "<3>:", "ata<4>:", "<5>:", "Bdistachyon<6>:", "<7>:","Sitalica<8>:", "<9>:", 
              "Sviridis<10>:", "<11>:", "sbi<12>:", "<13>:", "Eguineensis<14>:", "<15>:", 
              "Spolyrhiza<16>:", "<17>:","Vvinifera<18>:", "<19>:", "Egrandis<20>:", "<21>:", "gar<22>:", "<23>:", "ath<24>:", "<25>:", 
              "Bstricta<26>:", "<27>:", "pvu<28>:", "<29>:", "Mtruncatula<30>:", "<31>:", "Ptrichocarpa<32>:", "<33>:", 
              "Stuberosum<34>:", "<35>:", "sly<36>:", "Atrichopoda<38>:"),
    "Names"=c("osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Commelinids", "Eguineensis", "Monocots", "Spolyrhiza", "Mesangiospermae",
              "Vvinifera", "Rosids", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa", "Eudicots", "Stuberosum", "Solanaceae", "sly", 
              "Atrichopoda")),
  "Mesangiospermae" = data.frame(
    "Nodes"=c("osa<0>:", "<1>:", "Hvulgare<2>:", "<3>:", "ata<4>:", "<5>:", "Bdistachyon<6>:", "<7>:","Sitalica<8>:", "<9>:", 
              "Sviridis<10>:", "<11>:", "sbi<12>:", "<13>:", "Eguineensis<14>:", "<15>:", 
              "Spolyrhiza<16>:","Vvinifera<18>:", "<19>:", "Egrandis<20>:", "<21>:", "gar<22>:", "<23>:", "ath<24>:", "<25>:", 
              "Bstricta<26>:", "<27>:", "pvu<28>:", "<29>:", "Mtruncatula<30>:", "<31>:", "Ptrichocarpa<32>:", "<33>:", 
              "Stuberosum<34>:", "<35>:", "sly<36>:"),
    "Names"=c("osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Commelinids", "Eguineensis", "Monocots", "Spolyrhiza",
              "Vvinifera", "Rosids", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa", "Eudicots", "Stuberosum", "Solanaceae", "sly")),
  "Monocots" = data.frame(
    "Nodes"=c("osa<0>:", "<1>:", "Hvulgare<2>:", "<3>:", "ata<4>:", "<5>:", "Bdistachyon<6>:", "<7>:","Sitalica<8>:", "<9>:", 
              "Sviridis<10>:", "<11>:", "sbi<12>:", "<13>:", "Eguineensis<14>:", "Spolyrhiza<16>:"),
    "Names"=c("osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Commelinids", "Eguineensis", "Spolyrhiza")),
  "Commelinids" = data.frame(
    "Nodes"=c("osa<0>:", "<1>:", "Hvulgare<2>:", "<3>:", "ata<4>:", "<5>:", "Bdistachyon<6>:", "<7>:","Sitalica<8>:", "<9>:", 
              "Sviridis<10>:", "<11>:", "sbi<12>:", "Eguineensis<14>:"),
    "Names"=c("osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon", "Poaceae",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi", "Eguineensis")),
  "Poaceae" = data.frame(
    "Nodes"=c("osa<0>:", "<1>:", "Hvulgare<2>:", "<3>:", "ata<4>:", "<5>:", "Bdistachyon<6>:","Sitalica<8>:", "<9>:", 
              "Sviridis<10>:", "<11>:", "sbi<12>:"),
    "Names"=c("osa", "BOP", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon",
              "Sitalica", "Setaria", "Sviridis", "PACMAD", "sbi")),
  "BOP" = data.frame("Nodes"=c("osa<0>:", "Hvulgare<2>:", "<3>:", "ata<4>:", "<5>:", "Bdistachyon<6>:"),
                     "Names"=c("osa", "Hvulgare", "Hordeeae", "ata", "Pooideae", "Bdistachyon")),
  "Pooideae" = data.frame("Nodes"=c("Hvulgare<0>:", "<1>:", "ata<2>:", "Bdistachyon<4>:"),
                          "Names"=c("Hvulgare", "Hordeeae", "ata", "Bdistachyon")),
  "Hordeeae" = data.frame("Nodes"=c("Hvulgare<0>:", "ata<2>:"), "Names"=c("Hvulgare", "ata")),
  "PACMAD" = data.frame("Nodes"=c("Sitalica<0>:", "<1>:", "Sviridis<2>:", "sbi<4>:"),
                        "Names"=c("Sitalica", "Setaria", "Sviridis", "sbi")),
  "Paniceae" = data.frame("Nodes"=c("Sitalica<0>:", "Sviridis<2>:"), "Names"=c("Sitalica", "Sviridis")),
  "Setaria" = data.frame("Nodes"=c("Sitalica<0>:", "Sviridis<2>:"), "Names"=c("Sitalica", "Sviridis")),
  "Eudicots" = data.frame(
    "Nodes"=c("Vvinifera<0>:", "<1>:", "Egrandis<2>:", "<3>:", "gar<4>:", "<5>:", "ath<6>:", "<7>:", 
              "Bstricta<8>:", "<9>:", "pvu<10>:", "<11>:", "Mtruncatula<12>:", "<13>:", "Ptrichocarpa<14>:", 
              "Stuberosum<16>:", "<17>:", "sly<18>:"),
    "Names"=c("Vvinifera", "Rosids", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa", "Stuberosum", "Solanaceae", "sly")),
  "Solanaceae" = data.frame("Nodes"=c("Stuberosum<0>:", "sly<2>:"), "Names"=c("Stuberosum", "sly")),
  "Rosids" = data.frame(
    "Nodes"=c("Vvinifera<0>:", "Egrandis<2>:", "<3>:", "gar<4>:", "<5>:", "ath<6>:", "<7>:", 
              "Bstricta<8>:", "<9>:", "pvu<10>:", "<11>:", "Mtruncatula<12>:", "<13>:", "Ptrichocarpa<14>:"),
    "Names"=c("Vvinifera", "Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "Eurosids", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa")),
  "Eurosids" = data.frame(
    "Nodes"=c("Egrandis<0>:", "<1>:", "gar<2>:", "<3>:", "ath<4>:", "<5>:", 
              "Bstricta<6>:", "pvu<8>:", "<9>:", "Mtruncatula<10>:", "<11>:", "Ptrichocarpa<12>:"),
    "Names"=c("Egrandis", "Malvids", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta", "pvu", 
              "Fabaceae", "Mtruncatula", "Fabids", "Ptrichocarpa")),
  "Malvids" = data.frame(
    "Nodes"=c("Egrandis<0>:", "gar<2>:", "<3>:", "ath<4>:", "<5>:", "Bstricta<6>:"),
    "Names"=c("Egrandis", "gar", "subMalvids", "ath", "Brassicaceae", "Bstricta")),
  "subMalvids" = data.frame("Nodes"=c("gar<0>:", "ath<2>:", "<3>:", "Bstricta<4>:"), "Names"=c("gar", "ath", "Brassicaceae", "Bstricta")),
  "Brassicales" = data.frame("Nodes"=c("ath<0>:", "Bstricta<2>:"), "Names"=c("ath", "Bstricta")),
  "Brassicaceae" = data.frame("Nodes"=c("ath<0>:", "Bstricta<2>:"), "Names"=c("ath", "Bstricta")),
  "Fabids" = data.frame("Nodes"=c("pvu<0>:", "<1>:", "Mtruncatula<2>:", "Ptrichocarpa<4>:"),
                        "Names"=c("pvu", "Fabaceae", "Mtruncatula", "Ptrichocarpa")),
  "Fabaceae" = data.frame("Nodes"=c("pvu<0>:", "Mtruncatula<2>:"), "Names"=c("pvu", "Mtruncatula"))
)

ec.PS.all.add <- ec.PS.all.minus <- list()
ec.PS.rapid.add <- ec.ps.rapid.minus <- list()

for (tmp.PS in ec.PS.names){
  tmp.exp.con <- read.table(paste0(link.data.path, "cafe_result/",
                                   tmp.PS, "/", tmp.PS, "_all_fams.txt"), sep = "\t", 
                            row.names = 1, stringsAsFactors = F, comment.char = "")
  tmp.exp.con$OGs.add <- apply(tmp.exp.con, 1, function(x){
    paste0(unlist(strsplit( x[1], ","))[grepl(pattern = "\\+", x = unlist(strsplit(x[1], ",")))], collapse = ",")
  })
  tmp.exp.con$OGs.minus <- apply(tmp.exp.con, 1, function(x){
    paste0(unlist(strsplit( x[1], ","))[grepl(pattern = "\\-", x = unlist(strsplit(x[1], ",")))], collapse = ",")
  })
  tmp.exp.con$OGs.add <- gsub("\\[.+?\\]", "", tmp.exp.con$OGs.add)
  tmp.exp.con$OGs.minus <- gsub("\\[.+?\\]", "", tmp.exp.con$OGs.minus)
  tmp.exp.con$Rapid <- apply(tmp.exp.con, 1, function(x){
    paste0(unlist(strsplit( x[1], ","))[grepl(pattern = "\\*", x = unlist(strsplit( x[1], ",")))], collapse = ",")
  })
  tmp.exp.con$Rapid.add <- apply(tmp.exp.con, 1, function(x){
    paste0(unlist(strsplit( x[4], ","))[grepl(pattern = "\\+", x = unlist(strsplit( x[4], ",")))], collapse = ",")
  })
  tmp.exp.con$Rapid.minus <- apply(tmp.exp.con, 1, function(x){
    paste0(unlist(strsplit( x[4], ","))[grepl(pattern = "\\-", x = unlist(strsplit( x[4], ",")))], collapse = ",")
  })
  tmp.exp.con$Rapid.add <- gsub("\\[.+?\\]", "", tmp.exp.con$Rapid.add)
  tmp.exp.con$Rapid.minus <- gsub("\\[.+?\\]", "", tmp.exp.con$Rapid.minus)
  tmp.exp.con <- tmp.exp.con[3:nrow(tmp.exp.con), ]
  ##
  tmp.node.names <- name.change.list[[tmp.PS]]$Nodes
  cat(tmp.PS, length(tmp.node.names), 
      sum(tmp.node.names%in%rownames(tmp.exp.con)) == length(tmp.node.names), "\n")
  if (all(tmp.node.names%in%rownames(tmp.exp.con))){
    tmp.exp.con <- tmp.exp.con[intersect(tmp.node.names, rownames(tmp.exp.con)), ]
    rownames(tmp.exp.con) <- name.change.list[[tmp.PS]]$Names
  }
  if(any(tmp.exp.con$OGs.add != "")){
    tmp_all_add <- tmp.exp.con[tmp.exp.con$OGs.add != "", ]
    # tmp.count <- sapply(strsplit(tmp_all_add$OGs.add, ","), length)
    # ec.PS.count[rownames(tmp_all_add), tmp.PS] <- tmp.count
    for (ii in rownames(tmp_all_add)){
      if(ii %in% names(ec.PS.all.add)){
        ec.PS.all.add[[ii]] <- unique(c(ec.PS.all.add[[ii]], 
                                        unlist(strsplit(tmp_all_add[ii, "OGs.add"], split = ","))))
      } else {
        ec.PS.all.add[[ii]] <- unlist(strsplit(tmp_all_add[ii, "OGs.add"], 
                                               split = ","))
      }
    }
  }
  ##
  if (any(tmp.exp.con$OGs.minus != "")){
    tmp.all.minus <- tmp.exp.con[tmp.exp.con$OGs.minus != "", ]
    # tmp.count <- sapply(strsplit(tmp.all.minus$OGs.minus, ","), length)
    # ec.PS.count[rownames(tmp.all.minus), tmp.PS] <- tmp.count
    for (ii in rownames(tmp.all.minus)){
      if (ii %in% names(ec.PS.all.minus)){
        ec.PS.all.minus[[ii]] <- unique(c(ec.PS.all.minus[[ii]], 
                                          unlist(strsplit(tmp.all.minus[ii, "OGs.minus"], split = ","))))
      } else {
        ec.PS.all.minus[[ii]] <- unlist(strsplit(tmp.all.minus[ii, "OGs.minus"], 
                                                 split = ","))
      }
    }
  }
  ##
  if (any(tmp.exp.con$Rapid.add != "")){
    tmp.rapid.add <- tmp.exp.con[tmp.exp.con$Rapid.add != "", ]
    tmp.count <- sapply(strsplit(tmp.rapid.add$Rapid.add, ","), length)
    ec.PS.count[rownames(tmp.rapid.add), tmp.PS] <- tmp.count
    for (ii in rownames(tmp.rapid.add)){
      if (ii %in% names(ec.PS.rapid.add)){
        ec.PS.rapid.add[[ii]] <- unique(c(ec.PS.rapid.add[[ii]], 
                                          unlist(strsplit(tmp.rapid.add[ii, "Rapid.add"], split = ","))))
      } else {
        ec.PS.rapid.add[[ii]] <- unlist(strsplit(tmp.rapid.add[ii, "Rapid.add"], 
                                                 split = ","))
      }
    }
  }
  ##
  if (any(tmp.exp.con$Rapid.minus!="")){
    tmp.rapid.minus <- tmp.exp.con[tmp.exp.con$Rapid.minus!="", ]
    # tmp.rapid.minustmp.count <- sapply(strsplit(tmp.rapid.minus$Rapid.minus, ","), length)
    # ec.PS.count[rownames(tmp.rapid.minus), tmp.PS] <- tmp.count
    for (ii in rownames(tmp.rapid.minus)){
      if (ii %in% names(ec.ps.rapid.minus)){
        ec.ps.rapid.minus[[ii]] <- unique(c(ec.ps.rapid.minus[[ii]], 
                                            unlist(strsplit(tmp.rapid.minus[ii, "Rapid.minus"], split = ","))))
      } else {
        ec.ps.rapid.minus[[ii]] <- unlist(strsplit(tmp.rapid.minus[ii, "Rapid.minus"], 
                                                   split = ","))
      }
    }
  }
}


sp.evolution.path <- list(
  "ath"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae", 
          "Eudicots", "Rosids", "Eurosids", "Malvids", "subMalvids", "Brassicales", "Brassicaceae", "ath"),
  "gar"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae", 
          "Eudicots", "Rosids", "Eurosids", "Malvids", "subMalvids", "gar"),
  "pvu"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae", 
          "Eudicots", "Rosids", "Eurosids", "Fabids", "Fabaceae", "pvu"),
  "sly"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae", 
          "Eudicots", "Solanaceae", "sly"),
  "sbi"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae",
          "Monocots", "Commelinids", "Poaceae", "PACMAD", "sbi"),
  "ata"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae",
          "Monocots", "Commelinids", "Poaceae", "BOP", "Pooideae", "Hordeeae", "ata"),
  "osa"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae",
          "Monocots", "Commelinids", "Poaceae", "BOP", "osa"),
  "ppa"=c("Chloroplastida", "Embryophyta", "Stomatophyta", "ppa")
  )

## species path
sp.all.add <- sp.all.minus <- sp.rapid.add <- sp.rapid.minus <- list()
for(ii in names(sp.evolution.path)){
  ##
  tmp.overlap <- intersect(sp.evolution.path[[ii]], names(ec.PS.all.add))
  sp.all.add[[ii]] <- unique(unlist(ec.PS.all.add[tmp.overlap]))
  ##
  tmp.overlap <- intersect(sp.evolution.path[[ii]], names(ec.PS.all.minus))
  sp.all.minus[[ii]] <- unique(unlist(ec.PS.all.minus[tmp.overlap]))
  ##
  tmp.overlap <- intersect(sp.evolution.path[[ii]], names(ec.PS.rapid.add))
  sp.rapid.add[[ii]] <- unique(unlist(ec.PS.rapid.add[tmp.overlap]))
  ##
  tmp.overlap <- intersect(sp.evolution.path[[ii]], names(ec.ps.rapid.minus))
  sp.rapid.minus[[ii]] <- unique(unlist(ec.ps.rapid.minus[tmp.overlap]))
}


# merging node
node.order <- c("Ath","Bstricta","Gar","Egrandis","Pvu","Mtruncatula",
                "Ptrichocarpa","Vvinifera","Stuberosum","Sly","Sitalica",
                "Sviridis","Sbi","Ata","Hvulgare","Bdistachyon","Osa",
                "Eguineensis","Spolyrhiza","Atrichopoda","Paab","Gibi",
                "Smoellendorffii","Ppa","Mpolymorpha","Creinhardtii",
                "MpusillaCCMP1545","Embryophyta","Stomatophyta","Tracheophyta",
                "Spermatophyta","Acrogymnospermae","Angiosperm","Mesangiospermae",
                "Monocots","Commelinids","Poaceae","Pooideae","Hordeeae","PACMAD",
                "Setaria","Eudicots","Rosids","Eurosids","SubMalvids","Brassicaceae",
                "Fabaceae","Solanaceae")
tmp.ec.PS.all.add <- ec.PS.all.add
tmp.ec.PS.rapid.add <- ec.PS.rapid.add
names(tmp.ec.PS.all.add) <- FirstUp(names(tmp.ec.PS.all.add))
names(tmp.ec.PS.rapid.add) <- FirstUp(names(tmp.ec.PS.rapid.add))
all(node.order%in%names(tmp.ec.PS.all.add)) & all(names(tmp.ec.PS.all.add)%in%node.order)
all(node.order%in%names(tmp.ec.PS.rapid.add)) & all(names(tmp.ec.PS.rapid.add)%in%node.order)
tmp.ec.PS.all.add <- tmp.ec.PS.all.add[node.order]
tmp.ec.PS.rapid.add <- tmp.ec.PS.rapid.add[node.order]

all.add.OGs <- unique(unlist(tmp.ec.PS.all.add))
node.num <- length(tmp.ec.PS.all.add)
out.data <- matrix(nrow = length(all.add.OGs), 
                   ncol = node.num*2)
rownames(out.data) <- all.add.OGs
colnames(out.data) <- rep(names(tmp.ec.PS.all.add), each=2)
colnames(out.data)[(1:node.num)*2] <- "NULL"
dim(out.data)

for (ii in all.add.OGs){
  for (jj in 1:node.num){
    if (ii%in%tmp.ec.PS.all.add[[jj]]){
      out.data[ii, jj*2-1] <- "T"
    }
    if(ii%in%tmp.ec.PS.rapid.add[[jj]]){
      out.data[ii, jj*2] <- "T"
    }
  }
}
# write.table(out.data, file = "Table5.txt", quote = F, sep = "\t", row.names = T)


## Figure5A ####
### overlapped gene families
OGs.count <- table(unlist(sp.rapid.add[species7.names]))
overlap.OGs <- names(OGs.count)[OGs.count == 7]

m6Aratio.overlap.OGs <- vector()
for (tmp.spp in species7.names){
  for (tmp.OGs in overlap.OGs){
    cat(tmp.spp, tmp.OGs, "\n")
    tmp.PS.genes <- unlist(strsplit(OGs.34spp[tmp.OGs, tmp.spp], split = ", "))
    # tmp.PS.exp.genes <- sum(tmp.PS.genes%in%unlist(exp.list[tmp.spp]))
    tmp.PS.m6A.genes <- sum(tmp.PS.genes%in%unlist(m6A.list[tmp.spp]))
    m6Aratio.overlap.OGs <- c(m6Aratio.overlap.OGs, 
                              tmp.PS.m6A.genes/length(tmp.PS.genes)) ### 2021-08-01 
  }
}

m6Aratio.overlap.OGs.df <- data.frame(
  "Species"= rep(species7.names, each=length(overlap.OGs)),
  "Values"= m6Aratio.overlap.OGs, 
  "OGs"=rep(overlap.OGs, 7))

m6Aratio.overlap.OGs.df$Species <- factor(m6Aratio.overlap.OGs.df$Species, 
                                          levels = intersect(c("ath", "gar", "pvu", "sly", "sbi", 
                                                               "ata", "osa", "ppa"), 
                                                             unique(m6Aratio.overlap.OGs.df$Species)))
m6Aratio.overlap.OGs.data <- reshape2::dcast(m6Aratio.overlap.OGs.df, formula = OGs~Species, value.var="Values")
rownames(m6Aratio.overlap.OGs.data) <- m6Aratio.overlap.OGs.data[,1]
m6Aratio.overlap.OGs.data <- m6Aratio.overlap.OGs.data[, -1]
sum(is.na(m6Aratio.overlap.OGs.data))
m6Aratio.overlap.OGs.data[is.na(m6Aratio.overlap.OGs.data)] <- 0
m6Aratio.overlap.OGs.data <- m6Aratio.overlap.OGs.data[, species7.names]
m6Aratio.overlap.cor <- cor(m6Aratio.overlap.OGs.data)
m6Aratio.overlap.cor[m6Aratio.overlap.cor == 1] <- 0

pdf("figures/Figure5A.pdf", height = 4, width = 4)
corrplot(corr = m6Aratio.overlap.cor, cl.lim = c(0,1), cl.length = 6, #add=TRUE, 
         col = colorRampPalette(c("black", "green",
                                  "white", "#FB8D42", 
                                  "#7E0527"))(100), 
         type="upper", method="color",
         order="original", addCoef.col = "white", #"#EE7E32", "#E30713"
         tl.pos="td")
dev.off()


## Figure 5B ####
ref.values <- c(12568/27420, 15286/33268, 11777/27433, 
               11000/34075, 13031/34129, 13521/39614,
               15303/42189)*100
ref.value.list <- list()
ref.value.list[['ath']] <- c(12568, 27420-12568)
ref.value.list[['gar']] <- c(15286, 33268-15286)
ref.value.list[['pvu']] <- c(11777, 27433-11777)
ref.value.list[['sly']] <- c(11000, 34075-11000)
ref.value.list[['sbi']] <- c(13031, 34129-13031)
ref.value.list[['ata']] <- c(13521, 39614-13521)
ref.value.list[['osa']] <- c(15303, 42189-15303)


m6A.rapid.ratio <- vector()
for(tmp.sp in species7.names){
  tmp.OGs <- OGs.34spp[sp.rapid.add[[tmp.sp]], tmp.sp]
  tmp.OGs.genes <- unlist(strsplit(tmp.OGs, split = ", "))
  tmp.m6A <- sum(tmp.OGs.genes%in%m6A.list[[tmp.sp]])
  tmp.m6A.ratio <- tmp.m6A/length(tmp.OGs.genes)*100              ## 2021-08-02 
  # cat(tmp.sp, tmp.m6A, length(tmp.OGs.genes)-tmp.m6A, "\n") 
  cat(tmp.sp, chisq.test(matrix(c(tmp.m6A,
                                       length(tmp.OGs.genes)-tmp.m6A,
                                       ref.value.list[[tmp.sp]]),
                                     nrow = 2))$p.value, "\n")
  m6A.rapid.ratio <- c(m6A.rapid.ratio, tmp.m6A.ratio)
}

df.m6A.com <- data.frame("species"=rep(FirstUp(species7.names), 4), 
                         "value"=c(ref.values, 100 - ref.values, m6A.rapid.ratio, 
                                   100 - m6A.rapid.ratio),
                         "Type"=rep(c("ref", "rapid"), each=14),
                         "Type2"=rep(c("ref-m6A", "ref-Non-m6A",
                                       "rapid-m6A", "rapid-Non-m6A"), each=7))
df.m6A.com$species <- factor(df.m6A.com$species, levels = unique(df.m6A.com$species))
df.m6A.com$Type2 <- factor(df.m6A.com$Type2, levels = rev(unique(df.m6A.com$Type2)))

pdf("figures/Figure5B.pdf", height = 3, width = 7)
ggplot(data = df.m6A.com, aes(x=Type, y=value, fill=Type2))+
  geom_bar(stat='identity', width=0.7) +
  theme_bw(base_size = 14) +
  BaseTheme() +
  xlab("") + ylab(comlab.name) +
  scale_fill_manual("Type", values = c("#CCCCCC","#557773", "#CCCCCC", "#A68171")) + 
  facet_grid(.~species) + coord_cartesian(ylim = c(0, 50))
dev.off()


## Fgiure 5C and Supplementary Fig. 9 ####
Node1 <- c("Embryophyta", "Stomatophyta")
Node1.spp <- c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa", "ppa")
Node2 <- c("Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae")
Node2.spp <- c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa")
Node3 <- c("Eudicots")
Node3.spp <- c("ath", "gar", "pvu", "sly")
Node4 <- c("Rosids", "Eurosids")
Node4.spp <- c("ath", "gar", "pvu")
Node5 <- c("Malvids", "subMalvids")
Node5.spp <- c("ath", "gar")
Node6 <- c("Brassicales", "Brassicaceae", "ath")
Node6.spp <- c("ath")
Node7 <- c("gar")
Node7.spp <- c("gar")
Node8 <- c("Fabids", "Fabaceae", "pvu")
Node8.spp <- c("pvu")
Node9 <- c("Solanaceae", "sly")
Node9.spp <- c("sly")
Node10 <- c("Monocots", "Commelinids", "Poaceae")
Node10.spp <- c("sbi", "ata", "osa")
Node11 <- c("PACMAD", "sbi")
Node11.spp <- c("sbi")
# Node12 <- c("BOP")
# Node12.spp <- c("ata", "osa")
Node13 <- c("Pooideae", "Hordeeae", "ata")
Node13.spp <- c("ata")
Node14 <- c("osa")
Node14.spp <- c("osa")
Node15 <- c("ppa")
Node15.spp <- c("ppa")

# Supplementary Figure 17 ####
node.each.name <- node.all.names <- node.all.m6Aratio <- node.all.numbers <- vector()
node.m6A.ratio <- node.mean.gene <- vector()
for(nodenum in c(1:11, 13:15)){
  tmp.node <- eval(parse(text = paste0("Node", nodenum)))
  tmp.node.sp <- eval(parse(text = paste0("Node", nodenum, ".spp")))
  raw_num <- raw.exp <- raw.m6A <- raw.mean.num <- vector()
  for(tmp.PS in tmp.node){
    if(tmp.PS %in% names(ec.PS.rapid.add)){
      for(tmp.spp in tmp.node.sp){
        # cat(length(ec.PS.rapid.add[[tmp.PS]]), "\n")
        tmp.gene.list <- strsplit(OGs.34spp[ec.PS.rapid.add[[tmp.PS]], tmp.spp], split = ", ")
        tmp.PS.genes <- unlist(tmp.gene.list)
        raw_num <- c(raw_num, length(tmp.PS.genes))
        raw.m6A <- c(raw.m6A, sum(tmp.PS.genes%in%m6A.list[[tmp.spp]]))
        raw.mean.num <- c(raw.mean.num, sapply(tmp.gene.list, length))
        ## detailed information
        node.each.name <- c(node.each.name, 
                            rep(paste0("Node", nodenum), length(tmp.gene.list)))
        node.all.names <- c(node.all.names, 
                            rep(tmp.PS, length(tmp.gene.list)))
        node.all.m6Aratio <- c(node.all.m6Aratio, 
                               sapply(tmp.gene.list, function(x){
                                 sum(x%in%m6A.list[[tmp.spp]])/length(x)*100 # 2021-08-01 
                                 })) 
        node.all.numbers <- c(node.all.numbers, sapply(tmp.gene.list, length))
      }      
    }
  }
  node.m6A.ratio <- c(node.m6A.ratio, sum(raw.m6A)/sum(raw_num)*100) # 2021-08-01 
  node.mean.gene <- c(node.mean.gene, mean(unlist(raw.mean.num)))
  cat(paste0("Node", nodenum), tmp.node, round(sum(raw.m6A)/sum(raw_num)*100, digits = 2), 
      round(mean(unlist(raw.mean.num)), digits = 2), "\n")
}

## Node1 Node2 Node3 Node4 Node5 Node10
ratio2count.df <- data.frame("Nodes" = node.each.name, "PSs" = node.all.names,
                             "ratios"= node.all.m6Aratio, "count" = node.all.numbers)
ratio2count.df <- ratio2count.df[ratio2count.df$count >= 10, ]
ratio2count.df <- ratio2count.df[ratio2count.df$ratios != 0, ]
ratio2count.df <- ratio2count.df[ratio2count.df$Nodes%in% c("Node1", "Node2", 
                                                            "Node3", "Node4",
                                                            "Node5", "Node10"), ]
ratio2count.df$Nodes <- factor(ratio2count.df$Nodes, levels = unique(ratio2count.df$Nodes))

pdf("figures/FigureS17.pdf", height = 3, width = 12, useDingbats = FALSE)
ggplot(data = ratio2count.df, mapping = aes(x = count, y = ratios)) +
  geom_point() + facet_grid(.~Nodes, scales = "free") + #space = "free"
  labs(x = '', y = '', title = '') +
  theme_bw(base_size = 14) +
  BaseTheme() +
  geom_smooth(method="lm", se=F, color="#DA2416")
dev.off()

for (tmp.PSs in unique(ratio2count.df$Nodes)){
  ratio2count.filter <- ratio2count.df[ratio2count.df$Nodes == tmp.PSs,]
  ratio2count.cor <- cor.test(ratio2count.filter$ratios, ratio2count.filter$count)
  cat(tmp.PSs, ratio2count.cor$estimate, ratio2count.cor$p.value, "\n")
}


# Figure 5D  ####
ratio2ratio <- list()
for (tmp.sp in species8.names){
  tmp.df <- t(apply(data.frame(sp.rapid.add[[tmp.sp]]), 1, function(x){
    gene.list <- unlist(strsplit(OGs.34spp[x[1], tmp.sp], ", "))
    gene.num <- length(gene.list)
    exp.num <- sum(gene.list%in%exp.list[[tmp.sp]])
    m6A.num <- sum(gene.list%in%m6A.list[[tmp.sp]])
    return(c(tmp.sp, x[1], gene.num, exp.num, m6A.num, m6A.num/gene.num*100))   ## 2021-07-28
  }))
  tmp.df <- tmp.df[which(as.numeric(tmp.df[,6]) > 0 & as.numeric(tmp.df[,3]) > 10), ]
  tmp.cor.result <- cor.test(as.numeric(tmp.df[,3]), as.numeric(tmp.df[,6]))
  cat(tmp.sp, tmp.cor.result$estimate, tmp.cor.result$p.value, "\n")
  ratio2ratio[[tmp.sp]] <- tmp.df 
}

ratio2ratio.df <- as.data.frame(do.call(rbind, ratio2ratio), stringsAsFactors = F)
colnames(ratio2ratio.df) <- c("Species", "OGs", "GeneNumber", "ExpNumber", "m6ANumber", "m6ARatio")
ratio2ratio.df[,3] <- as.numeric(ratio2ratio.df[,3])
ratio2ratio.df[,6] <- as.numeric(ratio2ratio.df[,6])
# write.table(ratio2ratio.df[,c(1,2,3,6)], "Figure5C_data.txt", quote = F, row.names = F)
ratio2ratio.df$Species <- factor(ratio2ratio.df$Species, levels = species8.names)

pdf("figures/Figure5D.pdf", height = 3, width = 16, useDingbats = F)
ggplot(data = ratio2ratio.df, mapping = aes(x = GeneNumber, y = m6ARatio)) +
  geom_point() + facet_grid(.~Species, scales = "free") +
  labs(x = '', y = '', title = '') +
  theme_bw(base_size = 14) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        # axis.text.x = element_text(colour = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_smooth(method="lm", se=F, color="#A3A2A2", span = 2) # lm
dev.off()



# Figure5E Supplementary Fig. 18 Kinases and Transciption factor #####
df.remove <- "Kinases"  #Kinases TFs
FP.raw.df <- gene.attributes.df[, !colnames(gene.attributes.df)%in%df.remove]
colnames(FP.raw.df)[colnames(FP.raw.df)%in%c("TFs", "Kinases")] <- "FPs"
FP.raw.df <- FP.raw.df[FP.raw.df$FPs!="", ]

FP.raw.df <- FP.raw.df[FP.raw.df$Species%in%species7.names, ]
FP.df <- dcast(FP.raw.df, formula = FPs ~ Species, 
               value.var = "Gene", fun.aggregate=function(x) paste(x, collapse = ","))

rownames(FP.df) <- FP.df[,1]
FP.df <- FP.df[, -1]
FP.df <- FP.df[, species7.names]
FP.df.number <- apply(FP.df, 1, function(x){length(x[x!=""])})
table(FP.df.number)
FP.df <- FP.df[FP.df.number==7, ]

FP.number.mat <- t(apply(FP.df, 1, function(x){
  tmp.vec <- vector()
  for(i in 1:length(x)){
    tmp.sp <- colnames(FP.df)[i]
    tmp.genes <- unlist(strsplit(x[i], ","))
    tmp.vec <- c(tmp.vec, length(tmp.genes))
  }
  tmp.vec
}))
colnames(FP.number.mat) <- colnames(FP.df)
FP.index <- apply(FP.number.mat, 1, min) > 5
table(FP.index)
### TF 44  Kinase 38
FP.df <- FP.df[FP.index, ]
FP.number.mat <- FP.number.mat[FP.index, ]

FP.ratio.mat <- t(apply(FP.df, 1, function(x){
  tmp.vec <- vector()
  for(i in 1:length(x)){
    tmp.sp <- colnames(FP.df)[i]
    tmp.genes <- unlist(strsplit(x[i], ","))
    tmp.genes.m6A <- tmp.genes[tmp.genes%in%m6A.list[[tmp.sp]]]
    if(length(tmp.genes)>0){
        tmp.vec <- c(tmp.vec, length(tmp.genes.m6A)/length(tmp.genes)*100) # 2021-08-01
    }else{
      tmp.vec <- c(tmp.vec, NA)
    }
  }
  tmp.vec
}))
colnames(FP.ratio.mat) <- colnames(FP.df)

clu.met <- "ward.D"
FP.heatmap <- pheatmap(FP.ratio.mat, silent = T, clustering_method = clu.met,
                       color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),
                       cluster_cols = F, border_color = NA)
FP.names <- FP.heatmap$tree_row$labels[FP.heatmap$tree_row$order]
FP.cluster <- cutree(FP.heatmap$tree_row,k=2)[FP.names]
cluster1 <- names(FP.cluster)[FP.cluster == 1]
cluster2 <- names(FP.cluster)[FP.cluster == 2]
cat(length(cluster1), length(cluster2), "\n")

if (df.remove == "Kinases"){
  rev.level <- rev(FP.names[c(20:1, 21:44)])
}else{ 
  rev.level <- rev(FP.names[c(26:19, 27:38, 1:18)])
}

# family members
FP.count <- melt(FP.number.mat[FP.names, ])
wilcox.test(FP.count[FP.count[,1]%in%cluster1, "value"], 
       FP.count[FP.count[,1]%in%cluster2, "value"])
FP.count[,1] <- factor(FP.count[,1], levels = rev.level)
FP.count.plot <- ggplot(FP.count, aes(x=Var1, y=value)) +
  geom_boxplot(outlier.colour = NA) + 
  theme_bw(base_size = 14) +
  BaseTheme() +
  xlab("") + ylab("Count") +
  coord_flip(ylim = c(0, 240))

# cv values of family members
FP.cv <- data.frame("Name" = FP.names, "cv" = apply(FP.number.mat[FP.names, ], 1, cv))
# 
FP.cv$Name <- factor(FP.cv$Name, levels = rev.level)
FP.cv.plot <- ggplot(data = FP.cv, mapping = aes(x = Name, y = cv)) + 
  geom_bar(stat = 'identity') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y="CV") + coord_flip()

# cv values of m6A ratio
FP.m6A.cv <- data.frame("Name" = FP.names, "cv" = apply(FP.ratio.mat[FP.names, ], 1, cv))
FP.m6A.cv$Name <- factor(FP.m6A.cv$Name, levels = rev.level)
FP.m6A.cv.plot <- ggplot(data = FP.m6A.cv, mapping = aes(x = Name, y = cv)) + 
  geom_bar(stat = 'identity') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y="CV") + coord_flip()
cor.test(FP.cv$cv, FP.m6A.cv$cv, method = "spearman")

# pdf("figures/Figure5E.pdf", height = 11, width = 22) #12 16
# print(as.ggplot(pheatmap::pheatmap(FP.ratio.mat, clustering_method = clu.met,
#                              # clustering_distance_rows = "maximum",
#                              silent = T, cluster_rows = T,
#                              color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),
#                              cluster_cols = F, #annotation_row = annotation_row,
#                              border_color = NA)) +
#   FP.count.plot + FP.cv.plot + FP.m6A.cv.plot +
#   plot_layout(nrow = 1, guides = "collect", widths = c(2,0.9,0.9,0.9)))
# dev.off()

pdf("figures/FigureS18.pdf", height = 9, width = 18) #12 16
print(as.ggplot(pheatmap::pheatmap(FP.ratio.mat, clustering_method = clu.met,
                             # clustering_distance_rows = "maximum",
                             silent = T, cluster_rows = T,
                             color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),
                             cluster_cols = F, #annotation_row = annotation_row,
                             border_color = NA)) +
  FP.count.plot + FP.cv.plot + FP.m6A.cv.plot +
  plot_layout(nrow = 1, guides = "collect", widths = c(2,0.9,0.9,0.9)))
dev.off()


# Figure5F R genes ####
input.df <- read.table(paste0(link.data.path, "figure_data/Fig5F.csv"), sep = ",", header = T)
input.df$m6ARatio <- input.df$m6ARatio*100
pdf("figures/Figure5F.pdf", width = 4.3, height = 4, useDingbats = F)
CCvaluePlot(input.df, "count", "m6ARatio", c(400, 35, 400, 32), 
            "Number of R gene family members", comlab.name)
dev.off()

# Supplementary Figure20 R genes ####
TNL.df <- read.table(paste0(link.data.path, "figure_data/Sup_R_TNL.csv"), sep = ",", header = T)
CNL.df <- read.table(paste0(link.data.path, "figure_data/Sup_R_CNL.csv"), sep = ",", header = T)
TNL.df$m6ARatio <- TNL.df$m6ARatio*100
CNL.df$m6ARatio <- CNL.df$m6ARatio*100
TNL.CC.plot <- CCvaluePlot(TNL.df, "count", "m6ARatio", c(60, 35, 60, 32), 
                           "Number of TNL genes", comlab.name)
CNL.CC.plot <- CCvaluePlot(CNL.df, "count", "m6ARatio", c(400, 35, 400, 32), 
                           "Number of CNL genes", comlab.name)
pdf("figures/FigureS20.pdf", width = 10, height = 5, useDingbats = F)
TNL.CC.plot + CNL.CC.plot
dev.off()
