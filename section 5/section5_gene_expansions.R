## Significantly expanded gene families #####

## current wroking directory
setwd("~/a2z/m6A_13spp/")

## Options Settings
options(stringsAsFactors = F)
options(scipen=200)

## loading libraries
library(raster)
library(ggtree)
library(pheatmap)
library(raster)
library(ggplotify)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## reload data
load("notebooks/RData/01_m6A_peaks_genes.RData")
load("notebooks/RData/02_detailed_gene_attributes.RData")


## Computational Analysis of gene Family Evolution
SPECIES8 <- c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa", "ppa")
ec_ps_names <- c("Chlorophyta", "Chloroplastida", "Embryophyta", "Stomatophyta", 
                 "Tracheophyta", "Spermatophyta", "Acrogymnospermae",
                 "Angiosperm", "Mesangiospermae", "Monocots", "Commelinids", 
                 "Poaceae", "BOP", "Pooideae", "Hordeeae", "PACMAD", 
                 "Paniceae", "Setaria","Eudicots", "Solanaceae", "Rosids", 
                 "Eurosids", "Malvids", "subMalvids", "Brassicales", 
                 "Brassicaceae", "Fabids", "Fabaceae")
spp_names <- colnames(OGs_34spp)[1:34]

ec_ps_count <- matrix(nrow = length(c(ec_ps_names, spp_names)) , 
                      ncol = length(ec_ps_names))
colnames(ec_ps_count) <- ec_ps_names
rownames(ec_ps_count) <- c(ec_ps_names, spp_names)

for(tmp_ps in ec_ps_names){
  tmp_exp_con <- read.table(paste0("associated_data/cafe_result/",
                                   tmp_ps, "/", tmp_ps, "_all_fams.txt"), sep = "\t", 
                            row.names = 1, stringsAsFactors = F, comment.char = "")
  cat("\n### ", tmp_ps, " ###\n", tmp_exp_con[1,1], "\n")
}

name_change <- list(
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

ec_ps_all_add <- ec_ps_all_minus <- ec_ps_rapid_add <- ec_ps_rapid_minus <- ec_spp_output <- list()

for(tmp_ps in ec_ps_names){
  tmp_exp_con <- read.table(paste0("associated_data/cafe_result/",
                                   tmp_ps, "/", tmp_ps, "_all_fams.txt"), sep = "\t", 
                            row.names = 1, stringsAsFactors = F, comment.char = "")
  tmp_exp_con$OGs_add <- apply(tmp_exp_con, 1, function(x){
    paste0(unlist(strsplit( x[1], ","))[grepl(pattern = "\\+", x = unlist(strsplit(x[1], ",")))], collapse = ",")
  })
  tmp_exp_con$OGs_minus <- apply(tmp_exp_con, 1, function(x){
    paste0(unlist(strsplit( x[1], ","))[grepl(pattern = "\\-", x = unlist(strsplit(x[1], ",")))], collapse = ",")
  })
  tmp_exp_con$OGs_add <- gsub("\\[.+?\\]", "", tmp_exp_con$OGs_add)
  tmp_exp_con$OGs_minus <- gsub("\\[.+?\\]", "", tmp_exp_con$OGs_minus)
  tmp_exp_con$Rapid <- apply(tmp_exp_con, 1, function(x){
    paste0(unlist(strsplit( x[1], ","))[grepl(pattern = "\\*", x = unlist(strsplit( x[1], ",")))], collapse = ",")
  })
  tmp_exp_con$Rapid_add <- apply(tmp_exp_con, 1, function(x){
    paste0(unlist(strsplit( x[4], ","))[grepl(pattern = "\\+", x = unlist(strsplit( x[4], ",")))], collapse = ",")
  })
  tmp_exp_con$Rapid_minus <- apply(tmp_exp_con, 1, function(x){
    paste0(unlist(strsplit( x[4], ","))[grepl(pattern = "\\-", x = unlist(strsplit( x[4], ",")))], collapse = ",")
  })
  tmp_exp_con$Rapid_add <- gsub("\\[.+?\\]", "", tmp_exp_con$Rapid_add)
  tmp_exp_con$Rapid_minus <- gsub("\\[.+?\\]", "", tmp_exp_con$Rapid_minus)
  tmp_exp_con <- tmp_exp_con[3:nrow(tmp_exp_con), ]
  ##
  tmp_node_names <- name_change[[tmp_ps]]$Nodes
  cat(tmp_ps, length(tmp_node_names), sum(tmp_node_names%in%rownames(tmp_exp_con))==length(tmp_node_names), "\n")
  if(all(tmp_node_names%in%rownames(tmp_exp_con))){
    tmp_exp_con <- tmp_exp_con[intersect(tmp_node_names, rownames(tmp_exp_con)), ]
    rownames(tmp_exp_con) <- name_change[[tmp_ps]]$Names
  }
  ## ec_ps_all_add <- ec_ps_all_minus <-  ec_ps_rapid_add <- ec_ps_rapid_minus <- ec_spp_output 
  if(any(tmp_exp_con$OGs_add!="")){
    tmp_all_add <- tmp_exp_con[tmp_exp_con$OGs_add!="", ]
    # tmp_count <- sapply(strsplit(tmp_all_add$OGs_add, ","), length)
    # ec_ps_count[rownames(tmp_all_add), tmp_ps] <- tmp_count
    for(ii in rownames(tmp_all_add)){
      if(ii %in% names(ec_ps_all_add)){
        ec_ps_all_add[[ii]] <- unique(c(ec_ps_all_add[[ii]], unlist(strsplit(tmp_all_add[ii, "OGs_add"], split = ","))))
      } else {
        ec_ps_all_add[[ii]] <- unlist(strsplit(tmp_all_add[ii, "OGs_add"], split = ","))
      }
    }
  }
  ##
  if(any(tmp_exp_con$OGs_minus!="")){
    tmp_all_minus <- tmp_exp_con[tmp_exp_con$OGs_minus!="", ]
    # tmp_count <- sapply(strsplit(tmp_all_minus$OGs_minus, ","), length)
    # ec_ps_count[rownames(tmp_all_minus), tmp_ps] <- tmp_count
    for(ii in rownames(tmp_all_minus)){
      if(ii %in% names(ec_ps_all_minus)){
        ec_ps_all_minus[[ii]] <- unique(c(ec_ps_all_minus[[ii]], unlist(strsplit(tmp_all_minus[ii, "OGs_minus"], split = ","))))
      } else {
        ec_ps_all_minus[[ii]] <- unlist(strsplit(tmp_all_minus[ii, "OGs_minus"], split = ","))
      }
    }
  }
  ##
  if(any(tmp_exp_con$Rapid_add!="")){
    tmp_rapid_add <- tmp_exp_con[tmp_exp_con$Rapid_add!="", ]
    tmp_count <- sapply(strsplit(tmp_rapid_add$Rapid_add, ","), length)
    ec_ps_count[rownames(tmp_rapid_add), tmp_ps] <- tmp_count
    for(ii in rownames(tmp_rapid_add)){
      if(ii %in% names(ec_ps_rapid_add)){
        ec_ps_rapid_add[[ii]] <- unique(c(ec_ps_rapid_add[[ii]], unlist(strsplit(tmp_rapid_add[ii, "Rapid_add"], split = ","))))
      } else {
        ec_ps_rapid_add[[ii]] <- unlist(strsplit(tmp_rapid_add[ii, "Rapid_add"], split = ","))
      }
    }
  }
  ##
  if(any(tmp_exp_con$Rapid_minus!="")){
    tmp_rapid_minus <- tmp_exp_con[tmp_exp_con$Rapid_minus!="", ]
    # tmp_rapid_minustmp_count <- sapply(strsplit(tmp_rapid_minus$Rapid_minus, ","), length)
    # ec_ps_count[rownames(tmp_rapid_minus), tmp_ps] <- tmp_count
    for(ii in rownames(tmp_rapid_minus)){
      if(ii %in% names(ec_ps_rapid_minus)){
        ec_ps_rapid_minus[[ii]] <- unique(c(ec_ps_rapid_minus[[ii]], unlist(strsplit(tmp_rapid_minus[ii, "Rapid_minus"], split = ","))))
      } else {
        ec_ps_rapid_minus[[ii]] <- unlist(strsplit(tmp_rapid_minus[ii, "Rapid_minus"], split = ","))
      }
    }
  }
}


species_path <- list(
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
species_all_add <- species_all_minus <- species_rapid_add <- species_rapid_minus <- list()
for(ii in names(species_path)){
  ##
  tmp_overlap <- intersect(species_path[[ii]], names(ec_ps_all_add))
  species_all_add[[ii]] <- unique(unlist(ec_ps_all_add[tmp_overlap]))
  ##
  tmp_overlap <- intersect(species_path[[ii]], names(ec_ps_all_minus))
  species_all_minus[[ii]] <- unique(unlist(ec_ps_all_minus[tmp_overlap]))
  ##
  tmp_overlap <- intersect(species_path[[ii]], names(ec_ps_rapid_add))
  species_rapid_add[[ii]] <- unique(unlist(ec_ps_rapid_add[tmp_overlap]))
  ##
  tmp_overlap <- intersect(species_path[[ii]], names(ec_ps_rapid_minus))
  species_rapid_minus[[ii]] <- unique(unlist(ec_ps_rapid_minus[tmp_overlap]))
}

## merging node
node_order <- c("Ath","Bstricta","Gar","Egrandis","Pvu","Mtruncatula",
                "Ptrichocarpa","Vvinifera","Stuberosum","Sly","Sitalica",
                "Sviridis","Sbi","Ata","Hvulgare","Bdistachyon","Osa",
                "Eguineensis","Spolyrhiza","Atrichopoda","Paab","Gibi",
                "Smoellendorffii","Ppa","Mpolymorpha","Creinhardtii",
                "MpusillaCCMP1545","Embryophyta","Stomatophyta","Tracheophyta",
                "Spermatophyta","Acrogymnospermae","Angiosperm","Mesangiospermae",
                "Monocots","Commelinids","Poaceae","Pooideae","Hordeeae","PACMAD",
                "Setaria","Eudicots","Rosids","Eurosids","SubMalvids","Brassicaceae",
                "Fabaceae","Solanaceae")
tmp_ec_ps_all_add <- ec_ps_all_add
tmp_ec_ps_rapid_add <- ec_ps_rapid_add
names(tmp_ec_ps_all_add) <- firstup(names(tmp_ec_ps_all_add))
names(tmp_ec_ps_rapid_add) <- firstup(names(tmp_ec_ps_rapid_add))
all(node_order%in%names(tmp_ec_ps_all_add))&all(names(tmp_ec_ps_all_add)%in%node_order)
all(node_order%in%names(tmp_ec_ps_rapid_add))&all(names(tmp_ec_ps_rapid_add)%in%node_order)
tmp_ec_ps_all_add <- tmp_ec_ps_all_add[node_order]
tmp_ec_ps_rapid_add <- tmp_ec_ps_rapid_add[node_order]

all_add_OGs <- unique(unlist(tmp_ec_ps_all_add))
node_num <- length(tmp_ec_ps_all_add)
out_data <- matrix(nrow = length(all_add_OGs), 
                   ncol = node_num*2)
rownames(out_data) <- all_add_OGs
colnames(out_data) <- rep(names(tmp_ec_ps_all_add), each=2)
colnames(out_data)[(1:node_num)*2] <- "NULL"

dim(out_data)
for(ii in all_add_OGs){
  for(jj in 1:node_num){
    if(ii%in%tmp_ec_ps_all_add[[jj]]){
      out_data[ii, jj*2-1] <- "T"
    }
    if(ii%in%tmp_ec_ps_rapid_add[[jj]]){
      out_data[ii, jj*2] <- "T"
    }
  }
}


significant_expansion <- vector()
for(ii in names(tmp_ec_ps_all_add)){
  tmp_vec <- rep("No", length(tmp_ec_ps_all_add[[ii]]))
  tmp_vec[tmp_ec_ps_all_add[[ii]]%in%tmp_ec_ps_rapid_add[[ii]]] <- "Yes"
  significant_expansion <- c(significant_expansion, tmp_vec)
}

ec_ps_add_output <- data.frame("Expanded OGs"=unlist(tmp_ec_ps_all_add),
                               "Node"=rep(names(tmp_ec_ps_all_add),
                                          sapply(tmp_ec_ps_all_add, length)),
                               "Significant expansion"=significant_expansion)

SPECIES27 <- SPECIES34[!SPECIES34%in%c("Thassleriana", "ghi", "gma", "Pvirgatum",
                                       "zma", "tdi", "tae")]
ec_ps_add_output <- cbind(ec_ps_add_output, OGs_34spp_count[unlist(tmp_ec_ps_all_add), SPECIES27])
colnames(ec_ps_add_output) <- firstup(colnames(ec_ps_add_output))

ec_ps_add_output[ec_ps_add_output$Node=="SubMalvids", "Node"] <- "BM"

write.table(out_data, file = "Table5.txt", quote = F, sep = "\t", row.names = T)


## each species
ps_add_minus <- ec_ps_count[,1:4]
colnames(ps_add_minus) <- c("all_add", "all_minus", "rapid_add", "rapid_minus")
tmp_vec <- sapply(ec_ps_all_add, length)
tmp_overlap <- intersect(rownames(ps_add_minus), names(tmp_vec))
ps_add_minus[tmp_overlap, "all_add"] <- tmp_vec[tmp_overlap]
tmp_vec <- sapply(ec_ps_all_minus, length)
tmp_overlap <- intersect(rownames(ps_add_minus), names(tmp_vec))
ps_add_minus[tmp_overlap, "all_minus"] <- tmp_vec[tmp_overlap]
tmp_vec <- sapply(ec_ps_rapid_add, length)
tmp_overlap <- intersect(rownames(ps_add_minus), names(tmp_vec))
ps_add_minus[tmp_overlap, "rapid_add"] <- tmp_vec[tmp_overlap]
tmp_vec <- sapply(ec_ps_rapid_minus, length)
tmp_overlap <- intersect(rownames(ps_add_minus), names(tmp_vec))
ps_add_minus[tmp_overlap, "rapid_minus"] <- tmp_vec[tmp_overlap]
ps_add_minus <- ps_add_minus[,c(1,3,2,4)]


## Figure5A ####
### overlapped gene families
overlap_OGs <- intersect(species_rapid_add[["ath"]], 
                         intersect(species_rapid_add[["gar"]], 
                                   intersect(species_rapid_add[["pvu"]],
                                             intersect(species_rapid_add[["sly"]],
                                                       intersect(species_rapid_add[["sbi"]],
                                                                 intersect(species_rapid_add[["ata"]],
                                                                           species_rapid_add[["osa"]]))))))
m6Aratio_overlap_OGs <- vector()
for(tmp_spp in c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa")){ #, "ppa"
  for(tmp_ogs in overlap_OGs){
    cat(tmp_spp, tmp_ogs, "\n")
    tmp_ps_genes <- unlist(strsplit(OGs_34spp[tmp_ogs, tmp_spp], split = ", "))
    tmp_ps_expgenes <- sum(tmp_ps_genes%in%unlist(exp_list[tmp_spp]))
    tmp_ps_m6Agenes <- sum(tmp_ps_genes%in%unlist(m6A_list[tmp_spp]))
    m6Aratio_overlap_OGs <- c(m6Aratio_overlap_OGs, 
                              tmp_ps_m6Agenes/tmp_ps_expgenes) ### 2021-08-02 length(tmp_ps_genes))
  }
}

m6Aratio_overlap_OGs_df <- data.frame(
  "Species"= rep(c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa"), #, "ppa"
                 each=length(overlap_OGs)),
  "Values"= m6Aratio_overlap_OGs, 
  "OGs"=rep(overlap_OGs, 7))

m6Aratio_overlap_OGs_df$Species <- factor(m6Aratio_overlap_OGs_df$Species, 
                                          levels = intersect(c("ath", "gar", "pvu", "sly", "sbi", 
                                                               "ata", "osa", "ppa"), 
                                                             unique(m6Aratio_overlap_OGs_df$Species)))

ggplot(m6Aratio_overlap_OGs_df,aes(OGs,Values,group=Species,colour=Species))+
  geom_line(size=1,linetype=1)+
  scale_color_manual(values = c(rep("red", 4), rep("green", 3))) +
  geom_point()

m6Aratio_overlap_OGs_data <- reshape2::dcast(m6Aratio_overlap_OGs_df, formula = OGs~Species, value.var="Values")
rownames(m6Aratio_overlap_OGs_data) <- m6Aratio_overlap_OGs_data[,1]
m6Aratio_overlap_OGs_data <- m6Aratio_overlap_OGs_data[, -1]
m6Aratio_overlap_OGs_data[is.na(m6Aratio_overlap_OGs_data)] <- 0
m6Aratio_overlap_OGs_data <- m6Aratio_overlap_OGs_data[, c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa")]
m6Aratio_overlap_cor <- cor(m6Aratio_overlap_OGs_data)
m6Aratio_overlap_cor[m6Aratio_overlap_cor==1] <- 0

# pdf("Figure26.pdf", height = 4, width = 4)
corrplot(corr = m6Aratio_overlap_cor, cl.lim = c(0,1), cl.length = 6, #add=TRUE, 
         col = colorRampPalette(c("black","green",
                                  "white","#FB8D42", 
                                  "#7E0527"))(100), 
         type="lower", method="color",
         order="original", addCoef.col = "white", #"#EE7E32", "#E30713"
         tl.pos="td")
# dev.off()


## Figure 5B ####
ref_value <- c(12224/18556, 15075/26298, 11935/19807, 
               11000/18528, 13062/20212, 13823/18904,
               15325/24399)*100
ref_value_list <- list()
ref_value_list[['ath']] <- c(12224, 18556-12224)
ref_value_list[['gar']] <- c(15075, 26298-15075)
ref_value_list[['pvu']] <- c(11935, 19807-11935)
ref_value_list[['sly']] <- c(11000, 18528-11000)
ref_value_list[['sbi']] <- c(13062, 20212-13062)
ref_value_list[['ata']] <- c(13823, 18904-13823)
ref_value_list[['osa']] <- c(15325, 24399-15325)
# ref_value_list[['ppa']] <- c(16483, 19251-16483)

m6A_rapid_ratio <- vector()
for(tmp_species in SPECIES7){
  tmp_OGs <- OGs_34spp[species_rapid_add[[tmp_species]], tmp_species]
  tmp_OG_genes <- unlist(strsplit(tmp_OGs, split = ", "))
  tmp_exp <- sum(tmp_OG_genes%in%exp_list[[tmp_species]])
  tmp_m6A <- sum(tmp_OG_genes%in%m6A_list[[tmp_species]])
  tmp_m6A_ratio <- tmp_m6A/tmp_exp*100  ## 2021-08-02 length(tmp_OG_genes)
  cat(tmp_species, tmp_m6A, tmp_exp-tmp_m6A, "\n") ### length(tmp_OG_genes)
  cat(tmp_species, chisq.test(matrix(c(tmp_m6A,
                                       tmp_exp-tmp_m6A, ### length(tmp_OG_genes)
                                       ref_value_list[[tmp_species]]),
                                 nrow = 2))$p.value, "\n")
  m6A_rapid_ratio <- c(m6A_rapid_ratio, tmp_m6A_ratio)
}

df_m6A_com <- data.frame("species"=rep(firstup(SPECIES7), 4), 
                         "value"=c(ref_value, 100-ref_value, m6A_rapid_ratio, 100-m6A_rapid_ratio),
                         "Type"=rep(c("ref", "rapid"), each=14),
                         "Type2"=rep(c("ref-m6A", "ref-Non-m6A","rapid-m6A", "rapid-Non-m6A"), each=7))
df_m6A_com$species <- factor(df_m6A_com$species, levels = unique(df_m6A_com$species))
df_m6A_com$Type2 <- factor(df_m6A_com$Type2, levels = rev(unique(df_m6A_com$Type2)))

# pdf("Figure24.pdf", height = 3, width = 7)
ggplot(data = df_m6A_com, aes(x=Type, y=value, fill=Type2))+
  geom_bar(stat='identity', width=0.7) +  #
  theme_bw(base_size = 14) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab(expression("m"^"6"*"A level (%)")) +
  scale_fill_manual("Type", values = c("#CCCCCC","#557773", "#CCCCCC", "#A68171")) + 
  facet_grid(.~species) + coord_cartesian(ylim = c(0, 80))
# dev.off()



## Fgiure 5C and Supplementary Fig. 9 ####
### eight species with six nodes
species_path <- list(
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

Node1 <- c("Embryophyta", "Stomatophyta")
Node1_spp <- c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa", "ppa")
Node2 <- c("Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae")
Node2_spp <- c("ath", "gar", "pvu", "sly", "sbi", "ata", "osa")
Node3 <- c("Eudicots")
Node3_spp <- c("ath", "gar", "pvu", "sly")
Node4 <- c("Rosids", "Eurosids")
Node4_spp <- c("ath", "gar", "pvu")
Node5 <- c("Malvids", "subMalvids")
Node5_spp <- c("ath", "gar")
Node6 <- c("Brassicales", "Brassicaceae", "ath")
Node6_spp <- c("ath")
Node7 <- c("gar")
Node7_spp <- c("gar")
Node8 <- c("Fabids", "Fabaceae", "pvu")
Node8_spp <- c("pvu")
Node9 <- c("Solanaceae", "sly")
Node9_spp <- c("sly")
Node10 <- c("Monocots", "Commelinids", "Poaceae")
Node10_spp <- c("sbi", "ata", "osa")
Node11 <- c("PACMAD", "sbi")
Node11_spp <- c("sbi")
# Node12 <- c("BOP")
# Node12_spp <- c("ata", "osa")
Node13 <- c("Pooideae", "Hordeeae", "ata")
Node13_spp <- c("ata")
Node14 <- c("osa")
Node14_spp <- c("osa")
Node15 <- c("ppa")
Node15_spp <- c("ppa")

node_each_name <- node_all_names <- node_all_m6Aratio <- node_all_numbers <- vector()
node_m6Aratio <- node_mean_gene <- vector()
for(nodenum in c(1:11, 13:15)){
  tmp_node <- eval(parse(text = paste0("Node", nodenum)))
  tmp_node_spp <- eval(parse(text = paste0("Node", nodenum, "_spp")))
  raw_num <- raw_exp <- raw_m6A <- raw_mean_num <- vector()
  for(tmp_ps in tmp_node){
    if(tmp_ps %in% names(ec_ps_rapid_add)){
      for(tmp_spp in tmp_node_spp){
        # cat(length(ec_ps_rapid_add[[tmp_ps]]), "\n")
        tmp_gene_list <- strsplit(OGs_34spp[ec_ps_rapid_add[[tmp_ps]], tmp_spp], split = ", ")
        tmp_ps_genes <- unlist(tmp_gene_list)
        raw_num <- c(raw_num, length(tmp_ps_genes)) #sum(tmp_ps_genes%in%exp_list[[tmp_spp]]))
        raw_m6A <- c(raw_m6A, sum(tmp_ps_genes%in%m6A_list[[tmp_spp]]))
        raw_exp <- c(raw_exp, sum(tmp_ps_genes%in%exp_list[[tmp_spp]]))
        raw_mean_num <- c(raw_mean_num, sapply(tmp_gene_list, length))
        ## detailed information
        node_each_name <- c(node_each_name, rep(paste0("Node", nodenum), length(tmp_gene_list)))
        node_all_names <- c(node_all_names, rep(tmp_ps, length(tmp_gene_list)))
        node_all_m6Aratio <- c(node_all_m6Aratio, sapply(tmp_gene_list, function(x){
          sum(x%in%m6A_list[[tmp_spp]])/sum(x%in%exp_list[[tmp_spp]]) ## 2021-07-31 length(x)*100
        }))
        node_all_numbers <- c(node_all_numbers, sapply(tmp_gene_list, length))
      }      
    }
  }
  node_m6Aratio <- c(node_m6Aratio, sum(raw_m6A)/sum(raw_exp)*100) ## 2021-08-02 sum(raw_num)
  node_mean_gene <- c(node_mean_gene, mean(unlist(raw_mean_num)))
  cat(tmp_node, round(sum(raw_m6A)/sum(raw_exp)*100,digits = 2),   ## 2021-08-02 sum(raw_num)
      round(mean(unlist(raw_mean_num)), digits = 2), "\n")
}

ratio2count_df <- data.frame("Nodes" = node_each_name,
                             "PSs"=node_all_names,
                             "ratios"=node_all_m6Aratio,
                             "count"=node_all_numbers)

## Node1 Node2 Node3 Node4 Node5 Node10
ratio2count_df <- ratio2count_df[ratio2count_df$count>=10, ]
ratio2count_df <- ratio2count_df[ratio2count_df$ratios!=0, ]
ratio2count_df <- ratio2count_df[ratio2count_df$Nodes%in%
                                   c("Node1", "Node2", "Node3", "Node4", 
                                     "Node5", "Node10"), ]
ratio2count_df$Nodes <- factor(ratio2count_df$Nodes, 
                               levels = unique(ratio2count_df$Nodes))

# pdf("Figure_internode.pdf", height = 3, width = 12)
ggplot(data = ratio2count_df, mapping = aes(x = count, y = ratios)) +
  geom_point() + facet_grid(.~Nodes, scales = "free") + #space = "free"
  labs(x = '', y = '', title = '') +
  theme_bw(base_size = 14) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        # axis.text.x = element_text(colour = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_smooth(method="lm", se=F, color="#DA2416")
# dev.off()

for(tmp_PSs in unique(ratio2count_df$Nodes)){
  ratio2count_filter <- ratio2count_df[ratio2count_df$Nodes==tmp_PSs,]
  ratio2count_cor <- cor.test(ratio2count_filter$ratios, 
                              ratio2count_filter$count, 
                              method = "pearson")
  cat(tmp_PSs, ratio2count_cor$estimate,
      ratio2count_cor$p.value, "\n")
}



## Figure 5D  ####
ratio2ratio <- list()
for(tmp_species in SPECIES8){
  tmp_df <- t(apply(data.frame(species_rapid_add[[tmp_species]]), 1, function(x){
    gene_list <- unlist(strsplit(OGs_34spp[x[1], tmp_species], ", "))
    gene_num <- length(gene_list)
    exp_num <- sum(gene_list%in%exp_list[[tmp_species]])
    m6A_num <- sum(gene_list%in%m6A_list[[tmp_species]])
    return(c(tmp_species, x[1], gene_num, exp_num, 
             m6A_num, m6A_num/exp_num*100)) ## 2021-07-28 m6A_num/gene_num*100
  }))
  tmp_df <- tmp_df[which(as.numeric(tmp_df[,6])>0&as.numeric(tmp_df[,3])>10), ]
  tmp_cor_result <- cor.test(as.numeric(tmp_df[,3]), as.numeric(tmp_df[,6]))
  cat(tmp_species, tmp_cor_result$estimate, tmp_cor_result$p.value, "\n")
  ratio2ratio[[tmp_species]] <- tmp_df 
}

ratio2ratio_df <- as.data.frame(do.call(rbind, ratio2ratio), stringsAsFactors = F)
colnames(ratio2ratio_df) <- c("Species", "OGs", "GeneNumber", "ExpNumber", "m6ANumber", "m6ARatio")
ratio2ratio_df[,3] <- as.numeric(ratio2ratio_df[,3])
ratio2ratio_df[,6] <- as.numeric(ratio2ratio_df[,6])
# write.table(ratio2ratio_df[,c(1,2,3,6)], "Figure5C_data.txt", quote = F, row.names = F)

ratio2ratio_df$Species <- factor(ratio2ratio_df$Species, levels = SPECIES8)
# pdf("Figure25_raw.pdf", height = 3, width = 10)
ggplot(data = ratio2ratio_df, mapping = aes(x = GeneNumber, y = m6ARatio)) +
  geom_point() + facet_grid(.~Species, scales = "free", 
                            space = "free") +
  labs(x = '', y = '', title = '') +
  theme_bw(base_size = 14) +
  theme(axis.line = element_line(colour = "black", size = 0.5),
        # axis.text.x = element_text(colour = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_smooth(method="lm", se=F, color="#A3A2A2", span = 2) # lm
# dev.off()

for(ii in SPECIES8){
  tmp_df <- ratio2ratio_df[ratio2ratio_df$Species==ii, ]
  fit_stat <- summary(lm(as.numeric(tmp_df$GeneNumber)~as.numeric(tmp_df$m6ARatio)))
  R2_adj <- fit_stat$adj.r.squared
  p_value <- fit_stat$coefficients[2,4]
  cat(ii, R2_adj, p_value, cor(as.numeric(tmp_df$GeneNumber),as.numeric(tmp_df$m6ARatio)),
      cor.test(as.numeric(tmp_df$GeneNumber),as.numeric(tmp_df$m6ARatio))$p.value, "\n")
}


# Figure5E Supplementary Fig. 10 Kinases and Transciption factor #####

df_remove <- "TFs"  #Kinases TFs
FP_raw_df <- gene_attributes_df[, !colnames(gene_attributes_df)%in%df_remove]
colnames(FP_raw_df)[colnames(FP_raw_df)%in%c("TFs", "Kinases")] <- "FPs"
FP_raw_df <- FP_raw_df[FP_raw_df$FPs!="", ] #&FP_raw_df$Exp!=""

FP_raw_df <- FP_raw_df[FP_raw_df$Species%in%SPECIES7, ]
FP_df <- dcast(FP_raw_df, formula = FPs~Species, 
               value.var = "Gene", fun.aggregate=function(x) paste(x, collapse = ","))

rownames(FP_df) <- FP_df[,1]
FP_df <- FP_df[, -1]
FP_df <- FP_df[, SPECIES7]
FP_df_number <- apply(FP_df, 1, function(x){length(x[x!=""])})
table(FP_df_number)
FP_df <- FP_df[FP_df_number==7, ]

FP_number_mat <- t(apply(FP_df, 1, function(x){
  tmp_vec <- vector()
  for(i in 1:length(x)){
    tmp_species <- colnames(FP_df)[i]
    tmp_genes <- unlist(strsplit(x[i], ","))
    # tmp_genes_exp <- tmp_genes[tmp_genes%in%exp_list[[tmp_species]]]
    tmp_vec <- c(tmp_vec, length(tmp_genes))
  }
  tmp_vec
}))
colnames(FP_number_mat) <- colnames(FP_df)
FP_index <- apply(FP_number_mat, 1, min) > 5
table(FP_index)
### TF 44  Kinase 38
FP_df <- FP_df[FP_index, ]
FP_number_mat <- FP_number_mat[FP_index, ]

FP_ratio_mat <- t(apply(FP_df, 1, function(x){
  tmp_vec <- vector()
  for(i in 1:length(x)){
    tmp_species <- colnames(FP_df)[i]
    tmp_genes <- unlist(strsplit(x[i], ","))
    tmp_genes_exp <- tmp_genes[tmp_genes%in%exp_list[[tmp_species]]]
    tmp_genes_m6A <- tmp_genes[tmp_genes%in%m6A_list[[tmp_species]]]
    if(length(tmp_genes)>0){
        tmp_vec <- c(tmp_vec, length(tmp_genes_m6A)/length(tmp_genes_exp)*100) ## 2021-08-02 length(tmp_genes)
    }else{
      tmp_vec <- c(tmp_vec, NA)
    }
  }
  tmp_vec
}))

colnames(FP_ratio_mat) <- colnames(FP_df)

if(df_remove == "Kinases"){clu_met <- "ward.D2"}else{clu_met <- "ward.D"}
FP_heatmap <- pheatmap(FP_ratio_mat, silent = T, clustering_method = clu_met,
                       color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),
                       cluster_cols = F, border_color = NA)

FP_names <- FP_heatmap$tree_row$labels[FP_heatmap$tree_row$order]
FP_cluster <- cutree(FP_heatmap$tree_row,k=2)[FP_names]
cluster1 <- names(FP_cluster)[FP_cluster==1]
cluster2 <- names(FP_cluster)[FP_cluster==2]
cat(length(cluster1), length(cluster2), "\n")

### cv values of family members
FP_cv <- data.frame("Name" = FP_names,
                    "cv" = apply(FP_number_mat[FP_names, ], 1, cv))
wilcox.test(FP_cv[FP_cv$Name%in%cluster1, "cv"], 
            FP_cv[FP_cv$Name%in%cluster2, "cv"])

if(df_remove == "Kinases"){rev_level <- rev(FP_names)}else{rev_level <- FP_names}
FP_cv$Name <- factor(FP_cv$Name, levels = rev_level)
FP_cv_plot <- ggplot(data = FP_cv, mapping = aes(x = Name, y = cv)) + 
  geom_bar(stat = 'identity') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y="CV") + coord_flip()

### cv values of m6A ratio
FP_m6A_cv <- data.frame("Name" = FP_names,
                    "cv" = apply(FP_ratio_mat[FP_names, ], 1, cv))
wilcox.test(FP_m6A_cv[FP_m6A_cv$Name%in%cluster1, "cv"], 
            FP_m6A_cv[FP_m6A_cv$Name%in%cluster2, "cv"])


if(df_remove == "Kinases"){rev_level <- rev(FP_names)}else{rev_level <- FP_names}
FP_m6A_cv$Name <- factor(FP_m6A_cv$Name, levels = rev_level)
FP_m6A_cv_plot <- ggplot(data = FP_m6A_cv, mapping = aes(x = Name, y = cv)) + 
  geom_bar(stat = 'identity') +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y="CV") + coord_flip()

cor.test(FP_cv$cv, FP_m6A_cv$cv)

##
FP_count <- melt(FP_number_mat[FP_names, ])
wilcox.test(FP_count[FP_count[,1]%in%cluster1, "value"], 
       FP_count[FP_count[,1]%in%cluster2, "value"])

FP_count[,1] <- factor(FP_count[,1], levels = rev_level)
FP_count_plot <- ggplot(FP_count, aes(x=Var1, y=value)) +
  geom_boxplot(outlier.color = "black", outlier.size = 0.8) + 
  theme_bw(base_size = 14) +
  theme(axis.text.x= element_text(angle=15,hjust = 1,vjust = 1, colour="black",size=10),
        axis.line = element_line(colour = "black",size=0.5),
        legend.text=element_text(colour="black", size=16),
        legend.title=element_text(colour="black", size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab("Count") + coord_flip()

# if(df_remove == "Kinases"){fig_num <- 27}else{fig_num <- 28}
if(df_remove == "Kinases"){width_num <- 16}else{width_num <- 22}
# pdf(paste0("Figure",fig_num, ".pdf"), height = 8, width = width_num) #12 16
as.ggplot(pheatmap::pheatmap(FP_ratio_mat, clustering_method = clu_met,
                             # clustering_distance_rows = "maximum",
                             silent = T, cluster_rows = T,
                             color = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(100),
                             cluster_cols = F, #annotation_row = annotation_row,
                             border_color = NA)) + 
  FP_m6A_cv_plot + FP_cv_plot + FP_count_plot + 
  plot_layout(nrow = 1, guides = "collect", widths = c(2,0.9,0.9,0.9))
# dev.off()
