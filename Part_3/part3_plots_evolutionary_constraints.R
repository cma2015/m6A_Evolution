# Evolutionary dynamics of m6A modification between orthologous genes 
#

# current wroking directory
rm(list=ls())
setwd("~/a2z/m6A_13spp/Part_3/")
options(stringsAsFactors = F)

# loading libraries ####
library(ggplot2)
library(ggpubr)
library(patchwork)
library(progress)
library(circlize)
library(RColorBrewer)
library(seqinr)
library(foreach)
library(doParallel)

## reload data
load("../Part_1/data1_m6A_peaks_genes.RData")
load("../Part_1/data2_gene_attributes.RData")

link.data.path <- "../associated_data/"
comlab.name <- "m6A methylation ratio (%)"
system("echo 30 > proc")
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


# Figure3A Species divergenece time and m6A divergence ####
input.3A.data <- read.table(paste0(link.data.path, "figure_data/Fig3A.csv"), 
                            sep = ",", header = T)
input.3A.data[input.3A.data$divergence>1, "divergence"] <- 1
input.3A.data[input.3A.data$type == "Dicots vs Bryophyta", "type"] <- "Bryophyta vs Dicots"
input.3A.data[input.3A.data$type == "Monocots vs Bryophyta", "type"] <- "Bryophyta vs Monocots"
input.3A.data$type <- factor(input.3A.data$type, 
                      levels = c("Monocots vs Monocots", "Dicots vs Dicots", 
                                 "Dicots vs Monocots", "Bryophyta vs Monocots", 
                                 "Bryophyta vs Dicots"))

pdf("figures/Figure3A.pdf", width = 7.2, height = 4.6, useDingbats = F)
ggplot(input.3A.data, mapping = aes(x=time, y=divergence, color=type)) +
  geom_point(size=3.2) + 
  scale_color_manual(values = c("#278D6C", "#C75D19", "#A67BAE",
                                "#BB3372", "#27579F")) +
  theme_bw(base_size = 14)+ 
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Species divergenece time (MYA)") + ylab(comlab.name) + 
  ylim(c(0.2, 1.05)) +
  geom_smooth(method="lm", se=F, color="#A3A2A2") #method="lm", se=T, color="#DA2416"
dev.off()


# Figure3C: plot single-copy OGs in 7 diploid ####
kaks.comb.list <- list(
  c("ath", "gar"), c("ath", "pvu"), c("ath", "sly"),
  c("ath", "sbi"), c("ath", "ata"), c("ath", "osa"),
  c("gar", "pvu"), c("sly", "gar"), c("sbi", "gar"),
  c("gar", "ata"), c("osa", "gar"), c("sly", "pvu"),
  c("sbi", "pvu"), c("pvu", "ata"), c("osa", "pvu"),
  c("sly", "sbi"), c("sly", "ata"), c("sly", "osa"),
  c("sbi", "ata"), c("sbi", "osa"), c("osa", "ata")
)
kaks.comb.names <- sapply(kaks.comb.list, function(x){paste0(x, collapse = "-")})

# pb <- progress_bar$new(total = length(kaks.comb.list), clear = FALSE)
# pb$tick()
plot_list1 <- plot_list2 <- plot_list3 <- list()
evo.cor.values <- uniden.values <- vector()

for (tmp.two.spp in kaks.comb.list){
  tmp.two.spp.name <- paste0(tmp.two.spp, collapse = "-")
  PS.index.names <- c("Chloroplastida", "Embryophyta", 
                      "Spermatophyta", "Angiosperm") # c("Chloroplastida", "Embryophyta", "Stomatophyta", "Tracheophyta", "Spermatophyta", "Angiosperm", "Mesangiospermae")
  tmp.spp.ortho <- read.table(paste0(link.data.path, 
                                     "Results_Jan16_34spp/Orthologues/Orthologues_",
                                     tmp.two.spp[1], "/", tmp.two.spp[1], 
                                     "__v__", tmp.two.spp[2],".tsv"), 
                              header = T, sep = "\t")
  one.divergence <- one.kaks <- list()
  for (PS.index in PS.index.names) {
    tmp.l.name <- paste0(tmp.two.spp.name, "_", PS.index)
    PS.tmp.two.spp <- OGs.13spp[OGs.13spp$PSnames == PS.index, tmp.two.spp]
    speciesone <- unlist(strsplit(PS.tmp.two.spp[,1], split = ", "))
    speciestwo <- unlist(strsplit(PS.tmp.two.spp[,2], split = ", "))
    species.gene.index <- apply(tmp.spp.ortho, 1, function(x){
      tmp.index1 <- all(unlist(strsplit(x[2], split = ", ")) %in% speciesone)
      tmp.index2 <- all(unlist(strsplit(x[3], split = ", ")) %in% speciestwo)
      return(tmp.index1 & tmp.index2)
    })
    tmp.ortho.filter <- tmp.spp.ortho[species.gene.index, 2:3]
    PS.tmp.type.index <- apply(tmp.ortho.filter, 1, function(x){
      tmp1 <- unlist(strsplit(x[1], split = ", "))
      tmp2 <- unlist(strsplit(x[2], split = ", "))
      if (length(tmp1)>1){a <- "m"} else {a <- "1"}
      if (length(tmp2)>1){b <- "m"} else {b <- "1"}
      return(paste0(c(a, b), collapse = ":"))
    })
    tmp.ortho.filter <- tmp.ortho.filter[PS.tmp.type.index=="1:1", ]
    PS.tmp.exp <- apply(tmp.ortho.filter, 1, function(x){
      if(x[1]%in%exp.list[[tmp.two.spp[1]]]){a <- "1"}else{a <- "0"}
      if(x[2]%in%exp.list[[tmp.two.spp[2]]]){b <- "1"}else{b <- "0"}
      return(paste0(c(a, b), collapse = ":"))
    })
    tmp.ortho.filter <- tmp.ortho.filter[PS.tmp.exp=="1:1", ]
    PS.tmp.m6A <- apply(tmp.ortho.filter, 1, function(x){
      if(x[1]%in%m6A.list[[tmp.two.spp[1]]]){a <- "1"}else{a <- "0"}
      if(x[2]%in%m6A.list[[tmp.two.spp[2]]]){b <- "1"}else{b <- "0"}
      return(paste0(c(a, b), collapse = ":"))
    })
    m6A.divergence.index <- PS.tmp.m6A%in%c("1:0", "0:1")
    uniden.part <- sum(m6A.divergence.index)/length(m6A.divergence.index)*100
    iden.part <- 100-uniden.part
    uniden.values <- c(uniden.values, uniden.part)
    one.divergence[[PS.index]] <- data.frame("Species"=tmp.two.spp.name, "PSs"=PS.index,
                                                "Type"=c("Identical pattern", "Diverged pattern"),
                                                "count"=c(sum(!m6A.divergence.index), sum(m6A.divergence.index)),
                                                "value"=c(iden.part, uniden.part))
    if (sum(m6A.divergence.index) > 0 ){
      m6A.divergence.df <- tmp.ortho.filter[m6A.divergence.index, ]
      ## ka ks
      if("run" == "norun"){
        system(paste0("mkdir -p ", link.data.path, "KaKs_m6A/", tmp.l.name))
        write.table(m6A.divergence.df, 
                    paste0(link.data.path, "KaKs_m6A/", tmp.l.name, "/gene_pairs.txt"), 
                    quote = F, sep = "\t", col.names = F, row.names = F)
        write(paste0("cat ", link.data.path, "pep_cds/cds/", tmp.two.spp[1],
                     ".fa ", link.data.path, "pep_cds/cds/", tmp.two.spp[2],
                     ".fa >", link.data.path, "KaKs_m6A/", tmp.l.name, 
                     "/spp_merge_cds.fa"),
              file = "KaKs_m6A_divergence.sh", append = T)
        write(paste0("cat ", link.data.path, "pep_cds/pep/", tmp.two.spp[1],
                     ".fa ", link.data.path, "pep_cds/pep/", tmp.two.spp[2],
                     ".fa >", link.data.path, "KaKs_m6A/", tmp.l.name, 
                     "/spp_merge_pep.fa"),
              file = "KaKs_m6A_divergence.sh", append = T)
        write(paste0("../tools/ParaAT2.0/ParaAT.pl -h ", link.data.path, 
                     "KaKs_m6A/", tmp.l.name,
                     "/gene_pairs.txt -a ", link.data.path, "KaKs_m6A/",
                     tmp.l.name, "/spp_merge_pep.fa ",
                     " -n ", link.data.path, "KaKs_m6A/", tmp.l.name,
                     "/spp_merge_cds.fa -p proc -o ", link.data.path, "KaKs_m6A/",
                     tmp.l.name, "/output -f axt"),
              file = "KaKs_m6A_divergence.sh", append = T)
        write(paste0("ls ", link.data.path, "KaKs_m6A/", tmp.l.name,
                     "/output | grep '.axt' | parallel -I% ../tools/KaKs_Calculator ",
                     "-i ", link.data.path, "KaKs_m6A/", tmp.l.name,
                     "/output/%  -o ", link.data.path, "KaKs_m6A/", tmp.l.name,
                     "/output/%.kaks -m MA"), file = "KaKs_m6A_divergence.sh", append = T)
        write(paste0("cat ", link.data.path, "KaKs_m6A/", tmp.l.name,
                     "/output/*.kaks | awk '{if(!(NR>1&&$1~/Sequence/)){print $0}}'",
                     " >", link.data.path, "KaKs_m6A/", tmp.l.name,
                     "/kaks_", tmp.l.name, "_out.txt"),
              file = "KaKs_m6A_divergence.sh", append = T)
      } else {
        tmp.kaks <- read.table(paste0(link.data.path, "KaKs_m6A/", tmp.l.name, 
                                      "/kaks_", tmp.l.name, "_out.txt"),
                               sep = "\t", header = T, stringsAsFactors = F)
        one.kaks[[PS.index]] <-  data.frame("Species"=tmp.two.spp.name, 
                                              "PSs"=PS.index,
                                              "value"=tmp.kaks$Ka.Ks)
        evo.cor.values <- c(evo.cor.values, median(one.kaks[[PS.index]]$value))
      }
    }
  }
  ## statistics
  test.chi <- function(inputdata, var1, var2){
    chisq.test(matrix(c(inputdata[[var1]][1, 4], inputdata[[var1]][2, 4], 
                        inputdata[[var2]][1, 4], inputdata[[var2]][2, 4]),
                      nrow = 2, byrow = T))$p.value
  }
  test.t <- function(inputdata, var1, var2){
    t.test(inputdata[[var1]]$value, inputdata[[var2]]$value, alternative = "less")$p.value
  }
  write(paste( tmp.two.spp.name, "1diver",
    test.chi(inputdata = one.divergence, var1 = PS.index.names[1], var2 = PS.index.names[2]),
    test.chi(inputdata = one.divergence, var1 = PS.index.names[2], var2 = PS.index.names[3]),
    test.chi(inputdata = one.divergence, var1 = PS.index.names[3], var2 = PS.index.names[4]),
    test.chi(inputdata = one.divergence, var1 = PS.index.names[1], var2 = PS.index.names[3]),
    test.chi(inputdata = one.divergence, var1 = PS.index.names[1], var2 = PS.index.names[4]),
    test.chi(inputdata = one.divergence, var1 = PS.index.names[2], var2 = PS.index.names[4]), "\n", sep = "\t\t"
    ), 
    file = "figures/figure3C_pvalue.txt", append = T)
  write(paste( tmp.two.spp.name, "2ka/ks",
    test.t(inputdata = one.kaks, var1 = PS.index.names[1], var2 = PS.index.names[2]),
    test.t(inputdata = one.kaks, var1 = PS.index.names[2], var2 = PS.index.names[3]),
    test.t(inputdata = one.kaks, var1 = PS.index.names[3], var2 = PS.index.names[4]),
    test.t(inputdata = one.kaks, var1 = PS.index.names[1], var2 = PS.index.names[3]),
    test.t(inputdata = one.kaks, var1 = PS.index.names[1], var2 = PS.index.names[4]),
    test.t(inputdata = one.kaks, var1 = PS.index.names[2], var2 = PS.index.names[4]), "\n",  sep = "\t\t"
    ),
    file = "figures/figure3C_pvalue.txt", append = T)
  ## percentage barplot
  one.divergence.df <- do.call(rbind, one.divergence)
  one.divergence.df <- one.divergence.df[one.divergence.df$Type=="Diverged pattern", ]
  one.divergence.df <- one.divergence.df[one.divergence.df$PSs%in%PS.index.names, ]
  one.divergence.df$PSs <- factor(one.divergence.df$PSs, levels = PS.index.names)
  div.val <- max(one.divergence.df[, "value"])
  plot_list1[[tmp.two.spp.name]] <- ggplot(data = one.divergence.df, aes(x=PSs, y=value, fill=PSs))+
    geom_bar(stat = "identity", width=0.6) +
    theme_bw(base_size = 16) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none") +
    BaseTheme() +
    xlab("") + ylab("") +
    coord_cartesian(ylim = c(0, div.val+8)) +
    scale_fill_manual(values = c("#036F6F", "#B06A6A",
                                  "#9B70A7", "#C88F17")) 
  # + scale_y_continuous(position = "right")
  ### Ka/Ks boxplot
  one.kaks.df <- do.call(rbind, one.kaks)
  limval <- max(c(boxplot(one.kaks.df[one.kaks.df$PSs=="Angiosperm", ]$value, plot=F)[[1]][5,1] + 0.08,
                  boxplot(one.kaks.df[one.kaks.df$PSs=="Spermatophyta", ]$value, plot=F)[[1]][5,1] + 0.08))
  if (tmp.two.spp.name == "ath-ata"){limval <- 0.45}
  if (tmp.two.spp.name == "ath-sbi"){limval <- 0.48}
  one.kaks.df <- one.kaks.df[one.kaks.df$PSs%in%PS.index.names, ]
  one.kaks.df$PSs <- factor(one.kaks.df$PSs, levels = PS.index.names)
  plot_list2[[tmp.two.spp.name]] <- ggplot(data = one.kaks.df, aes(x=PSs, y=value, fill=PSs))+
    geom_boxplot(outlier.color = NA) +
    theme_bw(base_size = 16) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    BaseTheme() +
    xlab("") + ylab("") +
    coord_cartesian(ylim = c(0, limval)) +
    scale_fill_manual(values = c("#036F6F", "#B06A6A", "#9B70A7", "#C88F17"))
}

# Percentage of m6A divergence 
# Ka/Ks
pdf("figures/Figure3C.pdf", width = 20.5, height = 15) # remove plot_list2[["ath-gar"]] in figure
(plot_list2[["ath-gar"]] | plot_list1[["ath-gar"]] | plot_list1[["ath-pvu"]] | plot_list1[["ath-sly"]] | plot_list1[["ath-sbi"]] | plot_list1[["ath-ata"]] | plot_list1[["ath-osa"]]) /
  (plot_list2[["ath-gar"]] | plot_list2[["ath-gar"]] | plot_list1[["gar-pvu"]] | plot_list1[["sly-gar"]] | plot_list1[["sbi-gar"]] | plot_list1[["gar-ata"]] | plot_list1[["osa-gar"]]) /
  (plot_list2[["ath-pvu"]] | plot_list2[["gar-pvu"]] | plot_list2[["ath-gar"]] | plot_list1[["sly-pvu"]] | plot_list1[["sbi-pvu"]] | plot_list1[["pvu-ata"]] | plot_list1[["osa-pvu"]]) /
  (plot_list2[["ath-sly"]] | plot_list2[["sly-gar"]] | plot_list2[["sly-pvu"]] | plot_list2[["ath-gar"]] | plot_list1[["sly-sbi"]] | plot_list1[["sly-ata"]] | plot_list1[["sly-osa"]]) /
  (plot_list2[["ath-sbi"]] | plot_list2[["sbi-gar"]] | plot_list2[["sbi-pvu"]] | plot_list2[["sly-sbi"]] | plot_list2[["ath-gar"]] | plot_list1[["sbi-ata"]] | plot_list1[["sbi-osa"]]) /
  (plot_list2[["ath-ata"]] | plot_list2[["gar-ata"]] | plot_list2[["pvu-ata"]] | plot_list2[["sly-ata"]] | plot_list2[["sbi-ata"]] | plot_list2[["ath-gar"]] | plot_list1[["osa-ata"]]) /
  (plot_list2[["ath-osa"]] | plot_list2[["osa-gar"]] | plot_list2[["osa-pvu"]] | plot_list2[["sly-osa"]] | plot_list2[["sbi-osa"]] | plot_list2[["osa-ata"]] | plot_list2[["ath-gar"]]) +
  plot_layout(guides = 'collect')
dev.off()


# Figure 3D ####
evo.div.df <- data.frame(
  "m6A_divergence" = uniden.values,
  "evolutionary.rate" = evo.cor.values,
  "Type" = rep(PS.index.names, 21)
)

evo.div.df <- evo.div.df[evo.div.df$Type%in%PS.index.names,]
cor.test(evo.div.df$m6A_divergence, evo.div.df$evolutionary.rate)

fit.stat2 <- summary(lm(as.numeric(evo.div.df[["m6A_divergence"]])~as.numeric(evo.div.df[["evolutionary.rate"]])))
cat(fit.stat2$adj.r.squared, fit.stat2$coefficients[2,4], "\n")

evo.div.df$Type <- factor(evo.div.df$Type, levels = PS.index.names)

pdf("figures/Figure3D.pdf", width = 4, height = 3.5, useDingbats = F)
ggplot(data = evo.div.df, mapping = aes(x=as.numeric(m6A_divergence), 
                                        y=as.numeric(evolutionary.rate))) +
  geom_point(size=2) +
  scale_color_manual(values = c("#036F6F", "#B06A6A", "#9B70A7", "#C88F17")) +
  labs(x="m6A_divergence", y=expression("evolutionary.rate")) +
  theme_bw(base_size = 12)+ 
  BaseTheme() + 
  # expand_limits(x=0, y=0) +
  labs(x="Percentage of m6A divergence (%)", y = "Ka/Ks") +
  geom_smooth(method="lm", se=F, color="#A3A2A2", size=1)  #+
dev.off()


# Supplementary Fig. 10 ####
cl <- makeCluster(7)
registerDoParallel(cl)
system.time(
  one2one.output <- foreach(tmp.two.spp = kaks.comb.list,
                            .packages = c("ggplot2", "ggpubr", "circlize",
                                          "RColorBrewer", "seqinr")) %dopar%
    {
      seq.align.path <- paste0(link.data.path, "seq_align/")
      system(paste0("mkdir -p ", seq.align.path))
      needle.path <- "~/miniconda3/envs/metaPlants/bin/_needle"
      tmp.two.spp.name <- paste0(tmp.two.spp, collapse = "-")
      PS.index.names <- c("Chloroplastida", "Embryophyta", 
                          "Spermatophyta", "Angiosperm")
      tmp.spp.ortho <- read.table(paste0(link.data.path, 
                                         "Results_Jan16_34spp/Orthologues/Orthologues_",
                                         tmp.two.spp[1], "/", tmp.two.spp[1], 
                                         "__v__", tmp.two.spp[2],".tsv"), 
                                  header = T, sep = "\t")
      sp1.cdna <- read.fasta(paste0(link.data.path, "pep_cds/cDNA/", 
                                    tmp.two.spp[1], ".fa"),
                             seqtype = "DNA", as.string = T, forceDNAtolower =F)
      sp2.cdna <- read.fasta(paste0(link.data.path, "pep_cds/cDNA/", 
                                    tmp.two.spp[2], ".fa"),
                             seqtype = "DNA", as.string = T, forceDNAtolower =F)
      output.list <- list()
      for (PS.index in PS.index.names) {
        tmp.l.name <- paste0(tmp.two.spp.name, "_", PS.index)
        PS.tmp.two.spp <- OGs.13spp[OGs.13spp$PSnames == PS.index, tmp.two.spp]
        speciesone <- unlist(strsplit(PS.tmp.two.spp[,1], split = ", "))
        speciestwo <- unlist(strsplit(PS.tmp.two.spp[,2], split = ", "))
        species.gene.index <- apply(tmp.spp.ortho, 1, function(x){
          tmp.index1 <- all(unlist(strsplit(x[2], split = ", ")) %in% speciesone)
          tmp.index2 <- all(unlist(strsplit(x[3], split = ", ")) %in% speciestwo)
          return(tmp.index1 & tmp.index2)
        })
        tmp.ortho.filter <- tmp.spp.ortho[species.gene.index, 2:3]
        PS.tmp.type.index <- apply(tmp.ortho.filter, 1, function(x){
          tmp1 <- unlist(strsplit(x[1], split = ", "))
          tmp2 <- unlist(strsplit(x[2], split = ", "))
          if (length(tmp1)>1){a <- "m"} else {a <- "1"}
          if (length(tmp2)>1){b <- "m"} else {b <- "1"}
          return(paste0(c(a, b), collapse = ":"))
        })
        tmp.ortho.filter <- tmp.ortho.filter[PS.tmp.type.index=="1:1", ]
        PS.tmp.exp <- apply(tmp.ortho.filter, 1, function(x){
          if(x[1]%in%exp.list[[tmp.two.spp[1]]]){a <- "1"}else{a <- "0"}
          if(x[2]%in%exp.list[[tmp.two.spp[2]]]){b <- "1"}else{b <- "0"}
          return(paste0(c(a, b), collapse = ":"))
        })
        tmp.ortho.filter <- tmp.ortho.filter[PS.tmp.exp=="1:1", ]
        PS.tmp.m6A <- apply(tmp.ortho.filter, 1, function(x){
          if(x[1]%in%m6A.list[[tmp.two.spp[1]]]){a <- "1"}else{a <- "0"}
          if(x[2]%in%m6A.list[[tmp.two.spp[2]]]){b <- "1"}else{b <- "0"}
          return(paste0(c(a, b), collapse = ":"))
        })
        m6A.divergence.index <- PS.tmp.m6A%in%c("1:0", "0:1")
        if (sum(m6A.divergence.index) > 0 ){
          m6A.divergence.df <- tmp.ortho.filter[m6A.divergence.index, ]
          identity.score <- apply(m6A.divergence.df, 1, function(x){
            # write.fasta(sequences = sp1.cdna[[x[1]]], names = x[1], file.out = paste0(seq.align.path, x[1],".fa"))
            # write.fasta(sequences = sp2.cdna[[x[2]]], names = x[2], file.out = paste0(seq.align.path, x[2],".fa"))
            # system(command = paste0(needle.path, " -gapopen=10 -gapextend=0.5 -asequence ",
            #                         paste0(seq.align.path, x[1],".fa"), " -bsequence ",
            #                         paste0(seq.align.path, x[2],".fa"), " ",
            #                         paste0(seq.align.path, x[1],"_", x[2],".needle")))
            system(command = paste0("grep 'Identity:' ", 
                                    paste0(seq.align.path, x[1], "_", x[2],".needle | grep -oP '[0-9]+\\.[0-9]+%'")),
                   intern = T)
          })
          output.list[[tmp.l.name]] <- gsub("%", "", identity.score)
        }
      }
      output.list
    }
)
stopImplicitCluster()
stopCluster(cl)

one2one.output2 <- list()
for(i in 1:21){
  tmp.two.spp_name <- strsplit(names(one2one.output[[i]])[1], "_")[[1]][1]
  tmp.two.spp.PS <- as.character(sapply(names(one2one.output[[i]]), function(x){strsplit(x, "_")[[1]][2]}))
  one2one.output2[[tmp.two.spp_name]] <- list()
  for(jj in tmp.two.spp.PS){
    tmp_data <- as.numeric(one2one.output[[i]][[paste0(tmp.two.spp_name, "_", jj)]])
    if(any(is.na(tmp_data))){cat(i, tmp.two.spp_name, jj, "\n")}
    one2one.output2[[tmp.two.spp_name]][[jj]] <- 100-tmp_data
  }
}

plot_list3 <- list()
for (i in kaks.comb.names){
  tmp_identity_list <- one2one.output2[[i]]
  type_counts <- sapply(tmp_identity_list, length)
  tmp_rep <- rep(names(type_counts), type_counts)
  tmp_value <- unlist(tmp_identity_list)
  tmp_df_data <- data.frame("type"=tmp_rep, "value"=tmp_value)
  tmp_df_data$type <- factor(tmp_df_data$type, levels = unique(tmp_df_data$type))
  liminval <- min(tmp_value)
  limaxval <- boxplot(tmp_value, plot=F)[[1]][5,1] + 13
  PS.index.names <- c("Chloroplastida", "Embryophyta", "Spermatophyta", "Angiosperm")
  plot_list3[[i]] <- ggplot(data = tmp_df_data, aes(x=type, y=value, fill=type))+
    geom_boxplot(outlier.colour = NA)+
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom", 
          legend.direction = "horizontal") +
    BaseTheme() +
    coord_cartesian(ylim = c(liminval, limaxval)) + 
    labs(title = i) +
    xlab("") + ylab("") + 
    scale_fill_manual(values = c("#036F6F", "#B06A6A", "#9B70A7", "#C88F17"))
  test.t <- function(inputdata, var1, var2){
    t.test(inputdata[[var1]], inputdata[[var2]], alternative = "less")$p.value
  }
  write(paste(i,
      test.t(inputdata = tmp_identity_list, var1 = PS.index.names[1], var2 = PS.index.names[2]),
      test.t(inputdata = tmp_identity_list, var1 = PS.index.names[2], var2 = PS.index.names[3]),
      test.t(inputdata = tmp_identity_list, var1 = PS.index.names[3], var2 = PS.index.names[4]),
      test.t(inputdata = tmp_identity_list, var1 = PS.index.names[1], var2 = PS.index.names[3]),
      test.t(inputdata = tmp_identity_list, var1 = PS.index.names[1], var2 = PS.index.names[4]),
      test.t(inputdata = tmp_identity_list, var1 = PS.index.names[2], var2 = PS.index.names[4]), "\n",  sep = "\t\t"
  ),
  file = "figures/figureS10_pvalue.txt", append = T)
}

# Frequency of\nsequence mutation (%)
pdf("figures/FigureS10 sequence mutations.pdf", width = 20, height = 13)
(plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] ) /
(plot_list3[["ath-pvu"]] | plot_list3[["gar-pvu"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] ) /
(plot_list3[["ath-sly"]] | plot_list3[["sly-gar"]] | plot_list3[["sly-pvu"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] ) /
(plot_list3[["ath-sbi"]] | plot_list3[["sbi-gar"]] | plot_list3[["sbi-pvu"]] | plot_list3[["sly-sbi"]] | plot_list3[["ath-gar"]] | plot_list3[["ath-gar"]] ) /
(plot_list3[["ath-ata"]] | plot_list3[["gar-ata"]] | plot_list3[["pvu-ata"]] | plot_list3[["sly-ata"]] | plot_list3[["sbi-ata"]] | plot_list3[["ath-gar"]] ) /
(plot_list3[["ath-osa"]] | plot_list3[["osa-gar"]] | plot_list3[["osa-pvu"]] | plot_list3[["sly-osa"]] | plot_list3[["sbi-osa"]] | plot_list3[["osa-ata"]] ) + plot_layout(guides = 'collect')
dev.off()
