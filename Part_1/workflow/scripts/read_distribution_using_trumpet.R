
#"ghi", "pvu", "gma", "sbi", "zma"

library(foreach)
library(doParallel)




cl <- makeCluster(3)
registerDoParallel(cl, cores=3)

system.time(
  plots_list1 <- foreach(spp=c("ata","tdi","zma")) %dopar%
    {
      library(ggplot2)
      library(GenomicFeatures)
      library(Guitar) ## version 1.5.0
      library(reshape2)
      library(scales)
      
      get.sampleid <- function(IP_BAM, Input_BAM, BAM_SUFFIX) {
        IP_name <- paste0("IP:",gsub(BAM_SUFFIX, "", basename(IP_BAM))) #paste0("IP", seq_along(IP_BAM))
        Input_name <- paste0("Input:",gsub(BAM_SUFFIX, "", basename(Input_BAM))) #paste0("Input", seq_along(Input_BAM))
        cat(IP_name, "\n")
        cat(Input_name, "\n")
        referIP_name <- NULL
        referInput_name <- NULL
        
        sample_name <- list()
        sample_name$IP_name <- IP_name
        sample_name$Input_name <- Input_name
        sample_name$referIP_name <- referIP_name
        sample_name$referInput_name <- referInput_name
        return(sample_name)
      }
      
      
      get_readscount <- function(IP_BAM, Input_BAM, GENE_ANNO_GTF,BAM_SUFFIX) {
        txdb <- makeTxDbFromGFF(GENE_ANNO_GTF, format = "gtf")
        gc <- makeGuitarCoordsFromTxDb(txdb, noBins = 20, minimalNcRNALength = 0)
        gc_info <- mcols(gc)
        index <- c(IP_BAM, Input_BAM)
        index <- as.character(index)
        file <- index
        sample_name <- get.sampleid(IP_BAM, Input_BAM, BAM_SUFFIX)
        sample_name <- unlist(sample_name)
        sample_name <- as.vector(sample_name)
        transform_table <- cbind(file, sample_name)
        transform_table <- as.data.frame(transform_table)
        colnames(transform_table) <- c("files", "sample ID")
        # get reads count
        result2 <- gc_info
        noFiles <- length(file)
        total_reads <- vector(length = noFiles)
        exon_reads <- vector(length = noFiles)
        intron_reads <- vector(length = noFiles)
        no_genic <- vector(length = noFiles)
        percent_exon <- vector(length = noFiles)
        percent_intron <- vector(length = noFiles)
        percent_nogenic <- vector(length = noFiles)
        UTR5_reads <- vector(length = noFiles)
        CDS_reads <- vector(length = noFiles)
        UTR3_reads <- vector(length = noFiles)
        percent_UTR5 <- vector(length = noFiles)
        percent_CDS <- vector(length = noFiles)
        percent_UTR3 <- vector(length = noFiles)
        
        for (i in seq_len(noFiles)) {
          print(paste("working on the ", i, "-th bam file ...", sep = ""))
          library(GenomicAlignments)
          bam <- readGAlignments(file[i])
          total_reads[i] <- paste0(round(length(bam)/10^6, 2), "M")
          bin_count <- countOverlaps(gc, bam)
          bin_count <- data.frame(bin_count)
          names(bin_count) <- sample_name[i]
          result2 <- data.frame(result2, bin_count)
          exon <- exonsBy(txdb, by = "tx")
          exon_count <- countOverlaps(bam, exon)
          exon_reads[i] <- paste0(round(sum(exon_count>0)/10^6, 2), 
                                  "M")
          percent_exon[i] <- paste0(round((sum(exon_count>0)/(length(bam)))*100,2), "%")
          percent_exon[i] <- paste0("(", percent_exon[i], ")")
          
          intron <- intronsByTranscript(txdb)
          intron_count <- countOverlaps(bam, intron)
          intron_reads[i] <- paste0(round(sum(intron_count > 0)/10^6, 
                                          2), "M")
          percent_intron[i] <- paste0(round((sum(intron_count > 0)/(length(bam)))*100,2), "%")
          percent_intron[i] <- paste0("(", percent_intron[i], ")")
          
          no_genic[i] <- paste0(round((length(bam)-sum(exon_count>0)-sum(intron_count>0))/10^6,2), "M")
          percent_nogenic[i] <- paste0(round(((length(bam)-sum(exon_count>0)-sum(intron_count>0))/length(bam))*100,2), "%")
          percent_nogenic[i] <- paste0("(", percent_nogenic[i], ")")
          
          utr5 <- fiveUTRsByTranscript(txdb)
          utr5_count <- countOverlaps(bam, utr5)
          UTR5_reads[i] <- paste0(round(sum(utr5_count > 0)/10^6, 2), 
                                  "M")
          cds <- cdsBy(txdb, by = "tx")
          cds_count <- countOverlaps(bam, cds)
          CDS_reads[i] <- paste0(round(sum(cds_count > 0)/10^6, 2), "M")
          utr3 <- threeUTRsByTranscript(txdb)
          utr3_count <- countOverlaps(bam, utr3)
          UTR3_reads[i] <- paste0(round(sum(utr3_count > 0)/10^6, 2), 
                                  "M")
          sum_component <- sum(utr5_count > 0) + sum(cds_count > 0) + 
            sum(utr3_count > 0)
          percent_UTR5[i] <- paste0(round((sum(utr5_count > 0)/(sum_component)) * 
                                            100, 2), "%")
          percent_UTR5[i] <- paste0("(", percent_UTR5[i], ")")
          percent_CDS[i] <- paste0(round((sum(cds_count > 0)/(sum_component)) * 
                                           100, 2), "%")
          percent_CDS[i] <- paste0("(", percent_CDS[i], ")")
          percent_UTR3[i] <- paste0(round((100 - round((sum(utr5_count > 0)/(sum_component)) * 
                                                         100, 2) - round((sum(cds_count > 0)/(sum_component)) * 
                                                                           100, 2)),2), "%")
          percent_UTR3[i] <- paste0("(", percent_UTR3[i], ")")
        }
        reads_exon <- paste(exon_reads, percent_exon)
        reads_intron <- paste(intron_reads, percent_intron)
        reads_nogenic <- paste(no_genic, percent_nogenic)
        reads_UTR5 <- paste(UTR5_reads, percent_UTR5)
        reads_CDS <- paste(CDS_reads, percent_CDS)
        reads_UTR3 <- paste(UTR3_reads, percent_UTR3)
        read_alignment_summary <- cbind(sample_name, total_reads,reads_exon, reads_intron, reads_nogenic, reads_UTR5, reads_CDS, reads_UTR3)
        read_alignment_summary <- as.data.frame(read_alignment_summary)
        colnames(read_alignment_summary) <- c("sample ID","total reads#", "exon reads", "intron reads", "no genic reads", "5'UTR reads", "CDS reads", "3'UTR reads")
        # remove DNA regions
        i <- which(result2$comp == "Front")  # remove front DNA
        result2 <- result2[-i, ]
        i <- which(result2$comp == "Back")  # remove tail DNA
        result <- result2[-i, ]
        t <- data.frame(result)
        t1 <- t[(t$category) == "mRNA", ]
        t2 <- t1[(t1$comp == "UTR5"), ]
        t3 <- t1[(t1$comp == "CDS"), ]
        t3$pos <- t3$pos + 1
        t4 <- t1[(t1$comp == "UTR3"), ]
        t4$pos <- t4$pos + 2
        t0 <- rbind(t2, t3, t4)
        s <- data.frame()
        s <- aggregate(cbind(t0[, 6]) ~ pos + txid, t0, mean)
        for (i in (length(t0) - noFiles + 2):length(t0)) {
          w <- aggregate(cbind(t0[, i]) ~ pos + txid, t0, mean)
          w <- w[, -(1:2)]
          s <- cbind(s, w)
        }
        s <- as.data.frame(s)
        colnames(s) <- c("pos", "txid", sample_name)
        read_result <- list()
        read_result[[1]] <- txdb
        read_result[[2]] <- list(s, read_alignment_summary, transform_table)
        return(read_result)
      }
      
      normalize_sample<-function(s1){
        row.sum<-rowSums(s1)
        z<-which(row.sum>10)
        s2<-s1[z,]
        row.mean<-rowMeans(s2)
        for(i in seq_len(length(row.mean))){
          if(row.mean[i]<2)
            row.mean[i]<-2
        }
        s3<-apply(s2,2,function(x,a)x/a, a=row.mean)
        s3<-as.matrix(s3)
        return(s3)
      }
      
      read_cover <- function(read_count, ind, bam_count) {
        cat("Sample Reads:", bam_count, "\n")
        r <- rowSums(read_count)
        read_count <- cbind(read_count, r)
        read_count <- as.data.frame(read_count)
        read_count <- read_count[order(read_count[, (length(ind) + 1)]),                            ]
        read_count <- read_count[, -(length(ind) + 1)]
        s3 <- normalize_sample(read_count)
        b <- vector(mode = "numeric", length = 0)
        c <- vector(mode = "numeric", length = 0)
        d <- vector(mode = "numeric", length = 0)
        e <- vector(mode = "numeric", length = 0)
        for (i in seq_len(length(ind))) {
          b[i] <- quantile(s3[, i], 0.25)
          c[i] <- quantile(s3[, i], 0.5)
          d[i] <- quantile(s3[, i], 0.75)
          e[i] <- sum(read_count[, i])/bam_count
        }
        qv <- cbind(b, c, d, e)
        pos <- seq(0.025, 2.975, 0.05)
        dt <- cbind(pos, qv)
        dt <- as.data.frame(dt)
        colnames(dt) <- c("pos", "25%", "50%", "75%","100%")
        df <- melt(dt, id = "pos")
        df <- data.frame(df)
        return(df)
      }
      
      singleBAMreads <- function(bam, se, len, n) {
        w <- vector(mode = "numeric", length = 0)
        for (i in seq_len(len)) {
          se <- seq(i, n, len)
          w <- cbind(w, bam[se])
        }
        return(w)
      }
      
      unified_sample<-function(group,se,ind,len,n)
      {
        m<-matrix(nrow = length(se),ncol = length(ind))
        v<-matrix(data=0,nrow=length(se),ncol = length(ind))
        group<-as.matrix(group)
        for(i in seq_len(ncol(group))){
          m<-singleBAMreads(group[,i],se,len,n)
          v<-m+v
        }
        v<-v/length(ncol(group))  
        return(v)
      }
      
      read_distribute <- function(result, IP_BAM, Input_BAM, BAM_SUFFIX) {
        txdb <- result[[1]] #makeTxDbFromGFF(GENE_ANNO_GTF, format = "gtf")
        ## get component
        utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
        cds <- cdsBy(txdb, by = "tx",  use.names=TRUE)
        utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)
        ## get componet name
        utr5_name <- names(utr5)
        utr5_name <- as.character(utr5_name)
        cds_name <- names(cds)
        cds_name <- as.character(cds_name)
        utr3_name <- names(utr3)
        utr3_name <- as.character(utr3_name)
        ## Select tx 
        s <- result[[2]][[1]]
        max_transcript_name <- as.character(unique(s$txid))
        ## Select component
        select_utr5_name <- as.numeric(match(utr5_name, max_transcript_name)) 
        select_utr5_name <- select_utr5_name[-which(is.na(select_utr5_name))] 
        select_utr5 <- utr5[select_utr5_name]
        select_cds_name <- as.numeric(match( cds_name, max_transcript_name))
        select_cds_name <- select_cds_name[-which(is.na(select_cds_name))] 
        select_cds <- cds[select_cds_name]
        select_utr3_name <- as.numeric(match(utr3_name, max_transcript_name))
        select_utr3_name <- select_utr3_name[-which(is.na(select_utr3_name))]
        select_utr3 <- utr3[select_utr3_name]
        
        select_name <- c(names(select_utr5), names(select_cds), names(select_utr3))
        last_select_name <- rownames(as.data.frame( which(table(select_name)==3))) 
        last_cds_num <- match(last_select_name, names(select_cds))
        last_cds <- select_cds[last_cds_num ]                     
        last_utr5_num <- match(last_select_name, names(select_utr5))
        last_utr5 <- select_utr5[last_utr5_num]
        last_utr3_num <- match(last_select_name, names(select_utr3))
        last_utr3 <- select_utr3[last_utr3_num]
        ## component size
        sum_comb <- function(i, component){
          sum_width <- sum(width(component[[i]]))
        }
        utr5_size <- unlist(lapply(1:length(last_utr5), sum_comb, component=last_utr5)) 
        cds_size <- unlist(lapply(1:length(last_cds), sum_comb, component=last_cds)) 
        utr3_size <- unlist(lapply(1:length(last_utr3), sum_comb, component=last_utr3)) 
        ## component size factor
        utr5.SF <- round((median(utr5_size)/median(cds_size)), 2)
        utr3.SF <- round((median(utr3_size)/median(cds_size)), 2) 
        ind <- unique(s$pos)
        len <- length(ind)
        n <- nrow(s)
        se <- seq(1, n, len)
        sa <- s[, -(1:2)]
        sample_name <- get.sampleid(IP_BAM, Input_BAM, BAM_SUFFIX)
        IP_groupname <- sample_name[[1]]
        Input_groupname <- sample_name[[2]]
        
        group <- sa[, seq_len(length(IP_groupname) + length(Input_groupname))]
        p <- data.frame()
        q <- data.frame()
        group <- as.matrix(group)
        for (i in seq_len(ncol(group))) {
          bam_count <- as.numeric(system(paste0("samtools view -c -F 260 ", c(IP_BAM, Input_BAM)[i]), intern = T))
          d <- singleBAMreads(bam = group[, i], se = se, len = len, n = n)
          p <- read_cover(read_count = d, ind = ind, bam_count = bam_count)
          q <- rbind(q, p)
        }
        Group_IP <- q[seq_len(length(IP_groupname) * length(unique(q$pos)) * 
                                (nrow(q[q$pos == (unique(q$pos)[1]), ]))/ncol(group)), ]
        Group_Input <- q[-(seq_len(nrow(Group_IP))), ]
        fr_num <- length(unique(q$pos)) * (nrow(q[q$pos == (unique(q$pos)[1]), 
                                                  ]))/ncol(group)
        id_name1 <- c(IP_groupname, Input_groupname)
        ID1 <- rep(id_name1, rep(fr_num, length(id_name1)))
        df1 <- rbind(Group_IP, Group_Input)
        df1 <- cbind(df1, ID1)
        df1 <- as.data.frame(df1)
        
        colnames(df1)<-c("pos","Quantile","value","Samples")
        pos <- df1$pos
        value <- df1$value
        Quantile <- df1$Quantile
        pos_utr5 <- unique(pos)[1:20]
        pos_cds <- unique(pos)[21:40]
        pos_utr3 <- unique(pos)[41:60]
        rescale_utr5 <- rescale(pos_utr5, c(1-utr5.SF,1), from = c(0,1))
        rescale_utr3 <- rescale(pos_utr3, to=c(2,2+utr3.SF), from = c(2,3))
        pos <- rep( c( rescale_utr5,  pos_cds, rescale_utr3), 4*length(unique(df1$Samples)))
        df1 <- cbind(pos, df1[,-1])
        max_value <- max(df1[df1$Quantile=="100%",]$value)
        p1 <- ggplot(df1[df1$Quantile=="100%",], aes(pos, value, colour = Samples)) + 
          geom_line() + 
          # facet_grid(ID1~.) +  
          geom_vline(xintercept = 1:2, linetype = "dotted", colour="red") + 
          annotate("rect", xmin = pos[1], xmax = pos[length(unique(pos))/3], ymin = -0.03*max_value, ymax = -0.02*max_value, alpha = 0.99, colour = "black") + 
          annotate("rect", xmin = pos[length(unique(pos))/3], xmax = pos[(length(unique(pos))*2/3+1)]-0.007, ymin = -0.04*max_value, ymax = -0.01*max_value, alpha = 0.99, colour = "black") + 
          annotate("rect", xmin = pos[(length(unique(pos))*2/3+1)]-0.007, xmax = pos[length(unique(pos))], ymin = -0.03*max_value, ymax = -0.02*max_value, alpha = 0.99, colour = "black") + 
          annotate("text", x = (pos[length(unique(pos))/6]+0.005), y= 0, label = "5'UTR", size = 3) + 
          annotate("text", x = 1.5, y = 0, label = "CDS", size = 3) + 
          annotate("text", x = (pos[(length(unique(pos))*2/3+1)]+pos[length(unique(pos))])/2, y = 0, label = "3'UTR", size = 3) + 
          theme_bw(base_size = 16)+
          theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=16),
                axis.text = element_text(colour = "black"),
                panel.border = element_rect(colour = "black"),
                axis.text.x = element_blank(), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                axis.ticks.x = element_blank(),
                title = element_text(size = 16),
                plot.title = element_text(hjust = 0.5),
                legend.key.height=unit(0.5,'cm'),
                legend.key.width=unit(0.25,'cm'),
                legend.text=element_text(size=9),
                legend.title=element_text(size=9))+
          labs(title = "", x = "", y = "Normalized Reads coverage")
        p1
      }
      
      
      firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
      }
      
      spp_nn <- firstup(spp)
      bam_suffix <- "_Aligned.sortedByCoord.out.bam"
      GENE_ANNO_GTF <- paste0("~/a2z/metaPlants/workflow/genome/", 
                              spp_nn, "/Annotation/", spp_nn, ".exons.gtf")
      if(!file.exists(GENE_ANNO_GTF)){
        system(paste0("~/miniconda3/envs/metaPlants/bin/python ", 
                      "~/a2z/metaPlants/workflow/scripts/gff3_gtf.py ", 
                      gsub("gtf$", "gff3", GENE_ANNO_GTF)))
        gene_count_n <- system(paste0("awk '$3==\"gene\"' ", 
                                      basename(GENE_ANNO_GTF), " | wc -l"), intern = T)
        system(paste0("mv ", basename(GENE_ANNO_GTF), 
                      " ", GENE_ANNO_GTF))
        cat( spp, gene_count_n, "\n")
      }
      GENE_ANNO_GTF_TMP <- paste0("~/a2z/metaPlants/workflow/genome/", 
                                  spp_nn, "/Annotation/", spp_nn, ".exons_tmp.gtf")
      IP_BAM <- paste0("data_12ss/", spp, "/02STAR/", spp, "_IP", 1:3, bam_suffix)
      Input_BAM <- paste0("data_12ss/", spp, "/02STAR/", spp, "_Input", 1:3, bam_suffix)
      cat(spp, file.exists(GENE_ANNO_GTF_TMP), file.exists(IP_BAM[1]), "\n")
      if(file.exists(IP_BAM[1])){
        result <- get_readscount(IP_BAM = IP_BAM, Input_BAM = Input_BAM,
                                 GENE_ANNO_GTF = GENE_ANNO_GTF_TMP, BAM_SUFFIX=bam_suffix)
        plot_result <- read_distribute(result = result, IP_BAM = IP_BAM,
                                       Input_BAM = Input_BAM, BAM_SUFFIX=bam_suffix)
        ##
        outfile_name <- paste0("Read_distribution_", spp, ".pdf")
        pdf(outfile_name, width = 10, height = 7)
        plot_result
        dev.off()
      }
      plot_result
    }
)
stopImplicitCluster()
stopCluster(cl)

outfile_name <- paste0("Read_distribution1.pdf")
pdf(outfile_name, width = 10, height = 7)
plots_list1[[1]]
plots_list1[[2]]
plots_list1[[3]]
dev.off()

gc <- makeGuitarCoordsFromTxDb(makeTxDbFromGFF("~/a2z/metaPlants/workflow/genome/Zma/Annotation/Zma.exons_tmp.gtf", format = "gtf"),
                               noBins = 20, minimalNcRNALength = 0)


