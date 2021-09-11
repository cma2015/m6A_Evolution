## plot for overlap between replication
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
  rp[1] <- substitute(expression(Cor == valueA),
                      list(valueA = format(r2,dig = 3)))[2]
  plot(x, y, pch=19, cex=0.3, las=1, xlab = Xlab,ylab = Ylab,col = pointCol,...)
  legend(legendPos, inset=-0.01,legend = rp, bty = 'n', 
         cex=legendCex)
}

## 
suppressMessages(library(DiffBind))
suppressMessages(library(rtracklayer))
suppressMessages(library(venn))
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)

args <- commandArgs(trailingOnly = TRUE)

tamoxifen <- dba(sampleSheet=args[1])  
#tamoxifen$samples

# tamoxifen_count <- dba.count(tamoxifen, score = "DBA_SCORE_RPKM_FOLD", bUseSummarizeOverlaps=TRUE)
# tamoxifen_count$peaks[[2]][1:20,]

sampName <- args[2]

peak_granges <- lapply(tamoxifen$samples$Peaks, import)
peak_grangeslist <- GRangesList(peak_granges)

#"venn_plot"
if(length(peak_grangeslist)==3){
  Rep1 <- peak_grangeslist[[1]]
  Rep2 <- peak_grangeslist[[2]]
  Rep3 <- peak_grangeslist[[3]]
  
  peak_coverage <- coverage(peak_grangeslist)
  covered_ranges <- slice(peak_coverage, lower=2, rangesOnly=T)
  covered_granges <- GRanges(covered_ranges)
  
  cenover <- slice(peak_coverage, lower=3, rangesOnly=T)
  cenover <- GRanges(cenover)
  
  rep1 <- length(intersect(Rep1, covered_granges))-length(cenover)
  rep2 <- length(intersect(Rep2, covered_granges))-length(cenover)
  rep3 <- length(intersect(Rep3, covered_granges))-length(cenover)
  
  venn3 <- length(intersect(Rep2,Rep3))-length(cenover)
  venn5 <- length(intersect(Rep1,Rep3))-length(cenover)
  venn6 <- length(intersect(Rep1,Rep2))-length(cenover)
  venn7 <- length(cenover)
  venn1 <- length(Rep3)-venn7-venn3-venn5
  venn2 <- length(Rep2)-venn7-venn3-venn6
  venn4 <- length(Rep1)-venn7-venn5-venn6
  
  pdf(args[3] , width = 4, height=4)
  venn(3, snames = c(paste0("Rep1(", length(Rep1), ")"), paste0("Rep2(", length(Rep2), ")"), paste0("Rep3(", length(Rep3), ")")),
       count=c(0,venn1,venn2,venn3,venn4,venn5,venn6,venn7), zcolor=rainbow(3), cexil = 0.9, cexsn = 1)
  dev.off()
  
  pdf(args[4] , width = 5, height=8)
  par(mfrow=c(3,2))
  for(cycVec in list(c(1,2), c(1,3), c(2,3))){
    peak_coverage_12 <- coverage(peak_grangeslist[cycVec])
    covered_ranges_12 <-  slice(peak_coverage_12, lower=1, rangesOnly=T)
    covered_granges_12 <- GRanges(covered_ranges_12)
    cat(sampName, "_", cycVec[1], "_", cycVec[2], "DBA_count_start!!\n")
    tamoxifen_count_12 <- dba.count(tamoxifen, score = "DBA_SCORE_RPKM_FOLD", peaks = covered_granges_12, 
                                    filter = NULL, bUseSummarizeOverlaps = T)
    CorrelationScatterPlot(log2(tamoxifen_count_12$peaks[[cycVec[1]]][,5]+1), 
                         log2(tamoxifen_count_12$peaks[[cycVec[2]]][,5]+1), Xlab="", Ylab="")
    mtext(side=1,text=paste0(sampName, "_IP_Rep", cycVec[1]),line=2, cex = 0.6)
    mtext(side=1,text="log2(RPKM+1)",line=3, cex = 0.6)  
    mtext(side=2,text=paste0(sampName, "_IP_Rep", cycVec[2]),line=2, cex = 0.6)
    mtext(side=2,text="log2(RPKM+1)",line=3, cex = 0.6)
    CorrelationScatterPlot(log2(tamoxifen_count_12$peaks[[cycVec[1]]][,7]+1), 
                         log2(tamoxifen_count_12$peaks[[cycVec[2]]][,7]+1), Xlab="", Ylab="")
    #if(sum(cycVec==c(1,2))==2){
    mtext(side=1,text=paste0(sampName, "_Input_Rep", cycVec[1]),line=2, cex = 0.6)
    mtext(side=1,text="log2(RPKM+1)",line=3, cex = 0.6)
    mtext(side=2,text=paste0(sampName, "_Input_Rep", cycVec[2]),line=2, cex = 0.6)
    mtext(side=2,text="log2(RPKM+1)",line=3, cex = 0.6) 
    #}else{
    #  mtext(side=1,text="log2(RPKM)",line=3, cex = 0.6)        
    #}
  }
  
  dev.off()
}

if(length(peak_grangeslist)==2){
  Rep1 <- peak_grangeslist[[1]]
  Rep2 <- peak_grangeslist[[2]]

  peak_coverage <- coverage(peak_grangeslist)
  covered_ranges <- slice(peak_coverage, lower=2, rangesOnly=T)
  covered_granges <- GRanges(covered_ranges)
  
  venn1 <- length(Rep1)-length(covered_granges)
  venn2 <- length(Rep2)-length(covered_granges)
  venn3 <- length(covered_granges)

  pdf(args[3] , width = 4, height=4)
  venn(2, snames = c(paste0("Rep1(", length(Rep1), ")"), paste0("Rep2(", length(Rep2), ")")),
       count=c(0,venn2,venn1,venn3), zcolor=rainbow(2), cexil = 0.9, cexsn = 1)
  dev.off()
  
  pdf(args[4] , width = 7, height=4)
  par(mfrow=c(1,2))
  for(cycVec in list(c(1,2))){
    peak_coverage_12 <- coverage(peak_grangeslist[cycVec])
    covered_ranges_12 <-  slice(peak_coverage_12, lower=1, rangesOnly=T)
    covered_granges_12 <- GRanges(covered_ranges_12)
    cat(sampName, "_", cycVec[1], "_", cycVec[2], "DBA_count_start!!\n")
    tamoxifen_count_12 <- dba.count(tamoxifen, score = "DBA_SCORE_RPKM_FOLD", peaks = covered_granges_12, 
                                    filter = NULL, bUseSummarizeOverlaps = T)
    CorrelationScatterPlot(log2(tamoxifen_count_12$peaks[[cycVec[1]]][,5]+1), 
                         log2(tamoxifen_count_12$peaks[[cycVec[2]]][,5]+1), Xlab="", Ylab="")
    mtext(side=1,text=paste0(sampName, "_IP_Rep", cycVec[1]),line=2, cex = 0.6)
    mtext(side=1,text="log2(RPKM+1)",line=3, cex = 0.6)  
    mtext(side=2,text=paste0(sampName, "_IP_Rep", cycVec[2]),line=2, cex = 0.6)
    mtext(side=2,text="log2(RPKM+1)",line=3, cex = 0.6)
    CorrelationScatterPlot(log2(tamoxifen_count_12$peaks[[cycVec[1]]][,7]+1), 
                         log2(tamoxifen_count_12$peaks[[cycVec[2]]][,7]+1), Xlab="", Ylab="")

    mtext(side=1,text=paste0(sampName, "_Input_Rep", cycVec[1]),line=2, cex = 0.6)
    mtext(side=1,text="log2(RPKM+1)",line=3, cex = 0.6)
    mtext(side=2,text=paste0(sampName, "_Input_Rep", cycVec[2]),line=2, cex = 0.6)
    mtext(side=2,text="log2(RPKM+1)",line=3, cex = 0.6) 
  }
  dev.off()
}