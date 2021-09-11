library(DiffBind);
library(ChIPpeakAnno);
library(rtracklayer);


peak_files <- list()
args <- commandArgs(trailingOnly = TRUE)

if(length(args)==4){
  peak_files[[1]] <- args[1]
  peak_files[[2]] <- args[2]
  peak_files[[3]] <- args[3]
  
  peak_granges <- lapply(peak_files, import)
  peak_grangeslist <- GRangesList(peak_granges)
  peak_coverage <- coverage(peak_grangeslist)
  covered_ranges <- slice(peak_coverage, lower=2, rangesOnly=T)
  covered_granges <- GRanges(covered_ranges)
  
  export(covered_granges, args[4], "bed")
}

if(length(args)==3){
  peak_files[[1]] <- args[1]
  peak_files[[2]] <- args[2]

  peak_granges <- lapply(peak_files, import)
  peak_grangeslist <- GRangesList(peak_granges)
  peak_coverage <- coverage(peak_grangeslist)
  covered_ranges <- slice(peak_coverage, lower=2, rangesOnly=T)
  covered_granges <- GRanges(covered_ranges)
  
  export(covered_granges, args[3], "bed")
}
