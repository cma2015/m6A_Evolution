  ##
suppressMessages(library(DiffBind))
suppressMessages(library(rtracklayer))
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
# options(scipen = 200)

args <- commandArgs(trailingOnly = TRUE)

peakMerge <- import(args[2])
tamoxifen <- dba(sampleSheet=args[1])  
system.time(
  txfCount <- dba.count(tamoxifen, peaks = peakMerge, bUseSummarizeOverlaps = T, minOverlap = 1)
)

if(length(tamoxifen$peaks) == 3){
  outDf <- data.frame("Chr"= txfCount$peaks[[1]][,1], "Start"= txfCount$peaks[[1]][,2]-1, "End" = txfCount$peaks[[1]][,3],
                    "IP1"= round(as.numeric(txfCount$peaks[[1]][,5]), 2), "IP2"= round(as.numeric(txfCount$peaks[[2]][,5]), 2), 
                    "IP3"= round(as.numeric(txfCount$peaks[[3]][,5]), 2), "Input1"= round(as.numeric(txfCount$peaks[[1]][,7]), 2), 
                    "Input2"= round(as.numeric(txfCount$peaks[[2]][,7]), 2), "Input3"= round(as.numeric(txfCount$peaks[[3]][,7]), 2)
                    )

}
if(length(tamoxifen$peaks) == 2){
  outDf <- data.frame("Chr"= txfCount$peaks[[1]][,1], "Start"= txfCount$peaks[[1]][,2]-1, "End" = txfCount$peaks[[1]][,3],
                      "IP1"= round(as.numeric(txfCount$peaks[[1]][,5]), 2), "IP2"= round(as.numeric(txfCount$peaks[[2]][,5]), 2),
                      "Input1"= round(as.numeric(txfCount$peaks[[1]][,7]), 2), "Input2"= round(as.numeric(txfCount$peaks[[2]][,7]), 2)
                      )
  
}

write.table(outDf, args[3], sep = "\t", quote = F, row.names = F)
