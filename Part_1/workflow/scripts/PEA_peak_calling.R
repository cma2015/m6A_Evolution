###########CMR (Chemical Modification of RNA) Calling
##########author: Jingjing Zhai, Chuang Ma
##########Contact: zhaijingjing603@gmail.com

.libPaths('/home/zhangt/R/x86_64-pc-linux-gnu-library/3.6')

# Import required libraries
library(getopt)
library(data.table)

options(stringAsfactors = F, useFancyQuotes = F)

option_specification = matrix(c('IPBED', 'i', 2, 'character',
                                'INPUTBED', 'n', 2, 'character',
                                'IPBG', 'p', 2, 'character',
                                'INPUTBG', 'q', 2, 'character',
                                'GENOMESIZE', 'g', 2, 'character',
                                'OUTFILE', 'o', 2, 'character'
),byrow=TRUE, ncol=4);
# Parse options
opt = getopt(option_specification);


.runMACS2 <- function(IPBAM, inputBAM, broad = FALSE, 
                      expName = NULL, paired = F,MACS2Dir, ...){
  
  tt <- list(...)
  if(length(tt) != 0){
    callPara <- tt[[1]]
  }else{
    callPara <- NULL
  }
  
  if(is.null(expName)){
    expName <- "example"
  }
  
  if(!broad){
    if(!paired){
      macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                            " -f BAM -n ", expName, " ", callPara)
    }else{
      macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                            " -f BAMPE -n ", expName, " ", callPara)
    }
    
  }else{
    
    if(!paired){
      macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                            " --broad ", callPara)
    }else{
      macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                            " -f BAMPE --broad ", callPara)
    }
    
  }
  
  system(command = macsCommand)
  
  peakDir <- paste0(expName, "_peaks.narrowPeak")
  peaks <- as.matrix(read.table(file = peakDir, sep = "\t", header = F, quote = ""))
  
  resList <- list(macsCommand = macsCommand, peaks = peaks)
  resList
}



.bam2bed <- function(BAM){
  
  #cat("Note: this function is used for converting bam format to bed format using bedtools!")
  BAMFile <- unlist(strsplit(x = BAM, split = "/", fixed = T))
  BAMFile <- BAMFile[length(BAMFile)]
  BAMName <- unlist(strsplit(x = BAMFile, split = ".", fixed = T))[1]
  bedCommand <- paste0("bedtools bamtobed -i ", BAM, " > ", BAMName, ".bed")
  system(bedCommand)
  resBed <- paste0(BAMName, ".bed")
  resBed
}



########################################Bisulfite-seq######################
##########################meRanTK#########################################

.m5Call <- function(meRanDir = "/home/software/meRanTK-1.2.0/", fq, paired = FALSE, refGenome, cpus = 1, 
                    method = c("meRanGh", "meRanGs"), GTF = NULL,  
                    resDic = NULL, maxDup = 50, conversionRate = 0.99,
                    minBaseQ = 30, ...){
  
  if(length(method) != 1){
    method <- method[1]
    cat("Note: methods are not defined, the first meRanGh will be used!")
  }
  
  if(method == "meRanGs" & is.null(GTF)){
    stop("The GTF/GTF file is not provided, which is required for meRanGs!")
  }
  
  if(is.null(resDic)){
    resDic <- paste0(getwd(), "/")
  }
  
  
  ############Building index
  idxDir <- paste0(resDic, method, "IDX")
  if(method == "meRanGh"){
    idxCommand <- paste0(meRanDir, method, " mkbsidx -t ", cpus, " -fa ", refGenome, 
                         " -id ", idxDir)
  }else{
    idxCommand <- paste0(meRanDir, method, " mkbsidx -t ", cpus, " -fa ", refGenome, 
                         " -GTF ", GTF, " -GTFtagEPT Parent -GTFtagEPG gene",  
                         " -id ", idxDir)
  }
  
  system(command = idxCommand)
  
  #########Alignment##########################
  SRAID <- unlist(strsplit(x = fq[1], "/"))
  SRAID <- SRAID[length(SRAID)]
  SRAID <- unlist(strsplit(x = SRAID, split = ".", fixed = TRUE))[1]
  alignDir <- paste0(resDic, method)
  dir.create(path = alignDir, showWarnings = FALSE)
  
  if(!paired){
    if(method == "meRanGh"){
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq, 
                             " -id ", idxDir, " -GTF ", GTF, " -bg -o ", alignDir,
                             " -S ", method, "_", SRAID, ".sam -MM -un")
    }else{
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq, 
                             " -id ", idxDir, " -GTF ", GTF, " -bg -o ", alignDir,
                             " -S ", method, "_", SRAID, ".sam -MM -un --star_genomeLoad")
    }
  }else{
    if(method == "meRanGh"){
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq[1], 
                             " ", fq[2], " -id ", idxDir, " -GTF ", GTF, " -bg -o ",
                             alignDir, " -S ", method, "_", SRAID, ".sam -MM -un")
    }else{
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq[1], 
                             " ", fq[2], " -id ", idxDir, " -GTF ", GTF, " -bg -o ",
                             alignDir, " -S ", method, "_", SRAID, ".sam -MM -un --star_genomeLoad")
    }
  }
  system(command = alignCommand)
  
  ########m5C Calling#############
  callCommand <- paste0(meRanDir, "meRanCall -p ", cpus, " -s ", alignDir, "/", 
                        method, "_", SRAID, "_sorted.bam -f ", refGenome, " -gref -o ",
                        resDic, method, "_meRanCall_m5C.txt", " -md ", maxDup,
                        " -cr ", conversionRate, " -mBQ ", minBaseQ, " -bed63 ")
  system(callCommand)
  m5CDir <- paste0(resDic, method, "_meRanCall_m5C.txt")
  m5CresMat <- as.matrix(read.table(file = m5CDir, sep = "\t", header = F, 
                                    quote = ""))
  
  resList <- list(idxCommand = idxCommand, alignCommand = alignCommand, 
                  callCommand = callCommand, m5C = m5CresMat)
  
  resList
}

#' @export
extractCov <- function(BAM, refGenome, method = c("bedtools", "Rsamtools")){
  
  if(length(method) > 1){
    method <- method[1]
  }
  
  if(method == "bedtools"){
    # command <- paste0("genomeCoverageBed -split -dz -ibam ", BAM, " -g ", 
    #                   refGenome, " > ", paste0(BAM, ".dz"))
    # bed file
    command <- paste0("genomeCoverageBed -split -dz -i ", BAM, " -g ", 
                      refGenome, " > ", paste0(BAM, ".dz"))
    if(!file.exists(paste0(BAM, ".dz"))){
      cat(command, "\n")
      system(command = command)
    }
  }else{
    refGenome <- read.fasta(file = refGenome, as.string = T)
    Chr <- names(refGenome)
    resMat <- NULL
    for(i in 1:length(Chr)){
      seqID <- Chr[i]
      resPos <- 1:(nchar(refGenome[[i]]))
      USiteIRange <- IRanges(resPos, resPos)
      names(USiteIRange) <- seqID
      tmpList <- list()
      tmpList[[seqID]] <- USiteIRange
      which <- RangesList(unlist(tmpList))
      p1 <- ScanBamParam(which = which, what = scanBamWhat())
      res1 <- scanBam(BAM, param = p1)
      number <- lapply(res1, function(x) nrow(as.data.frame(x)))
      tmpMat <- data.frame(Chr = rep(seqID, length(number)), Position = 1:(nchar(refGenome[[i]])),
                           number = unlist(number))
      resMat <- rbind(resMat, res1)
    }
    write.table(resMat, file = paste0(BAM, ".dz"), sep = "\t",
                quote = F, row.names = F, col.names = F)
  }
  res <- paste0(BAM, ".dz")
  res
}


.findContinuous <- function(inputVec){
  Breaks <- c(0, which(diff(inputVec) != 1), length(inputVec)) 
  res <- sapply(seq(length(Breaks) - 1), 
                function(i) inputVec[(Breaks[i] + 1):Breaks[i+1]]) 
  res
}



.intervalCov <- function(dzFile){
  
  baseCov <- fread(file = dzFile, sep = "\t", header = F, stringsAsFactors = F)
  baseCov <- data.frame(baseCov, stringsAsFactors = F)
  winDow <- rep(NA, nrow(baseCov))
  baseCov <- cbind(baseCov, winDow)
  colnames(baseCov) <- c("Chr", "Position", "readsNumber", "window")
  baseCov$window <- ceiling(baseCov$Position/25)
  
  curRes <- by(data = baseCov[,c(3,4)], 
               INDICES = factor(baseCov$Chr, levels = unique(baseCov$Chr)),
               FUN = function(x) by(x$readsNumber, factor(x$window, 
                                                          levels = unique(x$window)),
                                    mean))
  curRes <- lapply(curRes, as.numeric)
  Chr <- unique(baseCov$Chr)
  resList <- vector("list", length = length(Chr))
  names(resList) <- Chr
  for(i in 1:length(resList)){
    curChr <- Chr[i]
    curMat <- subset(x = baseCov, baseCov$Chr == curChr)
    curMat <- subset(curMat, !duplicated(curMat$window))
    curWave <- rep(NA, nrow(curMat))
    curMat <- cbind(curMat, curWave)
    curMat$curWave <- curRes[[i]]
    resList[[i]] <- curMat
  }
  resList
}
.getPvalue <- function(inputVec, mappedInput, mappedRIP){
  testMat <- matrix(c(as.numeric(inputVec[10]),
                      mappedInput,
                      as.numeric(inputVec[11]),
                      mappedRIP), nrow = 2, ncol = 2)
  p.value <- fisher.test(x = testMat)$p
  ratio <- log2(((as.numeric(inputVec[10]) + 1)*mappedRIP)/((as.numeric(inputVec[11]) + 1)*mappedInput))
  res <- c(ratio, p.value)
  res
}

#' @export
SlidingWindow <- function(input, RIP, mappedInput = NULL, 
                          mappedRIP = NULL, level = 0.05, 
                          ratio = -1, cpus = 1, readsCount = 10,
                          concatenate = 4, ...){
  cat("Note: parameter mappedInput and mappedRIP represent the number of
      reads mapped to reference genome in input and RIP experiments, respectively!\n")
  input.dz <- extractCov(BAM = input, ...)
  RIP.dz <- extractCov(BAM = RIP, ...)
  input <- .intervalCov(dzFile = input.dz)
  RIP <- .intervalCov(dzFile = RIP.dz)
  
  if(!all(names(input) == names(RIP))){
    cat("Note: The chromosomes in the input and RIP are not consistent!\n", 
        "the interactions will be used!")
    interNames <- intersect(names(input), names(RIP))
    input <- input[interNames]
    RIP <- RIP[interNames]
  }
  
  resList <- vector("list", length = length(RIP))
  names(resList) <- names(RIP)
  #i=1
  for(i in 1:length(input)){
    curMat <- merge(input[[i]], RIP[[i]], by = 'window', all=TRUE)
    
    curMat$curWave.x[is.na(curMat$curWave.x)] <- 0
    curMat$curWave.y[is.na(curMat$curWave.y)] <- 0
    curMat[, 'windowave.input'] <- round(curMat$curWave.x)
    curMat[, 'windowave.RIP'] <- round(curMat$curWave.y)
    curMat <- curMat[curMat$windowave.RIP != 0, ]
    curPvalue <- rep(NA, nrow(curMat))
    curRatio <- rep(NA, nrow(curMat))
    curFDR <- rep(NA, nrow(curMat))
    curMat <- cbind(curMat, curPvalue, curRatio, curFDR)
    cat("Chromosome: ", names(input)[i], "...... cpus: ", cpus,"\n")
    
    if(cpus == 1){
      system.time(pvalue <- t(apply(curMat, 1, .getPvalue, mappedInput = mappedInput,
                        mappedRIP = mappedRIP)))
    }else{
      library(foreach)
      library(doParallel)
      cl <- makeCluster(cpus)
      registerDoParallel(cl)
      
      system.time( pvalue <- foreach(x = 1:nrow(curMat), .combine = "rbind") %dopar%{
        testMat <- matrix(c(as.numeric(curMat[x, 10]),
                            mappedInput,
                            as.numeric(curMat[x, 11]),
                            mappedRIP), nrow = 2, ncol = 2)
        p.value <- fisher.test(x = testMat)$p
        ratio <- log2(((as.numeric(curMat[x, 10]) + 1)*mappedRIP)/((as.numeric(curMat[x, 11]) + 1)*mappedInput))
        res <- c(ratio, p.value)
        res
      })
      stopImplicitCluster()
      stopCluster(cl)
      # cat(nrow(curMat), nrow(pvalue), "\n")
    }
    curMat$curPvalue <- as.numeric(pvalue[,2])
    curMat$curRatio <- as.numeric(pvalue[,1])
    curMat$curFDR <- p.adjust(curMat$curPvalue, "fdr")
    resList[[i]] <- curMat
  }
  resMat <- do.call(what = rbind, args = resList)
  resMat <- subset(resMat, resMat$curRatio < ratio & resMat$curFDR < level & 
                     resMat$curWave.y >= readsCount)
  resMat <- resMat[which(!is.na(resMat$Chr.y)), ]
  tt <- by(resMat$window, factor(x = resMat$Chr.y, levels = unique(resMat$Chr.y)),
           .findContinuous)
  resPeaks <- NULL
  for(i in 1:length(tt)){
    curList <- tt[[i]]
    curLen <- unlist(lapply(curList, length))
    curList <- curList[which(curLen >= concatenate)]
    if(length(curList) == 0){
      next
    }
    Start <- unlist(lapply(curList, function(x) (x[1]-1)*25+1))
    End <- unlist(lapply(curList, function(x) (x[length(x)]*25)))
    curMat <- subset(resMat, resMat$Chr.y == names(tt)[i])
    curFDR <- lapply(curList, function(x)  curMat$curFDR[match(x, curMat$window)])
    meanFDR <- unlist(lapply(curFDR, mean))
    maxFDR <- unlist(lapply(curFDR, max))
    minFDR <- unlist(lapply(curFDR, min))
    curRatio <- lapply(curList, function(x)  curMat$curRatio[match(x, curMat$window)])
    meanRatio <- unlist(lapply(curRatio, mean))
    maxRatio <- unlist(lapply(curRatio, max))
    minRatio <- unlist(lapply(curRatio, min))
    windowNumber <- unlist(lapply(curList, length))
    curPeaks <- cbind(names(tt)[i], Start, End, windowNumber,
                      meanFDR, maxFDR, minFDR,
                      meanRatio, maxRatio, minRatio)
    resPeaks <- rbind(resPeaks, curPeaks)
  }
  colnames(resPeaks) <- c("Chromosome", "Start(1-based)", "End", "Bin number",
                          "Mean FDR", "Max FDR", "Minimum FDR",
                          "Mean Ratio", "Max Ratio", "Minimum Ratio")
  resPeaks
}

#' @export
peakCalling <- function(IPBAM, inputBAM, GTF, expName = NULL,
                        method = c("SlidingWindow","exomePeak", "MetPeak", 
                                   "MACS2", "BayesPeak"), 
                        paired = F, ...){
  
  
  
  if(length(method) != 1){
    method <- method[1]
    cat("Note: multiple peak calling methods are provided, the first one
        will be used!")
  }
  
  if(is.null(expName)){
    expName <- "example"
  }
  
  if(method == "exomePeak"){
    cat("Start peak calling using exomePeak...", "\n")
    results <- exomepeak(GENE_ANNO_GTF = GTF, IP_BAM = IPBAM,
                         INPUT_BAM = inputBAM, ...)
    peaks <- read.table(paste0(getwd(), '/exomePeak_output/peak.bed'), 
                        sep = "\t", quote = "", stringsAsFactors = F)
    #peaks <- results$all_peaks
    
  }else if(method == "BayesPeak"){
    cat("Start peak calling using BayesPeak...", "\n")
    IPBed <- .bam2bed(BAM = IPBAM)
    inputBed <- .bam2bed(BAM = inputBAM)
    res <- bayespeak(treatment = IPBed, control = inputBed, ...)
    peaks <- res$peaks
  }else if(method == "MACS2"){
    cat("Start peak calling using MACS2...", "\n")
    MACS2Dir <- system("which macs2", intern = TRUE)
    resList <- .runMACS2(IPBAM = IPBAM, inputBAM = inputBAM, 
                         paired = paired, MACS2Dir = MACS2Dir, ...)
    peaks <- resList$peaks
    
  }else if(method == "MetPeak"){
    cat("Start peak calling using MetPeak...", "\n")
    results <- metpeak(GENE_ANNO_GTF = GTF,IP_BAM = IPBAM, INPUT_BAM = inputBAM,
                       EXPERIMENT_NAME = expName, ...)
    resDir <- paste0(getwd(), "/", expName, "/peak.bed")
    peaks <- as.matrix(read.table(file = resDir, sep = "\t", quote = "", header = F))
  }else{
    cat("Start peak calling using Fisher exact test-based sliding window...", "\n")
    peaks <- SlidingWindow(input = inputBAM, RIP = IPBAM, ...)
  }
  resList <- list(peaks = peaks, method = method)
  }



.sub <- function(number, resPos, res1){
  curNumber <- as.numeric(resPos[number])
  mapped <- as.data.frame(res1[[number]])
  if(nrow(mapped) == 0){
    ratio <- 0
  }else{
    ratio <- length(which(mapped[,5] == curNumber))/nrow(mapped)
  }
}

.pseudoURatio <- function(RNAseq, inputBAM){
  seqID <- attr(RNAseq, 'name')
  resPos <- words.pos("[Tt]", text = RNAseq)
  names(resPos) <- seqID
  USiteIRange <- IRanges(resPos, resPos)
  names(USiteIRange) <- seqID
  tmpList <- list()
  tmpList[[seqID]] <- USiteIRange
  which <- RangesList(unlist(tmpList))
  p1 <- ScanBamParam(which = which, what = scanBamWhat())
  res1 <- scanBam(inputBAM, param = p1)
  resRatio <- apply(matrix(1:length(res1), nrow = length(res1), ncol = 1), 1, .sub, 
                    resPos = resPos, res1 = res1)
  
  resMat <- matrix(NA, nrow = length(resRatio), ncol = 3)
  resMat[,1] <- seqID
  resMat[,2] <- resPos
  resMat[,3] <- resRatio
  colnames(resMat) <- c("Transcript", "U.Position", "PseudoU.Ratio")
  resMat
}


#' @export
pseudoURatio <- function(refGenome, inputBAM, cpus = 1){
  refGenome <- read.fasta(refGenome, as.string = T)
  resList <- vector(mode = "list", length = length(refGenome))
  names(resList) <- names(refGenome)
  if(cpus == 1){
    for(i in 1:length(resList)){
      resList[[i]] <- .pseudoURatio(RNAseq = refGenome[[i]],
                                    inputBAM = inputBAM)
    }
  }else{
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".pseudoURatio", namespace = "PEA")
    sfExport(".sub", namespace = "PEA")
    sfLibrary("Rsamtools", character.only = TRUE)
    sfLibrary("seqinr", character.only = TRUE)
    sfLibrary("data.table", character.only = TRUE)
    resList <- sfLapply(refGenome, .pseudoURatio, inputBAM = inputBAM)
    sfStop()
  }
  resMat <- do.call(rbind, resList)
  resMat
}


#' @export
CMRCalling <- function(CMR = c("m6A", "m6Am", "m5C", "hm5C", "pseudoU", "m1A"),
                       cpus = 1, IPBAM = NULL, inputBAM = NULL, 
                       GTF = NULL,  paired = FALSE, ...){
  
  if(length(CMR) != 1){
    stop("Please specify a type of CMR!")
  }
  
  if(is.element(CMR, c("m6A", "m6Am", "m1A"))){
    resMat <- peakCalling(IPBAM = IPBAM, inputBAM = inputBAM, cpus = cpus, GTF = GTF, ...)
    resMat <- resMat$peaks
  }else if(is.element(CMR, c("m5C", "hm5C"))){
    resMat <- .m5Call(paired = paired, cpus = cpus, GTF = GTF, ...)
    #resMat <- resMat$peaks
  }else{
    resMat <- pseudoURatio(refGenome = refGenome, inputBAM = inputBAM, cpus = cpus)
    # resMat <- resMat$peaks
    
  }
  
  resMat
}

RIP.bed <- opt$IPBED
input.bed <- opt$INPUTBED
IPVex <- as.numeric(opt$IPBG) #as.numeric(system(paste0("cat ", RIP.bed, " | wc -l "), intern = T))
inputVec <- as.numeric(opt$INPUTBG) #as.numeric(system(paste0("cat ", input.bed, " | wc -l "), intern = T))
cat(IPVex, inputVec, opt$GENOMESIZE)

cmrMat <- CMRCalling(CMR = 'm6A', cpus = 15, IPBAM = RIP.bed, inputBAM = input.bed, 
                     method = 'SlidingWindow', mappedInput = inputVec,
                     mappedRIP = IPVex, refGenome = opt$GENOMESIZE)

cmrMat[,2] <- format(as.numeric(cmrMat[,2]) , scientific=F)
cmrMat[,3] <- format(as.numeric(cmrMat[,3]) , scientific=F)

write.table(cmrMat, file = opt$OUTFILE,
            sep = '\t', quote = F, row.names = F, col.names = F)
