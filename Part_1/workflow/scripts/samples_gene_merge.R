
library(reshape2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

filePath <- args[1]
detailName <- args[2]
outFile <- args[3]
fileList <- list.files(filePath)
fileList <- fileList[grepl("gene.txt", fileList)]
fileIndex <- grepl(detailName, fileList)
colNum <- 9
fileList <- fileList[fileIndex]

for(i in 1:length(fileList)){
  fileName <- fileList[i]
  fileInput <- read.delim2(paste0(filePath, "/", fileName), sep = "\t", 
                          header = T, stringsAsFactors = F)
  fileInput$Rep <- rep(paste0("R", i), nrow(fileInput))
  fileInput <- fileInput[as.numeric(fileInput[, colNum])>0, ]
  nameTable <- table(fileInput[,1])[table(fileInput[,1])>1]
  nameRepIndex <- vector()
  if(length(nameTable)>0){
    for(nameRep in names(nameTable)){
      nameRepIndex <- c(nameRepIndex, which(fileInput[,1]%in%nameRep)[-1])
    }
    fileInput <- fileInput[-nameRepIndex, ]
  }
  cat(i, nrow(fileInput), "\n")
  if(i == 1){
    tmpDf <- fileInput
  }else{
    tmpDf <- rbind(tmpDf, fileInput)
  }
}
tmpDfll <- tmpDf[, c(1,10,colNum)]

dataMerge <- tmpDfll%>% 
  group_by(Rep, Gene.ID) %>%
  mutate(id = row_number()) %>% 
  dcast(Gene.ID+id ~ Rep, value.var = colnames(tmpDfll)[3]) %>%
  select(-id)

rownames(dataMerge) <- dataMerge[,1]
dataMerge <- dataMerge[, -1]
dataMerge[is.na(dataMerge)] <- 0

dataMerge$expMean <- apply(dataMerge, 1, function(x){mean(as.numeric(x))})

write.table(dataMerge, outFile, sep = "\t", quote = F, col.names = F)
  