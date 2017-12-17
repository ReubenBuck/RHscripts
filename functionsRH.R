
library(illuminaio)

### Used to read in data
### Simply provide all the file names and it will create a convinant list obkect that is pretty efficien
readIDATFiles <- function(idatFiles = NA){
  idatList = NULL
  for(i in 1:length(idatFiles)){
    idatFile <- idatFiles[i]
    if(!grepl("\\.idat", idatFile))
      next()
    idatData <- readIDAT(file = idatFile)
    idatData0 <- data.frame(midBlock = idatData$MidBlock,idatData$Quants)
    barcodeID <- paste(idatData$Barcode, idatData$Unknowns$MostlyA, sep = "_")
    rownames(idatData0) <- NULL
    if( is.null(names(idatList)[length(idatList)]) ){
      idatList <- c(idatList, list(NULL))
      names(idatList)[length(idatList)] <- barcodeID
    }
    if(grepl("Grn", idatFile)){
      idatList[[barcodeID]]$Grn <- idatData0
    }else if(grepl("Red", idatFile)){
      idatList[[barcodeID]]$Red <- idatData0
    }else{
      stop("Idat file with no colour in file name")
    }
  }
  return(idatList)
}


# used to pull snps out of the list object
snpRawData <- function(idatData = NA, probeColor = NA, statistic = NA, snpRows, sentrixPositions){
  sapply(lapply(idatData, "[[", probeColor), "[[", statistic)[snpRows, sentrixPositions]
}


snpRawIntensity <- function(idatData = NA, snpRows, sentrixPositions){
  snpGrn <- snpRawData(idatData, probeColor = "Grn", statistic = "Mean", snpRows, sentrixPositions)
  snpRed <- snpRawData(idatData, probeColor = "Red", statistic = "Mean", snpRows, sentrixPositions)
  return(sqrt((snpGrn^2 + snpRed^2)))
}


snpRawGeno <- function(idatData = NA, snpRows, sentrixPositions){
  snpGrn <- snpRawData(idatData, probeColor = "Grn", statistic = "Mean", snpRows, sentrixPositions)
  snpRed <- snpRawData(idatData, probeColor = "Red", statistic = "Mean", snpRows, sentrixPositions)
  (snpGrn - snpRed)/pmax(snpGrn, snpRed)
}

plotRawSnp <- function(snpI, snpG, snpRow = NA, log = FALSE, ...){
  if(is.na(snpRow)){
    stop("Must pick a specific snp to plot")
  }
  snpG = snpG[snpRow,]
  snpI = snpI[snpRow,]
  if(log){
    snpI = log10(snpI)
    ylab = "probe intensity (log 10)"
  }else{
    ylab = "probe intensity"
  }
  plot(snpG, snpI, xaxt = "n", xlab = "genotype", ylab = ylab, xlim = c(-1,1), ..., pch = 16)
  axis(side = 1, at = c(-1,0,1), labels = c("AA", "AB", "BB"))
}


# clustering function
clustFun <- function(snpIntensity){
  kClust <- kmeans(snpIntensity, centers = 2)
  kClust$cluster[kClust$cluster == which(kClust$centers == min(kClust$centers))] <- 0
  kClust$cluster[kClust$cluster == which(kClust$centers == max(kClust$centers))] <- 1
  return(list(cluster = kClust$cluster, totss = kClust$totss, tot.withinss = kClust$tot.withinss))
}

# provide a n snp by m clones matrix
clusterStats <- function(x, type){
  if(type == "clone"){
    clust <- apply(t(x),1,clustFun)
  }else if(type == "snp"){
    clust <- apply(x,1,clustFun)
  }else{
    stop("type must be either clone or snp")
  }
  cluster <- lapply(clust, "[[", "cluster")
  if(type == "clone"){
    cluster <- do.call(cbind,lapply(cluster,matrix,nrow=nrow(cloneInteFilt),byrow=FALSE))
  }else if(type == "snp"){
    cluster <- do.call(rbind,lapply(cluster,matrix,ncol=ncol(cloneInteFilt),byrow=TRUE))
  }
  colnames(cluster) <- colnames(x)
  rownames(cluster) <- rownames(x)
  tss <- sapply(clust, "[[", "totss")
  wss <- sapply(clust, "[[", "tot.withinss")
  bss <-  tss-wss
  return(list(cluster = cluster, tss = tss, wss = wss, var.explained = bss/tss))
}



# convert allele matrix to a list
matrix2list <- function(mat){
  lapply(as.list(as.data.frame(t(mat))), matrix, nrow = 1)
}


# count allels
countAB <-function(x,y){
  sum(x == 1 & y == 1)
}

countA <- function(x,y){
  sum(x == 1 & y == 0)
}

countB <- function(x,y){
  sum(x == 0 & y == 1)
}

countO <- function(x,y){
  sum(x == 0 & y == 0)
}

alleleCount <- function(lMat){
  if(length(lMat) > 1000){
    stop("input is too large, must have length <= 120")
  }
  ab <- outer(lMat, lMat, Vectorize(countAB))
  a <-outer(lMat, lMat, Vectorize(countA))
  b <- outer(lMat, lMat, Vectorize(countB))
  o <- outer(lMat, lMat, Vectorize(countO))
  return(list(AB = ab, A = a, B = b, O = o))
}



# Get lod scores

lod <- function(AB, A, B, O){
  Tall = AB + A + B + O
  Ra = (AB + A)/Tall
  Rb = (AB + B)/Tall
  theta = (A + B)/(Tall * (Ra + Rb - (2*Ra*Rb)))
  Rab = ((AB/Tall) - (theta * Ra * Rb))/(1-theta)
  # id the breakage frequency is equal to 1, then the proportion of a and b together should be 0
  Rab[theta == 1] = 0
  L1 = ((1-theta)*Rab) + (theta * Ra * Rb)
  L2 = theta * Ra * (1 - Rb)
  L3 = theta * (1-Ra) * Rb
  L4 = ((1-theta) * (1 - Rab)) + (theta * (1-Ra) * (1-Rb))
  Lo <- (L1^AB) * (L2^A) * (L3^B) * (L4^O)
  L1a = ((1-1)*Rab) + (1 * Ra * Rb)
  L2a = 1 * Ra * (1 - Rb)
  L3a = 1 * (1-Ra) * Rb
  L4a = ((1-1) * (1 - Rab)) + (1 * (1-Ra) * (1-Rb))
  La <- (L1a^AB) * (L2a^A) * (L3a^B) * (L4a^O)
  return(log10(Lo/La))
}

markerDist <- function(AB, A, B, O){
  Tall = AB + A + B + O
  Ra = (AB + A)/Tall
  Rb = (AB + B)/Tall
  theta = (A + B)/(Tall * (Ra + Rb - (2*Ra*Rb)))
  distance <- -log(1-theta)
}

breakage <- function(AB, A, B, O){
  Tall = AB + A + B + O
  Ra = (A + AB)/Tall
  Rb = (B + AB)/Tall
  theta = (A + B)/(Tall * (Ra + Rb - (2*Ra*Rb)))
  
  return(theta)
}

# extractr linkage groups from a semetrical lod score matrix
linkGroups <- function(lodMat, threshold){
  pass <- lodMat >= threshold
  pass[lower.tri(pass,diag = TRUE)] <- FALSE
  g <- data.frame(from = rownames(pass)[which(pass)%%nrow(pass)], 
                  to = colnames(pass)[ceiling(which(pass)/nrow(pass))]
  )
  
  id <- rep(0, nrow(g))
  idCount <- 0
  for(i in 1:length(id)){
    if(any(as.matrix(g[id > 0,]) == g$from[i]) |  any(as.matrix(g[id > 0,]) == g$to[i])   ){
      idMatch <- id[id > 0][g$from[id > 0] == g$from[i] | 
                              g$to[id > 0] == g$from[i] | 
                              g$from[id > 0] == g$to[i] |
                              g$to[id > 0] == g$to[i]]
      id[i] <- idMatch[1]
      id[id %in% idMatch] <- id[i]
    }else{
      idCount = idCount + 1
      id[i] <- idCount
    }
  }
  return(lapply(split(as.matrix(g), f = as.factor(id)), unique))
}
