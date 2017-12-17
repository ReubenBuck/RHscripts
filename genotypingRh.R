## now know that lists are the best way to read in data.



idatFiles <- list.files("~/Desktop/Illumina RH data/Lyons Feline 19aug2011/UCD Feline Leslie Lyons idats 15aug2011/UCD Feline Leslie Lyons idats 15aug2011", 
                       full.names = TRUE,recursive = TRUE)



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

snpData <- function(idatData = NA, probeColor = NA, statistic = NA, snpRows, sentrixPositions){
  sapply(lapply(idatData, "[[", probeColor), "[[", statistic)[snpRows, sentrixPositions]
}




dat <- readIDATFiles(idatFiles = idatFiles)

snpChoice <- sample(x = 1:63310, size = 1000, replace = FALSE)

aGrn <- snpData(idatData = dat, probeColor = "Grn", statistic = "Mean", 
             snpRows = snpChoice)

aRed <- snpData(idatData = dat, probeColor = "Red", statistic = "Mean", 
                snpRows = snpChoice)


plot(log10(aGrn[1,]), log10(aRed[1,]))



i = 999
plot(log10(aGrn[i,]), log10(aRed[i,]))



a1 <- snpData(idatData = dat, probeColor = "Grn", statistic = "Mean")


dim(a)

head(a)


image(log10(t(a)))


barplot(colSums((a1)))




boxplot((a), log = "y", outline = FALSE)




barplot(apply(a1, 2, quantile, probs = .1)[1:100])
barplot(apply(a1, 2, quantile, probs = .2)[1:100])
barplot(apply(a1, 2, quantile, probs = .3)[1:100])
barplot(apply(a1, 2, quantile, probs = .4)[1:100])
barplot(apply(a1, 2, quantile, probs = .5)[1:100])
barplot(apply(a1, 2, quantile, probs = .6)[1:100])
barplot(apply(a1, 2, quantile, probs = .7)[1:100])
barplot(apply(a1, 2, quantile, probs = .8)[1:100])
barplot(apply(a1, 2, quantile, probs = .9)[1:100])
barplot(apply(a1, 2, quantile, probs = .95)[1:100])
barplot(apply(a1, 2, quantile, probs = .99)[1:100])






quants <- apply(a1, 2, quantile, probs = seq(0,1,.01))


image(scale(t(quants[,1:100]))[order(quants[10,1:100]),], col = heat.colors(100))

barplot(quants[90,1:100][order(quants[20,1:100])])



qOrder <- order(rowSums(quants[,1:100]))


scQuants <- scale(t(quants[,1:100]))


image(scQuants[order(rowSums(scQuants)),], col = heat.colors(100))

# is there a way to smooth over all of the orderings

# basically 


# look at the plot and see where the instability begins 
# this is where the signals begin 
# gives us a sense of how much cat DNA made it into the clone

# we are also assuming that signal intensity is constant



# to convert the signal, we want to measure distance from origin
# that way we can merge red and green



aGrn <- snpData(idatData = dat, probeColor = "Grn", statistic = "Mean")

aRed <- snpData(idatData = dat, probeColor = "Red", statistic = "Mean")


aDist <- sqrt((aGrn^2 + aRed^2))


i = 4000
plot((aGrn[i,] - aRed[i,])/pmax(aRed[i,], aGrn[i,]), (aDist[i,]), xlim = c(-1,1))






quants <- apply(aDist, 2, quantile, probs = seq(.001,1,.001))

qOrder <- order(rowSums(quants[,1:100]))
scQuants <- scale(t(quants[,1:100]))


image(scQuants[order(rowSums(scQuants)),], col = heat.colors(100))


image(t(log10(quants)[,1:100])[order(rowSums(scQuants)),], col = heat.colors(100))


# we ask at what quantile does the dat reach its peak

# this is how we can estimate how much of the cat got into the hamster
# if our signal shifts at the 60% mark, this sugests that 40% of the cat genome got into the hamster
# for what cat sites are we getting hits?

# there's variation in amoutn of 

# are normalisation factors constant?

# measuring the relative distances in the curves, maybe we can figure out what exactly is happening 
# what we would like is where the signal changes
layout(1:2)
plot((quants[,26]), type = "l")
lines((quants[,21]), type = "l", col = 2)

plot(log10((quants[,22])), type = "l")
lines(log10(quants[,21]), type = "l", col = 2)


sm <- predict(smooth.spline(((quants[,26])),all.knots = TRUE), deriv = 1)
plot(sm, type = "l", ylim = c(0, 300))

sm <- predict(smooth.spline(((quants[,26])),all.knots = TRUE), deriv = 2)
lines(sm, type = "l", ylim = c(0, .02), col = 2)




         




                               