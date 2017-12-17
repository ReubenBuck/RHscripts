rm(list = ls())

library(illuminaio)
source("~/Desktop/RHscripts/functionsRH.R")


options(stringsAsFactors = FALSE)


sampMap <- read.delim("~/Desktop/Illumina RH data/Lyons Feline 19aug2011/Sample_Map/Sample_Map.txt")




idatFiles <- list.files("~/Desktop/Illumina RH data/Lyons Feline 19aug2011/UCD Feline Leslie Lyons idats 15aug2011/UCD Feline Leslie Lyons idats 15aug2011", 
                        full.names = TRUE,recursive = TRUE)
idatList <- readIDATFiles(idatFiles)

aInte <- log10(snpRawIntensity(idatData = idatList))
aGeno <- snpRawGeno(idatData = idatList)


cloneName <- sampMap$SentrixPosition[grep("C[0-9]+.", sampMap$Name)]
negControl <- sampMap$SentrixPosition[sampMap$Name == "A23"]
posControl <- sampMap$SentrixPosition[c(grep("12682", sampMap$Name), grep("Aby", sampMap$Name)) ]




cloneInte <- aInte[,cloneName]
cloneGeno <- aGeno[,cloneName]
posInte <- aInte[,posControl]
posGeno <- aGeno[,posControl]
negInte <- aInte[,negControl]
negGeno <- aGeno[,negControl]


negOut <- which(negInte >= boxplot.stats(negInte)$stats[5])

posOut <- NULL
for(i in 1:ncol(posInte)){
  posOut <- c(posOut, list(which(posInte[,i] <= boxplot.stats(posInte[,i])$stats[1])))
}
names(posOut) <- colnames(posInte)
venn::venn(posOut)
posOut <- unique(unlist(posOut))

# now we can start to make example plots of what the data actually looks like 

allInte <- cbind(negInte,posInte,cloneInte)
allGeno <- cbind(negGeno, posGeno, cloneGeno)
allColours <- c("blue", rep("red", ncol(posInte)),
                rep(scales::alpha("black", .2), ncol(cloneInte)))

dev.off()
plotRawSnp(allInte, allGeno, snpRow = (1:nrow(allInte))[-unique(c(posOut, negOut))][50011], log = FALSE,col = allColours)
plotRawSnp(allInte, allGeno, snpRow = negOut[1], log = FALSE,col = allColours)


# filter out clones 
# maybe we could use tukey outlier coeffcient 
removeRow <- unique(c(row(allInte)[is.infinite(allInte)], negOut, posOut))

cloneInteFilt <- cloneInte[-removeRow,]
cloneGenoFilt <- cloneGeno[-removeRow,]


dim(cloneInteFilt)

## After filtering, I measure lod scores, breakage, and distance across a range of min thresholds

s <- sample(x = 1:nrow(cloneType), size = 100)

lodList <- NULL
distList <- NULL
breakList <- NULL
countsList <- NULL
for(i in seq(3,4,.05)){
  
  cloneType <- cloneInteFilt
  cloneType[cloneInteFilt < i] <- 0
  cloneType[cloneType > 0] <- 1
  cloneTypeSamp <- cloneType[s,]
  
  
  lMat <- matrix2list(cloneTypeSamp)
  counts <- alleleCount(lMat)
  lods <- lod(counts$AB, counts$A, counts$B, counts$O)
  mDist <- markerDist(counts$AB, counts$A, counts$B, counts$O)
  breaks <- breakage(counts$AB, counts$A, counts$B, counts$O)
  
  countsList <- c(countsList, list(counts[upper.tri(counts)]))
  distList <- c(distList, list(mDist[upper.tri(mDist)]))
  lodList <- c(lodList, list(lods[upper.tri(lods)]))
  breakList <- c(breakList, list(breaks[upper.tri(breaks)]))
}



hist(breakList[[1]] , breaks = seq(0,100,.03), xlim = c(0,1.5), col = scales::alpha(1, .2), border = NA,
     main = "", xlab = "breakage frequency")
hist(breakList[[6]], breaks = seq(0,100,.03), xlim = c(0,1.5), add = TRUE, col = scales::alpha(2, .2), border = NA)
hist(breakList[[14]], breaks = seq(0,100,.03), xlim = c(0,1.5), add = TRUE, col = scales::alpha(3, .2), border = NA)
hist(breakList[[21]], breaks = seq(0,100,.03), xlim = c(0,1.5), add = TRUE, col = scales::alpha(4, .2), border = NA)


hist(breakList[[1]][lodList[[1]] > 5] , breaks = seq(0,100,.03), xlim = c(0,1.5), col = scales::alpha(1, .2), border = NA,
     main = "", xlab = "breakage frequency", ylim = c(0,10))
hist(breakList[[6]][lodList[[6]] > 5], breaks = seq(0,100,.03), xlim = c(0,1.5), add = TRUE, col = scales::alpha(2, .2), border = NA)
hist(breakList[[14]][lodList[[14]] > 5], breaks = seq(0,100,.03), xlim = c(0,1.5), add = TRUE, col = scales::alpha(3, .2), border = NA)
hist(breakList[[21]][lodList[[21]] > 5], breaks = seq(0,100,.03), xlim = c(0,1.5), add = TRUE, col = scales::alpha(4, .2), border = NA)

hist(distList[[1]] , breaks = seq(0,100,.1), xlim = c(0,8), col = scales::alpha(1, .2), border = NA,
     main = "", xlab = "cR(15000)")
hist(distList[[6]], breaks = seq(0,100,.1), xlim = c(0,8), add = TRUE, col = scales::alpha(2, .2), border = NA)
hist(distList[[14]], breaks = seq(0,100,.1), xlim = c(0,8), add = TRUE, col = scales::alpha(3, .2), border = NA)
hist(distList[[21]], breaks = seq(0,100,.1), xlim = c(0,8), add = TRUE, col = scales::alpha(4, .2), border = NA)

hist(distList[[1]][lodList[[1]] > 5] , breaks = seq(0,100,.1), xlim = c(0,.7), col = scales::alpha(1, .2), border = NA,
     main = "", xlab = "cR(15000)", ylim= c(0,5))
hist(distList[[6]][lodList[[6]] > 3], breaks = seq(0,100,.01), xlim = c(0,8), add = TRUE, col = scales::alpha(2, .2), border = NA)
hist(distList[[14]][lodList[[14]] > 3], breaks = seq(0,100,.01), xlim = c(0,8), add = TRUE, col = scales::alpha(3, .2), border = NA)
hist(distList[[21]][lodList[[21]] > 3], breaks = seq(0,100,.01), xlim = c(0,8), add = TRUE, col = scales::alpha(4, .2), border = NA)

# if we were to sample a series of points along a string,
# what kind of distances would we expect?
# an important paramater is the distance between markers
# this small dataset works by calculating known distances between a set of markers,
# distances are converted to breakage frequencies which are then used to recalculate distances
# the important step here is that we remove distances that would now be observable 
# based on what we would calcualte form breakage frequencies alone

# the average proximity between ajacent markers
simMarkerCoeff <- .001
n = 100
simDist0 <- dist(sample(1:nrow(cloneInteFilt), n, replace = FALSE) + runif(n = n, min = -.5, max = -.5)) * simMarkerCoeff
simBreakage <- 1 - exp(-simDist0)
simDist1 <- -log(1 - simBreakage[simBreakage < .5])

hist(simBreakage, breaks = 20)
hist(simDist0, breaks = 20)
hist(simDist1, breaks = 100)
1 - (sum(is.infinite(simDist1))/((n^2 - n) / 2))
abline(h = sum(!is.infinite(simDist1))/20)
# now we have an idea of how the distance distribution is expected to behave

# using x clones can we calculate the relevant lod socres for certain breakage frequencies.
# what is a score that is likely to indicate linkage
# in addition no matter what breakage frequency we use, there seems to be a limit on the depth we can calculate between two markers
# this probably has something to do with taking the log of certain numbers



# now we know what expected distributions can look like
# we have to think about sources of bias in our calculation. 
# which clones are giving us scores greater than one and what are their inputs
# my suspicion is that AB is lower than what would be expected by chance for tow unlinked markers


i == 8

x = which(breakList[[i]] > 1)
inputs <- lapply(countsList[[i]], "[[", x[1002])
inputs

breakage(inputs$AB, inputs$A, inputs$B, inputs$O)

AB = inputs$AB; A = inputs$A; B = inputs$B; O = inputs$O
Tall = AB + A + B + O
Tall
Ra = (A + AB)/Tall
Ra
Rb = (B + AB)/Tall
Rb
theta = (A + B)/(Tall * (Ra + Rb - (2*Ra*Rb)))
theta
Tall * (Ra*Rb)
Tall * Ra
Tall * Rb
# for one of the examples both A and B are being observed too many times
# they have retention frequenies at about 50%

# based on their retention frequencies, we would expect to see them together more based on random chance
# why are these retention frequency issues causing prolems

# is this because of clones or snps
# what is the liklihood this is caused by too many false positives

# potentially we are miscalling about 2 or 3

# however I do have a way to measure the expected ratios acorss our data

# it seems that breakage frequency estimates are not vary robust
# an estimate error of 1 or 2 clones can give you the wrong breakage frequency quite eaily



## I need to justify why these scores work, even though I know 



### An alternative genotyping procedure
cloneCluster <- clusterStats(cloneInteFilt, type = "clone")
hist(rowSums(cloneCluster$cluster ))
hist(colSums(cloneCluster$cluster ))
plot(sort(rowSums(cloneCluster$cluster)))


# do determine whihc cluster levels should be used
layout(matrix(1:9, nrow = 3, byrow = TRUE))
plot(colSums(cloneInteFilt * (cloneCluster$cluster))/colSums(cloneCluster$cluster),
     cloneCluster$var.explained , 
     xlim = c(2.3,4.7), col = 2, ylim = c(0,1), xlab = "mean probe intensity of active cluster",
     ylab = "var explained")
q = identify(colSums(cloneInteFilt * (cloneCluster$cluster))/colSums(cloneCluster$cluster),
             cloneCluster$var.explained, n = 8 )
for(i in q){
  plot(density(cloneInteFilt[,i]), main = i, xlim = c(2.3,4.7))
  abline(v = min(cloneInteFilt[,i][cloneCluster$cluster[,i] == 1]) )
}

# chose a threshold to remove clones
cloneKeep = which(cloneCluster$var.explained > .8)

# filter out clones below threshold
cloneInteFilt <- cloneInteFilt[,cloneKeep]
cloneGenoFilt <- cloneGenoFilt[,cloneKeep]


# re analyse clone clustering
cloneCluster <- clusterStats(cloneInteFilt, type = "clone")

# ctreate snp clustering on filtered clones
snpCluster <- clusterStats(cloneInteFilt, type = "snp")


layout(matrix(1:9, nrow = 3, byrow = TRUE))
smoothScatter(snpCluster$var.explained, (rowSums(snpCluster$cluster)/ncol(snpCluster$cluster)),
              ylab = "retention frequency", xlab = "var explained", main = "snps")
a <- identify(snpCluster$var.explained, (rowSums(snpCluster$cluster)/ncol(snpCluster$cluster)), n = 7)
plot.new()
legend("center", legend = c("both", "snp only", "clone only", "neither"), 
       fill  = c(5,3,2,1), bty = "n", title = "active cluster", cex = 2)
for(i in a){
  plotRawSnp(snpI = cloneInteFilt, snpG = cloneGenoFilt, snpRow = i, col =(cloneCluster$cluster[i,] + 2) * (snpCluster$cluster[i,] + 1) -1, main = i)
}


# from here i can begin to choose some of my clustering stats 

# .8 and .4

snpKeep <- snpCluster$cluster[snpCluster$var.explained > .8 & rowSums(snpCluster$cluster)/ncol(snpCluster$cluster) < .4, ]


s <- sample(1:nrow(snpKeep), size = 400)
cloneTypeSamp <- snpKeep[s,]


lMat <- matrix2list(cloneTypeSamp)
counts <- alleleCount(lMat)
lods <- lod(counts$AB, counts$A, counts$B, counts$O)
mDist <- markerDist(counts$AB, counts$A, counts$B, counts$O)
breaks <- breakage(counts$AB, counts$A, counts$B, counts$O)

lods <- lods[upper.tri(lods)]
mDist <- mDist[upper.tri(mDist)] * 100
breaks <- breaks[upper.tri(breaks)]


layout(matrix(1:4, nrow = 2))

hist(breaks, breaks = 100)
hist(lods, breaks = 100)
hist(mDist, breaks = 100)

plot(mDist, lods, col = scales::alpha(1, .1), pch = 16, xlim = c(0,800))
plot(breaks, lods, col = scales::alpha(1, .1), pch = 16)
smoothScatter(breaks, lods)
plot(breaks, mDist) # of course one is just taking the log of the other

hist(mDist[lods > 5], breaks = 100, xlim = c(0, 100))

# we know the threshold tells us what the liklihood is that the two markers are linked
# where linked means they are on the same fragment

# I'm pretty happy with this so far, but is there some way we can be sure of our calls, without going through with building a map?



# markers that are all within 80 cR(15000)

# based on
# how many individual markers are we seeing?
# build a matrix and extract marker pairs that pass lod thresholds 


# so maybe we can look at linkage groups
# the question is, how do we know if they are good.
# we know a certain dist score or breakage frequency corresponds to a specific lod score
# we can compare our sample to that
# we could use the estimted distance to 




# looks like there are a bunch of markers that are linked


# this needs to be fixed



# alternatively we could simulate the process of breaking the genome and putting it into clones,
# we can use that to caluclate lod scores
# There probably needs to be some liklihood estimate around the break frequencies
# That is brobably what the lod score is



# to run a simulation, we can test different breakage frequencies 
# the proportion of the time that two adjacent markers will be seen in the same clone

# if we sample markers and put them into clones
# there is a particular sampling distribution for each marker
# esentially that distribution could be permuted acording to the breakage frequency between adjacent markers
# it could be done in a way that after 5 permutations the breakage frequency approciamtes to 1
# at first lets try a rearrangement 

# the hamming distance is a way to measure the difference between the permutations
# we could do a bunch of purmtations and then sort them according to hamming distance to get a sense of marker arrangement
# is there some way in which this can produce an undesirable result?
# I think the best way is to generate each new permutation by performing a series of know permutations on the last ordering

# we could pick a certain number to be zeros and a certain number to be ones 
# sample these from a distribution that corresponds to our observed retention frequency

# first work out the requried hamming distance 


A <- sample(0:1, size = 130, replace = TRUE, prob = c(.8,.2))
B <- A
n = 45
b <- sample(1:130,n)
B[b] <- (B[b] + sample(0:1,size = n, replace = TRUE, prob = c(.8,.2)))%%2

A <- matrix(sample(0:1, size = 130, replace = TRUE, prob = c(.8,.2)), nrow = 1)
n = 1
for(i in 1:9000){
  B <- A[nrow(A),]
  B[sample(which(A[nrow(A),] == 0), size = n)] <- 1
  B[sample(which(A[nrow(A),] == 1), size = n)] <- 0
  #B[b] <- (B[b] + sample(1:0,size = n, replace = TRUE, prob = c(1-(sum(B)/length(B)) ,sum(B)/length(B))))%%2
  #B[b] <- (B[b] + sample(1:0,size = n, replace = TRUE, prob = c(.2,.8)))%%2
  A <- rbind(A,B)
}
plot(rowSums(A))

# how do we keep the frequencies at about 20%
# currently it very quickly converged to .5
# if its 0 the probability of it changing should be

count <- alleleCount(matrix2list(A[sample(1:nrow(A), 100),]))
broke <- breakage(count$AB, count$A, count$B, count$O)
hist(broke[upper.tri(broke)])

md <- markerDist(count$AB, count$A, count$B, count$O)

ll <- lod(count$AB, count$A, count$B, count$O)
ll <- ll[upper.tri(ll)]
broke <- broke[upper.tri(broke)]
smoothScatter(broke, ll)
# need to randomise the number of markers that are vissible










