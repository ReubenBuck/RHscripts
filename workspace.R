
rm(list = ls())

library(illuminaio)
source("~/Desktop/RHscripts/functionsRH.R")


options(stringsAsFactors = FALSE)


sampMap <- read.delim("~/Desktop/Illumina RH data/Lyons Feline 19aug2011/Sample_Map/Sample_Map.txt")




idatFiles <- list.files("~/Desktop/Illumina RH data/Lyons Feline 19aug2011/UCD Feline Leslie Lyons idats 15aug2011/UCD Feline Leslie Lyons idats 15aug2011", 
                        full.names = TRUE,recursive = TRUE)
idatList <- readIDATFiles(idatFiles)

aInte <- snpRawIntensity(idatData = idatList)
aGeno <- snpRawGeno(idatData = idatList)

layout(1)
plotRawSnp(log10(aInte), aGeno, snpRow = 100, log = FALSE)


clones <- c(1:131,133:137, 139:143)

# the classic way of genotyping/removing problem snps
plot(density(log10(aInte[,"6159231072_R06C01"])))
lines(density(log10((aInte[,c(1:131,133:137, 139:143)]))), lty = 2)
lines(density(log10(aInte[,"6159231062_R06C02"])), col= 2)
lines(density(log10(aInte[,"6159231072_R06C02"])), col= 2)
lines(density(log10(aInte[,"6159231077_R01C01"])), col= 2)



plot(density(log10(aInte[,"6159231072_R06C01"]), bw = .01), lwd = 2, col = 4)
lines(density(log10((aInte[,clones])), bw = .01), lty = 1, col = 1, lwd = 2)
lines(density(log10(aInte[,"6159231062_R06C02"]), bw = .01), col= 2, lwd = 2)
lines(density(log10(aInte[,"6159231072_R06C02"]), bw = .01), col= 2, lwd = 2)
lines(density(log10(aInte[,"6159231077_R01C01"]), bw = .01), col= 2, lwd = 2)
for(i in clones){
lines(density(log10((aInte[,i])), bw = .01), lty = 1, col = scales::alpha(colour = 1, alpha = .05))
}

samp <- sample(x = 1:nrow(aInte), size = 10000,replace = FALSE)
plot(ecdf(log10(aInte[samp,"6159231072_R06C01"])), lwd = 3, col = 4, main = "", xlab = "probe intensity", ylab = "cumulative density")
lines(ecdf(log10((aInte[samp,clones]))), lty = 1, col = 1, lwd = 2)
lines(ecdf(log10(aInte[samp,"6159231062_R06C02"])), col= 3, lwd = 2)
lines(ecdf(log10(aInte[samp,"6159231072_R06C02"])), col= 3, lwd = 2)
lines(ecdf(log10(aInte[samp,"6159231077_R01C01"])), col= 3, lwd = 2)
for(i in clones){
  lines(ecdf(log10((aInte[samp,i]))), lty = 1, col = scales::alpha(colour = 1, alpha = .05))
}



# how do we get rid of the cat ones
length(aInte[log10(aInte[,"6159231072_R06C01"]) > 2.8,"6159231072_R06C01"])
abline(v = 2.8)


# filter out high scoring negative controls and low scoring positive controls
# how best to do that?
# is there some distribution these are expected to follow to identify out liers
boxplot(log10(aInte[,"6159231072_R06C01"]))


x <- log10(aInte[,"6159231072_R06C01"])


maxNeg <- boxplot.stats(x)$stats[5]
abline(v = maxNeg)

outs <- which(x > maxNeg)


# remove snps with a score to high for negative control and snps with a score too low for each positive contril


# they plotted both of these out 
# picked all the snps that were nonoverlapping 
# and then made their calls 

# there ma be some variation in SNP integration


# Other papers used quite simple genotyping methods
# we could be a little more sophisticated
# also what is the cost/benfit of such an approach



# is the difference in signal intensity constant?

i = 100

sel <- aInte[,1:i]

selSort <- apply(sel, 2, sort)
selSort <- selSort[as.integer(seq(from = 1, to = nrow(selSort), length.out = 1000)),]
selSort <- log10(selSort)

percentCat <- NULL
for(i in 1:ncol(selSort)){
  percentCat <- c(percentCat, which(selSort[,i] > 3.5)[1])
}

selSort <- selSort[,order(percentCat)]

selDist <- NULL

eGrid <- expand.grid(1:i,1:i)
selDist <- abs(selSort[,eGrid$Var1] - selSort[,eGrid$Var2])



layout(1:2)
par(mar= c(1,1,1,1))
image((t(selDist)), col = heat.colors(100))
image(t(selSort), col = heat.colors(100))

plot(apply(selDist, 1, sd))
boxplot(t(selDist))

#selDist <- apply(selSort,2, dist)

# is there sidderence in clone score between activated snps and non activated snps



# is this some sort of balancing thing 
# what is the probability 























hist((rowMeans(log10(aInte[,1:100]))), breaks = 100)
hist((apply(log10(aInte[,1:100]), MARGIN = 1, FUN = median)), breaks = 100)



quants <- apply(aDist, 2, quantile, probs = seq(.001,1,.001))

qOrder <- order(rowSums(quants[,1:100]))
scQuants <- scale(t(quants[,1:100]))


image(scQuants[order(rowSums(scQuants)),], col = heat.colors(100))


image(t(log10(quants)[,1:100])[order(rowSums(scQuants)),], col = heat.colors(100))



aInteClones <- aInte[,sampMap$SentrixPosition[grepl("C[0-9][0-9][0-9]",sampMap$Name)]]
aGenoClones <- aGeno[,sampMap$SentrixPosition[grepl("C[0-9][0-9][0-9]",sampMap$Name)]]

# we should try clustering our samples at k = 2

plotRawSnp(snpI = aInteClones, snpG = aGenoClones, snpRow = 10)

stripchart(data.frame(t(log10(aInteClones[1:50,]))), cex = .1, yaxt = "n");axis(side = 2, at = 1:50)


clustFun <- function(snpIntensity){
  kClust <- kmeans(snpIntensity, centers = 2)
  kClust$cluster[kClust$cluster == which(kClust$centers == min(kClust$centers))] <- 0
  kClust$cluster[kClust$cluster == which(kClust$centers == max(kClust$centers))] <- 1
  return(list(cluster = kClust$cluster + 1, totss = kClust$totss))
}

cols <- unlist(as.list(t(apply(aInteClones[1:50,], 1,clustFun ))))


pdf(file = "Desktop/cutoffRH.pdf")
layout(1:2, heights = c(1,3))
par(mar = c(0,5,4,3))
plot(density(log10(aInteClones[1:length(aInteClones)])), xlim = c(2.2,4.8), xaxt = "n")
abline(v = 3.35)


snpChoice <- sort(sample(x = 1:63000, 50))

k <- t(apply(log10(aInteClones[snpChoice,] + 1), 1,clustFun ))
cols <- sapply(k, "[[", "cluster")

par(mar = c(5,5,0,3))
plot(y = rep(1:length(snpChoice),ncol(aInteClones)),
     x = unlist(as.list((log10(aInteClones[snpChoice,])))),
     #col = unlist(as.list(t(cols))), 
     col = scales::alpha("black", .2), pch = 16,
     cex = .5 , xlim = c(2.2,4.8),
     xlab = "Intensity (log10)", ylab = "SNP"
     );abline(v = 3.35)

dev.off()
# how robust is map construction

# we can measure clustering scores too
# within distance and between distance

# the idea is that the distributions should overlap


# if we were to get rid of problem snps
# we could follow old methods 
# anything in our hamsterdata that overlaps with a cat score we remove. 
# then we kame a series of 

k <- t(apply(log10(aInteClones + 1), 1,clustFun ))
hist(sapply(k, "[[", "totss"), breaks = 100)

# we can probably use totss to measure accuracy 

# how do we find whihc clones need to be adjusted

# and which way they need to be adjusted

# the approximate frequency a

# how often is decleared as 2 and how often is it decleared as one


# ths might tell us which way it needs to be pushed

# need to report its distance from center 1


clustFun <- function(snpIntensity){
  kClust <- kmeans(snpIntensity, centers = 2)
  kClust$cluster[kClust$cluster == which(kClust$centers == min(kClust$centers))] <- 0
  kClust$cluster[kClust$cluster == which(kClust$centers == max(kClust$centers))] <- 1
  
  kClust$centers <- sort(kClust$centers)
  
  return(list(cluster = kClust$cluster + 1, 
              totss = kClust$totss, 
              clust1Dist = snpIntensity - kClust$centers[1],
              clust2Dist = snpIntensity - kClust$centers[2],
              snpRank = rank(snpIntensity)))
}

k <- t(apply(log10(aInteClones + 1), 1,clustFun ))
hist(sapply(k, "[[", "totss"), breaks = 100)

clust1 <- t(sapply(k, "[[", "clust1Dist"))
clust2 <- t(sapply(k, "[[", "clust2Dist"))
snpRank <- t()


boxplot(clust1[,1],clust2[,1], clust1[,2], clust2[,2])


# we can determine their bias and nudge them according to a step size

# it might be better to use rank
# we expact the average rank to be middle









snpGroups <- apply(k, 1, table)
cloneGroups <- apply(k,2,table)


plot(density(snpGroups[1,]/ncol(aInteClones)), xlim = c(0,1))
lines(density(snpGroups[2,]/ncol(aInteClones)), col = 2)

plot(density(cloneGroups[1,]/nrow(aInteClones)), xlim = c(0,1))
lines(density(cloneGroups[2,]/nrow(aInteClones)), col = 2)


# are there clones with a higher rate than expected?


#cluster snps by clones

k2 <- apply(aInteClones,2,clustFun)


sum((k2 != k)[k2 == 1])


sum((k2 != k)[k2 == 2])


df <- data.frame(k = unlist(as.list(k)), k2 = unlist(as.list(k2)))

tDF <- table(df)

fisher.test(tDF)


plot(colSums(k != k2))
plot(density(rowSums(k != k2)))








# how many times is each snp seen?
# we probably need to do a filtering at this stage 




# Here we want to figure out linkage stats from the called genotypes

# there is also retentrion frequencies of our snps/clones

# the data will be a i by j matrix with i snps and j colones
# for presence of a snp within a clone it gets a 1, for absence it gets a zero

# From this table we would like to get pairwise linkage statistics 
# the apriwise linkage stats will be the lod scores
# we should set a distace threshold
# a resolution step
# is it usful to use the local maxima 

# using the hits package we could probably sort into a graph with single linkage
# try and work out wich markers belong to whcih group

# and then sort them into linkage groups

# can use matrix multiplication between rows to get pairwise frequencies
# Estimated 

# create the two matrixceis

set.seed(1)
mat <- matrix(sample(x = c(0,1), size = 10*100, replace = TRUE), nrow = 100, ncol = 10, 
              dimnames = list(snp = paste("snp", 1:100, sep = ""), clone = paste("c", 1:10, sep = "")))

lMat <- as.list(as.data.frame(t(mat)))
lMat <- lapply(lMat, matrix, nrow = 1)
lMat <- as.matrix(lMat)

outer(as.array(lMat), as.array(lMat), FUN="%*%")

recom <- mat[1,] %*% mat[2,]
nonRecom <- sum(mat[1,] + mat[2,] == 1)

recom <- list(mat[1,], mat[2,]) %*% list(mat[2,], mat[3,])


library(tensorA)
matT <- to.tensor(mat)



m<-1:4
x<-list(m,m,m,m,m)

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
  



mat <- matrix(sample(x = c(0,1), size = 130*60000, replace = TRUE, prob = c(.8,.2)), nrow = 60000, ncol = 130, 
              dimnames = list(snp = paste("snp", 1:60000, sep = ""), clone = paste("c", 1:130, sep = "")))

set.seed(1)
samp <- sample(x = 1:length(lMat), size = 100)

lMat <- matrix2list(mat)
counts <- alleleCount(lMat[samp])
lods <- lod(counts$AB, counts$A, counts$B, counts$O)
g1 <- linkGroups(lodMat, 3)

# next function is to get distances


# now we can calculate lod too




# how am i gettting snps mapping to more than one ID 


hist(sapply(g1, length))


unlist(g1)[duplicated(unlist(g1))]


# no we know which elements we can workout if there is any shared nodes
# need to assign them to groups
# one way is with a for loop

# would we like a list?
# each list contains the nodes that belong to each group

a <- matrix(1:100, nrow = 20)

x = 35
a[x%%nrow(a), ceiling(x/nrow(a))]


lodMad[is.na(lodMad)]

ab1 <- ab[is.na(lodMad)]
a1 <- a[is.na(lodMad)]
b1 <- b[is.na(lodMad)]
o1 <- o[is.na(lodMad)]

lodMad[ab == 6 & a == 20 & b == 24 & o == 80]
lodMad[ab == 6 & a == 24 & b == 20 & o == 80]

# that is peculiar, some number combinations end up with NAs.

# there is still a lot of pairs to consider.

# what if we take a random sample?
# surely we shoudl be able to get linkage groups out of those

# if we set some sort of seed



lodMad <- lodMad[upper.tri(lodMad)]
hist(lodMad)

# with random data we see no linked loci 
# interesting

# the next thing is to pick a threshold and build linkage groups from that







