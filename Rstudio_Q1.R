#set working directory
setwd("C:/Users/jade-/Downloads/LIFE4141_coursework_resources-20240202")
options(warn=1) #for warning message generation

#load packages
library(adegenet)
library(adegraphics)
library(vcfR)
library(pegas)
library(StAMPP)
library(ade4)
library(MASS)

#modified function - allows for tetraploids to be converted to genlight objects
vcfR2genlight.tetra <- function(x, n.cores = 1)
{
  bi <- is.biallelic(x)
  if (sum(!bi > 0)) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two allels"))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object"))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}

#PCA calculation on genlight objects
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit = alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar<- glVar(x, alleleAsUnit = alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  
#convert to full data. dividing by ploidy data keeps the NA values
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x) 
  
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  
#center and scale
  mx <- scale(mx,
             center = if(center) vecMeans else F, 
              scale = if (scale) vecVar else F)
             
  allProd <- tcrossprod(mx) / nInd(x)
  
#perform eigen analysis
  eigRes <- eigen(allProd, symmetric = TRUE, only.values = FALSE)
  rank = sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[,1:rank, drop = FALSE]
  
  if(is.null(nf)){
    barplot(eigRes$values, main = "Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n=1))
  }
  
#rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x))
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  
  if(loadings){
    if(scale){
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow = nLoc(x), ncol=nf)
    
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center){
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.a(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadingss <- res$loadings + matrix(temp) %*% eigRes$vectors[k,1:nf, drop = FALSE]
    }
    res$loadings <- res$loadings / nInd(x)
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  
#format output
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownaes(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
     colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x), alles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
   }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  } 
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}
  
#import SNP data from VCF file
vcf <- read.vcfR("LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz")

#convert vcf to genlight object
aa.genlight <- vcfR2genlight.tetra(vcf) #using modified function
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") #add real SNP names
pop(aa.genlight) <- substr(indNames(aa.genlight),1,3) #add pop names (first 3 chars)

#check on made objects
aa.genlight
indNames(aa.genlight)
ploidy(aa.genlight)

#run PCA
pca.1 <- glPcaFast(aa.genlight, nf=300) #use modified function

#plot
scatter(pca.1, posi="bottomright") #plot scatter w individual labels
loadingplot(pca.1)

#proportion of variance by the first 3 axis
pca.1$eig[1]/sum(pca.1$eig)
pca.1$eig[2]/sum(pca.1$eig)
pca.1$eig[3]/sum(pca.1$eig)

#colour populations by using a palette
col <- funky(10)
s.class(pca.1$scores, pop(aa.genlight), xax = 1, yax = 2, col = transp(col,.6), ellipseSize = 0, starSize = 0, ppoints.cex = 4, paxes.draw = T, pgrid.draw = F)

#save the figures into a single pdf file
pdf("PCA_all_SNPs_ax12_1k_less.pdf", width = 14, height = 7)
g1 <- s.class(pca.1$scores, pop(aa.genlight), xax = 1, yax = 2, col = transp(col,.6),  ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE),  optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)

#formats the pdf with the use of adegraphics package
ADEgS(c(g1, g2), layout = c(1,2))
dev.off()
