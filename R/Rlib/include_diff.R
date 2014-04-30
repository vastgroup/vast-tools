#!/usr/bin/Rscript
#

#  This function takes a qual and returns c(post_alpha, post_beta)
#  Increments by prior alpha and prior distribution beta, uniform by default
parseQual <- function(qual, prior_alpha=1, prior_beta=1) {
  res <- as.numeric(unlist(strsplit(unlist(strsplit(qual, "@"))[2], ",")))
  res[1] <- res[1] + prior_alpha
  res[2] <- res[2] + prior_beta
  res
}

# calculate the probability that the first dist is > than second
#  P(psi1 > psi2) when alpha=0; more generally we are determining the probability
#  that Psi1 is greater than Psi2 by alpha.. eg. P((psi1 - psi2) > alpha) = 
## IMPORTANT: run this function as sample(firstDist, length(firstDist)) AND
##								           sample(secondDist, length(secondDist))
## UNLESS you have paired data, then don't sample.
pDiff  <- function(firstDist, secondDist, alpha=0.15) {
	N <- length(firstDist)	
	pass <- length( which( (firstDist - secondDist) > alpha  ) )
	pass / N
}

#  This function finds the maximum difference between psi1 and psi2 given
#  at least some acceptable probability for difference.  Defaults to 0.8
#  distributions where no diff value exists with a probability > 0.8 are given
#  maxDiff of 0.
maxDiff <- function(firstDist, secondDist, acceptProb=0.9) {
	alphaSet <- seq(0,1,0.01)  #make this global?
	probs <- unlist(lapply(alphaSet, function(x) { pDiff(firstDist, secondDist, x) }))
    ind <- max(c(which(probs > acceptProb), 1))
    alphaSet[ind]
}

shuffle <- function(x) {
	sample(x, length(x))	
}
