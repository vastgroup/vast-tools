#!/usr/bin/Rscript

# calculate the probability that the first dist is > than second
#  P(psi1 > psi2) when alpha=0; more generally we determining the probability
#  that Psi1 is greater than Psi2 by alpha.. eg. P((psi1 - psi2) > alpha) = 
## IMPORTANT: run this function as sample(firstDist, length(firstDist)) AND
##								   sample(secondDist, length(secondDist))
pDiff  <- function(firstDist, secondDist, alpha=0.15) {
	N <- length(firstDist)	
	err <- length( which( (firstDist - secondDist) > alpha  ) )
	err / N
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
