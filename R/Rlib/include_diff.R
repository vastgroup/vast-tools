#!/usr/bin/Rscript
#
# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca

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

#
#  return the beta variance
betaVar <- function(alpha, beta) {
  var <- alpha*beta / (
		 ((alpha + beta) ** 2) * (alpha + beta + 1)
  )
  var
}

# 
betaCI <- function(betaDist, percentile = c(0.05, 0.95)) {
  quantile(betaDist, p=percentile, na.rm = T)
}

# Extention of betaCI function that includes the sampling step
betaCISample <- function(alpha, beta, n = 5000) {
  if (is.na(alpha) || is.na(beta)) {
    sample <- NA 
  } else {
    sample <- rbeta(n, alpha, beta)
  }
  return(betaCI(sample))
}


### MAKE VISUAL OUTPUT
plotDiff <- function(eventName, inpOne, inpTwo, maxD, medOne, medTwo, sampOneName, sampTwoName, rever ) {
#  dput(inpOne)
#  dput(inpTwo)
#  dput(maxD)
#  dput(medOne)
#  dput(medTwo)
#  dput(sampOneName)
#  dput(sampTwoName)

  if(rever) {   #write this better. ;-)
    curCol <- cbb[3:2]
  } else {
    curCol <- cbb[2:3]
  }

  distPlot <- ggplot(melt(as.data.frame(
         do.call(cbind,list(inpOne, inpTwo))
         )), aes(fill=variable, x=value))+
#         geom_vline(x=medOne, col=cbb[2])+
#         geom_vline(x=medTwo, col=cbb[3])+
         geom_histogram(aes(y=..density..),alpha=0.5, col="grey", position="identity")+
         theme_bw()+xlim(c(0,1))+xlab(expression(hat(Psi)))+
         scale_fill_manual(values=curCol, labels=c(sampOneName, sampTwoName), name="Samples")

  probPlot <- ggplot(as.data.frame(cbind(seq(0,1,0.01),
            unlist(lapply(alphaList, function(x) {
               pDiff(inpOne, inpTwo, x)
            })))), aes(x=V1, y=V2))+
            geom_line()+theme_bw()+
            geom_vline(x=maxD, lty="dashed", col=cbb[7])+
            ylab(expression(P((hat(Psi)[1]-hat(Psi)[2]) > x)))+
            xlab(expression(x))+ylim(c(0,1))+
            annotate("text",x=(maxD+0.08), y=0.05, label=maxD, col=cbb[7])

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2, widths = unit(c(5, 4), "null"), heights = unit(c(1, 5), "null"))))
  grid.text(eventName,gp=gpar(font=2), draw=T, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  print(distPlot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(probPlot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

}

# Shuffle...
shuffle <- function(x) {
	sample(x, length(x))	
}
