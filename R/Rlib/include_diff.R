#!/usr/bin/env Rscript
#
# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca

#  This function takes a qual and returns c(post_alpha, post_beta)
#  Increments by prior alpha and prior distribution beta, uniform by default
parseQual <- function(qual, prior_alpha=1, prior_beta=1) {
    if(is.na(qual) || !grepl("@", qual)) { return(c(prior_alpha, prior_beta)) }  ## for INT NA Columns
    res <- unlist(strsplit(unlist(strsplit(as.character(qual), "@"))[2], ","))
    if(is.na(res[1]) || is.na(res[2])) { return(c(prior_alpha, prior_beta)) }
    if(is.null(res[1]) || is.null(res[2])) { return(c(prior_alpha, prior_beta)) }
    if(res[1] == "NA" || res[2] == "NA") { return(c(prior_alpha, prior_beta)) }
    res <- as.numeric(res)
    if(is.nan(res[1]) || is.nan(res[2])) { return(c(prior_alpha, prior_beta)) }
    if(is.infinite(res[1]) || is.infinite(res[2])) { return(c(prior_alpha, prior_beta)) }
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
plotDiff <- function(inpOne, inpTwo, expOne, expTwo, maxD, medOne, medTwo, sampOneName, sampTwoName, rever ) {

  if(rever) {   #write this better. ;-)
    curCol <- cbb[3:2]
  } else {
    curCol <- cbb[2:3]
  }

  if(length(expOne) == 0 || length(expTwo) == 0) { return(NULL) }
  one <- data.frame(x=expOne, y=-0.5)
  two <- data.frame(x=expTwo, y=-0.5)

  distPlot <- ggplot(melt(as.data.frame(
         do.call(cbind,list(inpOne, inpTwo))
         ), measure.vars=c("V1","V2")), aes(fill=variable, x=value))+
         geom_histogram(aes(y=..density..), binwidth=0.03333,alpha=0.5, col="grey", position="identity")+
         theme_bw()+xlim(c(0,1))+xlab(expression(hat(Psi)))+
         scale_fill_manual(values=curCol, labels=c(sampOneName, sampTwoName), name="Samples")+
         geom_point(data=one, mapping=aes(x=x, y=y), col=cbb[2], fill=cbb[2], alpha=0.85, inherit.aes = FALSE)+
         geom_point(data=two, mapping=aes(x=x, y=y), col=cbb[3], fill=cbb[3], alpha=0.85, inherit.aes = FALSE)

  probPlot <- ggplot(as.data.frame(cbind(seq(0,1,0.01),
            unlist(lapply(alphaList, function(x) {
               pDiff(inpOne, inpTwo, x)
            })))), aes(x=V1, y=V2))+
            geom_line()+theme_bw()+
            geom_vline(xintercept=maxD, lty="dashed", col=cbb[7])+
            ylab(expression(P((hat(Psi)[1]-hat(Psi)[2]) > x)))+
            xlab(expression(x))+ylim(c(0,1))+
            annotate("text", x=(maxD+0.08), y=0.05, label=maxD, col=cbb[7])

  return(list(distPlot, probPlot))
}

plotPrint <- function(plotTitle, plotCoord, plotList) {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 2, widths = unit(c(5, 4), "null"), heights = unit(c(1, 0.5, 5), "null"))))
    grid.text(as.character(plotTitle), check.overlap=TRUE, gp=gpar(font=2), draw=T, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    grid.text(as.character(plotCoord), check.overlap=TRUE, gp=gpar(font=1), draw=T, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    print(plotList[[1]], vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(plotList[[2]], vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
    popViewport(1)
}

# Shuffle...
shuffle <- function(x) {
    sample(x, length(x))	
}
