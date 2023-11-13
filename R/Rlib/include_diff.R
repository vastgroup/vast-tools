#!/usr/bin/env Rscript
#
# Author: Tim Sterne-Weiler, Ulrich Braunschweig 2014-2023
# u.braunschweig@utoronto.ca

## Check that columns of INCLUSION... table are what we think they are
checkHeader <- function(x, replicateA, replicateB) {
    reps <- c(replicateA, replicateB)
    sampInd <- unlist(sapply(reps, FUN=function(y) {which(x == y)}))
    if (length(sampInd) != length(reps)) {
        stop("[vast diff error]: Not all sampleNames found in input header!")
    }
    namesOK <- all(paste0(reps, "-Q") == x[sampInd + 1])
    if (!namesOK |  !all(x[1:3] == c("GENE", "EVENT", "COORD"))) {
        stop("[vast diff error]: Input table does not have expected format!")
    }
}

##  This function takes a qual and returns c(post_alpha, post_beta)
##  Increments by prior alpha and prior distribution beta, uniform by default
parseQual <- function(qual, prior_alpha=1, prior_beta=1) {
    res1 <- as.numeric(sub("[^@]*@([^,]*),.*", "\\1", qual))
    res2 <- as.numeric(sub("[^@]*@[^,]*,(.*)", "\\1", qual))

    irreg <- is.na(res1) | is.null(res1) | is.infinite(res1) |
        is.na(res2) | is.null(res2) | is.infinite(res2)
    res1[irreg] <- 0
    res2[irreg] <- 0
    res1 <- res1 + prior_alpha
    res2 <- res2 + prior_beta
    cbind(res1, res2)
}

## Calculate the probability that the first dist is > than second
##  P(psi1 > psi2) when alpha=0; more generally we are determining the probability
##  that Psi1 is greater than Psi2 by alpha.. eg. P((psi1 - psi2) > alpha) = 
## IMPORTANT: run this function as sample(firstDist, length(firstDist)) AND
##	      sample(secondDist, length(secondDist))
## UNLESS you have paired data, then don't sample.
pDiff  <- function(firstDist, secondDist, alpha=0.15) {
    N <- length(firstDist)	
    pass <- length( which( (firstDist - secondDist) > alpha  ) )
    pass / N
}

##  This function finds the maximum difference between psi1 and psi2 given
##  at least some acceptable probability for difference.  Defaults to 0.8
##  distributions where no diff value exists with a probability > 0.8 are given
##  maxDiff of 0.
maxDiff <- function(firstDist, secondDist, acceptProb=0.9) {
    alphaSet <- seq(0,1,0.01)  #make this global?
    probs <- unlist(lapply(alphaSet, function(x) { pDiff(firstDist, secondDist, x) }))
    ind <- max(c(which(probs > acceptProb), 1))
    alphaSet[ind]
}

##  Return the beta variance
betaVar <- function(alpha, beta) {
    var <- alpha*beta / (
	 ((alpha + beta) ** 2) * (alpha + beta + 1)
    )
    var
}

## Return a confidence interval 
betaCI <- function(betaDist, percentile = c(0.05, 0.95)) {
    quantile(betaDist, p=percentile, na.rm = T)
}

## Extension of betaCI function that includes the sampling step
betaCISample <- function(alpha, beta, n = 5000) {
    if (is.na(alpha) || is.na(beta)) {
      sample <- NA 
    } else {
      sample <- rbeta(n, alpha, beta)
    }
    return(betaCI(sample))
}


## MAKE VISUAL OUTPUT
plotDiff <- function(inpOne, inpTwo, expOne, expTwo, maxD, medOne, medTwo, sampOneName, sampTwoName, rever) {

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

diffBeta <- function(i, lines, opt,
                     shapeFirst, shapeSecond,
                     totalFirst, totalSecond,
                     shapeFirstAve, shapeSecondAve,
                     expFirst, expSecond,
                     repA.qualInd, repB.qualInd) {
## Main diff functionality; fit beta distributions to sample groups for one event and compare
## Is applied to each line of the current nLines of the INCLUSION... table

    ## if no data, next; # adapted from @lpantano's fork  7/22/2015
    if ((sum(totalFirst[i,] > (opt$minReads + opt$alpha + opt$beta)) < opt$minSamples) ||
        (sum(totalSecond[i,] > (opt$minReads + opt$alpha + opt$beta)) < opt$minSamples) ) {
              return(NULL)
    }

    ## Sample Posterior Distributions
    psiFirst <- lapply(1:(dim(shapeFirst)[2]), function(x) {
      ## sample here from rbeta(N, alpha, beta) if > -e
      if (totalFirst[i,x] < opt$minReads) { return(NULL) }
      rbeta(opt$size, shape1=shapeFirst[i,x,1], shape2=shapeFirst[i,x,2])
    })

    psiSecond <- lapply(1:(dim(shapeSecond)[2]), function(x) {
      ## sample here from rbeta(N, alpha, beta) if > -e
      if (totalSecond[i,x] < opt$minReads) { return(NULL) }
      rbeta(opt$size, shape1=shapeSecond[i,x,1], shape2=shapeSecond[i,x,2])
    })

    if (opt$paired) { 
        ## make sure both samples have a non-NULL replicate
        pairedNull <- sapply(psiFirst, is.null) & sapply(psiSecond, is.null)
        psiFirst[pairedNull] <- NULL
        psiSecond[pairedNull] <- NULL
    }

    ## Create non-parametric Joint Distributions
    psiFirstComb  <- do.call(c, psiFirst)
    psiSecondComb <- do.call(c, psiSecond)

    if (length(psiFirstComb) <= 0 || length(psiSecondComb) <= 0 ) { return(NULL) }

    ## if they aren't paired, then shuffle the joint distributions...
    if (!opt$paired ) {
      paramFirst <- try (suppressWarnings(
                              fitdistr(psiFirstComb,
                                      "beta",
                                      list( shape1=shapeFirstAve[i,1], shape2=shapeFirstAve[i,2])
                              )$estimate ), TRUE )
      paramSecond <- try (suppressWarnings(
                              fitdistr(psiSecondComb,
                                      "beta",
                                      list( shape1=shapeSecondAve[i,1], shape2=shapeSecondAve[i,2])
                              )$estimate ), TRUE )
      ## if optimization fails its because the distribution is too narrow
      ## in which case our starting shapes should already be good enough
      if (class(paramFirst) != "try-error") {
        psiFirstComb <- rbeta(opt$size, shape1=paramFirst[1], shape2=paramFirst[2])
      }
      if (class(paramSecond) != "try-error") {
        psiSecondComb <- rbeta(opt$size, shape1=paramSecond[1], shape2=paramSecond[2])
      }
    }

    ## get empirical posterior median of psi
    medOne <- median(psiFirstComb)
    medTwo <- median(psiSecondComb)

    ## look for a max difference given prob cutoff...
    if (medOne > medTwo) {
      max <- maxDiff(psiFirstComb, psiSecondComb, opt$prob)
    } else {
      max <- maxDiff(psiSecondComb, psiFirstComb, opt$prob)
    }

    filtOut <- sprintf("%s\t%s\t%f\t%f\t%f\t%s", lines[[i]][1], lines[[i]][2], medOne, medTwo, medOne - medTwo, round(max,2))

    ## check for significant difference
    if (max < opt$minDiff) {
      ## non-sig, return null plots and text output
      return(list(NULL, NULL, NULL, NULL, filtOut))
    } else { 
      ## SIGNIFICANT from here on out
      sigInd <- i

      if (opt$noPDF) {
        eventTitle <- NULL
        eventCoord <- NULL
        retPlot    <- NULL
      } else {
        eventTitle <- paste(c("Gene: ", lines[[i]][1], "  Event: ", lines[[i]][2]), collapse="")
        eventCoord <- paste(c("Coordinates: ", lines[[i]][3]), collapse="")

        ## Print visual output to pdf;
        if (medOne > medTwo) {
            retPlot <- plotDiff(psiFirstComb, psiSecondComb,
                                expFirst[i,], expSecond[i,], max, medOne, medTwo, sampOneName, sampTwoName, FALSE)
        } else {
            retPlot <- plotDiff(psiSecondComb, psiFirstComb,
                                expFirst[i,], expSecond[i,], max, medTwo, medOne, sampTwoName, sampOneName, TRUE)
        }
      }
      ## sig event return
      return(list(retPlot, eventTitle, eventCoord, sigInd, filtOut))  #return of mclapply function
    }
}
