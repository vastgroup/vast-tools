#!/usr/bin/env Rscript
#
# Author: Tim Sterne-Weiler, Ulrich Braunschweig 2014-2024
# u.braunschweig@utoronto.ca

## Check that columns of INCLUSION... table are what we think they are
checkHeader <- function(x, replicateA, replicateB) {
    reps <- c(replicateA, replicateB)
    sampInd <- unlist(sapply(reps, FUN=function(y) {which(x == y)}))
    if (length(sampInd) != length(reps)) {
        stop("[vast diff error]: Not all replicate names found in input header!\n")
    }
    namesOK <- all(paste0(reps, "-Q") == x[sampInd + 1])
    if (!namesOK |  !all(x[1:3] == c("GENE", "EVENT", "COORD"))) {
        stop("[vast diff error]: Input table does not have expected format!\n")
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
## P(psi1 > psi2) when alpha=0; more generally we are determining the probability
## that psi1 is greater than psi2 by alpha.. eg. P((psi1 - psi2) > alpha) = 
## IMPORTANT: run this function as sample(firstDist, length(firstDist)) AND
##	      sample(secondDist, length(secondDist))
## UNLESS you have paired data, then don't sample.
pDiff  <- function(firstDist, secondDist, alpha=0.15) {
    N <- length(firstDist)
    pass <- length(which((firstDist - secondDist) > alpha))
    pass / N
}

##  This function finds the maximum difference between psi1 and psi2 given
##  at least some acceptable probability for difference.  Defaults to 0.8
##  distributions where no diff value exists with a probability > 0.8 are given
##  maxDiff of 0.
maxDiff <- function(firstDist, secondDist, acceptProb=0.9, alphaSet) {
    probs <- unlist(lapply(alphaSet, pDiff, 
        firstDist=firstDist, secondDist=secondDist))
    ind <- max(c(which(probs > acceptProb), 1))
    alphaSet[ind]
}

##  Return the beta variance
betaVar <- function(alpha, beta) {
    var <- alpha * beta / (
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

plotDiff <- function(inpOne, inpTwo, expOne, expTwo, maxD, medOne, medTwo, 
    sampOneName, sampTwoName, alphaSet, rever) {
### Make visual output

  if(rever) {   #write this better. ;-)
    curCol <- cbb[3:2]
  } else {
    curCol <- cbb[2:3]
  }

  if (length(expOne) == 0 || length(expTwo) == 0) {return(NULL)}
  one <- data.frame(x=expOne, y=-0.5)
  two <- data.frame(x=expTwo, y=-0.5)

  distPlot <- ggplot(melt(as.data.frame(do.call(cbind, list(inpOne, inpTwo))),
                  measure.vars=c("V1","V2")), 
      aes(fill=variable, x=value)) +
      geom_histogram(aes(y=after_stat(density)), binwidth=0.03333, alpha=0.5, col="grey", position="identity") +
      theme_bw() + scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) + xlab(expression(hat(Psi))) +
      scale_fill_manual(values=curCol, labels=c(sampOneName, sampTwoName), name="Samples") +
      geom_point(data=one, mapping=aes(x=x, y=y), col=cbb[2], fill=cbb[2], alpha=0.50, inherit.aes = FALSE) +
      geom_point(data=two, mapping=aes(x=x, y=y), col=cbb[3], fill=cbb[3], alpha=0.50, inherit.aes = FALSE)

  probPlot <- ggplot(as.data.frame(cbind(seq(0,1,0.01),
            unlist(lapply(alphaSet, function(x) {
               pDiff(inpOne, inpTwo, x)
            })))), aes(x=V1, y=V2)) +
            geom_line() + theme_bw() +
            geom_vline(xintercept=maxD, lty="dashed", col=cbb[7]) +
            ylab(expression(P((hat(Psi)[1] - hat(Psi)[2]) > x))) +
            xlab(expression(x)) + scale_y_continuous(limits = c(0, 1), oob = scales::oob_keep) +
            annotate("text", x=(maxD + 0.08), y=0.05, label=maxD, col=cbb[7])

  return(list(distPlot, probPlot))#
}

plotPrint <- function(plotTitle, plotCoord, plotList) {
### Print saved plotting information to output file
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
                     repA.qualInd, repB.qualInd,
                     okFirst, okSecond, skip, 
                     alphaSet) {
## Main diff functionality; fit beta distributions to sample groups for one event and compare
## Is applied to each line of the current nLines of the INCLUSION... table

    ## Sample Posterior Distributions
    psiFirst <- lapply(1:(dim(shapeFirst)[2]), function(x) {
        ## sample here from rbeta(N, alpha, beta) if > -e
        if (!okFirst[i,x]) {return(NULL)}
        rbeta(opt$size, shape1=shapeFirst[i,x,1], shape2=shapeFirst[i,x,2])
    })

    psiSecond <- lapply(1:(dim(shapeSecond)[2]), function(x) {
        ## sample here from rbeta(N, alpha, beta) if > -e
        if (!okSecond[i,x]) {return(NULL)}
        rbeta(opt$size, shape1=shapeSecond[i,x,1], shape2=shapeSecond[i,x,2])
    })

    ## Create non-parametric Joint Distributions
    psiFirstComb  <- do.call(c, psiFirst)
    psiSecondComb <- do.call(c, psiSecond)

    ## Try to to fit a beta distribution on the replicates' parameter estimates.
    if (!opt$paired) {
        if (length(psiFirstComb) > 0) {
            paramFirst <- try(suppressWarnings(
                                    fitdistr(psiFirstComb,
                                            "beta",
                                            list(shape1=shapeFirstAve[i,1], shape2=shapeFirstAve[i,2])
                                   )$estimate ), TRUE )
        }
        if (length(psiSecondComb)) {
            paramSecond <- try(suppressWarnings(
                                    fitdistr(psiSecondComb,
                                            "beta",
                                            list(shape1=shapeSecondAve[i,1], shape2=shapeSecondAve[i,2])
                                    )$estimate ), TRUE )
        }
        ## If optimization fails it's because the distribution is too narrow
        ## in which case our starting shapes should already be good enough.
        ## Shuffle since not paired, and downsample to size of smaller number
        ## but don't downsize if the other sample type is absent anyway.
        minSample <- min(c(length(psiFirstComb), length(psiSecondComb)))
      
        if (length(psiFirstComb) > 0) {
            if (!inherits(paramFirst, "try-error")) {
                psiFirstComb <- rbeta(opt$size, shape1=paramFirst[1], shape2=paramFirst[2])
            } else {
                if (length(psiSecondComb) > 0) {
                    psiFirstComb <- sample(psiFirstComb, minSample)
                }
            }
        }
        if (length(psiSecondComb) > 0) {
            if (!inherits(paramSecond, "try-error")) {
                psiSecondComb <- rbeta(opt$size, shape1=paramSecond[1], shape2=paramSecond[2])
            } else {
                if (length(psiFirstComb) > 0) {
                    psiSecondComb <- sample(psiSecondComb, minSample)
                }
            }
        }
    }

    ## get empirical posterior median of psi
    medOne <- median(psiFirstComb)
    medTwo <- median(psiSecondComb)

    ## look for a max difference given prob cutoff...
    if (skip[i]) {
        if (is.null(medOne)) {medOne <- NA}
        if (is.null(medTwo)) {medTwo <- NA}
        max <- NA
    } else {
        if (medOne > medTwo) {
            max <- maxDiff(psiFirstComb, psiSecondComb, opt$prob, alphaSet)
        } else {
            max <- maxDiff(psiSecondComb, psiFirstComb, opt$prob, alphaSet)
        }
    }

    filtOut <- sprintf("%s\t%s\t%f\t%f\t%f\t%s", 
                   lines[[i]][1], lines[[i]][2], medOne, medTwo, medOne-medTwo, round(max,2))

    ## check for significant difference
    if (opt$noPDF || (is.na(max) || max < opt$minDiff)) {
        ## non-sig, return null plots and text output
        return(list(NULL, NULL, NULL, NULL, filtOut))
    } else { 
        ## significant; plot
        sigInd <- i
        eventTitle <- paste(c("Gene: ", lines[[i]][1], "  Event: ", lines[[i]][2]), collapse="")
        eventCoord <- paste(c("Coordinates: ", lines[[i]][3]), collapse="")

        ## Store visual output to be printed to pdf
        if (medOne > medTwo) {
            retPlot <- plotDiff(psiFirstComb, psiSecondComb,
                                expFirst[i,][okFirst[i,]], expSecond[i,][okSecond[i,]], 
                                max, medOne, medTwo, sampOneName, sampTwoName, alphaSet, FALSE)
        } else {
            retPlot <- plotDiff(psiSecondComb, psiFirstComb,
                                expFirst[i,][okFirst[i,]], expSecond[i,][okSecond[i,]],
                                max, medTwo, medOne, sampTwoName, sampOneName, alphaSet, TRUE)
        }
      
        ## sig event return
        return(list(retPlot, eventTitle, eventCoord, sigInd, filtOut))
    }
}

writeLog <- function(vastPath, opt) {
### Add a line with call information to the general vast-tools log file

    logName <- file.path(opt$output, "VTS_LOG_commands.txt")
    vastVer <- scan(file.path(vastPath, "VERSION"), what="character", quiet = T)

    ## Get date, time and vast-tools version
    systime <- Sys.time()
    sdate <- paste(
	      sub("([0-9]{4})-.+", "\\1", systime),
	      sub("[0-9]{4}-[0-9]{2}-([0-9]{2}).+", "\\1", systime),
	      sub("[0-9]{4}-([0-9]{2})-.+", "\\1", systime),
	      sep="-")
    stime <- sub("[^ ]+ (.{5}):[0-9]{2}.*", "\\1", systime)

    msg1 <- paste0(
      "[VAST-TOOLS v", vastVer, ", ",
      sdate, " (", stime, ")]"
    )

    msg2 <- "vast-tools diff"

    ## Format input options
    msg3 <- paste0(
        "-a ", opt$replicateA, " ",
        "-b ", opt$replicateB, " ",
        "--sampleNameA ", opt$sampleNameA, " ",
        "--sampleNameB ", opt$sampleNameB, " ",
        "-i ", opt$input, " ",
        "-n ", opt$nLines, " ",
        "-p ", opt$paired, " ",
        "-d ", opt$baseName, " ",
        "-o ", opt$output, " ",
        "-r ", opt$prob, " ",
        "-m ", opt$minDiff, " ",
        "-e ", opt$minReads, " ",
        "-S ", opt$minSamples, " ",
        "--alpha ", opt$alpha, " ",
        "--beta ", opt$beta, " ",
        "-s ", opt$size, " ",
        "-z ", opt$seed
    )

    msg <- paste(msg1, msg2, msg3, paste = " ")

    try_error <- try(logFile <- file(logName, "a"), silent = TRUE)
    if (inherits(try_error, "try-error")) {
        cat("[vast diff warning]: Could not open log file\n")
    } else {
        write(msg, logFile)
        close(logFile)
    }
}