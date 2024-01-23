#!/usr/bin/env Rscript

# Author: Tim Sterne-Weiler, Ulrich Braunschweig 2014-2024
# u.braunschweig@utoronto.ca

# Copyright (C) 2014 Tim Sterne-Weiler, Ulrich Braunschweig
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


argv <- commandArgs(trailingOnly = F)
scriptPath <- dirname(sub("--file=", "", argv[grep("--file", argv)]))

# Source Rlib.
source(file.path(scriptPath, "Rlib", "include.R"))
source(file.path(scriptPath, "Rlib", "include_diff.R"))

# custom install from include.R
loadPackages(c("optparse"), local.lib=file.path(scriptPath, "Rlib"))

argv <- commandArgs(TRUE)

# optparse..
option.list <- list(
    make_option(c("-a", "--replicateA"), type = "character", default = NULL, metavar = "SampleA,SampleB,SampleC",
        help = "Required, 1:n sample names separated by , [mandatory!]"),
    make_option(c("-b", "--replicateB"), type = "character", default = NULL, metavar = "SampleA,SampleB,SampleC",
        help = "Required, 1:n sample names separated by , [mandatory!]\n
	
[input options]"),
    make_option(c("--sampleNameA"), type = "character", default = NULL, metavar = "string",
        help = "Name of the replicate set A, [default is first element of --replicateA]"),
    make_option(c("--sampleNameB"), type = "character", default = NULL, metavar = "string",
        help = "Name of the replicate set B, [default is first element of --replicateB]"),
    make_option(c("-i", "--input"), type = "character", default = "INCLUSION_LEVELS",
        help = "Exact or Partial match to PSI table in output directory [default %default]"),
    make_option(c("-n", "--nLines"), type = "integer", default = "5000",
        help = "Number of lines to read/process in parallel at a time... 
                lower number = less memory = greater overhead [default %default]"),
    make_option(c("-p", "--paired"), action="store_true", default = FALSE,
        help = "Samples are paired, -a pairOneA,pairTwoA,.. -b pairOneB,pairTwoB,.. [default %default]\n

[output options]"),
    make_option(c("-d", "--baseName"), type = "character", default = "input.DIFF", metavar="FILE",
        help = "Base name for output files [default %default]"),
    make_option("--noPDF", action="store_true", default = FALSE,
        help = "Do not save PDF output [default %default]"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
        help = "Vast-tools output directory, [default vast_out]\n

[statistical options]"),
    make_option(c("-r", "--prob"), type = "numeric", default = 0.95,
        help = "Probability threshold for P( (psi1 - psi2) > x ) > threshold [default %default]"),
    make_option(c("-m", "--minDiff"), type = "numeric", default = 0.1,
        help = "Threshold for min diff where P( (psi1 - psi2) > threshold ) > --prob [default %default]"),
    make_option(c("-e", "--minReads"), type = "numeric", default = 10,
        help = "Threshold for min reads in a sample (use this flag unless you believe the prior) [default %default]"),
    make_option(c("-S","--minSamples"), type = "integer", default = 1,
        help = "Threshold for min samples with min reads in a group  [default %default]"),
    make_option(c("--alpha"), type = "numeric", default = 1,
        help = "First shape parameter for the Beta prior distribution P(psi), Uniform by default [default %default]"),
    make_option(c("--beta"), type = "numeric", default = 1,
        help = "Second shape parameter for the Beta prior distribution P(psi), Uniform by default [default %default]"),
    make_option(c("-s", "--size"), type = "integer", default = 1000,
        help = "Size of the posterior emperical distribution over psi, lower = faster 
                but less precise [default %default]\n

[general options]"),
    make_option(c("-c", "--cores"), type = "integer", default = 1, metavar="int",
        help="Number of cores to use for plot processing.. [default %default]"),
    make_option(c("-z", "--seed"), type = "integer", default = 10, metavar="int",
        help="Seed the RNG for a deterministic result.. [default %default]"),
    make_option(c("-v", "--verbose"), type = "logical", default = TRUE, metavar=NULL,
        help="Enable verbose [default %default]")
)

parser <- OptionParser(option_list = option.list,
    usage       = "vast-tools diff -a SampleA,..,SampleD -b SampleF,..,SampleG [options]",
	description = "Score differential alternative splicing between two groups of samples (set A - set B)"
)

optpar <- parse_args(parser, argv, positional_arguments = TRUE)
opt <- optpar$options

if (length(commandArgs(TRUE)) == 2 & commandArgs(TRUE)[1] == "-o") {
    print_help(parser)
    stop("No input")
}

## input checks
if (opt$prob < 0 | opt$prob >= 1) {
    stop ("--prob must be a positive number < 1")
}

if (opt$minDiff < 0 | opt$minDiff >= 1) {
    stop ("--minDiff must be a positive number < 1")
}	


loadPackages(c("MASS", "reshape2", "ggplot2", "grid", "parallel"), 
    local.lib=file.path(scriptPath, "Rlib"))

## move to output directory
setwd(opt$output)

## seed RNG
set.seed(opt$seed)

## try and find the input file if they aren't exact
if (!file.exists(opt$input)) {
    potentialFiles <- Sys.glob(paste(c("*", opt$input, "*"), collapse=""))
    if (length(potentialFiles) >= 1) { 
        # now sort and take the one with the 'biggest' number of samples
        potentialFiles_sort <- rev(sort(potentialFiles))
        opt$input <- potentialFiles_sort[1]
    } else {
        # Still can't find input after searching...
	    print_help(parser)
        stop("[vast diff error]: No input file given!")
    }
}

## Setting input files
inputFile <- file(opt$input, 'r')

firstRepSet  <- unlist(strsplit(as.character(opt$replicateA), ","))
secondRepSet <- unlist(strsplit(as.character(opt$replicateB), ","))
firstRepN  <- length(firstRepSet)
secondRepN <- length(secondRepSet)

if (firstRepN <= 0 || secondRepN <= 0) { 
    print_help(parser) 
    stop("[vast diff error]: No replicate sample names given! -a sampA,sampB -b sampC,sampD")
}

## Set number of replicates
if (opt$paired && firstRepN != secondRepN) {
    stop("Paired samples but unequal numbers")
}

## Make sure there are sample names
if (is.null(opt$sampleNameA)) {
    opt$sampleNameA <- firstRepSet[1]
}
if (is.null(opt$sampleNameB)) {
    opt$sampleNameB <- secondRepSet[1]
}
## Set output sample names for plot
sampOneName <- substr(opt$sampleNameA, 1, 9)
sampTwoName <- substr(opt$sampleNameB, 1, 9)


## Get header
head_n <- unlist(strsplit(readLines(inputFile, n=1), "\t"))

## Check if header is correct
checkHeader(head_n, firstRepSet, secondRepSet)
if (opt$verbose) {cat("[vast diff]: Input table format OK.\n")}

## Indexes of samples of interest
repAind <- which(head_n %in% firstRepSet)
repBind <- which(head_n %in% secondRepSet)

## Indexes of Quals
repA.qualInd <- repAind + 1
repB.qualInd <- repBind + 1


## CONST
## dPSI levels at which to test for significant difference
alphaSet <- seq(0, 1, 0.01)


## TMP OUT
if (opt$baseName == "input.DIFF") {
    pdfname <- sub("\\.[^.]*(\\.gz)?$", ".DIFF_plots.pdf", basename(opt$input))
    outname <- sub("\\.[^.]*(\\.gz)?$", ".DIFF.txt", basename(opt$input))
} else {
    pdfname <- paste0(opt$baseName, ".pdf")
    outname <- paste0(opt$baseName, ".tab")
}

outhandle <- file(outname, "w")

if (!opt$noPDF) {
    pdf(pdfname, width=7, height=3.5, family="sans", compress=FALSE)
}
writeLines(sprintf("GENE\tEVENT\t%s\t%s\tE[dPsi]\tMV[dPsi]_at_%s", opt$sampleNameA, opt$sampleNameB, opt$prob), 
  outhandle)


### BEGIN READ INPUT ###
## Iterate through input, 'nLines' at a time to reduce overhead/memory
while (length(lines <- readLines(inputFile, n=opt$nLines)) > 0) {
    lines <- strsplit(lines, split="\t")
    ## use parallel computing to store plots in plotListed
    ## then print them to the pdf afterwards before next chunk of nLines from file.
    
    ## Get inclusion/exclusion read numbers (i.e., beta shape parameters) from quality columns
    shapeFirst  <- array(dim=c(length(lines), length(repAind), 2))
    shapeSecond <- array(dim=c(length(lines), length(repBind), 2))
    for (i in 1:length(repA.qualInd)) {
        shapeFirst[,i,] <- parseQual(sapply(lines, "[[", repA.qualInd[i]),
                                     prior_alpha=opt$alpha, prior_beta=opt$beta)
    }
    for (i in 1:length(repB.qualInd)) {
        shapeSecond[,i,] <- parseQual(sapply(lines, "[[", repB.qualInd[i]),
                                      prior_alpha=opt$alpha, prior_beta=opt$beta)
    }

    totalFirst  <- apply(shapeFirst,  MAR=2, rowSums)
    totalSecond <- apply(shapeSecond, MAR=2, rowSums)

    shapeFirstAve  <- apply(shapeFirst,  MAR=3, rowMeans)
    shapeSecondAve <- apply(shapeSecond, MAR=3, rowMeans)


    ## Expected PSI for each replicate
    expFirst  <- sapply(1:(dim(shapeFirst)[2]), FUN=function(x)  {shapeFirst[,x,1] / totalFirst[,x]})
    expSecond <- sapply(1:(dim(shapeSecond)[2]), FUN=function(x) {shapeSecond[,x,1] / totalSecond[,x]})

    ## Figure out which lines to skip
    okFirst  <- totalFirst  > opt$minReads + opt$alpha + opt$beta - 1
    okSecond <- totalSecond > opt$minReads + opt$alpha + opt$beta - 1

    if (opt$paired) {
        skip <- rowSums(okFirst & okSecond) < opt$minSamples
        ## Flag matched samples in both if absent in one type
        pairedBad <- !okFirst | !okSecond
        okFirst[pairedBad]  <- FALSE
        okSecond[pairedBad] <- FALSE
    } else {
        skipFirst  <- rowSums(okFirst)  < opt$minSamples
        skipSecond <- rowSums(okSecond) < opt$minSamples
        skip <- skipFirst | skipSecond
    }

    ## Simulate distributions and score overlap, create plot data
    plotListed <- mclapply(1:length(lines), diffBeta, lines=lines, opt=opt,
                        shapeFirst, shapeSecond,
                        totalFirst, totalSecond,
                        shapeFirstAve, shapeSecondAve,
                        expFirst, expSecond,
                        repA.qualInd, repB.qualInd,
                        okFirst, okSecond, skip,
                        alphaSet,
                        mc.cores=opt$cores, mc.preschedule=TRUE, mc.cleanup=TRUE)
    

    for (it in 1:length(lines)) {
        write(plotListed[[it]][[5]], outhandle)
        
        ## PRINT LIST OF PLOTS
        if (opt$noPDF || is.null(plotListed[[it]][[1]])) {next}
        plotPrint(plotListed[[it]][[2]], plotListed[[it]][[3]], plotListed[[it]][[1]])
    }

} # End While

if (!opt$noPDF) {garbage <- dev.off()}

writeLog(file.path(scriptPath, ".."), opt)

flush(outhandle)
close(outhandle)
if (opt$verbose) {cat("[vast diff]: DONE\n\n")}
q(status=0)
