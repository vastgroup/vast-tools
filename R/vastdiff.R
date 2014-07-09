#!/usr/bin/Rscript --vanilla

# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca

# Copyright (C) 2014 Tim Sterne-Weiler
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
scriptPath <- dirname(sub("--file=","",argv[grep("--file",argv)]))

# Source Rlib.
source(paste(c(scriptPath,"/Rlib/include.R"), collapse=""))
source(paste(c(scriptPath,"/Rlib/include_diff.R"), collapse=""))

# custom install from include.R
loadPackages(c("optparse", "RColorBrewer", "reshape2", "ggplot2", "grid", "parallel"), local.lib=paste(c(scriptPath,"/Rlib"), collapse=""))

argv <- commandArgs(TRUE)

# optparse..
option.list <- list(
    make_option(c("-a", "--replicateA"), type = "character", default = NULL, metavar = "SampleA@SampleB@SampleC",
        help = "Required, 1:n sample names separated by @ [mandatory!]"),
    make_option(c("-b", "--replicateB"), type = "character", default = NULL, metavar = "SampleA@SampleB@SampleC",
        help = "Required, 1:n sample names separated by @ [mandatory!]\n

[input options]"),
    make_option(c("--sampleNameA"), type = "character", default = NULL, metavar = "string",
        help = "Name of the replicate set A, [default is first element of --replicateA]"),
    make_option(c("--sampleNameB"), type = "character", default = NULL, metavar = "string",
        help = "Name of the replicate set B, [default is first element of --replicateB]"),
    make_option(c("-i", "--input"), type = "character", default = "INCLUSION_LEVELS",
        help = "Exact or Partial match to PSI table in output directory [default %default]"),
    make_option(c("-n", "--nLines"), type = "integer", default = "100",
        help = "Number of lines to read/process in parallel at a time... lower number = less memory = greater overhead [default %default]"),
    make_option(c("-p", "--paired"), type = "logical", default = FALSE,
        help = "Samples are paired, -a pairOneA@pairTwoA@.. -b pairOneB@pairTwoB [default %default]\n

[output options]"),
    make_option(c("-f", "--filter"), type = "logical", default = TRUE,
        help = "Filter output for differential events only [default %default]"),
    make_option(c("-d", "--pdf"), type = "character", default = "plotDiff_out", metavar="FILE",
        help = "Plot visual output (pdf) for differential events into FILE [default %default]"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
        help = "Output directory, [default vast_out]\n

[statistical options]"),
    make_option(c("-r", "--prob"), type = "numeric", default = 0.95,
        help = "Probability threshold for P( (psi1 - psi2) > x ) > threshold [default %default]"),
    make_option(c("-m", "--minDiff"), type = "numeric", default = 0.1,
        help = "Threshold for min diff where P( (psi1 - psi2) > threshold ) > --prob [default %default]"),
    make_option(c("-e", "--minReads"), type = "numeric", default = 5,
        help = "Threshold for min reads in a sample (use this flag unless you believe the prior) [default %default]"),
    make_option(c("--alpha"), type = "numeric", default = 1,
        help = "First shape parameter for the Beta prior distribution P(psi), Uniform by default [default %default]"),
    make_option(c("--beta"), type = "numeric", default = 1,
        help = "Second shape parameter for the Beta prior distribution P(psi), Uniform by default [default %default]"),
    make_option(c("-s", "--size"), type = "integer", default = 5000,
        help = "Size of the posterior emperical distribution over psi, lower = faster... [default %default]\n

[general options]"),
    make_option(c("-c", "--cores"), type = "integer", default = 1, metavar="int",
        help="Number of cores to use for plot processing.. [default %default]"),
    make_option(c("-v", "--verbose"), type = "logical", default = TRUE, metavar=NULL,
        help="Enable verbose [default %default]")
)

parser <- OptionParser(option_list = option.list,
            usage = "vast-tools diff -a SampleA@..@SampleD -b SampleF@..@SampleG [options]\n\nAuthor: Tim Sterne-Weiler\nQuestions? OR Bug Reports: tim.sterne.weiler@utoronto.ca")
optpar <- parse_args(parser, argv, positional_arguments = TRUE)
opt <- optpar$options

## move to output directory
setwd(opt$output)

## try and find the input file if they aren't exact
if(!file.exists(opt$input)) {
  potentialFiles <- Sys.glob( paste(c("*",opt$input,"*"), collapse="") )
  if( length( potentialFiles ) >= 1) { 
	 # now sort and take the one with the 'biggest' number of samples
    potentialFiles_sort <- rev( sort( potentialFiles ) )
    opt$input <- potentialFiles_sort[1]
  } else {
    # Still can't find input after searching...
	 print_help(parser)
    stop("[vast diff error]: No input file given!")
  }
}

## Setting input files.
inputFile <- file( opt$input, 'r' )

## Make plot directory
#if(opt$pdf) {
  #dir.create(paste(c(opt$output, "/diff_out"), collapse=""))
#}

#-replicatesA=name1@name2@name3 -replicatesB=name4@name5
firstRepSet <- unlist(strsplit( as.character(opt$replicateA) , "@" ))
secondRepSet <- unlist(strsplit( as.character(opt$replicateB), "@" ))

if( length(firstRepSet) <= 0 || 
	 length(secondRepSet) <= 0) { 
  print_help(parser) 
  stop("[vast diff error]: No replicate sample names given!!! -a sampA@sampB -b sampC@sampD")
}

# Set number of replicates
firstRepN <- length(firstRepSet)
secondRepN <- length(secondRepSet)

# Make sure there are sample names
if(is.null( opt$sampleNameA ) ) {
  opt$sampleNameA <- firstRepSet[1]
}
if(is.null( opt$sampleNameB ) ) {
  opt$sampleNameB <- secondRepSet[1]
}
# Set output sample names for plot
sampOneName <- paste(c(substr(opt$sampleNameA, 1, 4), "(n=", as.character(firstRepN), ")"), collapse="")
sampTwoName <- paste(c(substr(opt$sampleNameB, 1, 4), "(n=", as.character(secondRepN), ")"), collapse="")


## INITIALIZE LISTS ##
shapeFirst <- vector("list", firstRepN)
shapeSecond <- vector("list", secondRepN)

psiFirst <- vector("list", firstRepN)
psiSecond <- vector("list", secondRepN)

# Get header
head <- readLines( inputFile, n=1 )
head_n <- unlist( strsplit( head, "\t" ) )

# if we are to filter to stdout, then print header
if( opt$filter ) {
  writeLines(head, stdout())
}

# check if header is correct..  TODO

# Indexes of samples of interest
repAind <- which( head_n %in% firstRepSet  )
repBind <- which( head_n %in% secondRepSet )

if(length(repAind) == 0 ||
	length(repBind) == 0) { 
   print_help(parser)
   stop("[vast diff error]: Incorrect sampleNames given, One or more do not exist!!!\n") 
}

# Indexes of Quals
repA.qualInd <- repAind + 1
repB.qualInd <- repBind + 1

# make sure this succeeded  TODO

# CONST
alphaList <- seq(0,1,0.01)

### TMP OUT
pdf(paste(c(opt$pdf, ".pdf"), collapse=""), width=7, height=3)

### BEGIN READ INPUT ###
# Iterate through input, 'nLines' at a time to reduce overhead/memory
while(length( lines <- readLines(inputFile, n=opt$nLines) ) > 0) { 

  # use parallel computing to store plots in plotListed
  # then print them to the pdf afterwards before next chunk of nLines from file.
  plotListed <- vector("list", length(lines))
  eventTitleListed <- vector("list", length(lines))

  plotListed <- mclapply(1:length(lines), function(i) {
 
    tabLine <- unlist( strsplit( lines[i], "\t" ) )
	 #writeLines(paste(tabLine[repA.qualInd], collapse="\t"), stderr());
	
	 # Posterior parameters... Prior given from command line --alpha, --beta
	 shapeFirst <- lapply( tabLine[repA.qualInd], function(x) { 
										parseQual(x, opt$alpha, opt$beta) 
								} )
	 shapeSecond <- lapply( tabLine[repB.qualInd], function(x) {
										parseQual(x, opt$alpha, opt$beta)
								} )

    totalFirst <- unlist(lapply( shapeFirst, function(x) { x[1] + x[2] }))
    totalSecond <- unlist(lapply( shapeSecond, function(x) { x[1] + x[2] }))

    # if no data, next;
    if( all(totalFirst < (opt$minReads + opt$alpha + opt$beta)) ||
		  all(totalSecond < (opt$minReads + opt$alpha + opt$beta)) ) {
		return(NULL)
	 }

	 # Sample Posterior Distributions
	 psiFirst <- lapply( shapeFirst, function(x) {
      #sample here from rbeta(N, alpha, beta)
      rbeta(opt$size, shape1=x[1], shape2=x[2])
    })
	 psiSecond <- lapply( shapeSecond, function(x) {
      #sample here from rbeta(N, alpha, beta)
      rbeta(opt$size, shape1=x[1], shape2=x[2])
    })

	 # Create non-parametric Joint Distributions
	 psiFirstComb <- do.call(c, psiFirst)
    psiSecondComb <- do.call(c, psiSecond)

#    print(length(psiFirstComb))

    # if they aren't paired, then shuffle the joint distributions...
    if( !opt$paired ) {
      psiFirstComb <- shuffle(psiFirstComb)
	   psiSecondComb <- shuffle(psiSecondComb)
    }

	 # get emperical posterior median of psi
    medOne <- median(psiFirstComb)
    medTwo <- median(psiSecondComb)

    # look for a max difference given prob cutoff...
    if(medOne > medTwo) {
      max <- maxDiff(psiFirstComb, psiSecondComb, opt$prob)
    } else {
      max <- maxDiff(psiSecondComb, psiFirstComb, opt$prob)
    }
#    writeLines(lines[i], stderr()) ### DEBUGGING
    # check for significant difference
    if(max < opt$minDiff) { return(NULL) } # or continue...

    # SIGNIFICANT from here on out:
	 if( opt$filter ) { 
		writeLines(lines[i], stdout())
    }

    eventTitle <- paste(c("Gene: ", tabLine[1], "     ", "Event: ", tabLine[2]), collapse="")
#    eventTitleListed[[i]] <- paste(c("Gene: ", tabLine[1], "     ", "Event: ", tabLine[2]), collapse="")

	 # Print visual output to pdf;
#    if( opt$pdf ) {
      if( medOne > medTwo ) {
        retPlot <- plotDiff(psiFirstComb, psiSecondComb, max, medOne, medTwo, sampOneName, sampTwoName , FALSE)
      } else {
        retPlot <- plotDiff(psiSecondComb, psiFirstComb, max, medTwo, medOne, sampTwoName, sampOneName , TRUE)
      }
#    }


	 return(list(retPlot, eventTitle))  #return of mclapply function
  }, mc.cores=opt$cores) #End For

  for(it in 1:length(lines)) {
  # PRINT LIST OF PLOTS.
    if(is.null(plotListed[[it]][[1]])) { next; }
    plotPrint(plotListed[[it]][[2]], plotListed[[it]][[1]])
  }

} #End While

garbage <- dev.off()

q(status=0)
