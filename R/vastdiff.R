#!/usr/bin/Rscript --vanilla

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
loadPackages(c("optparse", "RColorBrewer", "reshape2", "ggplot2", "grid"))

argv <- commandArgs(TRUE)

# optparse..
option.list <- list(
    make_option(c("-v", "--verbose"), type = "logical", default = TRUE,
        help="Enable verbose [%default]"),
    make_option(c("-i", "--input"), type = "character", default = "INCLUSION_LEVELS",
        help = "Exact or Partial match to PSI table in output directory [%default]"),
    make_option(c("-a", "--replicateA"), type = "character", default = NULL,
        help = "Required, SampleA@SampleB@SampleC etc.. [first %default]"),
    make_option(c("-b", "--replicateB"), type = "character", default = NULL,
        help = "Required, SampleA@SampleB@SampleC etc.. [first %default]"),
    make_option(c("-f", "--filter"), type = "logical", default = TRUE,
        help = "Filter output for differential events only [first %default]"),
    make_option(c("-p", "--plot"), type = "logical", default = FALSE,
        help = "Plot visual output for differential events [first %default]"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
        help = "Output directory, passed from vast [%default]")
)

print(option.list)

parser <- OptionParser(option_list = option.list,
            usage = "usage: %prog -a SampleA@..@SampleD -b SampleF@..@SampleG [options]")
optpar <- parse_args(parser, argv, positional_arguments = TRUE)
opt <- optpar$options

# move to output directory
setwd(opt$output)

# try and find the input file if they aren't exact
if(!file.exists(opt$input)) {
  potentialFiles <- Sys.glob( paste(c("*",opt$input,"*"), collapse="") )
  if( length( potentialFiles ) >= 1) { 
	 # now sort and take the one with the 'biggest' number of samples
    potentialFiles_sort <- rev( sort( potentialFiles ) )
    opt$input <- potentialFiles_sort[1]
  } else {
    # Still can't find input after searching...
    stop("[vast diff error]: No input file given!")
  }
}

# Setting input files.
inputFile <- file( opt$input, 'r' )

if(opt$plot) {
  dir.create(paste(c(opt$output, "/diff_out"), collapse=""))
  q()
}


print_help <- function() {
  text <- ""
  writeLines(text, stderr())
  write("Updated: 2014-04-09", stderr())
}


#-replicatesA=name1@name2@name3 -replicatesB=name4@name5
firstRepSet <- unlist(strsplit( as.character(opt$replicateA) , "@" ))
secondRepSet <- unlist(strsplit( as.character(opt$replicateB), "@" ))

firstRepN <- length(firstRepSet)
secondRepN <- length(secondRepSet)

## INITIALIZE LISTS ##
shapeFirst <- vector("list", firstRepN)
shapeSecond <- vector("list", secondRepN)

psiFirst <- vector("list", firstRepN)
psiSecond <- vector("list", secondRepN)

### READ INPUT ###

# Get header
head <- readLines( inputFile, n=1 )
head_n <- unlist( strsplit( head, "\t" ))

#print(head_n)
#print(firstRepSet)

# Indexes of samples of interest
repAind <- which( head_n %in% firstRepSet  )
repBind <- which( head_n %in% secondRepSet )

#print(repAind)

# Indexes of Quals
repA.qualInd <- repAind + 1
repB.qualInd <- repBind + 1

#write(repA.qualInd, stdout())
#q()

### BEGIN READ INPUT ###
# Iterate through input, 1000 lines at a time to reduce overhead/memory
while(length( lines <- readLines(inputFile, n=1000) ) > 0) { 
  for(i in 1:length(lines)) { 
    tabLine <- unlist( strsplit( lines[i], "\t" ) )
	 #writeLines(paste(tabLine[repA.qualInd], collapse="\t"), stderr());
		 First <- lapply(
								
							 )	  
  } 
}

q()
### END READ INPUT ###

#sample here from rbeta(N, alpha, beta)

psiFirstComb <- do.call(c, psiFirst)
psiSecondComb <- do.call(c, psiSecond)

shuffOne <- shuffle(psiFirstComb)  #unless paired=T
shuffTwo <- shuffle(psiSecondComb) # unless paired=T


#GET FROM INPUT
Sample_1_Name <- "repA"
Sample_2_Name <- "repB"

sampOneName <- paste(c(substr(Sample_1_Name, 1, 4), "(n=", as.character(firstRepN), ")"), collapse="")
sampTwoName <- paste(c(substr(Sample_2_Name, 1, 4), "(n=", as.character(secondRepN), ")"), collapse="")

# calculate the probability that the first dist is > than second
distPlot <- ggplot(melt(as.data.frame(
			do.call(cbind,list(psiFirstComb, psiSecondComb))
			)), aes(fill=variable, x=value))+
			geom_histogram(aes(y=..density..),alpha=0.5, col="grey", position="identity")+
			theme_bw()+xlim(c(0,1))+xlab(expression(hat(Psi)))+
			scale_fill_manual(values=cbb[2:3], labels=c(sampOneName, sampTwoName), name="Samples")


probPlot <- ggplot(as.data.frame(cbind(seq(0,1,0.01), 
				unlist(lapply(seq(0,1,0.01), function(x) { 
					pDiff(shuffTwo, shuffOne, x) 
				})))), aes(x=V1, y=V2))+
				geom_line()+theme_bw()+
				geom_vline(x=maxDiff(shuffTwo, shuffOne), lty="dashed")+
				ylab(expression(P((hat(Psi)[1]-hat(Psi)[2]) > x)))+
				xlab(expression(x))+
				geom_text(x=maxDiff(shuffTwo, shuffOne), y=-0.1, label=maxDiff(shuffTwo, shuffOne))

pdf(sprintf("%s_%d_mer.pdf", argv[3], n), width=7, height=6)
multiplot(distPlot, probPlot, cols=2)
dev.off()

q(status=0)
