#!/usr/bin/env Rscript

# Author: Kevin Ha, 2014
# k.ha@mail.utoronto.ca

# Copyright (C) 2014 Kevin Ha
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

MAX_ENTRIES <- 1000

args <- commandArgs(trailingOnly = F)
scriptPath <- dirname(sub("--file=","", args[grep("--file",args)]))
source(paste(c(scriptPath,"/Rlib/include.R"), collapse=""))

loadPackages(c("optparse"), local.lib=paste(c(scriptPath,"/Rlib"),
                                                       collapse=""))

#### Arguments #################################################################
# - Input file
# - Tissue group or Species

args <- commandArgs(TRUE)

desc <- "Script for generating PSI plots (scatterplot) across samples.

Input:
  PSI data - one AS event per row - using the standard PSI format
      e.g. GENE  EVENT  COORD  LENGTH FullCO  COMPLEX  Tissue1 Tissue1_Q ...
  Recommended to use only a subset of AS events instead of the full table
  otherwise the resulting PDF file will be very large. Use option -m/--max to limit
  the maximum number of plots to generate.

  PSI values that are \"NA\" or have \"NA\" quality scores will not be plotted
  (not point will be drawn).

  If no input file is provided, standard input will be used.

Output:
  A PDF file will be created with one PSI plot per page.

Customizing plots [optional]:
  The color and ordering of samples can be customized by supplying a plot
  configuration file (e.g. psiplotter.config). This file is tab-delimited and must be
  manually created. The format of psiplotter.config is the following (the header
  line is required):
  Order    SampleName    GroupName    RColorCode
  1        Ooctye        EarlyDev     blue
  2        Embr_2C       EarlyDev     red
  etc..

  Order 	: The ordering of the samples from left to right.
  SampleName 	: Name of the sample. MUST match sample name in input table*.
  GroupName	: Group name. Use for plotting the average PSI of samples belonging
    to the same group (need to use option -u/--group-means)
  RColorCode	: Any of the following R color specifications:
    1) color name (see colors())
    2) hex color code (#rrggbb)

  *The samples under SampleName MUST MATCH the names in the PSI input table.
  Only the samples listed in the config file will be represented in the
  resulting plots. Other samples in the input file but not in the config
  file will be ignored. This may be useful if you want to customize the
  type of samples in your plots.
"

option.list <- list(
  make_option(c("-v", "--verbose"), type = "logical", default = TRUE,
              meta="TRUE|FALSE",
              help="Enable verbose [%default]"),
  make_option(c("-c", "--config"), type = "character", default = NULL,
              help = "Plot configuration file. Used for customizing order and color
            [%default]"),
  make_option(c("-m", "--max"), type = "integer", default = MAX_ENTRIES,
              help = "Maximum number of AS events to plot [first %default]"),
  make_option(c("-l", "--gridLines"), type = "logical", default = TRUE,
              meta="TRUE|FALSE", dest = "gridLines",
              help = "Show grid lines [%default]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              meta="DIR",
              help = "Output directory where pdf will be saved
                    [default is same location as input data]"),
  make_option(c("-E", "--noErrorBar"), type = "logical", default = FALSE,
              meta="TRUE|FALSE", dest = "noErrorBar",
              help = "Do not plot 95% confidence interval as error bars [%default]"),
  make_option(c("-u", "--groupMeans"), type = "logical", default = FALSE,
              meta="TRUE|FALSE", dest = "plotGroupMeans",
              help = "Plot mean PSIs for groups defined in config file. Requires
              --config option. [%default]"),
  make_option(c("-W", "--width"), type = "numeric", default = NULL, dest = "width",
              help = "Width of graphics region in inches (similar to width in
              pdf()) [%default]"),
  make_option(c("-H", "--height"), type = "numeric", default = NULL,
              dest = "height",
              help = "Height of graphics region in inches (similar to height in
              pdf()) [%default]"),
  make_option(c("--gene"), type = "character", default = NULL,
              dest = "gene",
              help = "Filter events by the GENE column. Can be any 
              valid R regular expression. [%default]")
)
parser <- OptionParser(option_list = option.list,
                       desc = desc,
                       usage = "usage: %prog [options] INCLUSION_LEVELS.tab")
opt <- parse_args(parser, args = args, positional_arguments = TRUE)

# Load remaining packages
loadPackages(c("psiplot", "methods"), local.lib=paste(c(scriptPath,"/Rlib"),
                                                       collapse=""))

# Check for correct version of psiplot
v <- as.character(packageVersion('psiplot'))
# The minimum required version of psiplot
required <- '2.0.0'
if (compareVersion(v, required) == -1) {
  stop(paste("Your version of psiplot, the R package that this script uses, is",
             v, "and is out of date.\n",
             "Please update before running this script.",
             "See https://github.com/kcha/psiplot for additional details.")
  )
}

# Check options
if (length(opt$args) == 0) {
  print_help(parser)
  stop("Missing arguments")
}

using_stdin <- FALSE
file <- opt$args[1]
if (file == "-") {
    file <- file('stdin')
    using_stdin <- TRUE
} else if (!file.exists(file)) {
  stop(paste("Input PSI file", file, "doesn't exist!"))
}

config_file <- opt$options$config
if (!(is.null(config_file) || file.exists(config_file)))
  stop(paste("Tissue Group file", config_file, "doesn't exist!"))

verbPrint <- function(s) {
  if (opt$options$verbose) {
    write(s, stderr())
  }
}

verbPrint(paste("\nPSI Plotter"))
verbPrint(paste("\n// Input file:", ifelse(using_stdin, "STDIN", file)))
verbPrint(paste("// Tissue Group file:",
                ifelse(is.null(config_file), "Did not provide", config_file)))

all_events <- read.delim(file, stringsAsFactors=FALSE)
if(is.null(config_file)) {
  config <- NULL
  nsamples <- (ncol(all_events) - 6)/2
} else {
  config <- read.delim(config_file, stringsAsFactors=FALSE)
  nsamples <- length(which(colnames(all_events) %in% config$SampleName))
}

# Filter list if necessary
if (!is.null(opt$options$gene)) {
  all_events <- all_events[grep(opt$options$gene, all_events$GENE),]
  
  if (nrow(all_events) == 0) {
    stop("No matching events found.")
  } else {
    verbPrint(paste("// Filtered", nrow(all_events), "events that match pattern", 
                    opt$options$gene))
  }
}


# Perform some checks #########################################################
if (!grepl("^GENE", colnames(all_events)[1])) {
  stop("Invalid column names. Does your input file contain the correct header?")
}

if (nrow(all_events) > opt$options$max) {
  warning(paste("Too many entries in input file. Plotting only the first",
                opt$options$max, ". Try splitting your input file into smaller",
                "files and running them separately."))
}

#### Prepare plotting ##########################################################
verbPrint("// Plotting...")
if (!is.null(opt$options$config)) {
    verbPrint(paste("// Plot group means as horizontal lines:",
                opt$options$plotGroupMeans))
}

# Set output file
outfile <- paste0("PSI_plots.pdf")
if (!using_stdin) {
  outfile <- sub("\\.[^.]*(\\.gz)?$", ".PSI_plots.pdf", basename(file))
}

# Check if output directory was specified
if (is.null(opt$options$output)) {
  if (!using_stdin) {
    outfile <- file.path(dirname(file), outfile)
  }
} else {
  # Create directory if necessary
  if (!file.exists(opt$options$output))
    dir.create(opt$options$output, recursive = TRUE)
  outfile <- file.path(opt$options$output, outfile)
}

# Set width and height of graphics
W <- ifelse(is.null(opt$options$width), 3.5 + 0.1*max(0, nsamples - 8),
            opt$options$width)
H <- ifelse(is.null(opt$options$height), 3.2 + 0.05*max(0, nsamples - 8),
            opt$options$height)
verbPrint(paste("// Width = ", round(W, 2), "in, Height =", round(H, 2), "in"))

#### Plot ##########################################################
pdf(outfile, width = W, height = H)
par(mfrow = c(1,1), las = 2) #3 graphs per row; 2=label always perpendicular to the axis
nplot <- min(nrow(all_events), opt$options$max)
pb <- txtProgressBar(style = 3, file = stderr())
for (i in 1:nplot) {
  result <- plot_event(all_events[i,], config = config,
             errorbar = !opt$options$noErrorBar,
             groupmean = opt$options$plotGroupMeans,
             gridlines = opt$options$gridLines,
             cex.xaxis = 10, cex.yaxis = 10, cex.main = 8)
  setTxtProgressBar(pb, i/nplot)
}
setTxtProgressBar(pb, 1)
verbPrint("")
dev.off()

verbPrint("// Done!\n")
verbPrint(paste("//", nplot, "plots are saved in:", outfile))
####
