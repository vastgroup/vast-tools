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
source(file.path(scriptPath, "Rlib", "preprocess_sample_colors.R"))
source(paste(c(scriptPath,"/Rlib/include.R"), collapse=""))
source(paste(c(scriptPath,"/Rlib/include_diff.R"), collapse=""))

loadPackages(c("optparse"), local.lib=paste(c(scriptPath,"/Rlib"), collapse=""))

#### Arguments #################################################################
# - Input file
# - Tissue group or Species

args <- commandArgs(TRUE)

desc <- "Script for generating PSI plots (scatterplot) across samples.

Input:
  PSI data - one AS event per row - using the standard PSI format
      e.g. GENE  EVENT  COORD  LENGTH FullCO  COMPLEX  Tissue1 Tissue1_Q ... 
  Recommended to use only a subset of AS events instead of the full table
  otherwise the resulting PDF file will be very large. See options for
  customizing the maximum number of plots to generate.

  PSI values that are \"NA\" or have \"NA\" quality scores will not be plotted
  (not point will be drawn).

  If no input file is provided, standard input will be used.

Output:
  A PDF file will be created with one PSI plot per page.

Customizing plots [optional]:
  The color and ordering of samples can be customized by supplying a plot
  configuration file (psiplotter.config). This file is tab-delimited and must be
  manually created. The format of psiplotter.config is the following (the header
  line is required): 
  Order    SampleName    GroupName    RColorCode
  1        Ooctye        EarlyDev     blue
  2        Embr_2C       EarlyDev     red
  etc..
 
  Order 	: The ordering of the samples from left to right.
  SampleName 	: Name of the sample. MUST match sample name in input table.
  GroupName	: Group name. Use for plotting the average PSI of samples belonging
    to the same group (need to use option -u/--group-means)
  RColorCode	: Any of the three kinds of R color specifications:
    1) color name (as specified by colors())
    2) hex color code (#rrggbb)

  The samples under SampleName MUST MATCH the names in the PSI input table.
  Only the samples listed in the config file will be represented in the 
  resulting plots. Other samples in the PSI table but not in the config 
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
  make_option(c("-o", "--output"), type = "character", default = NULL,
              meta="DIR",
              help = "Output directory where pdf will be saved
                    [default is same location as input data]"),
  make_option(c("-E", "--noErrorBar"), type = "logical", default = FALSE,
              meta="TRUE|FALSE", dest = "noErrorBar",
              help = "Do not plot error bars [%default]"),
  make_option(c("-u", "--group-means"), type = "logical", default = FALSE,
              meta="TRUE|FALSE", dest = "plotGroupMeans",
              help = "Plot mean PSIs for groups defined in config file. Requires
              --config option. [%default]")
)
parser <- OptionParser(option_list = option.list,
                       desc = desc,
                       usage = "usage: %prog [options] INCLUSION_LEVELS.tab")
opt <- parse_args(parser, args = args, positional_arguments = TRUE)

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

verbPrint(paste("PSI Plotter"))
verbPrint(paste("\n// Input file:", ifelse(using_stdin, "STDIN", file)))
verbPrint(paste("// Tissue Group file:", 
                ifelse(is.null(config_file), "Did not provide", config_file)))

#### Format input data functions ##############################################
convert_psi <- function(t) {
  # Helper function to filter and return PSI values
  # PSIs are converted to NA if first coverage code is 'N'
  # e.g. PSI=100, Coverage=N,N,N,OK,S ---> PSI=NA
  #
  # Input: original PSI plus quality scores WITHOUT the first 7 columns
  
  stopifnot(ncol(t) %% 2 == 0)
  psi <- t
  
  for (i in seq(1, ncol(psi), 2)) {
    cov <- strsplit(as.character(psi[,i+1]), split = ",")
    cov <- sapply(cov, "[", 1)
    
    na <- which(cov == "N")
    if (length(na) > 0) {
      psi[na, i] <- NA
      psi[na, i+1] <- NA
    }
  }
  return(psi)
}

get_beta_ci <- function(q) {
  # Helper function to filter and return confidence intervals based on beta
  # distribution from Q scores
  # 
  # Input: original PSI plus quality scores WITHOUT the first 7 columns
  
  parameters <- sapply(q, function(x) parseQual(as.character(x)))
  ci <- lapply(1:ncol(parameters), function(j) betaCISample(parameters[1,j], 
                                                            parameters[2,j]))
  ci <- do.call("rbind", ci) * 100
  return(ci)
}

format_table <- function(m) {
  # Format table to keep only PSIs and convert exon metadata as rownames
  id <- paste(m$COMPLEX, m$GENE, m$COORD, m$LENGTH, sep="|")
  
  # Extract PSIs
  r <- convert_psi(m[,7:ncol(m)])
  rownames(r) <- id
  return(r)
}

all_events <- read.csv(file, sep="\t")

# Perform some checks #########################################################
if (!grepl("^GENE", colnames(all_events)[1])) {
  stop("Invalid column names. Does your input file contain the correct header?")
}

if (nrow(all_events) > opt$options$max) {
  warning(paste("Too many entries in input file. Plotting only the first",
                opt$options$max, ". Try splitting your input file into smaller",
                "files and running them separately."))
}

# Format input data ###########################################################
verbPrint("// Formatting input data for plotting")
all_events_formatted <- format_table(all_events)
# Call function to re-order columns of PSI data
#
# returns a list containing four elements:
#   data        - the PSI data with sample columsn re-ordered
#   col         - vector of colours that will be plotted
#   group.index - list of indices for each sample group (e.g. ESC, Neural, etc.)
#   group.col   - corresponding color for sample group
reordered <- preprocess_sample_colors(all_events_formatted, config_file)
verbPrint(paste("//", ncol(reordered$data), "out of", 
                ncol(all_events_formatted) / 2, "samples selected"))
PSIs <- as.matrix(reordered$data)
#ALLev <- row.names(PSIs)
samples <- colnames(PSIs)

#### Prepare plotting ##########################################################
verbPrint("// Plotting...")
if (!is.null(opt$options$config)) {
    verbPrint(paste("// Plot group means as horizontal lines:", 
                opt$options$plotGroupMeans))
}

# tissuegroups <- c("ESC", "Neural", "Muscle", "Tissues")

# assign list of colors
supercolors <- reordered$col

# Set output file
outfile <- "PSI_plots.pdf"
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

pdf(outfile, width = 8.5, height = 5.5)
par(mfrow = c(1,1), las = 2) #3 graphs per row; 2=label always perpendicular to the axis
nplot <- min(nrow(PSIs), opt$options$max)
for (i in 1:nplot) {
  # Set plot title
  event <- strsplit(rownames(PSIs)[i], split = "\\|")[[1]]
  title <- sprintf("%s (position = %s, length = %s, type = %s)", 
    event[2], event[3], event[4], event[1])


  # Set up plot
  plot(NA,
       main=title,
       ylab="PSI", xlab="", xaxt="n",
       ylim=c(1,100), xlim=c(1, ncol(PSIs)),
       cex.main=0.9, cex.axis=0.8)
  axis(1, at=seq(1, ncol(PSIs), by=1), labels = FALSE)
  text(seq(1, ncol(PSIs), by=1), 
       par("usr")[3] - 3.5, 
       labels = samples, 
       srt = 45, adj=c(1,1), xpd = TRUE,cex=0.6)
  
  
  # Draw error bars
  if (! opt$options$noErrorBar) {
    ci <- get_beta_ci(reordered$qual[i,])
    ci[which(is.na(reordered$data[i,])),] <- NA
    
    arrows(1:ncol(PSIs), ci[,1],
           1:ncol(PSIs), ci[,2],
           length = 0.025,
           angle = 90,
           code = 3, col = as.character(supercolors))
  }
  
  # Draw horizontal lines for groups
  if (!is.null(config_file) && opt$options$plotGroupMeans) {
    seen <- vector()
    groups <- names(reordered$group.index)
    for (t in 1:length(groups)) {
      abline(h=mean(PSIs[i, reordered$group.index[[groups[t]]] ], 
                    na.rm=TRUE), 
             col=reordered$group.col[groups[t]], lwd=0.5)
      seen <- append(seen, paste(groups[t]))
    }
    
    # plot legend for mean group values
    if (length(seen) > 0) {
      legend_position <- ifelse(PSIs[i,ncol(PSIs)] > 50, "bottomright", "topright")
      legend(legend_position, legend = seen, lty = 1, col = reordered$group.col,
             title = "Group Means", cex = 0.7, ncol = 2)  
    }
  }
  
  # Draw grid lines
  abline(v=1:ncol(PSIs), col="grey", lwd=0.3, lty=2)
  abline(h=seq(0,100,10), col="grey", lwd=0.3, lty=2)
  
  # Draw PSIs
  points(1:ncol(PSIs), as.numeric(PSIs[i,]), col=as.character(supercolors), 
         pch=20, cex = 1)
}
dev.off()

verbPrint("// Done!\n")
verbPrint(paste("//", nplot, "plots are saved in:", outfile))
####
