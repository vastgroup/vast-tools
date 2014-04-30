#!/usr/bin/env Rscript
#
# PSI Plotter script

version <- function() {
  return("0.4_vast")
}

print_help <- function() {
  text <- "**** PSI Plotter ****
Script for generating PSI plots across samples (aka Nuno plots).

Usage: ./PSI_Plotter.R PSI_Input.tab[.gz] [Tissue_Groups.txt]

Arguments:
  1) Input PSI data - one AS event per row - using the standard PSI format
    e.g. GENE  EVENT  COORD  LENGTH FullCO  COMPLEX  Tissue1_PSI Tissue1_Q ... 
  2) [optional] Tissue group file or species (currently supports Hsa or Mmu)
    Use for customizing the order and colors of the scatter plot.

Options:
  --version     Display the version    
  --help        This help message

Output:
  A PDF file will be created with one PSI plot per page.

Test run:
  ./PSI_Plotter.R test_data/INCLUSION_LEVELS-ALL3m-Mmu89-SELECTED.test.tab \\
    test_data/Tissues.Mmu.txt
"
  writeLines(text, stderr())
  write(paste("Version:", version()), stderr())
  write("Updated: 2014-04-30", stderr())
}

preprocess_sample_colors <- function(psi, database) {
   # Given a list of sample names from PSI data, generate a corresponding color 
   # code sequence and re-order the columns according to a specified order (see
   # below).
   # 
   # Color codes and sample ordering are taken from a "master" samples database
   # in the following format:
   # Order    SampleName    GroupName    RColorCode
   # 1        Ooctye        EarlyDev     36
   # 2        Embr_2C       EarlyDev     36
   # etc..
   #
   # The RColorCode value corresponds to the index of the vector produced by
   # colors(). For example, RColorCode = 36 corresponds to:
   # > cols <- colors()
   # > mycolour <- cols[36]
   #
   # Args:
   #    psi: a n x m data frame of PSI values where n is the number of AS events
   #          and m is the number of samples
   #    database: filename of the samples database in tab-delimited format
   #
   # Returns:
   #    A list containing the re-ordered data.frame ("data") and a sequence of
   #    corresponding color codes ("col")
   R <- list()

   if (is.null(tissueFile)) {
        mycols <- rep("black", ncol(psi))
        R <- list(data=psi, col=mycols, group.index=NULL, group.col=NULL)
   } else {
       db <- read.csv(database, sep="\t")
       
       # check input file
       if (ncol(db) < 4) {
         stop("The tissues database file is not formatted correctly. Please
          double check")
       }
       
       # check if all samples in input data is in the database
       unk.samples <- colnames(psi) %in% db$SampleName
       if (!all(unk.samples)) {
         s <- colnames(psi)[!unk.samples]
         stop(paste("The following samples are not in the tissues database:", 
                    paste(s, collapse=", ")))
       }
      
       # keep only tissue groups that are present in input data
       # (to take into account samples that might have been excluded)
       db <- db[db$SampleName %in% colnames(psi),]
       
       # Re-order the PSI table
       db <- db[order(db$Order),]
       db$Order <- 1:nrow(db)
       new.column.idx <- sapply(db$SampleName, 
                                function(x) which(colnames(psi) == x))
       
       psi.new <- psi[,new.column.idx]
       
       # Generate a corresponding color code sequence
       mycols <- colors()[db$RColorCode]
       names(mycols) <- db$SampleName
       
       # Store indices of columns for each group
       groups <- unique(db$GroupName)
       mygroups <- list()
       mygroupcol <- rep(NA, length(groups))
       for (i in 1:length(groups)) {
          mygroups[[i]] <- which(colnames(psi.new) %in%  
                                     db[db$GroupName == groups[i],"SampleName"])
          mygroupcol[i] <- colors()[unique(db[db$GroupName == groups[i], 
                                              "RColorCode"])]
       }
       names(mygroups) <- groups
       names(mygroupcol) <- groups
       R <- list(data=psi.new, col=mycols, 
                   group.index=mygroups, group.col=mygroupcol)
   }
   
   return(R)
}

#### Arguments #################################################################
# - Input file
# - Tissue group or Species

args <- commandArgs(TRUE)

if (length(args) < 1) {
  print_help()
  stop("Missing arguments")
}
if (args[1] == "--help") {
  print_help()
  stop("Terminating")
}
if (args[1] == "--version") {
  write(paste("Version:", version()), stderr())
  stop("Terminating")
}

file <- args[1]
if (!file.exists(file))
  stop(paste("Input PSI file", file, "doesn't exist!"))

tissueFile <- NULL
if (length(args) == 2) {
    tissueFile <- args[2]
    if (!file.exists(tissueFile))
      stop(paste("Tissue Group file", tissueFile, "doesn't exist!"))
}

write(paste("PSI Plotter - Version", version()), stderr())
write(paste("\n// Input file:", file), stderr())
write(paste("// Tissue Group file:", tissueFile), stderr())

#### Format input data #########################################################

all_events <- read.csv(file, sep="\t")

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
    if (length(na) > 0) 
      psi[na, i] <- NA
  }
  return(psi[, seq(1, ncol(psi), 2)])
}

format_table <- function(m) {
  # Format table to keep only PSIs and convert exon metadata as rownames
  id <- paste(m$COMPLEX, m$GENE, m$COORD, m$LENGTH, sep="=")
  
  # Extract PSIs
  psi <- convert_psi(m[,7:ncol(m)])
  rownames(psi) <- id
  return(psi)
}

if (!grepl("^GENE", colnames(all_events)[1])) {
  stop("Invalid column names. Does your input file contain the correct header?")
}
write("// Brewing some coffee...", stderr())

write("// Formatting input data for plotting...", stderr())
PSIs <- format_table(all_events)
# Call function to re-order columns of PSI data
#
# returns a list containing four elements:
#   data        - the PSI data with sample columsn re-ordered
#   col         - vector of colours that will be plotted
#   group.index - list of indices for each sample group (e.g. ESC, Neural, etc.)
#   group.col   - corresponding color for sample group
reordered.PSI <- preprocess_sample_colors(PSIs, tissueFile)
PSIs <- as.matrix(reordered.PSI$data)
ALLev <- row.names(PSIs)
samples <- colnames(PSIs)

write(paste("//", ncol(PSIs), "samples detected"), stderr())

#### Prepare plotting ##########################################################
write("// Plotting...", stderr())

# assign list of colors
supercolors <- reordered.PSI$col

# Set output file
outfile <- sub("\\.[^.]*(\\.gz)?$", ".PSI_plots.pdf", file)

pdf(outfile, width = 8.5, height = 5.5)
par(mfrow = c(1,1), las = 2) #3 graphs per row; 2=label always perpendicular to the axis
for (i in 1:nrow(PSIs)) {
  plot(as.numeric(PSIs[i,]),
       col=supercolors,
       pch=20,
       main=rownames(PSIs)[i],
       ylab="PSI", xlab="", xaxt="n",
       ylim=c(1,100),
       cex=0.8, cex.main=0.9, cex.axis=0.8)
  axis(1, at=seq(1, ncol(PSIs), by=1), labels = FALSE)
  text(seq(1, ncol(PSIs), by=1), 
       par("usr")[3] - 3.5, 
       labels = samples, 
       srt = 45, adj=c(1,1), xpd = TRUE,cex=0.5)
  
  if (!is.null(tissueFile)) {
      abline(h=mean(PSIs[i, reordered.PSI$group.index[["ESC"]] ], na.rm=TRUE), 
             col=reordered.PSI$group.col["ESC"], lwd=0.5)
      abline(h=mean(PSIs[i, reordered.PSI$group.index[["Neural"]] ], na.rm=TRUE),
             col=reordered.PSI$group.col["Neural"], lwd=0.5)
      abline(h=mean(PSIs[i, reordered.PSI$group.index[["Muscle"]] ], na.rm=TRUE),
             col=reordered.PSI$group.col["Muscle"], lwd=0.5)
      abline(h=mean(PSIs[i, reordered.PSI$group.index[["Tissues"]] ], na.rm=TRUE),
             col=reordered.PSI$group.col["Tissues"], lwd=0.5)
  }

  abline(v=1:ncol(PSIs), col="grey", lwd=0.3, lty=2)
  abline(h=seq(0,100,10), col="grey", lwd=0.3, lty=2)
}
dev.off()

write("// Done!\n", stderr())
write(paste("//", nrow(PSIs), "plots are saved in:", outfile), stderr())
####
