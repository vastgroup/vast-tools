#!/usr/bin/Rscript --vanilla

# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca

argv <- commandArgs(trailingOnly = F)
scriptPath <- dirname(sub("--file=","",argv[grep("--file",argv)]))

# Source Rlib.
source(paste(c(scriptPath,"/R/Rlib/include.R"), collapse=""))

# custom install from include.R
loadPackages(c("optparse", "RColorBrewer", "reshape2", "ggplot2", "grid", "parallel"), local.lib=paste(c(scriptPath,"/R/Rlib/"), collapse=""))

if(Sys.chmod(paste(c(scriptPath, "/vast-tools"), collapse=""), mode = "755")) {
  writeLines("Setting vast-tools permissions... success!", stderr());
}

if(system("which bowtie") > 0) {
  stop("Cannot find 'bowtie' in path!!!  Please install this properly (e.x. /usr/bin/ ) or supply -bowtieProg flag to vast-tools!");
} else {
  writeLines("Found bowtie... Everything looks --OK");
}

q(status=0)
