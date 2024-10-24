#!/usr/bin/env Rscript

# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca

# Updates: Manuel Irimia, 2015-present
# mirimia@gmail.com

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-p", "--prompt"), action="store_true", default=FALSE, type="logical", help="User prompt during installation [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="prompt",  type="logical", help="Quiet installation, no prompt"),
  make_option(c("-f", "--file"), action="store", default=".", type='character', help="From where to run install script")
)

opt <- parse_args(OptionParser(option_list=option_list))

scriptPath <- dirname( opt$file )

joinStr <- function(x,y) {
  return(paste(c(as.character(x), as.character(y)), collapse=""))
}

writeLines(joinStr("Using ", R.Version()$version.string))

# Source Rlib.
source(paste(c(scriptPath,"/R/Rlib/include.R"), collapse=""))

downloadDb <- function(speUrl, speFile) {
   system("mkdir -p VASTDB")
   if(system(joinStr("wget ", speUrl)) > 0) {
     stop(joinStr("Cannot download ", speUrl))
   }
   if(system(paste("tar xzvf ", speFile, " -C VASTDB")) > 0) {
     stop(joinStr("Cannot tar xzvf ", speFile))
   }
   if(system(joinStr("rm ", speFile)) > 0) {
     stop(joinStr("Cannot rm ", speFile))
   }
}

# Predefined database selection (without user prompt)
vastdbFiles <- c(
  "vastdb.hsa.23.06.20.tar.gz",
  "vastdb.hs2.23.06.20.tar.gz",
  "vastdb.mmu.23.06.20.tar.gz",
  "vastdb.mm2.23.06.20.tar.gz"
  # Add any additional databases here if needed
)

for (db in vastdbFiles) {
  Url <- paste("https://vastdb.crg.eu/libs/", db, sep = "")
  downloadDb(Url, db)
}

# custom install from include.R
loadPackages(c("MASS", "getopt", "optparse", "RColorBrewer", "reshape2", "ggplot2", "grid", "parallel", "devtools"), local.lib=paste(c(scriptPath,"/R/Rlib"), collapse=""))

# install github packages
if (!require('psiplot', character.only=T)) {
  install_github('kcha/psiplot', lib=paste(c(scriptPath,"/R/Rlib"), collapse=""))
  writeLines(sprintf("psiplot has been installed locally in R/Rlib!"), stderr())
}

if(Sys.chmod(paste(c(scriptPath, "/vast-tools"), collapse=""), mode = "755")) {
  writeLines("Setting vast-tools permissions... success!", stderr());
}

if(system("which bowtie") > 0) {
  stop("Cannot find 'bowtie' in path!!!  Please install this properly (e.x. /usr/bin/ ) or supply -bowtieProg flag to vast-tools!");
} else {
  writeLines("Found bowtie... Everything looks --OK");
}

q(status=0)
