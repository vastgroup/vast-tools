#!/usr/bin/env Rscript

# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca

# Updates: Manuel Irimia, 2015-present
# mirimia@gmail.com

argv <- commandArgs(trailingOnly = F)
scriptPath <- dirname(sub("--file=","",argv[grep("--file",argv)]))

joinStr <- function(x,y) {
  return(paste(c(as.character(x), as.character(y)), collapse=""))
}

writeLines(joinStr("Using ", R.Version()$version.string));

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

humanDbFile <- "vastdb.hsa.22.06.16.tar.gz"
mouseDbFile <- "vastdb.mmu.22.06.16.tar.gz"
chickenDbFile <- "vastdb.gga.13.11.15.tar.gz"
planariaDbFile <- "vastdb.sme.13.11.15.tar.gz"


humanUrl <- joinStr("http://vastdb.crg.eu/libs/", humanDbFile)
mouseUrl <- joinStr("http://vastdb.crg.eu/libs/", mouseDbFile)
chickenUrl <- joinStr("http://vastdb.crg.eu/libs/", chickenDbFile)
planariaUrl <- joinStr("http://vastdb.crg.eu/libs/", planariaDbFile)
#

writeLines("Looking for VAST Database [VASTDB]")
if(!file.exists("VASTDB")) {
  writeLines("Cannot find 'VASTDB'.. Do you want me to download it for you? [y/n]")
  auto <- readLines(file("stdin"),1)
  close(file("stdin"))
  if(as.character(auto) == 'y') {
    writeLines("OK I will try... Please choose database [hg19/mm9/galGal3/smed31/all] [h/m/g/p/a]")
    db <- readLines(file("stdin"),1)
    db <- as.character(db)
    close(file("stdin"))
    if(db == 'h' || db == 'a') {
      downloadDb(humanUrl, humanDbFile)
    }
    if(db == 'm' || db == 'a') {
      downloadDb(mouseUrl, mouseDbFile)
    }
    if(db == 'g' || db == 'a') {
      downloadDb(chickenUrl, chickenDbFile)
    }
    if(db == 'p' || db == 'a') {
      downloadDb(planariaUrl, planariaDbFile)
    }
  }
} else {
  writeLines("Found what appears to be VASTDB.. OK")
}

# custom install from include.R
loadPackages(c("MASS", "getopt", "optparse", "RColorBrewer", "reshape2", "ggplot2", "grid", "parallel", "devtools"), local.lib=paste(c(scriptPath,"/R/Rlib"), collapse=""))

# install github packages
if (!require('psiplot', character.only=T)) {
  with_lib(new=paste(c(scriptPath,"/R/Rlib"), collapse=""), install_github('kcha/psiplot'))
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
