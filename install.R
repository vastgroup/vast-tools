#!/usr/bin/env Rscript

# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca

# Updates: Manuel Irimia, 2015-present
# mirimia@gmail.com

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-p", "--prompt"), action="store_true", default=TRUE, type="logical", help="User prompt during installation [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="prompt",  type="logical", help="Quiet installation, no prompt"),
  make_option(c("-f", "--file"), action="store", default=".", type='character', help="From where to run install script")
)

opt <- parse_args(OptionParser(option_list=option_list))

scriptPath <- dirname( opt$file )

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

### To be updated with every release:
vastdbFiles<-c(
"vastdb.hsa.20.12.19.tar.gz",
"vastdb.hs2.20.12.19.tar.gz",
"vastdb.mmu.20.12.19.tar.gz",
"vastdb.mm2.20.12.19.tar.gz",
"vastdb.bta.20.12.19.tar.gz",
"vastdb.gg3.20.12.19.tar.gz",
"vastdb.gg4.06.04.20.tar.gz",
"vastdb.xt1.06.04.20.tar.gz",
"vastdb.dre.20.12.19.tar.gz",
"vastdb.bl1.20.12.19.tar.gz",
"vastdb.spu.20.12.19.tar.gz",
"vastdb.dme.20.12.19.tar.gz",
"vastdb.sma.20.12.19.tar.gz",
"vastdb.cel.20.12.19.tar.gz",
"vastdb.sme.20.12.19.tar.gz",
"vastdb.nve.20.12.19.tar.gz",
"vastdb.ath.20.12.19.tar.gz");


if ( opt$prompt ) {
  writeLines("Looking for VAST Database [VASTDB]")
  auto <- "invalid"
  if(!file.exists("VASTDB")) {
    cat("Cannot find 'VASTDB'.. Do you want to download it? [y/n]: ")
    auto <- readLines(file("stdin"),1)
    close(file("stdin"))
  }
  if(file.exists("VASTDB")) {
    cat("Found 'VASTDB'.. Do you want to download yet uninstalled datasets or update installed datasets? [y/n]: ")
    auto <- readLines(file("stdin"),1)
    close(file("stdin"))
  }
  
  if(as.character(auto) == 'y') {
    cat("Please choose one or more (comma-separated) from: 

    1 - Homo sapiens (hg19, Hsa). Scaffold annotation: Ensembl v60.
    2 - Homo sapiens (hg38, Hs2). Scaffold annotation: Ensembl v88.
    3 - Mus musculus (mm9, Mmu). Scaffold annotation: Ensembl v62.
    4 - Mus musculus (mm10, Mm2). Scaffold annotation: Ensembl v88.
    5 - Bos taurus (bosTau6, Bta). Scaffold annotation: Ensembl v76.
    6 - Gallus gallus (galGal3, Gg3). Scaffold annotation: Ensembl v65.
    7 - Gallus gallus (galGal4, Gg4). Scaffold annotation: Ensembl v83.
    8 - Xenopus tropicalis (xenTro3, Xt1). JGI_4.2. Scaffold annotation: Ensembl v84.
    9 - Danio rerio (danRer10, Dre). Zv10. Scaffold annotation: Ensembl v80.
   10 - Branchiostoma lanceolatum (braLan2, Bl1). Bl71nemr. Scaffold annotation: Ensembl Metazoa v46.
   11 - Strongylocentrotus purpuratus (strPur4, Spu). Spur3.1. Scaffold annotation: SpBase.
   12 - Drosophila melanogaster (dm6, Dme). Scaffold annotation: Ensembl Metazoa v26.
   13 - Strigamia maritima (strMar1, Sma). Smar1. Scaffold annotation: Ensembl Metazoa v26.
   14 - Caenorhabditis elegans (ce11, Cel). Scaffold annotation: Ensembl v87.
   15 - Schmidtea mediterranea (schMed31, Sme). v31. Scaffold annotation: Dresden transcriptome & transdecoder.
   16 - Nematostella vectensis (nemVec1, Nve). ASM20922v1/GCA_000209225.1.  Scaffold annotation: Ensembl Metazoa v36.
   17 - Arabidopsis thaliana (araTha10, Ath). TAIR10.  Scaffold annotation: Ensembl  Plants v31.

")
    db <- readLines(file("stdin"),1)
    db <- as.character(db)
    close(file("stdin"))
    
    dbs<-as.numeric(strsplit(db,"\\s*,\\s*",perl=TRUE)[[1]])
    for (db in dbs){
      if (db<1 || db>length(vastdbFiles)) {stop("Input ",db," is out of range.")}
      Url <- paste("http://vastdb.crg.eu/libs/", vastdbFiles[db], sep = "")
      downloadDb(Url, vastdbFiles[db])
    }
  }

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
