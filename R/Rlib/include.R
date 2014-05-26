#!/usr/bin/Rscript
#
# Author: Tim Sterne-Weiler, 2014
# tim.sterne.weiler@utoronto.ca
#
#  Install/loading of required Packages.


loadPackages <- function(toLoad, local.lib="Rlib/") {
  for(i in toLoad) {
    writeLines(sprintf("Trying to load required package: %s", toString(i)), stderr())
    if(!suppressWarnings(require(toString(i), character.only=T, quietly=T)) &&
		 !suppressWarnings(require(toString(i), character.only=T, lib.loc=local.lib, quietly=T))){
        writeLines(sprintf("%s did not load correctly! Now trying to install..", i), stderr()) 
        install.packages(i, , repos='http://cran.us.r-project.org', lib=local.lib)
        if(require(toString(i), character.only=T, lib.loc=local.lib)) {
          writeLines(sprintf("%s has been installed locally in R/Rlib!", i), stderr())
        } else {
          stop(sprintf("quitting!!! I could not install %s for you!", i)) 
		  }
        
    }
  }

}

#######
# Colorblind Palette!

cbb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
