#!/usr/bin/Rscript

##  Install/loading of required Packages.
toLoad <- c("getopt", "RColorBrewer", "reshape2", "ggplot2", "grid")

for(i in toLoad) {
  if(!require(toString(i), character.only=T)){
    print(sprintf("%s did not load correctly! Now trying to install..", i)) 
    install.packages(i, , repos='http://cran.us.r-project.org')
    if(require(toString(i), character.only=T)){
      print(sprintf("%s has been installed and loaded by vastdiff.R!", i)) 
    } else {
      stop(sprintf("vastdiff.R quitting!!! I could not install %s for you!", i)) 
    }   
  }
}
##

multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # Make the panel
    plotCols = cols                       # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1 
        if(is.null(plots[[i]])) { next }
        print(plots[[i]], vp = vplayout(curRow, curCol ))  
    }   

}
#######

cbb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
