# 2013-08-29

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

   if (is.null(database)) {
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
       #unk.samples <- colnames(psi) %in% db$SampleName
       #if (!all(unk.samples)) {
         #s <- colnames(psi)[!unk.samples]
         #stop(paste("The following samples are not in the tissues database:", 
                    #paste(s, collapse=", ")))
       #}
      
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
