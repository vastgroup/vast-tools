# 2013-08-29

preprocess_sample_colors <- function(data, database) {
   # Given a list of sample names from PSI + Quality Score data, generate a
   # corresponding color code sequence and re-order the columns according to a
   # specified order (see below).
   # 
   # Color codes and sample ordering are taken from a "master" samples database
   # in the following format:
   # Order    SampleName    GroupName    RColorCode
   # 1        Ooctye        EarlyDev     blue
   # 2        Embr_2C       EarlyDev     red
   # etc..
   #
   # RColorCode	: Any of the three kinds of R color specifications:
   #  1) color name (as specified by colors())
   #  2) hex color code (#rrggbb)
   #
   # Args:
   #    data: a n x m data frame of PSI values and quality scores where n is 
   #            the number of AS events and m/2 is the number of samples
   #    database: filename of the samples database in tab-delimited format
   #
   # Returns:
   #    A list containing the re-ordered data.frame of PSI ("data"), re-ordered 
   #    data.frame of quality scores ("qual"), and a sequence of corresponding 
   #    color codes ("col")

   R <- list()

   if (is.null(database)) {
     mycols <- rep("black", ncol(data) / 2)
     R <- list(data=data[, seq(1, ncol(data), 2)], 
               qual=data[, seq(2, ncol(data), 2)], 
               col=mycols, group.index=NULL, group.col=NULL)
   } else {
     db <- read.table(database, header = T, sep="\t", comment.char="")
     
     # check input file
     if (ncol(db) < 4) {
       stop("The tissues database file is not formatted correctly. Please
          double check")
     }
     
     # check if all samples in input data is in the database
     #unk.samples <- colnames(data) %in% db$SampleName
     #if (!all(unk.samples)) {
     #s <- colnames(data)[!unk.samples]
     #stop(paste("The following samples are not in the tissues database:", 
     #paste(s, collapse=", ")))
     #}
     
     # keep only tissue groups that are present in input data
     # (to take into account samples that might have been excluded)
     db <- db[db$SampleName %in% colnames(data),]
     
     # Re-order the PSI table
     db <- db[order(db$Order),]
     db$Order <- 1:nrow(db)
     new.column.idx <- sapply(db$SampleName, 
                              function(x) which(colnames(data) == x))
     
     data.new <- data[,new.column.idx]
     
     # Generate a corresponding color code sequence
     mycols <- db$RColorCode
     names(mycols) <- db$SampleName
     
     # Store indices of columns for each group
     groups <- unique(db$GroupName)
     mygroups <- list()
     mygroupcol <- rep(NA, length(groups))
     for (i in 1:length(groups)) {
       mygroups[[i]] <- which(colnames(data.new) %in%  
                                db[db$GroupName == groups[i],"SampleName"])
       mygroupcol[i] <- as.character(db[db$GroupName == groups[i], 
                                        "RColorCode"][1])
     }
     names(mygroups) <- groups
     names(mygroupcol) <- groups
     
     qual.new <- data[,new.column.idx + 1]
     R <- list(data=data.new, 
               qual=qual.new, 
               col=mycols, 
               group.index=mygroups, group.col=mygroupcol)
   }
   
   return(R)
}
