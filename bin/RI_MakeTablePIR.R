#!/usr/bin/Rscript

### Generate PIR (percent intron retention) tables from normalized count files containing
### EIJ1, EIJ2, EEJ, I read counts and a table giving sample names and files.
### Also required is a template containing the introns to be merged to.
### Note that thresholds on junction coverage/balance are currently hard-coded!
### Which samples are combined depends on which files *.cReadcount* are found in countDir
###
### Arguments:
###   sp           Absolute Path of the relevant species directory in vastDB, e.g. .../AS_PIPE_S/Hsa
###   countDir     (optional) Directory of read count tables; default: "RAW_READS/" in "sp"
###   outDir       (optional) Directory for output; default: "spli_out/" in "sp"
###   rmHigh       (optional) TRUE/1 or FALSE/0; should introns that have a PIR > 95 in all non-NA samples be 
###                set to NA? Default: TRUE
###   verb         (optional) TRUE/1 or FALSE/0; Print status messages; default: FALSE
###
### Returns a single table with two columns for each sample:
###   - PIR: filtered applying criteria
###   - Q:   quality, format: coverage,balance-p.val@alpha,beta where alpha and beta are reassigned pseudocounts for ret/const
###
### U. Braunschweig, The Donnelly Centre, University of Toronto - 05/2014
###
### Invoke like this from perl:
###   system "Rscript ~/AS_PIPE_S_dev/bin/RI_MakeTablePIR.R \"sp='/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/Hsa/'\" \"rmHigh=0\" \"verb=1\""


cArgs <- commandArgs(TRUE)
if (length(cArgs) < 1) {
    cat("Usage: Rscript RI.MakeTablePIR.R "sp=<dir>" ["countDir=<dir>" "outDir=<dir>" "rmHigh=0/1" "verb=0/1"]\n")
    stop("sp (species) not specified")
} 

for (i in 1:length(cArgs)) {
    eval(parse(text=cArgs[[i]]))
}

## Check input
if (!exists("sp") || !file.exists(sp)) {
    stop("sp (species) not found")
} else {
    sp <- paste(sub("/$?", "", sp), "/", sep="")  # make sure of trailing /
    dbDir <- sub("[^/]+$", "", sp)
    species <- sub("(.*/)?([^/]+)/$", "\\2", sp)
}
if (!exists("countDir"))     {countDir <- paste(sp, "RAW_READS/", sep="")}
if (!exists("outDir"))       {outDir   <- paste(sp, "spli_out/", sep="")}
if (!file.exists(outDir))    {dir.create(outDir, recursive=TRUE)}
if (!exists("rmHigh"))       {rmHigh   <- TRUE}
if (!exists("verb"))         {verb     <- FALSE}
templFile <- paste(sp, "TEMPLATES/", species, ".IR.Template.txt", sep="")
if (!file.exists(templFile)) {stop("Template file ", templFile, " not found")}
countDir  <- paste(sub("/$?", "", countDir), "/", sep="")
outDir    <- paste(sub("/$?", "", outDir), "/", sep="")
rmHigh    <- as.logical(rmHigh)
verb      <- as.logical(verb)



## Check which samples are there and load template
sampleFiles <- dir(countDir, pattern="*cReadcount")
if (length(sampleFiles) == 0) {
    stop("No IR samples found in ", countDir)
} else {
    if (verb) {cat("Merging IR of", length(sampleFiles), "samples...\n")}
}
samples <- data.frame(Sample = sub("\\.cReadcount.*", "", sampleFiles),
                      File   = sampleFiles,
                      stringsAsFactors=FALSE)
template <- read.delim(templFile)


## Hard-coded parameters
thresh.cov <- 10    # Threshold for coverage (# of reads)
thresh.bal <- 0.05  # Threshold for balance p.value of binom.test()
thresh.PIR <- 95    # Threshold for removing all-high introns (minimum PIR)


## prepare tables for coverage, balance and PIR
pir <- as.data.frame(lapply(1:(2*nrow(samples)), FUN=function(x) {rep(NA, nrow(template))}))
names(pir) <- paste(rep(samples$Sample, each=2), c("","-Q"), sep="")
#cov <- as.data.frame(lapply(1:nrow(samples), FUN=function(x) {rep(NA, nrow(template))}))
#names(cov) <- samples$Sample
#bal <- cov
#qal <- cov

for (i in 1:nrow(samples)) {
    ## Read data for one sample
    if (verb) {cat(samples$Sample[i], "\n")}
    dat <- read.delim(paste(countDir, samples$File[i], sep=""))
    if (names(dat)[1] != "Event") {
        dat <- read.delim(paste(countDir, samples$File[i], sep=""), header=F)
        names(dat) <- c("Event","EIJ1","EIJ2","EEJ","I")
    }
    dat <- dat[dat$Event %in% template$juncID,]
    dat <- dat[!(dat$Event %in% dat$Event[duplicated(dat$Event)]),]  # should never be present

    ## calculate PIR, coverage, balance
    pir.i <- 100 * (dat[,2] + dat[,3]) / (dat[,2] + dat[,3] + 2 * dat[,4])
    cov.i <- dat[,4] + apply(dat[,c(2,3,5)], MAR=1, FUN=median)
    tot.i <- dat[,2] + dat[,3] + dat[,4]
    alpha.i <- pir.i / 100 * tot.i
    beta.i  <- (1 - pir.i / 100) * tot.i
    
    xranges <- t(apply(dat[,c(2,3,5)], MAR=1, FUN=range))
    xranges <- apply(xranges, MAR=2, round)
    xranges[,2] <- xranges[,1] + xranges[,2]
    bal.i <- numeric(length=nrow(xranges))
    bal.i[xranges[,2] == 0] <- 1
    bal.i[xranges[,2] > 0] <- apply(xranges[xranges[,2] > 0,], MAR=1, FUN=function(x) {
        binom.test(x=x[1], n=x[2], p=1/3.5, alternative="less")$p.value         
    })

    ## make the 'quality' column: cov,bal@alpha,beta
    qal.i <- paste(round(cov.i, 1), ",", signif(bal.i, 3), "@", alpha.i, ",", beta.i, sep="")

    ## Remove values that do not pass criteria
    pir.i[which(cov.i <= thresh.cov | bal.i < thresh.bal)] <- NA

    ## bring into correct order and populate tables
    datMerge <- data.frame(dat$Event, datInd = 1:nrow(dat))
    datMerge <- merge(data.frame(template$juncID, intronInd=1:nrow(template)), datMerge, by=1, all.x=TRUE)
    datMerge <- datMerge[order(datMerge[,1]),3]

    pir[,2*i - i] <- pir.i[datMerge]
    pir[,2*i]     <- qal.i[datMerge]
#    cov[,i]       <- cov.i[datMerge]
#    bal[,i]       <- bal.i[datMerge]
}

## Remove values that do not pass criteria
#cleanpir <- pir
#for (i in 1:nrow(samples)) {cleanpir[which(cov[,i] <= thresh.cov | bal[,i] < thresh.bal), i] <- NA}
if (rmHigh) {
    minpir <- suppressWarnings(apply(as.data.frame(pir[,seq(from=1, by=2, length.out=nrow(samples))]), MAR=1, FUN=min, na.rm=T))
    pir[minpir > thresh.PIR, seq(from=1, by=2, length.out=nrow(samples))] <- NA
}


## Save table
pir <- data.frame(template[,1:6], pir)

write.table(pir, file=paste(outDir, "INCLUSION_LEVELS_IR-", species, nrow(samples), ".tab", sep=""),
            row.names=F, col.names=T, quote=F, sep='\t')

if (verb) {cat("... done.\n\n")}
