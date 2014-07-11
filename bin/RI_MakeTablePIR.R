#!/usr/bin/env Rscript

### Generate PIR (percent intron retention) tables from normalized count files containing
### EIJ1, EIJ2, EEJ, I read counts. All such files in the given collection will be merged.
### Also required is a template containing the introns to be merged to.
### Which samples are combined depends on which files *.IR are found in countDir
###
### Returns a single table with two columns for each sample:
###   - PIR: filtered applying criteria
###   - Q:   quality, format: coverage,balance-p.val@alpha,beta where alpha and beta are reassigned pseudocounts for ret/const
###
### U. Braunschweig, The Donnelly Centre, University of Toronto - 06/2014


argv <- commandArgs(trailingOnly = FALSE)
scriptPath <- dirname(sub("--file=","",argv[grep("--file",argv)]))

# Source Rlib. --TSW
source(paste(c(scriptPath,"/../R/Rlib/include.R"), collapse=""))
loadPackages(c("optparse"), local.lib=paste(c(scriptPath,"/../R/Rlib"), collapse=""))

argv <- commandArgs(trailingOnly = TRUE)

opt.list <- list(
    make_option(c("-s", "--species"),  action="store",
                help="Path of the vastdb branch that contains the current analysis, e.g. ~/vastdb/Hsa"),
    make_option(c("-q", "--quality"),  action="store",  # legacy qualities --UB
                help="Path of the file containing legacy quality scores for IR [no default]"),
    make_option(c("-c", "--countDir"), action="store",
                default="to_combine",
                help="Location of raw count tables [default: <species>/%default]"),
    make_option(c("-o", "--outDir"),   action="store",
                default="raw_incl",
                help="Location of output [default: <species>/%default]"),
    make_option(c("-r", "--rmHigh"),   action="store_true", default = FALSE,
                help="Remove values of events that are always above threshold?  [default: %default]"),
    make_option(c("-v", "--verbose"),  action="store", type = "integer",
                default = 1,
                help="Print status messages (0 - no/1 - yes)? [default: %default]"),
    make_option(c("-P", "--PIRthresh"),  action="store", default=95, type="numeric",
                help="Threshold for removal of events that are not below it in any sample [default: %default]"),
    make_option(c("-C", "--COVthresh"),  action="store", default=10, type="numeric",
                help="Threshold for junction read count [default: %default]"),
    make_option(c("-B", "--BALthresh"),  action="store", default=0.05, type="numeric",
                help="Threshold for p-value of balance binomial test [default: %default]")
    )
opt <- parse_args(OptionParser(option_list=opt.list), args=argv)


## Check input
if (is.null(opt$species))       stop("Species/collection is required")
if (!file.exists(opt$species))  stop("Species/collection ", opt$species, " not found")
opt$species <- paste(sub("/$?", "", opt$species), "/", sep="")  # make sure of trailing /
dbDir <- paste(dirname(opt$species), "/", sep="")
species <- basename(opt$species)

if (is.null(opt$quality))       stop("Quality file is required")  # Legacy qualities --UB
if (!file.exists(opt$quality))  stop("Quality file ", opt$quality, " not found")

if (opt$countDir == opt.list[[3]]@default) {opt$countDir <- paste(opt$species, opt$countDir, sep="")}
if (opt$outDir   == opt.list[[4]]@default) {opt$outDir   <- paste(opt$species, opt$outDir, sep="")}
if (!file.exists(opt$outDir))   {dir.create(opt$outDir, recursive=TRUE)}
countDir  <- paste(sub("/$?", "", opt$countDir), "/", sep="")
outDir    <- paste(sub("/$?", "", opt$outDir), "/", sep="")

templFile <- paste(opt$species, "TEMPLATES/", species, ".IR.Template.txt", sep="")
if (!file.exists(templFile))    {stop("Template file ", templFile, " not found")}
rmHigh    <- as.logical(opt$rmHigh)
verb      <- as.logical(opt$verbose)

if (is.na(rmHigh))                               {stop("Invalid value for rmHigh")}
if (is.na(opt$PIRthresh) || opt$PIRthresh > 100) {stop("Invalid value for PIRthresh")}
if (is.na(opt$COVthresh) || opt$COVthresh < 0)   {stop("Invalid value for COVthresh")}
if (is.na(opt$BALthresh) || opt$BALthresh > 1)   {stop("Invalid value for BALthresh")}


## Check which samples are there and load template
sampleFiles <- dir(countDir, pattern="*\\.IR$")
if (is.na(sampleFiles[1])) {
    stop("No IR samples found in ", countDir)
} else {
    if (verb) {cat("Merging IR of", length(sampleFiles), "sample(s)...\n")}
}
samples <- data.frame(Sample = sub("\\.IR", "", sampleFiles),
                      File   = sampleFiles,
                      stringsAsFactors=FALSE)
template <- read.delim(templFile)


## prepare tables for coverage, balance and PIR
pir <- as.data.frame(lapply(1:(2*nrow(samples)), FUN=function(x) {rep(NA, nrow(template))}))
names(pir) <- paste(rep(samples$Sample, each=2), c("","-Q"), sep="")

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

    ## check that filtered input data is not empty
    if (nrow(dat) == 0) {
        stop("Filtered input data contains 0 rows. Was the correct species being specified?")
    }

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

    ## Remove values that do not pass criteria  # not done any more --UB
    #pir.i[which(cov.i <= opt$COVthresh | bal.i < opt$BALthresh)] <- NA

    ## bring into correct order and populate tables
    datMerge <- data.frame(dat$Event, datInd = 1:nrow(dat))
    datMerge <- merge(data.frame(template$juncID, intronInd=1:nrow(template)), datMerge, by=1, all.x=TRUE)
    datMerge <- datMerge[order(datMerge[,1]),3]

    pir[,2*i - 1] <- pir.i[datMerge]
    pir[,2*i]     <- qal.i[datMerge]
}


## Add legacy quality scores for compatibility with downstream tools
legQual <- read.delim(opt$quality, as.is=TRUE, check.names = FALSE)
if (!all(names(legQual)[-1] %in% samples$Sample)) {stop("Samples in IR and IR quality file do not match")}

legQual <- legQual[legQual$EVENT %in% template$juncID,]
qualMerge <- data.frame(legQual$EVENT, qualInd=1:nrow(legQual))
qualMerge <- merge(data.frame(template$juncID, intronInd=1:nrow(template)), qualMerge, by=1, all.x=TRUE)
legQual <- legQual[qualMerge[order(qualMerge[,1]), 3],]  # reorder the same way as the template

for (i in 1:nrow(samples)) {
    legQual[is.na(legQual[,i + 1]),i + 1] <- "N,N,N,NA,0=0=0" # set qual scores when missing due to no reads
    pir[is.na(pir[,i * 2]),i * 2] <- "@NA,NA"                 # set pseudocounts -"-
    pir[,i * 2] <- paste(legQual[,1 + i], sub("[^,]+", "", pir[,i * 2]), sep="")
}

## Remove values from events that are never below PIRthresh
if (rmHigh) {
    minpir <- suppressWarnings(apply(as.data.frame(pir[,seq(from=1, by=2, length.out=nrow(samples))]), MAR=1, FUN=min, na.rm=T))
    pir[minpir > opt$PIRthresh, seq(from=1, by=2, length.out=nrow(samples))] <- NA
}


## Save table
pir <- data.frame(template[,1:6], pir)

write.table(pir, file=paste(outDir, "INCLUSION_LEVELS_IR-", species, nrow(samples), ".tab", sep=""),
            row.names=F, col.names=T, quote=F, sep='\t')

if (verb) {cat("... done.\n")}
