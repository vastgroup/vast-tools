#!/usr/bin/Rscript

### Generate PIR (percent intron retention) tables from normalized count files containing
### EIJ1, EIJ2, EEJ, I read counts. All such files in the given collection will be merged.
### Also required is a template containing the introns to be merged to.
### Which samples are combined depends on which files *.IR are found in countDir
###
### Returns a single table with two columns for each sample:
###   - PIR: filtered applying criteria
###   - Q:   quality, format: coverage,balance-p.val@alpha,beta where alpha and beta are reassigned pseudocounts for ret/const
###
### U. Braunschweig, The Donnelly Centre, University of Toronto - 05/2014


### FIXED --TSW
argv <- commandArgs(trailingOnly = F)
scriptPath <- dirname(sub("--file=","",argv[grep("--file",argv)]))

# Source Rlib. --TSW
source(paste(c(scriptPath,"/../R/Rlib/include.R"), collapse=""))

loadPackages("optparse")


### IS THIS CORRECT??
argv <- commandArgs(trailingOnly = TRUE)
if (any(grepl("^-s$", argv)) & any(grepl("--species=", argv))) stop("Species/collection specified multiple times")
if (any(grepl("^-s$", argv))) {
    scriptPath <- dirname(argv[grep("^-s$", argv) + 1])
} else {
    if (!any(grepl("--species=", argv))) stop("Species/collection is required")
    scriptPath <- dirname(sub("--species=","",argv[grep("--species",argv)]))
}


suppressPackageStartupMessages(require("optparse"))
opt.list <- list(
    make_option(c("-s", "--species"),  action="store",
                help="Path of the vastdb branch that contains the current analysis, e.g. ~/vastdb/Hsa"),
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
opt <- parse_args(OptionParser(option_list=opt.list))

## Check input
if (!file.exists(opt$species)) stop("Species/collection ", opt$species, " not found")
opt$species <- paste(sub("/$?", "", opt$species), "/", sep="")  # make sure of trailing /
dbDir <- paste(dirname(opt$species), "/", sep="")
species <- basename(opt$species)

if (opt$countDir == opt.list[[2]]@default) {opt$countDir <- paste(opt$species, opt$countDir, sep="")}
if (opt$outDir   == opt.list[[3]]@default) {opt$outDir   <- paste(opt$species, opt$outDir, sep="")}
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
samples <- data.frame(Sample = sub("\\.cReadcount.*", "", sampleFiles),
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
    pir.i[which(cov.i <= opt$COVthresh | bal.i < opt$BALthresh)] <- NA

    ## bring into correct order and populate tables
    datMerge <- data.frame(dat$Event, datInd = 1:nrow(dat))
    datMerge <- merge(data.frame(template$juncID, intronInd=1:nrow(template)), datMerge, by=1, all.x=TRUE)
    datMerge <- datMerge[order(datMerge[,1]),3]

    pir[,2*i - i] <- pir.i[datMerge]
    pir[,2*i]     <- qal.i[datMerge]
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

if (verb) {cat("... done.\n\n")}
