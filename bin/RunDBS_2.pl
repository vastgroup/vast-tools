#!/usr/bin/perl -w
# This pipeline takes PSI templates and adds PSIs from new samples.
use Cwd qw(abs_path);
use Getopt::Long;

# INITIALIZE
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

my $sp = "Hsa"; #species Hsa by default
my $dbDir;

my $verboseFlag = 1;
my $helpFlag = 0;

my $outdir;

GetOptions("help" => \$helpFlag, 
			  "dbDir=s" => \$dbDir,
			  "sp=s" => \$sp,
			  "verbose" => \$verboseFlag,
			  "outdir=s" => \$outdir);

if(!defined($dbDir)) {
  $dbDir = "$binPath/../$sp";
}
$dbDir = abs_path($dbDir);

chdir($outDir);

our $EXIT_STATUS = 0;

sub sysErrMsg {
  my $sysCommand = shift;
  not system($sysCommand) or die "[vast align error]: $sysCommand Failed in $0!";
}

sub errPrint {
  my $errMsg = shift;
  print STDERR "[vast align error]: $errMsg\n";
  $EXIT_STATUS = 1;
}

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast align]: $verbMsg\n";
  }
}

if ($helpFlag){
    print "\nCommand:\nRunDBS_2.pl Species (Hsa or Mmu)\n\n";
    print ">> It will add all the samples that are in the running folder (=SAMPLES)\n";
    print ">>>> Root names for each sample must be identical (they will if coming from RunDBS_1)\n";
    die ">> Not yet for IR\n\n";
}

mkdir("raw_incl") unless (-e "raw_incl"); # make new output directories.  --TSW
mkdir("raw_reads") unless (-e "raw_reads"); # ^

### Settings:
#$sp=$ARGV[0];
die "Needs species 3-letter key\n" if !defined($sp);  #ok for now, needs to be better. --TSW

@files=glob("spli_out/*exskX"); #gathers all exskX files (a priori, simple).
$N=$#files+1;

### Gets the PSIs for the events in the a posteriori pipeline
verbPrint "\nBuilding Table for COMBI (a posteriori pipeline)\n";
sysErrMsg "$binPath/Add_to_COMBI.pl -sp=$sp -dbDir=$dbDir";

### Gets the PSIs for the a priori, SIMPLE
verbPrint "\nBuilding Table for EXSK (a priori pipeline, single)\n";
sysErrMsg "$binPath/Add_to_APR.pl -sp=$sp -type=exskX -dbDir=$dbDir";

### Gets the PSIs for the a priori, COMPLEX
verbPrint "\nBuilding Table for MULTI (a priori pipeline, multiexon)\n";
sysErrMsg "$binPath/Add_to_APR.pl -sp=$sp -type=MULTI3X -dbDir=$dbDir";

### Gets the PSIs for the MIC pipeline
verbPrint "\nBuilding Table for MIC (microexons)\n";
sysErrMsg "$binPath/Add_to_MIC.pl -sp=$sp -dbDir=$dbDir";

### Adds those PSIs to the full database of PSIs (MERGE3m).
verbPrint "\nBuilding non-redundant PSI table (MERGE3m)\n";  # FIXED?? --TSW
sysErrMsg "$binPath/Add_to_MERGE3m.pl raw_incl/INCLUSION_LEVELS_EXSK-$sp$N-n.tab raw_incl/INCLUSION_LEVELS_MULTI-$sp$N-n.tab raw_incl/INCLUSION_LEVELS_COMBI-$sp$N-n.tab raw_incl/INCLUSION_LEVELS_MIC-$sp$N-n.tab";

### Gets PSIs for ALT5ss and adds them to the general database
verbPrint "\nBuilding Table for Alternative 5'ss choice events\n";
sysErrMsg "$binPath/Add_to_ALT5.pl -sp=$sp -dbDir=$dbDir";

### Gets PSIs for ALT3ss and adds them to the general database
verbPrint "\nBuilding Table for Alternative 3'ss choice events\n";
sysErrMsg "$binPath/Add_to_ALT3.pl -sp=$sp -dbDir=$dbDir";

