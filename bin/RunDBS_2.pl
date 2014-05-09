#!/usr/bin/perl -w
# This pipeline takes PSI templates and adds PSIs from new samples.
use strict;
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

my $globalLen = 50; # testing? not file specific any longer --TSW

my $outDir;

GetOptions("help" => \$helpFlag, 
			  "dbDir=s" => \$dbDir,
			  "sp=s" => \$sp,
			  "verbose" => \$verboseFlag,
			  "output=s" => \$outDir,
			  "o=s" => \$outDir);

if(!defined($dbDir)) {
  $dbDir = "$binPath/../$sp";
}
$dbDir = abs_path($dbDir);

chdir($outDir);

our $EXIT_STATUS = 0;

sub sysErrMsg {
  my $sysCommand = shift;
  not system($sysCommand) or die "[vast combine error]: $sysCommand Failed in $0!";
}

sub errPrint {
  my $errMsg = shift;
  print STDERR "[vast combine error]: $errMsg\n";
  $EXIT_STATUS = 1;
}

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine]: $verbMsg\n";
  }
}

if ($helpFlag){
    errPrint "Usage:

vast-tools combine -o OUTPUTDIR [options]

OPTIONS:
	-o OUTPUTDIR, --output OUTPUTDIR	:	output directory to combine samples from... [default vast_out]
	-dbDir DBDIR				:	Database directory
	-sp Hsa/Mmu				:	Species selection
	-v, --verbose				:	Verbose messages
	-h, --help				:	Print this message
";
  exit $EXIT_STATUS;
}

mkdir("raw_incl") unless (-e "raw_incl"); # make new output directories.  --TSW
mkdir("raw_reads") unless (-e "raw_reads"); # ^

### Settings:
#$sp=$ARGV[0];
die "Needs species 3-letter key\n" if !defined($sp);  #ok for now, needs to be better. --TSW

my @files=glob("spli_out/*exskX"); #gathers all exskX files (a priori, simple).
my $N=$#files+1;

### Gets the PSIs for the events in the a posteriori pipeline
verbPrint "Building Table for COMBI (a posteriori pipeline)\n";
my $combi = "$binPath/Add_to_COMBI.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

### Gets the PSIs for the a priori, SIMPLE
verbPrint "Building Table for EXSK (a priori pipeline, single)\n";
my $exsk = "$binPath/Add_to_APR.pl -sp=$sp -type=exskX -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

### Gets the PSIs for the a priori, COMPLEX
verbPrint "Building Table for MULTI (a priori pipeline, multiexon)\n";
my $multi = "$binPath/Add_to_APR.pl -sp=$sp -type=MULTI3X -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

### Gets the PSIs for the MIC pipeline
verbPrint "Building Table for MIC (microexons)\n";
my $mic = "$binPath/Add_to_MIC.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

### Chain the following commands together:
###     Get PSIs for the events in the a posteriori pipeline
###     Gets the PSIs for the a priori, SIMPLE
###     Gets the PSIs for the a priori, COMPLEX
###     Gets the PSIs for the MIC pipeline
my $cmd =  "\($combi\; $exsk\; $multi\; $mic\)";

### Adds those PSIs to the full database of PSIs (MERGE3m).
verbPrint "Building non-redundant PSI table (MERGE3m)\n";  # FIXED?? --TSW
my $merge = "$binPath/Add_to_MERGE3m.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

sysErrMsg "$cmd | $merge";

### Gets PSIs for ALT5ss and adds them to the general database
verbPrint "Building Table for Alternative 5'ss choice events\n";
sysErrMsg "$binPath/Add_to_ALT5.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

### Gets PSIs for ALT3ss and adds them to the general database
verbPrint "Building Table for Alternative 3'ss choice events\n";
sysErrMsg "$binPath/Add_to_ALT3.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";
