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

my $globalLen = 50; # testing? not file specific any longer --TSW

my $outDir;

GetOptions("help" => \$helpFlag, 
			  "dbDir=s" => \$dbDir,
			  "sp=s" => \$sp,
			  "verbose" => \$verboseFlag,
			  "outdir=s" => \$outDir);

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
verbPrint "Building Table for COMBI (a posteriori pipeline)\n";
sysErrMsg "$binPath/Add_to_COMBI.pl -sp=$sp -dbDir=$dbDir -len=$globalLen";

### Gets the PSIs for the a priori, SIMPLE
verbPrint "Building Table for EXSK (a priori pipeline, single)\n";
sysErrMsg "$binPath/Add_to_APR.pl -sp=$sp -type=exskX -dbDir=$dbDir -len=$globalLen";

### Gets the PSIs for the a priori, COMPLEX
verbPrint "Building Table for MULTI (a priori pipeline, multiexon)\n";
sysErrMsg "$binPath/Add_to_APR.pl -sp=$sp -type=MULTI3X -dbDir=$dbDir -len=$globalLen";

### Gets the PSIs for the MIC pipeline
verbPrint "Building Table for MIC (microexons)\n";
sysErrMsg "$binPath/Add_to_MIC.pl -sp=$sp -dbDir=$dbDir -len=$globalLen";

### Adds those PSIs to the full database of PSIs (MERGE3m).
verbPrint "Building non-redundant PSI table (MERGE3m)\n";  # FIXED?? --TSW
sysErrMsg "$binPath/Add_to_MERGE3m.pl raw_incl/INCLUSION_LEVELS_EXSK-$sp$N-n.tab raw_incl/INCLUSION_LEVELS_MULTI-$sp$N-n.tab raw_incl/INCLUSION_LEVELS_COMBI-$sp$N-n.tab raw_incl/INCLUSION_LEVELS_MIC-$sp$N-n.tab -sp=$sp -dbDir=$dbDir -len=$globalLen";

### Gets PSIs for ALT5ss and adds them to the general database
verbPrint "Building Table for Alternative 5'ss choice events\n";
sysErrMsg "$binPath/Add_to_ALT5.pl -sp=$sp -dbDir=$dbDir -len=$globalLen";

### Gets PSIs for ALT3ss and adds them to the general database
verbPrint "Building Table for Alternative 3'ss choice events\n";
sysErrMsg "$binPath/Add_to_ALT3.pl -sp=$sp -dbDir=$dbDir -len=$globalLen";

