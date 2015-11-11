#!/usr/bin/env perl
# This pipeline takes PSI templates and adds PSIs from new samples.
use warnings;
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
my $compress = 0;

my $noIRflag = 0; #don't use IR!
my $IR_version = 1; # either 1 or 2

my $cRPKMCounts = 0; # print a second cRPKM summary file containing read counts

GetOptions("help"  	 => \$helpFlag,
	   "dbDir=s"     => \$dbDir,
	   "sp=s"        => \$sp,
	   "verbose"     => \$verboseFlag,
	   "output=s"    => \$outDir,
	   "o=s"         => \$outDir,
           "z"           => \$compress,
	   "noIR"        => \$noIRflag,
	   "IR_version"  => \$IR_version,
           "C"           => \$cRPKMCounts);

our $EXIT_STATUS = 0;

sub sysErrMsg {
  my @sysCommand = (shift);
  not system(@sysCommand) or die "[vast combine error]: @sysCommand Failed in $0!";
}

sub errPrint {
  my $errMsg = shift;
  print STDERR "[vast combine error]: $errMsg\n";
  $EXIT_STATUS = 1;
}

sub errPrintDie {
  my $errMsg = shift;
  errPrint $errMsg;
  exit $EXIT_STATUS if ($EXIT_STATUS != 0);
}

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine]: $verbMsg\n";
  }
}

if ($helpFlag or !defined($ARGV[0])){
    print STDERR "
Usage: vast-tools combine -o OUTPUTDIR [options]

Combine multiple samples analyzed using \"vast-tools align\" into a single summary tables. 

OPTIONS:
	-o, --output 		Output directory to combine samples from (default vast_out)
	--dbDir DBDIR		Database directory
	-sp Hsa/Mmu/etc		Species selection
	-z			Compress all output files using gzip
	--noIR			Don't run intron retention pipeline (default off)
        --IR_version 1/2        Version of the IR analysis (default 1)
	-v, --verbose		Verbose messages
	-h, --help		Print this help message
	-C			Create a cRPKM plus read counts summary table. By default, a
    				table containing ONLY cRPKM is produced. This option is only
           			applicable when expression analysis is enabled.

*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)
					\n";

  exit $EXIT_STATUS;
}

errPrintDie "Need output directory" unless (defined $outDir);
errPrintDie "The output directory $outDir does not exist" unless (-e $outDir);
errPrintDie "IR version must be either 1 or 2." if ($IR_version != 1 && $IR_version != 2);

if(!defined($dbDir)) {
  $dbDir = "$binPath/../VASTDB";
}
$dbDir = abs_path($dbDir);
$dbDir .= "/$sp";
errPrintDie "The database directory $dbDir does not exist" unless (-e $dbDir);
verbPrint "Using VASTDB -> $dbDir";

chdir($outDir);

mkdir("raw_incl") unless (-e "raw_incl"); # make new output directories.  --TSW
mkdir("raw_reads") unless (-e "raw_reads"); # ^

### Settings:
errPrintDie "Needs species 3-letter key\n" if !defined($sp);  #ok for now, needs to be better. --TSW

my @files=glob("to_combine/*exskX"); #gathers all exskX files (a priori, simple).
my $N=$#files+1;

if ($N != 0) {
    ### Gets the PSIs for the events in the a posteriori pipeline
    verbPrint "Building Table for COMBI (a posteriori pipeline)\n";
    sysErrMsg "$binPath/Add_to_COMBI.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

    ### Gets the PSIs for the a priori, SIMPLE
    verbPrint "Building Table for EXSK (a priori pipeline, single)\n";
    sysErrMsg "$binPath/Add_to_APR.pl -sp=$sp -type=exskX -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

    ### Gets the PSIs for the a priori, COMPLEX
    verbPrint "Building Table for MULTI (a priori pipeline, multiexon)\n";
    sysErrMsg "$binPath/Add_to_APR.pl -sp=$sp -type=MULTI3X -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

    ### Gets the PSIs for the MIC pipeline
    verbPrint "Building Table for MIC (microexons)\n";
    sysErrMsg "$binPath/Add_to_MIC.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

    #my($verbRFlag) = ($verboseFlag) ? "T" : "F";

    my @irFiles = glob(abs_path("to_combine") . "/*.IR");
    $noIRflag = 1 if @irFiles == 0;

    unless($noIRflag) {
	# To define version [02/10/15]; minimize changes for users
        # $v => "" or "_v2" [v1/v2]
	my $v;
	if ($IR_version == 1){
	    $v="";
	}
	elsif ($IR_version == 2){
	    $v="_v2"; 
	}
	### Gets the PIRs for the Intron Retention pipeline
	verbPrint "Building quality score table for intron retention (version $IR_version)\n";
	sysErrMsg "$binPath/RI_MakeCoverageKey$v.pl -sp $sp -dbDir $dbDir " . abs_path("to_combine");
	verbPrint "Building Table for intron retention (version $IR_version)\n";
	sysErrMsg "$binPath/RI_MakeTablePIR.R --verbose $verboseFlag -s $dbDir --IR_version $IR_version" .
	    " -c " . abs_path("to_combine") .
	    " -q " . abs_path("to_combine") . "/Coverage_key$v-$sp$N.IRQ" .
	    " -o " . abs_path("raw_incl");
    }

    ### Adds those PSIs to the full database of PSIs (MERGE3m).
    # to be deprecated and replaced by Add_to_FULL (see below) --KH
    #verbPrint "Building non-redundant PSI table (MERGE3m)\n";
    #sysErrMsg "$binPath/Add_to_MERGE3m.pl " .
    #"raw_incl/INCLUSION_LEVELS_EXSK-$sp$N-n.tab " .
    #"raw_incl/INCLUSION_LEVELS_MULTI-$sp$N-n.tab " .
    #"raw_incl/INCLUSION_LEVELS_COMBI-$sp$N-n.tab " .
    #"raw_incl/INCLUSION_LEVELS_MIC-$sp$N-n.tab " .
    #"-sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

    ### Gets PSIs for ALT5ss and adds them to the general database
    verbPrint "Building Table for Alternative 5'ss choice events\n";
    sysErrMsg "$binPath/Add_to_ALT5.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

    ### Gets PSIs for ALT3ss and adds them to the general database
    verbPrint "Building Table for Alternative 3'ss choice events\n";
    sysErrMsg "$binPath/Add_to_ALT3.pl -sp=$sp -dbDir=$dbDir -len=$globalLen -verbose=$verboseFlag";

    ### Combine results into unified "FULL" table
    verbPrint "Combining results into a single table\n";
    my @input =    ("raw_incl/INCLUSION_LEVELS_EXSK-$sp$N-n.tab",
                    "raw_incl/INCLUSION_LEVELS_MULTI-$sp$N-n.tab",
                    "raw_incl/INCLUSION_LEVELS_COMBI-$sp$N-n.tab",
                    "raw_incl/INCLUSION_LEVELS_MIC-$sp$N-n.tab",
                    "raw_incl/INCLUSION_LEVELS_ALT3-$sp$N-n.tab",
                    "raw_incl/INCLUSION_LEVELS_ALT5-$sp$N-n.tab");

    unless($noIRflag) {
	push(@input, "raw_incl/INCLUSION_LEVELS_IR-$sp$N.tab");
    }
    
    my $finalOutput = "INCLUSION_LEVELS_FULL-$sp$N.tab";
    sysErrMsg "cat @input | $binPath/Add_to_FULL.pl -sp=$sp -dbDir=$dbDir " .
	"-len=$globalLen -verbose=$verboseFlag > $finalOutput";
    
    verbPrint "Final table saved as: " . abs_path($finalOutput) ."\n";
    
    if ($compress) {
      verbPrint "Compressing files\n";
      sysErrMsg "gzip -v raw_incl/*.tab raw_reads/*.tab $finalOutput";
      $finalOutput .= ".gz";
    }
}

### Combine cRPKM files, if present
my @rpkmFiles=glob("expr_out/*.cRPKM"); 
if (@rpkmFiles > 0) {
    verbPrint "Combining cRPKMs into a single table\n";
    my $cRPKMOutput = "cRPKM-$sp" . @rpkmFiles . ".tab";
    $cRPKMCounts = $cRPKMCounts ? "-C" : "";
    sysErrMsg "$binPath/MakeTableRPKMs.pl -sp=$sp -dbDir=$dbDir $cRPKMCounts";

    if ($compress) {
      verbPrint "Compressing files\n";
      sysErrMsg "gzip -v expr_out/*.cRPKM $cRPKMOutput";
      $cRPKMOutput .= ".gz";
    }

    verbPrint "Final cRPKM table saved as: " . abs_path($cRPKMOutput) . "\n";
}

if ($N + @rpkmFiles == 0) {
    verbPrint "Could not find any files to combine. If they are compressed, please decompress them first.\n";
}

verbPrint "Completed " . localtime;
