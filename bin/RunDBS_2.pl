#!/usr/bin/env perl
# This pipeline takes PSI templates and adds PSIs from new samples.
use warnings;
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use File::Copy 'move';

# INITIALIZE
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

my $sp;              #species Hsa no longer default
my $dbDir;

my $verboseFlag = 1;
my $helpFlag = 0;

my $globalLen = 50;  # testing? not file specific any longer --TSW

my $outDir;
my $compress = 0;

my $noIRflag = 0;    # don't use IR!
my $IR_version = 2;  # either 1 or 2

my $cRPKMCounts = 0; # print a second cRPKM summary file containing read counts

my $asmbly;       # for human and mouse: vts formats the output wrt. hg19/hg3, mm9/mm10 depending on user's choice of argument -a
 
GetOptions("help"  	 => \$helpFlag,
	   "dbDir=s"     => \$dbDir,
	   "sp=s"        => \$sp,
	   "a=s"         => \$asmbly,
	   "verbose"     => \$verboseFlag,
	   "output=s"    => \$outDir,
	   "o=s"         => \$outDir,
           "z"           => \$compress,
	   "noIR"        => \$noIRflag,
	   "IR_version=i" => \$IR_version,
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

if ($helpFlag or (!defined $sp)){
    print STDERR "
Usage: vast-tools combine -o OUTPUTDIR -sp [Hsa|Mmu|etc] [options]

Combine multiple samples analyzed using \"vast-tools align\" into a single summary tables. 

OPTIONS:
	-o, --output 		Output directory to combine samples from (default vast_out)
	-sp Hsa/Mmu/etc		Species selection (mandatory)
	-a			Genome assembly of the output coordinates (only for -sp Hsa or Mmu) 
				For -sp Hsa: hg19 or hg38, (default hg19)
				    - vast-tools works internally with hg19; 
                                      if you choose hg38, the output gets lifted-over to hg38
				For -sp Mmu: mm9 or mm10, (default mm9)
				    - vast-tools will work internally with mm9; 
                                      if you choose mm10, the output gets lifted-over to mm10
	--noIR			Don't run intron retention pipeline (default off)
        --IR_version 1/2        Version of the IR analysis (default 2)
	--dbDir DBDIR		Database directory
	-z			Compress all output files using gzip
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

# if species is not human nor mouse, we override $asmbly ignoring potential user input
if($sp ne "Hsa" && $sp ne "Mmu"){$asmbly="";}
# get assembly specification for human and mouse
if( $sp eq "Hsa" ){if(!defined($asmbly)){$asmbly="hg19";}; unless($asmbly =~ /(hg19|hg38)/){errPrintDie "Specified assmbly $asmbly either unknown or inapplicable for species $sp\n"}}
if( $sp eq "Mmu" ){if(!defined($asmbly)){$asmbly="mm9";};  unless($asmbly =~ /(mm9|mm10)/){errPrintDie "Specified assmbly $asmbly either unknown or inapplicable for species $sp\n"}}
# we add leading "-" for convenience during defining output file name later 
if($asmbly ne ""){$asmbly="-".$asmbly;}


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
    
    # To define version [02/10/15]; minimize changes for users
    # $v => "" or "_v2" [v1/v2]
    my $v;
    my @irFiles;
    if ($IR_version == 1){
	$v="";
	@irFiles = glob(abs_path("to_combine") . "/*.IR");
    }
    elsif ($IR_version == 2){
	$v="_v2";
	@irFiles = glob(abs_path("to_combine") . "/*.IR2");
    }
    
    $noIRflag = 1 if @irFiles == 0;

    unless($noIRflag) {
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
    
    my $finalOutput = "INCLUSION_LEVELS_FULL-$sp$N$asmbly.tab";
    sysErrMsg "cat @input | $binPath/Add_to_FULL.pl -sp=$sp -dbDir=$dbDir " .
	"-len=$globalLen -verbose=$verboseFlag > $finalOutput";
    
    # lift-over if necessary (hg19->hg38 or mm9->mm10)
    if( $asmbly=~/(hg38|mm10)/ ){
    	# select liftOvr dictionary
    	my $dictionary="lftOvr_dict_from_hg19_to_hg38.pdat"; if($asmbly=~/mm10/){$dictionary="lftOvr_dict_from_mm09_to_mm10.pdat";}
    	# do liftOvr
    	sysErrMsg "$binPath/LftOvr_INCLUSION_LEVELS_FULL.pl translate $finalOutput $dbDir/FILES/$dictionary ${finalOutput}.lifted";
    	# move files
    	move("${finalOutput}.lifted","${finalOutput}");
    }
    
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
