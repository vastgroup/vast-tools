#!/usr/bin/env perl

# Authors:
# Original Draft: Manuel Irimia, 2011-2014 
# 			      	mirimia@gmail.com
# Reworked: Tim Sterne-Weiler & Kevin Ha, 2014
# 				tim.sterne.weiler@utoronto.ca & k.ha@mail.utoronto.ca 
# Updates and improvements: Andre Gohr & Manuel Irimia, 2015-present
#                               andre.gohr@crg.eu & mirimia@gmail.com

use warnings;
use strict;
use Cwd qw(abs_path cwd);
use Getopt::Long;

# INITIALIZE PATH AND FLAGS--TSW
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

my $helpFlag = 0;
my $bowtie = "bowtie"; # by default;
my $species = "Hsa"; # by default;
my $dbDir; #default
my $pairedEnd = 0; # no by default
my $runExprFlag = 0; # no by default
my $onlyExprFlag = 0; # no by default
my $trim;
my $cores = 1; #default
my $readLength = ""; # default.
my $outdir;
my $noIRflag = 0;  # don't run intron retention (for speed..)
my $stringentIRflag = 0; # Run extra genome/eej subtraction step
my $IR_version = 2; # IR version [09/Nov/2015] [new default 01/04/16]
my $minReadNum;

my $legacyFlag = 0;
my $verboseFlag = 1;  # on for debugging 
my $keepFlag = 0;  # delete genome subtracted reads
my $tmpDir;

my $trimOnceFlag = 0;
my $trimStep = 25;
my $fastaOnly = 0; # Use this to trim to fasta not fastq;

my $useGenSub = 0;

my $trimLen; # This is an undocumented flag for Trim.pl (allows 48bp use)
my $bowtieV = 2; # This undocumented option for # of allowed mismatches..

my $trimmed = 0; # use pre-trimmed read set by Trim.pl
my $keep_trimmed= 0; # to keep the original file when pre-trimmed

my $ribofoot = 0; # flag for ribosome footprinting libraries

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(		  "bowtieProg=s" => \$bowtie,
			  "sp=s" => \$species,
			  "dbDir=s" => \$dbDir,
			  "c=i" => \$cores, 
			  "cores=i" => \$cores,
			  "expr" => \$runExprFlag,
			  "exprONLY" => \$onlyExprFlag,
			  "trim=s" => \$trim,
			  "help" => \$helpFlag,
			  "h" => \$helpFlag,
			  "legacy" => \$legacyFlag,
			  "verbose" => \$verboseFlag,
			  "v" => \$verboseFlag,
              		  "output=s" => \$outdir,
			  "o=s" => \$outdir,
			  "noIR" => \$noIRflag,
			  "stringentIR" => \$stringentIRflag,
			  "IR_version=i" => \$IR_version, 
			  "keep" => \$keepFlag,
			  "minReadDepth=i" => \$minReadNum, #to do
			  "tmpDir=s" => \$tmpDir,
			  "stepSize=i" => \$trimStep,
			  "trimOnce" => \$trimOnceFlag,
			  "findSubtracted" => \$useGenSub,
                          "trimLen=i" => \$trimLen,
                          "mismatchNum=i" => \$bowtieV,
                          "preTrimmed" => \$trimmed,
                          "useFastq" => \$fastaOnly,
			  "riboFoot" => \$ribofoot
			  );

our $EXIT_STATUS = 0;


sub extractReadLen {  # extracts automatically read length from fastq or fastq.gz file
	my $fn=$_[0]; # extracts the first 5000 reads and if all have the same length, returns this length
	              # if they don't have all the same length, returns -1
	              # From 24/12/16: they can have different lengths, but only those >= 50 are used
	
	my $fh;
	if(isZipped($fn)){
		open ( $fh, "-|", "gzip -dc $fn") or errPrintDie("$!");
	}else{
		open( $fh , $fn ) or errPrintDie("$!");
	}
	
	my $maxN=5000;
	my $c=0;
	my $check=1;
	my $readL;
	my %tally_readL;
	my $total=1;
	my $perc;
	while(<$fh>){if($check){if(substr($_,0,1) ne "@"){errPrintDie("Sequence data must be provided in FASTQ format but file $fn does not look like FASTQ format (first line does not start with @).");}$check=0;}
		chomp;
		$c++;
		if($c % 4==2){
			unless($readL){
			    $readL=length($_);
			    $tally_readL{$readL}=1;
			}
			else{ 
#			    if($readL != length($_)){$readL=-1;last;} # not need to 
			    $readL=length($_);
			    $tally_readL{$readL}++ if defined $tally_readL{$readL};
			    $tally_readL{$readL}=1 if !defined $tally_readL{$readL};
			    $total++;
			}
		}
		if($c/4 > $maxN){last;}
	}
	close($fh);

	### to get the most common length
	foreach my $temp_readL (sort {$a<=>$b} keys %tally_readL){
	    $perc = sprintf("%.2f",100*$tally_readL{$temp_readL}/$total);
	    $readL=$temp_readL;
	}

	return($readL,$perc);
}

sub sysErrMsg {
  my @sysCommand = @_;
  not system(@sysCommand) or die "[vast align error]: @sysCommand Failed in $0!";
}

sub errPrint {
  my $errMsg = shift;
  print STDERR "[vast align error]: $errMsg\n";
  $EXIT_STATUS++; 
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
    print STDERR "[vast align]: $verbMsg\n";
  }
}

sub isZipped {
  my $file = shift;
  return $file =~ /\.(gz)$/;
}

sub getPrefixCmd {
  my $file = shift;
  my $prefix = isZipped($file) ? "gzip -dc $file" : "cat $file";
  return $prefix;
}

my $inpType = !$fastaOnly ? "-f" : "-q"; 

# Check database directory
unless(defined($dbDir)) {
  $dbDir = "$binPath/../VASTDB";
}
$dbDir = abs_path($dbDir);
$dbDir .= "/$species";
errPrint "The database directory $dbDir does not exist" unless (-e $dbDir or $helpFlag);

if (!defined($ARGV[0]) or $helpFlag or $EXIT_STATUS){
    print "\nUsage: vast-tools align fastq_file_1 [fastq_file_2] [options]

Align a single RNA-Seq sample to VASTDB genome and junction libraries.
Length of reads must be at least 50 nt; for expression analysis, all reads
must be of same length.

OPTIONS:
	--sp Hsa/Mmu/Gga	Three letter code for the database (default Hsa)
	--dbDir db		Database directory (default VASTDB)
	--cores, -c i		Number of cores to use for bowtie (default 1)
	--output, -o		Output directory (default vast_out)
	--expr			For expression analyses: -expr 
				(PSIs plus cRPKM calculations) (default off)
	--exprONLY		For expression analyses: -exprONLY (only cRPKMs) 
				(default off)
	--bowtieProg prog	Default is to use the bowtie in PATH. Alternatively you can
				supply a specific bowtie program here (default `bowtie`)
	--noIR			Don't run intron retention pipeline 
				(substantially increases speed) (default off)
	--stringentIR		Don't run first filtering step of IR 
				(this will increase speed a little) (default off)
        --IR_version 1/2        Version of the Intron Retention analysis (default 2)
	--keep			Don't remove trimmed and genome-subtracted reads 
				after use. (default off)
	--findSubtracted	Set this flag to start alignment from genome-subtracted
				reads (default off). If enabled, must supply *-e.fq as input
	--trimOnce		Only use first 50bp of reads, if paired, only use 
					50 from fwd and 50 from rev (default off)
	--stepSize i		Trim 50bp every --stepSize (default is 25)
	--preTrimmed		If you are trying to use pre-trimmed fasta/q files 
					(only output from Trim.pl, default off)
	--useFastq		This option is only necessary if you have pre-trimmed reads 
					in fastq not fasta format (default off)
	-h, --help		Print this help message


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";

    exit $EXIT_STATUS;
}

# Command line flags here
if (defined $ARGV[1]) { $pairedEnd = 1; }

# Input sanity checks
errPrintDie "Needs species\n" if !$species;
errPrintDie "Input file " . $ARGV[0] . " does not exist!" if (! -e $ARGV[0]);
errPrintDie "Input file " . $ARGV[1] . " does not exist!" if ($pairedEnd and ! -e $ARGV[1]);
errPrintDie "Invalid number of cores. Must be at least 1." if ($cores !~ /^[1-9]\d*$/);
errPrintDie "Invalid step size." if ($trimStep !~ /^[1-9]\d*$/);
errPrintDie "IR version must be either 1 or 2." if ($IR_version != 1 && $IR_version != 2);

# FOR RIBOFOOT
if($ribofoot) {
  $trimOnceFlag = 1; # only trim once. no slide.
  $runExprFlag = 0; # no need for expression calculations.
  $readLength = 32;
  $trimLen = 32;   
#  $noIRflag = 1;  # temporary;
}

## Getting sample name and length:
my $fq1 = $ARGV[0];
unless(substr($fq1,0,1) eq "/" ){# file path is relative
	$fq1=cwd() . "/$fq1";    #  add to file path current working directory; necessary because later we change the working directory
}
my $fq2;
if($pairedEnd){
	$fq2 = $ARGV[1];
	unless(substr($fq2,0,1) eq "/" ){# file path is relative
		$fq2=cwd() . "/$fq2";    #  add to file path current working directory; necessary because later we change the working directory
	}
}

my $fq;     # takes the fastq file to be processed at each step

my $fileName1 = $fq1;
my $zipped = isZipped($fq1);
my $subtractedFq;

my($root, $length);

$fileName1 =~ s/^.*\///g; # strip path

my $genome_sub = 0;
my $length2="";
my($percF,$percF2);
if ($fileName1 =~ /\-e\.f/){ # it has to be a fastq file (not fasta)
    $genome_sub=1;
    ($root,$length)=$fileName1=~/(\S+?)\-(\d{1,4})\-e\.(fastq|fq|fasta|fa)(\.gz)?/;  #Fixed regex --TSW
    $fq=$&;
    $subtractedFq = $fq1;
    errPrint "Only for 50nt if genome subtracted\n" if $length!=50;
} else {
    if ($runExprFlag || $onlyExprFlag){ # only if GE is actives checks if readLength is provided
	if($ribofoot){  # length is set already in ribofoot mode
	    $length=$readLength;  # set to 32
	}else{
	    ($length,$percF)=extractReadLen($fq1); # it doesn't really matter any more (24/12/16)
	}
	$fileName1 =~ /(\S+)\.(fastq|fq|fasta|fa)(\.gz)?/;  # regex by --TSW
	$root = $1;
    }
    else { # anything is valid here
	($length,$percF)=extractReadLen($fq1);
	$fileName1 =~ /(\S+)\.(fastq|fq|fasta|fa)(\.gz)?/; 
	$root = $1;
    }
    if ($pairedEnd){
	($length2,$percF2)=extractReadLen($fq2);
    }
    $fq = $zipped ? "$root-50.fq.gz" : "$root-50.fq"; #only fastq files are allowed at this point; default trimmed length = 50
}
###

unless($fq2){verbPrint("Input RNA-seq file(s): $fq1");}else{verbPrint("Input RNA-seq file(s): $fq1 and $fq2");}

# if something went wrong with extraction of root of filenames
if($root eq ""){ errPrintDie("Could not extract the base name from the RNA-seq input files, which must look like *.(fastq|fastq.gz|fq|fq.gz|fasta|fasta.gz|fa|fa.gz)");}

unless($fq2){verbPrint("Most common read length detected for fq1: $length ($percF\%)");}
else{verbPrint("Most common read lengths detected for fq1 & fq2: $length ($percF\%) and $length2 ($percF2\%)");}

#if(($onlyExprFlag || $runExprFlag) && $length == -1){ # XXX reads are of variable length
#	verbPrint("Reads are of variable lengths in file $fq1.\nExpression analysis turned off as for this all reads must be of the same length."); 
#	if($onlyExprFlag){exit(1);}
#	$runExprFlag=0;
#}

verbPrint "Using VASTDB -> $dbDir";
# change directories
mkdir($outdir) unless (-e $outdir);
chdir($outdir) or errPrint "Unable to change directories into output" and die;
verbPrint "Setting output directory to $outdir";
mkdir("to_combine") unless (-e "to_combine");
mkdir("expr_out") if (($runExprFlag || $onlyExprFlag) && (! -e "expr_out"));

# set default tmpDir for sort;
verbPrint "Setting tmp directory..";
unless(defined($tmpDir)) {
  mkdir("tmp");
  $tmpDir = abs_path("tmp");  
} else {
  $tmpDir = abs_path($tmpDir);  # or try to find it
  unless(-e $tmpDir) {
    errPrint "$tmpDir does not exist!";
  }
}
unless($EXIT_STATUS > 0) {
  verbPrint "Set tmp directory to $tmpDir!";
}

#length options:
my ($le, $half_length);
my $difLE;

if ($length >= 50){
    $difLE = $length-50;
    $le = 50;
} 
elsif ($ribofoot) {
    $difLE = 0;
    $le = 32;
} 
elsif ($trimLen){
    $difLE = 0;
    $le = 50; # even if trimLen is shorter
} 
else {
    errPrint "Minimum reads length has to be 50nt\n";
}
#####

if ($EXIT_STATUS) {
    exit $EXIT_STATUS;
}

if (!$genome_sub and !$useGenSub){
 my $cmd;
#### Expression analysis (it maps only the first $le nucleotides of the read)
 if ($runExprFlag || $onlyExprFlag){
     verbPrint "Mapping reads against mRNA sequences";

#### Only the first $le nucleotides of the forward read --MI 25/12/14
#     if ($pairedEnd) {
#         $cmd = "$fq1 $fq2";  # altered this to cat both for/rev reads into bowtie.
#     } else {
         $cmd = "$fq1";
#     }

     $cmd = getPrefixCmd($cmd);
#    24/12/16 --MI
#    $cmd .= " | $bowtie -p $cores -m 1 -v $bowtieV -3 $difLE $dbDir/EXPRESSION/mRNA -"; 
     if (defined($trimLen)){
	 $cmd .= " | $binPath/Trim.pl --once --targetLen $trimLen -v | $bowtie -p $cores -m 1 -v $bowtieV $dbDir/EXPRESSION/mRNA -"; 
     }
     else {
	 $cmd .= " | $binPath/Trim.pl --once --targetLen 50 -v | $bowtie -p $cores -m 1 -v $bowtieV $dbDir/EXPRESSION/mRNA -"; 
     }
     
     verbPrint "Calculating cRPKMs\n";
     sysErrMsg "$cmd | $binPath/expr_RPKM.pl - $dbDir/EXPRESSION/$species"."_mRNA-$le.eff expr_out/$root > expr_out/$root\.cRPKM";
 }
 if ($onlyExprFlag){
     print STDERR "Expression analysis done\n";
     exit 0;
 }
}
###

#### Merge PE
# if ($pairedEnd){
#   verbPrint "Concatenating paired end reads";
     #sysErrMsg "cat $fq1 $fq2 > $fq";  # away with this as well? 
                                       # $fq is used in trimming below. but we
                                       # can pipe into it. KH
#   $fq = "$fq1 $fq2";
  #} else {
   $fq = $fq1; # Above is deprecated --TSW 7/14/14
#}

 
#### Trimming
#
#

$keep_trimmed=1 if $trimmed; #keeps the original file provided as pre-trimmed input

my $cmd = getPrefixCmd($fq);

unless($trimmed) {
    my $trimArgs = "--stepSize $trimStep -v"; # before, verbose by default (24/12/16)
    $trimArgs .= " --fasta" if(!$fastaOnly);
    $trimArgs .= " --once" if($trimOnceFlag);
    $trimArgs .= " --targetLen $trimLen" if(defined($trimLen));
    if($pairedEnd) {
	my $pairFq = isZipped($fq2) ? "<( gzip -dc $fq2 )" : $fq2;
	$trimArgs .= " --paired $pairFq";
    } 
    
    verbPrint "Trimming fastq sequences to $le nt sequences";
    ## Add min read depth?
    # Renamed fa/fq --MI [11/11/15]
    if ($fastaOnly){
	sysErrMsg("bash", "-c", "$cmd | $binPath/Trim.pl $trimArgs | gzip -c > $root-$le.fq.gz");
	$fq = "$root-$le.fq.gz"; # set new $fq with trimmed reads --KH
    }
    else { # default behaviour 
	sysErrMsg("bash", "-c", "$cmd | $binPath/Trim.pl $trimArgs | gzip -c > $root-$le.fa.gz");
	$fq = "$root-$le.fa.gz"; # set new $fq with trimmed reads --KH
    }
    $trimmed = 1;
}
####

 
#### Get effective reads (i.e. genome subtraction).
 $subtractedFq = "$root-$le-e.fa.gz" if !$useGenSub;
 $subtractedFq = "$root-$le-e.fq.gz" if $useGenSub;
 unless(-e $subtractedFq and $useGenSub) {
   verbPrint "Doing genome subtraction\n";
   # Force bash shell to support process substitution
   $cmd = getPrefixCmd($fq);
   $cmd .= " | $bowtie -p $cores $inpType -m 1 -v 2 --un >(gzip > $subtractedFq) --max /dev/null $dbDir/FILES/gDNA - /dev/null";
   sysErrMsg("bash", "-c", $cmd);
 } else {
   verbPrint "Found $subtractedFq. Skipping genome subtraction step...\n"; 
 }


####

if ($EXIT_STATUS) {
    exit $EXIT_STATUS;
}

#### Map to the EEJ:
my $runArgs = "-dbDir=$dbDir -sp=$species -readLen=$le -root=$root";
my $preCmd = getPrefixCmd($subtractedFq);
verbPrint "Mapping reads to the \"splice site-based\" (aka \"a posteriori\") EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie $inpType -p $cores -m 1 -v $bowtieV " .
                "$dbDir/FILES/$species"."_COMBI-M-$le - | " .
             "cut -f 1-4,8 - | sort -T $tmpDir -k 1,1 | " .
             "$binPath/Analyze_COMBI.pl deprecated " .
             "$dbDir/COMBI/$species/$species"."_COMBI-M-$le-gDNA.eff $runArgs";

verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") SIMPLE EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie $inpType -p $cores -m 1 -v $bowtieV " .
                "$dbDir/FILES/EXSK-$le - | " .
             "cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
             "$binPath/Analyze_EXSK.pl $runArgs";

verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") MULTI EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie $inpType -p $cores -m 1 -v $bowtieV " .
                "$dbDir/FILES/MULTI-$le - | " .
             "cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
             "$binPath/Analyze_MULTI.pl $runArgs";

verbPrint "Mapping reads to microexon EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie $inpType -p $cores -m 1 -v $bowtieV " .
                "$dbDir/FILES/$species"."_MIC-$le - | ".
            " cut -f 1-4,8 - | sort -T $tmpDir -k 1,1 | " .
            " $binPath/Analyze_MIC.pl $runArgs";

# Align to intron retention mapped reads here..
unless (($genome_sub and $useGenSub)  or $noIRflag) {
  verbPrint "Mapping reads to intron retention library (version $IR_version)...\n";

# To define version [02/10/15]; minimize changes for users
# $v => "" or "_v2" [v1/v2]
# $type => "new" or "ALL" [v1/v2]
  my $v;
  my $type;
  if ($IR_version == 1){
      $v="";
      $type="new";
  }
  elsif ($IR_version == 2){
      $v="_v2"; 
      $type="ALL";
  }
  
  $preCmd = getPrefixCmd($fq);
  sysErrMsg "$preCmd | $bowtie $inpType -p $cores -m 1 -v $bowtieV " .
              "$dbDir/FILES/$species.IntronJunctions.$type.$le.8 - | " .
              "cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
              "$binPath/MakeSummarySAM.pl | " .
              "$binPath/RI_summarize$v.pl - $runArgs";
  sysErrMsg "$preCmd | $bowtie $inpType -p $cores -m 1 -v $bowtieV " .
                  "$dbDir/FILES/$species.Introns.sample.200 - | " .
              "cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
              "$binPath/MakeSummarySAM.pl | " .
              "$binPath/RI_summarize_introns$v.pl - $runArgs";
} else {
  verbPrint "Skipping intron retention step...\n";
}

unless($keepFlag or $keep_trimmed) {
  verbPrint "Cleaning $fq files!";
  sysErrMsg "rm $fq";
}

unless($keepFlag) {
  verbPrint "Cleaning up $subtractedFq!";
  sysErrMsg "rm $subtractedFq";
}

unless($noIRflag || $IR_version == 2) {  # --UB
    my $juncAnnotationFile = "./to_combine/$root.IR.summary.txt";
    verbPrint "Cleaning up $juncAnnotationFile!";
    sysErrMsg "rm $juncAnnotationFile";
}


verbPrint "Completed " . localtime;
exit $EXIT_STATUS;
