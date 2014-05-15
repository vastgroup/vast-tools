#!/usr/bin/perl -w

use strict;
use File::Which;
use Cwd qw(abs_path);
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
my $readLength; 
my $outdir;

my $legacyFlag = 0;
my $verboseFlag = 1;  # on for debugging 

Getopt::Long::Configure("no_auto_abbrev");
GetOptions("bowtieProg=s" => \$bowtie,
			  "sp=s" => \$species,
			  "db=s" => \$dbDir,
			  "c=i" => \$cores,
			  "pe" => \$pairedEnd,
			  "expr" => \$runExprFlag,
			  "exprONLY" => \$onlyExprFlag,
			  "trim=s" => \$trim,
			  "help" => \$helpFlag,
			  "legacy" => \$legacyFlag,
			  "verbose" => \$verboseFlag,
			  "readLen=i" => \$readLength,
           "output=s" => \$outdir,
			  "o=s" => \$outdir);

our $EXIT_STATUS = 0;

sub sysErrMsg {
  my $sysCommand = shift;
  not system($sysCommand) or die "[vast align error]: $sysCommand Failed in $0!";
}

sub errPrint {
  my $errMsg = shift;
  print STDERR "[vast align error]: $errMsg\n";
  $EXIT_STATUS++; 
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

# Check database directory
unless(defined($dbDir)) {
  $dbDir = "$binPath/../VASTDB";
}
$dbDir = abs_path($dbDir);
$dbDir .= "/$species";
errPrint "The database directory $dbDir does not exist" unless (-e $dbDir);

if (!defined($ARGV[0]) or $helpFlag or $EXIT_STATUS){
    print "\nUsage:

vast-tools align fastq_file_1 [fastq_file_2] [options]

OPTIONS:
	-sp Mmu/Hsa		:	Three letter code for the database (default Hsa)
	-dbDir db		:	Database directory (default vastdb_curVer/Hsa)
	-pe			:	Paired end data? (defaults to off)
	-readLen i		:	Optional read length, otherwise fastq file naming convention enforced (see README)
	-c i			:	# of cores to use for bowtie and pigz (default 1)
	-trim once/twice	:	For trimming, it can be trimmed once (at 3') or twice (in an overlapping manner). (default is twice if length > 50)
	-output OUTPUT, -o	:	Output directory (default <current working directory>)
	-expr			:	For expression analyses: -expr (PSIs plus cRPKM calculations) (default off)
	-exprONLY		:	For expression analyses: -exprONLY (only cRPKMs) (default off)
	-bowtieProg path/bowtie	:	Default is to use the bowtie in PATH, instead you can specify here (default bowtie)

";#NOTE: Recommended to allow at least 15GB of RAM (~10GB are needed for mapping to the genome). For large files (~1 lane), >25GB

#";
  exit $EXIT_STATUS;
}

errPrint "Needs species\n" if !$species;

# Use pigz if installed  --KH
my $zip = which('pigz');
if ($zip eq '') {
    $zip = which('gzip');
    if ($zip eq '') {
        errPrint "couldn't find gzip or pigz\n";
    } 
} else {
    $zip .= " -fp $cores ";
}
verbPrint "Found $zip..." unless $zip eq '';


# Command line flags here
if($pairedEnd and !defined($ARGV[0]) and !defined($ARGV[1])) { $EXIT_STATUS = 1; }


## Getting sample name and length:
my $fq1 = abs_path($ARGV[0]);
my $fq2;
my $fq;     # takes the fastq file to be processed at each step

if(!defined($fq1)) {
  errPrint "No FASTQ file given!";
  $fq1 = "";
}

my $fileName1 = $fq1;
my $fileName2;
my $zipped = isZipped($fq1);
my $subtractedFq;

my($root, $length);

$fileName1 =~ s/^.*\///g; # strip path

my $genome_sub = 0;
if ($fileName1 =~ /\-e\.f/){
    $genome_sub=1;
    ($root,$length)=$fileName1=~/(\S+?)\-(\d{1,4})\-e\.(fastq|fq)(\.gz)?/;  #Fixed regex --TSW
    $fq=$&;
    $subtractedFq = $fq1;
    errPrint "Only for 50nt or 36nt if genome substracted\n" if $length!=36 && $length!=50;
} else {
    # allow readlength to be given by -readLen x --TSW
    if($readLength) {
         $length = $readLength;
         $fileName1 =~ /(\S+)\.(fastq|fq)(\.gz)?/; 
         $root = $1;
    } else { # default behavior by --MI
         ($root,$length)=$fileName1=~/(\S+?)\_?1?\-(\d{1,4})\.(fastq|fq)(\.gz)?/; #Fixed regex --TSW
			if(!defined($length) or $length eq "") { 
  				errPrint "You must either give read length as -readLen i, or rename your fq files name-len.fq";
			}
    }
    if ($pairedEnd){
		$fq2 = abs_path($ARGV[1]);
		$fileName2 = $fq2;
      $fileName2 =~ s/^.*\///g; # strip path
    }
    $fq = $zipped ? "$root-$length.fq.gz" : "$root-$length.fq";
}

#verbPrint "$fileName1\n$fq1\n;"; die ""; # for debugging.
###

# change directories
errPrint "The output directory \"$outdir\" does not exist" unless (-e $outdir);
chdir($outdir) or errPrint "Unable to change directories into output" and die;
verbPrint "Setting output directory to $outdir";
mkdir("spli_out") unless (-e "spli_out");
mkdir("expr_out") unless (-e "expr_out");

#length options:
my ($le, $half_length);
my $difLE;

if ($length >= 50){
    $difLE = $length-50;
    $le = 50;
} elsif ($length >= 36){
    $difLE = $length-36;
    $le=36;
} else {
    errPrint "Minimum reads length: 50nt (Human) and 36nt (Mouse)\n";
}
errPrint "Reads <50nt not available for Human\n" if $le==36 && $species eq "Hsa";
#####

if ($EXIT_STATUS) {
    exit $EXIT_STATUS;
}

if (!$genome_sub){
 my $cmd;
#### Expression analysis (it maps only the first $le nucleotides of the read)
 if ($runExprFlag || $onlyExprFlag){
     verbPrint "Mapping reads against mRNA sequences";

     if ($pairedEnd) {
         $cmd = "$fq1 $fq2";  # altered this to cat both for/rev reads into bowtie.
     } else {
         $cmd = "$fq1";
     }

     $cmd = getPrefixCmd($cmd);
     $cmd .= " | $bowtie -p $cores -m 1 -v 2 -3 $difLE $dbDir/EXPRESSION/mRNA - expr_out/$species"."_mRNA-$le-$root.out";

     sysErrMsg $cmd;
     verbPrint "Calculating cRPKMs\n";
     sysErrMsg "$binPath/expr_RPKM.pl expr_out/$species"."_mRNA-$le-$root.out $dbDir/EXPRESSION/$species"."_mRNA-$le.eff"; 
 }
 if ($onlyExprFlag){
     print STDERR "Expression analysis done\n";
     exit 0;
 }
###

#### Merge PE
 if ($pairedEnd){
     verbPrint "Concatenating paired end reads";
     sysErrMsg "cat $fq1 $fq2 > $fq";  # away with this as well? 
                                       # $fq is used in trimming below. but we
                                       # can pipe into it. KH
 } else {
   $fq = $fq1;
 }
 
#### Trimming
#
# TODO: substitute all of this with gawk 'NR%2==0{print substr($1,5,55)}NR%2==1' INPUT.fq... style 
 my $trimmed = 0;    # flag determining whether trimming occurred
 if ($difLE >= 10){
   if (!defined($trim) or $trim eq "twice"){
	  if ($length > ($le*2)+10){
	     $half_length = sprintf("%.0f", $length / 2);
         verbPrint "Trimming and splitting fastq sequences to $le nt sequences";
         sysErrMsg "cat $fq | $binPath/Trim-twiceOver.pl - $half_length | $binPath/Trim-twiceOver.pl - $le > $root-$le.fq";
	  } else {
	     verbPrint "Trimming and splitting fastq sequences to $le nt sequences";
	     sysErrMsg "$binPath/Trim-twiceOVER.pl $fq $le > $root-$le.fq";
	  }
   } elsif ($trim eq "once"){
	 verbPrint "Trimming fastq sequences to $le nt sequences";
	 sysErrMsg "$binPath/Trim-once.pl $fq $le > $root-$le.fq";
	# sysErrMsg "mv $root-$length-$le.fq $root-$le.fq"; #piping to stdout removes need for this --TSW
   }
   $fq = "$root-$le.fq"; # set new $fq with trimmed reads --KH
   $trimmed = 1;
 }
####
 
#### Get effective reads (i.e. genome substraction).
 verbPrint "Doing genome substraction\n";
 $subtractedFq = "$root-$le-e.fq.gz";
 # Updated genome subtraction command to handle trimmed or untrimmed input files
 $cmd = getPrefixCmd($fq);
 $cmd .= " | $bowtie -p $cores -m 1 -v 2 --un >(gzip > $subtractedFq) --max /dev/null $dbDir/FILES/gDNA - /dev/null";
 sysErrMsg $cmd;

 if ($trimmed) {
     verbPrint "Compressing trimmed reads";
     sysErrMsg "$zip $fq";
 }
 
####
}

#### Map to the EEJ:
my $preCmd = getPrefixCmd($subtractedFq);
verbPrint "Mapping reads to the \"splice site-based\" (aka \"a posteriori\") EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_COMBI-M-$le | cut -f 1-4,8 - | sort -Vu -k 1,1 - | $binPath/Analyze_COMBI.pl deprecated $dbDir/COMBI/$species/$species"."_COMBI-M-$le-gDNA.eff -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";

verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") SIMPLE EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie -p $cores -m 1 -v 2 $dbDir/FILES/EXSK-$le - | cut -f 1-4,8 | sort -Vu -k 1,1 - | $binPath/Analyze_EXSK.pl -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";  

verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") MULTI EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie -p $cores -m 1 -v 2 $dbDir/FILES/MULTI-$le - | cut -f 1-4,8 | sort -Vu -k 1,1 | $binPath/Analyze_MULTI.pl -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";

verbPrint "Mapping reads to microexon EEJ library and Analyzing...\n";
sysErrMsg "$preCmd | $bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_MIC-$le - | cut -f 1-4,8 - | sort -Vu -k 1,1 | $binPath/Analyze_MIC.pl -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";

# TODO Align to intron retention mapped reads here..
verbPrint "Mapping reads to intron retention library...\n";
$preCmd = getPrefixCmd("$fq.gz");
sysErrMsg "$preCmd | $bowtie -p $cores -m 1 -v 2 " .
                "$dbDir/FILES/$sp.IntronJunctions.new.$le.8 - | " .
            "cut -f 1-4,8 | sort -Vu -k 1,1 | " .
            "$binPath/MakeSummarySAM.pl | " .
            "$binPath/RI_summarize.pl";  
sysErrMsg "$preCmd | $bowtie -p $cores -m 1 -v 2 " .
                "$dbDir/FILES/$sp.Introns.sample.200 - | " .
            "cut -f 1-4,8 | sort -Vu -k 1,1 | " .
            "$binPath/MakeSummarySAM.pl | " .
            "$binPath/RI_summarize_introns.pl";

verbPrint "Completed " . localtime;
exit $EXIT_STATUS
