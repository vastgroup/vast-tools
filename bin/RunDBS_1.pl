#!/usr/bin/perl -w

use File::Which;
use Cwd qw(abs_path);
use Getopt::Long;

# INITIALIZE PATH AND FLAGS--TSW
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

my $bowtie = "bowtie"; # by default;
my $species = "Hsa"; # by default;
my $dbDir = "$binPath/../$species"; #default
my $pairedEnd = 0; # no by default
my $runExprFlag = 0; # no by default
my $onlyExprFlag = 0; # no by default
my $trim;
my $cores = 1; #default
my $readLength; 

my $legacyFlag = 0;
my $verboseFlag = 1;  # on for debugging 

GetOptions("bowtieProg=s" => \$bowtie,
			  "sp=s" => \$species,
			  "db=s" => \$dbDir,
			  "c" => \$cores,
			  "pe" => \$pairedEnd,
			  "expr" => \$runExprFlag,
			  "exprONLY" => \$onlyExprFlag,
			  "trim=s" => \$trim,
			  "help" => \$helpFlag,
			  "legacy" => \$legacyFlag,
			  "verbose" => \$verboseFlag,
			  "readLen=i" => \$readLength);

our $EXIT_STATUS = 0;

sub sysErrMsg {
  my $sysCommand = shift;
  not system($sysCommand) or die "[vastdb align error]: $sysCommand Failed in $0!";
}

sub errPrint {
  my $errMsg = shift;
  print STDERR "[vastdb align error]: $errMsg\n";
  $EXIT_STATUS = 1; 
}

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vastdb align]: $verbMsg\n";
  }
}

# Set up output directory structure
mkdir("align_out") unless (-e "align_out");
mkdir("expr_out") unless (-e "expr_out");
mkdir("align_out/$species") unless (-e "align_out/$species");
mkdir("expr_out/$species") unless (-e "align_out/$species");

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


### Getting sample name and length:
my $fq1 = $ARGV[0];
my $fq2;
my $fileName1 = $fq1;
my $fileName2;
my $zipped = ($fq1 =~ /\.gz$/) ? 1 : 0;

my($root, $length);

$fileName1 =~ s/^.*\///g; # strip path

if ($fileName1 =~ /\-e\.f/){
    $genome_sub=1;
    ($root,$length)=$fileName1=~/(\S+?)\-(\d{1,4})\-e\.(fastq|fq)(\.gz)?/;  #Fixed regex --TSW
    $fq=$&;
    errPrint "Only for 50nt or 36nt if genome substracted\n" if $length!=36 && $length!=50;
} else {
    # allow readlength to be given by -readLen x --TSW
    if($readLength) {
         $length = $readLength;
         $fileName1 =~ /(\S+)\.(fastq|fq)(\.gz)?/; 
         $root = $1;
    } else { # default behavior by --MI
         ($root,$length)=$fileName1=~/(\S+?)\_{0,1}1{0,1}\-(\d{1,4})\.(fastq|fq)(\.gz)?/; #Fixed regex --TSW
			if(!defined($length) or $length eq "") { 
  				errPrint "You must either give read length as -readLen i, or rename your fq files name-len.fq";
			}
    }
    if ($pairedEnd){
		$fq2 = $ARGV[1];
		$fileName2 = $fq2;
      $fileName2 =~ s/^.*\///g; # strip path
    }
}


#verbPrint "$fileName1\n$fq1\n;"; die ""; # for debugging.
###

#length options:
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

# move this.
if (!defined($ARGV[0]) or $helpFlag or $EXIT_STATUS){
    print "\nUsage:

vastdb align fastq_file_1 [fastq_file_2] [options]

OPTIONS:
	-sp Mmu/Hsa		:	Three letter code for the database (default Hsa)
	-dbDir db		:	Database directory (default vastdb_curVer/Hsa)
	-pe			:	Paired end data? (defaults to off)
	-readLen i		:	Optional read length, otherwise fastq file naming convention enforced (see README)
	-c i			:	# of cores to use for bowtie and pigz (default 1)
	-trim once/twice	:	For trimming, it can be trimmed once (at 3') or twice (in an overlapping manner). (default is twice if length > 50)
	-expr			:	For expression analyses: -expr (PSIs plus cRPKM calculations) (default off)
	-exprONLY		:	For expression analyses: -exprONLY (only cRPKMs) (default off)
	-bowtieProg path/bowtie	:	Default is to use the bowtie in PATH, instead you can specify here (default bowtie)

NOTE: Recommended to allow at least 15GB of RAM (~10GB are needed for mapping to the genome). For large files (~1 lane), >25GB

";
  exit $EXIT_STATUS;
}

die "Needs species\n" if !$species;


#sysErrMsg "gunzip $file" if $file=~/\.gz/;
#sysErrMsg "gunzip $file2" if $file2=~/\.gz/ && $pairedEnd;

if (!$genome_sub){
#### Expression analysis (it maps only the first $le nucleotides of the read)
 if ($runExprFlag || $onlyExprFlag){
     verbPrint "Mapping reads against mRNA sequences\n";

     my $cmd = "$bowtie -p $cores -m 1 -v 2 -3 $difLE $dbDir/EXPRESSION/mRNA - expr_out/$species"."_mRNA-$le-$root.out";

     if ($pairedEnd) {
         $cmd = "$fq1 $fq2 | $cmd";  # altered this to cat both for/rev reads into bowtie.
     } else {
         $cmd = "$fq1 | $cmd";
     }

     if ($zipped) {
			$cmd = "zcat $cmd";
	 # Is zcat a better choice?  --TSW
    #     $cmd = "gzip -dc $cmd";
     } else {
         $cmd = "cat $cmd";
     }

     sysErrMsg $cmd;
     verbPrint "Calculating cRPKMs\n";
     sysErrMsg "$binPath/expr_RPKM.pl $dbDir/EXPRESSION/$species"."mRNA-$le-$root.out expr_out/$species"."_mRNA-$le.eff"; 
 }
 if ($onlyExprFlag){
     #print "Compressing raw fastq files\n";
     #sysErrMsg "gzip $fq" if !$pairedEnd;
     #sysErrMsg "gzip $fq1 $fq2" if $pairedEnd;
     print STDERR "Expression analysis done\n";
     exit 0;
 }
###

#### Merge PE
 if ($pairedEnd){
     verbPrint "Concatenating paired end reads\n";
     sysErrMsg "cat $fq1 $fq2 > $fq";  # away with this as well? 
                                       # $fq is used in trimming below. but we
                                       # can pipe into it. KH
     #sysErrMsg "gzip $fq1 $fq2";
 } else {
   $fq = $fq1;
 }
 
#### Trimming
 if ($difLE >= 10){
   if (!defined($trim) or $trim eq "twice"){
	  if ($length > ($le*2)+10){
	     $half_length = sprintf("%.0f", $length / 2);
	     verbPrint "Trimming and splitting fastq sequences from $length to $half_length nt\n";
	     sysErrMsg "$binPath/Trim-twiceOVER.pl $fq $half_length > $root-$length-$half_length.fq";
	     verbPrint "Trimming and splitting fastq sequences to $le nt sequences\n";
	     sysErrMsg "$binPath/Trim-twiceOVER.pl $root-$length-$half_length.fq $le > $root-$le.fq";
	     sysErrMsg "rm $root-$length-$half_length.fq";
	     #sysErrMsg "mv $root-$length-$half_length-$le.fq $root-$le.fq";  #piping to stdout removes need for this --TSW
	  } else {
	     verbPrint "Trimming and splitting fastq sequences to $le nt sequences\n";
	     sysErrMsg "$binPath/Trim-twiceOVER.pl $fq $le > $root-$le.fq";
	  #   sysErrMsg "mv $root-$length-$le.fq $root-$le.fq"; #piping to stdout removes need for this --TSW
	  }
   } elsif ($trim eq "once"){
	 verbPrint "Trimming fastq sequences to $le nt sequences\n";
	 sysErrMsg "$binPath/Trim-once.pl $fq $le > $root-$le.fq";
	# sysErrMsg "mv $root-$length-$le.fq $root-$le.fq"; #piping to stdout removes need for this --TSW
   }
 }
 verbPrint "Compressing raw fastq file\n" unless $length==50 || $length==36;
 #sysErrMsg "gzip $fq" unless $length==50 || $length==36;
####
 
#### Get effective reads (i.e. genome substraction).
 verbPrint "Doing genome substraction\n";
 sysErrMsg "$bowtie -p $cores -m 1 -v 2 --un $root-$le-e.fq --max /dev/null $dbDir/FILES/gDNA $root-$le.fq /dev/null";
 verbPrint "Compressing trimmed reads\n";
 sysErrMsg "$zip $root-$le.fq";
####
}

#### Map to the EEJ:
verbPrint "Mapping reads to the \"splice site-based\" (aka \"a posteriori\") EEJ library and Analyzing...\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_COMBI-M-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - | $binPath/Analyze_COMBI.pl deprecated $dbDir/COMBI/$species/$species"."_COMBI-M-$le-gDNA.eff -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";
#sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_COMBI-M-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - > align_out/$species"."COMBI-M-$le-$root-e_s.out"; # DEPRECATED --TSW

verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") SIMPLE EEJ library and Analyzing...\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/EXSK-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - | $binPath/Analyze_EXSK.pl -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";  
#sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/EXSK-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - > align_out/$species"."EXSK-$le-$root-e_s.out"; # DEPRECATED --TSW

verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") MULTI EEJ library and Analyzing...\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/MULTI-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 | $binPath/Analyze_MULTI.pl -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";
#sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/MULTI-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 > align_out/$species"."MULTI-$le-$root-e_s.out"; # DEPRECATED --TSW

verbPrint "Mapping reads to microexon EEJ library and Analyzing...\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_MIC-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 | $binPath/Analyze_MIC.pl -dbDir=$dbDir -sp=$species -readLen=$le -root=$root";
#sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_MIC-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 > align_out/$species"."MIC-$le-$root-e.out"; # DEPRECATED --TSW

verbPrint "Compressing genome-substracted reads\n";
sysErrMsg "$zip $root-$le-e.fq";
####

## Analyze MIC
#verbPrint "Starting EEJ analyses:\n";
#verbPrint "Analyzing microexons\n";
#sysErrMsg "$binPath/Analyze_MIC.pl align_out/$species"."MIC-$le-$root-e.out -dbDir=$dbDir";

## Analyze MULTI and EXSK (A priori pipeline)
#print "Sorting a priori outputs\n";
#sysErrMsg "$binPath/sort_outs.pl align_out/$species"."MULTI-$le-$root-e.out";
#sysErrMsg "$binPath/sort_outs.pl align_out/$species"."EXSK-$le-$root-e.out";
#verbPrint "Analyzing a priori outputs\n";
#sysErrMsg "$binPath/Analyze_EXSK.pl align_out/$species"."EXSK-$le-$root-e_s.out -dbDir=$dbDir";
#sysErrMsg "$binPath/Analyze_MULTI.pl align_out/$species"."MULTI-$le-$root-e_s.out -dbDir=$dbDir";

## Analyze a posteriori pipeline
#print "Sorting a posteriori output\n";
#sysErrMsg "$binPath/sort_outs.pl align_out/$species"."COMBI-M-$le-$root-e.out";
#verbPrint "Analyzing a posteriori output for exon skippings\n";
#sysErrMsg "$binPath/Analyze_COMBI.pl align_out/$species"."COMBI-M-$le-$root-e_s.out $dbDir/COMBI/$species/$species"."_COMBI-M-$le-gDNA.eff";
##

sysErrMsg "$zip align_out/$species/*.out";  #should this just clean up instead?  --TSW

