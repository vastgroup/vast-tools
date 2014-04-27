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
my $trim = "once"; # 1 by default
my $cores = 1; #default
my $readLength; 

my $legacyFlag = 0;
my $verboseFlag = 0;

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

# Use pigz if installed  --KH
my $zip = which('pigz');
if ($zip eq '') {
    $zip = which('gzip');
    if ($zip eq '') {
        die "Error: need gzip or pigz\n";
    }
} else {
    print STDERR "[vastdb align] Found pigz...";
    $zip .= " -fp $cores ";
}

sub sysErrMsg {
  my $sysCommand = shift;
  not system($sysCommand) or die "[vastdb align error]: $sysCommand Failed in $0!";
}

# Set up output directory structure
unless($legacyFlag) {
  mkdir("align_out") unless (-e "align_out");
  mkdir("expr_out") unless (-e "expr_out");
  mkdir("align_out/$species") unless (-e "align_out/$species");
  mkdir("expr_out/$species") unless (-e "align_out/$species");
}

# Command line flags here
if($pairedEnd and !defined($ARGV[0]) and !defined($ARGV[1])) { $helpFlag = 1; }


if (!$ARGV[0] or $helpFlag){
    print "\nCommand: \nvastdb align fastq_file_1 [fastq_file_2] -sp Mmu/Hsa [-expr/-exprONLY -trim once/twice -pe -c 1 -bowtieProg bowtie]\n\n";
    
    print ">> fastq_file can be compressed or uncompressed. It MUST have this name format: Sample-readlength.fq\n";
    print ">>>> Only underscores (\"_\") are allowed in the Sample name. No dashes, dots, spaces, etc.\n";
    print ">>>> If running with the paired end option (\"\-PE\"), they must be called: Sample_1-length.fq Sample_2-length.fq\n";
    print ">> The process can be started with \"genome substracted\" samples if a Sample-lenght-e.fq is used\n";
    print ">> Species implemented: Human (hg19) or Mouse (mm9)\n";
    print ">> For expression analyses: -expr (PSIs plus cRPKM calculations) OR -exprONLY (only cRPKMs) OR none\n";
    print ">> For trimming, it can be trimmed once (at 3') or twice (in an overlapping manner).\n";
    print ">>>> If nothing is provided, the program will decide based on the length of the reads (Default is twice if length>50)\n";
    print ">> -c N: Number of threds used when running bowtie (default=1).\n";
    die ">> Recommended to allow at least 15GB of RAM (~10GB are needed for mapping to the genome). For large files (~1 lane), >25GB\n\n";
}

die "Needs species\n" if !$species;



### Getting sample name and length:
$file=$ARGV[0];
my $zipped = ($file =~ /\.gz$/) ? 1 : 0;
if ($file=~/\-e\.f/){
    $genome_sub=1;
    ($root,$length)=$file=~/(.+?)\-(.+?)\-e\.(fastq|fq)(\.gz)?/;
    $fq=$&;
    die "Only for 50nt or 36nt if genome substracted\n" if $length!=36 && $length!=50;
}
else {
    if ($pairedEnd){
	($root,$length)=$file=~/(.+?)\_1\-(.+?)\.(fastq|fq)(\.gz)?/;
	$fq1=$&;
	$fq2=$ARGV[1];
    #$fq2=~s/\.gz//;
	$fq="$root-$length.fq";
    }
    else {
	($root,$length)=$file=~/(.+?)\-(.+?)\.(fastq|fq)(\.gz)?/;
	$fq=$&;
    }
}
###

#length options:
if ($length>=50){
    $difLE=$length-50;
    $le=50;
}
elsif ($length>=36){
    $difLE=$length-36;
    $le=36;
}
else {
    die "Minimum reads length: 50nt (Human) and 36nt (Mouse)\n";
}
die "Reads <50nt not available for Human\n" if $le==36 && $species eq "Hsa";
#####

#sysErrMsg "gunzip $file" if $file=~/\.gz/;
#sysErrMsg "gunzip $file2" if $file2=~/\.gz/ && $pairedEnd;

if (!$genome_sub){
#### Expression analysis (it maps only the first $le nucleotides of the read)
 if ($runExprFlag || $onlyExprFlag){
     print STDERR "[vastdb align] Mapping reads against mRNA sequences\n";

     my $cmd = "$bowtie -p $cores -m 1 -v 2 -3 $difLE $dbDir/EXPRESSION/mRNA - expr_out/$species"."mRNA-$le-$root.out";

     if (!$pairedEnd) {
         $cmd = "$fq | $cmd";
     } else {
         $cmd = "$fq1 | $cmd";
     }

     if ($zipped) {
         $cmd = "gzip -dc $cmd";
     } else {
         $cmd = "cat $cmd";
     }

     sysErrMsg $cmd;
     print STDERR "[vastdb align] Calculating cRPKMs\n";
     not sysErrMsg "$binPath/expr_RPKM.pl $dbDir/EXPRESSION/$species"."mRNA-$le-$root.out expr_out/$species"."_mRNA-$le.eff"
					or die "expr_RPKM.pl failed!\n"; 
 }
 if ($onlyExprFlag){
     #print "Compressing raw fastq files\n";
     #sysErrMsg "gzip $fq" if !$pairedEnd;
     #sysErrMsg "gzip $fq1 $fq2" if $pairedEnd;
     die "Expression analysis done\n";
 }
###

#### Merge PE
 if ($pairedEnd){
     print STDERR "[vastdb align] Concatenating paired end reads\n";
     sysErrMsg "cat $fq1 $fq2 > $fq";
     #sysErrMsg "gzip $fq1 $fq2";
 }
 
#### Trimming
 if ($difLE>=10){
     if ($trim eq "twice" || !$trim){
	 if ($length > ($le*2)+10){
	     $half_length=sprintf("%.0f",$length/2);
	     print STDERR "[vastdb align] \nTrimming and splitting fastq sequences from $length to $half_length nt\n";
	     sysErrMsg "$binPath/Trim-twiceOVER.pl $fq $half_length";
	     print STDERR "[vastdb align] Trimming and splitting fastq sequences to $le nt sequences\n";
	     sysErrMsg "$binPath/Trim-twiceOVER.pl $root-$length-$half_length.fq $le";
	     sysErrMsg "rm $root-$length-$half_length.fq";
	     sysErrMsg "mv $root-$length-$half_length-$le.fq $root-$le.fq";
	 }
	 else {
	     print STDERR "[vastdb align] \nTrimming and splitting fastq sequences to $le nt sequences\n";
	     sysErrMsg "$binPath/Trim-twiceOVER.pl $fq $le";
	     sysErrMsg "mv $root-$length-$le.fq $root-$le.fq";
	 }
     }
     elsif ($trim eq "once"){
	 print STDERR "[vastdb align] \nTrimming fastq sequences to $le nt sequences\n";
	 sysErrMsg "$binPath/Trim-once.pl $fq $le";
	 sysErrMsg "mv $root-$length-$le.fq $root-$le.fq";
     }
 }
 print STDERR "[vastdb align] Compressing raw fastq file\n" unless $length==50 || $length==36;
 #sysErrMsg "gzip $fq" unless $length==50 || $length==36;
####
 
#### Get effective reads (i.e. genome substraction).
 print STDERR "[vastdb align] \nDoing genome substraction\n";
 sysErrMsg "$bowtie -p $cores -m 1 -v 2 --un $root-$le-e.fq --max /dev/null $dbDir/FILES/gDNA $root-$le.fq /dev/null";
 print STDERR "[vastdb align] Compressing trimmed reads\n";
 sysErrMsg "$zip $root-$le.fq";
####
}

#### Map to the EEJ:
print STDERR "[vastdb align] \nMapping reads to the \"splice site-based\" (aka \"a posteriori\") EEJ library\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_COMBI-M-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - > align_out/$species"."COMBI-M-$le-$root-e_s.out";
print STDERR "[vastdb align] Mapping reads to the \"transcript-based\" (aka \"a priori\") SIMPLE EEJ library\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/EXSK-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - > align_out/$species"."EXSK-$le-$root-e_s.out";  
print STDERR "[vastdb align] Mapping reads to the \"transcript-based\" (aka \"a priori\") MULTI EEJ library\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/MULTI-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 > align_out/$species"."MULTI-$le-$root-e.out";
print STDERR "[vastdb align] Mapping reads to microexon EEJ library\n";
sysErrMsg "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_MIC-$le $root-$le-e.fq $species"."MIC-$le-$root-e.bam";
print STDERR "[vastdb align] Compressing genome-substracted reads\n";
sysErrMsg "$zip $root-$le-e.fq";
####

## Analyze MIC
print STDERR "[vastdb align] \nStarting EEJ analyses:\n";
print STDERR "[vastdb align] Analyzing microexons\n";
sysErrMsg "$binPath/Analyze_MIC.pl align_out/$species"."MIC-$le-$root-e.out";

## Analyze MULTI and EXSK (A priori pipeline)
#print "Sorting a priori outputs\n";
#sysErrMsg "$binPath/sort_outs.pl align_out/$species"."MULTI-$le-$root-e.out";
#sysErrMsg "$binPath/sort_outs.pl align_out/$species"."EXSK-$le-$root-e.out";
print STDERR "[vastdb align] Analyzing a priori outputs\n";
sysErrMsg "$binPath/Analyze_EXSK.pl align_out/$species"."EXSK-$le-$root-e_s.out";
sysErrMsg "$binPath/Analyze_MULTI.pl align_out/$species"."MULTI-$le-$root-e_s.out";

## Analyze a posteriori pipeline
#print "Sorting a posteriori output\n";
#sysErrMsg "$binPath/sort_outs.pl align_out/$species"."COMBI-M-$le-$root-e.out";
print STDERR "[vastdb align] Analyzing a posteriori output for exon skippings\n";
sysErrMsg "$binPath/Analyze_COMBI.pl align_out/$species"."COMBI-M-$le-$root-e_s.out $dbDir/COMBI/$species/$species"."_COMBI-M-$le-gDNA.eff";
##

sysErrMsg "$zip align_out/$species/*.out";

