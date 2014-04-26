#!/usr/bin/perl -w

use File::Which;
use Cwd qw(abs_path);
use Getopt::Long;

# INITIALIZE
my $binPath = abs_path($0);
$binPath =~ s/\/$0$//;

my $bowtie = "bowtie"; # by default;
my $species = "Hsa"; # by default;
my $dbDir = "$binPath/../$species";
my $pairedEnd = 0; # no by default
my $runExprFlag = 0; # no by default
my $onlyExprFlag = 0; # no by default
my $trim;
my $cores = 1;
my $readLength;

my $legacyFlag = 0;
my $verboseFlag = 0;

GetOptions("bowtieProg=s" => \$bowtie,
			  "sp|SP=s" => \$species,
			  "db|DB=s" => \$dbDir,
			  "c" => \$cores,
			  "pe|PE" => \$pairedEnd,
			  "expr" => \$runExprFlag,
			  "exprONLY" => \$onlyExprFlag,
			  "trim=s" => \$trim,
			  "help" => \$helpFlag,
			  "legacy" => \$legacyFlag,
			  "verbose" => \$verboseFlag,
			  "readLen=i" => \$readLength);

# Use pigz if installed
my $zip = which('pigz');
if ($zip eq '') {
    $zip = which('gzip');
    if ($zip eq '') {
        die "Error: need gzip or pigz\n";
    }
} else {
    print "Found pigz...";
    $zip .= " -fp $cores ";
}

# Command line flags here

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
	$fq2=$file2;
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

#system "gunzip $file" if $file=~/\.gz/;
#system "gunzip $file2" if $file2=~/\.gz/ && $pairedEnd;

if (!$genome_sub){
#### Expression analysis (it maps only the first $le nucleotides of the read)
 if ($runExprFlag || $onlyExprFlag){
     print "Mapping reads against mRNA sequences\n";

     my $cmd = "$bowtie -p $cores -m 1 -v 2 -3 $difLE $dbDir/EXPRESSION/mRNA - $dbDir/EXPRESSION/OUTs/$species"."mRNA-$le-$root.out";

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

     system $cmd;
     print "Calculating cRPKMs\n";
     system "$binPath/expr_RPKM.pl $dbDir/EXPRESSION/$species"."mRNA-$le-$root.out $dbDir/EXPRESSION/OUTs/$species"."_mRNA-$le.eff"; 
 }
 if ($onlyExprFlag){
     #print "Compressing raw fastq files\n";
     #system "gzip $fq" if !$pairedEnd;
     #system "gzip $fq1 $fq2" if $pairedEnd;
     die "Expression analysis done\n";
 }
###

#### Merge PE
 if ($pairedEnd){
     print "Concatenating paired end reads\n";
     system "cat $fq1 $fq2 > $fq";
     #system "gzip $fq1 $fq2";
 }
 
#### Trimming
 if ($difLE>=10){
     if ($trim eq "twice" || !$trim){
	 if ($length > ($le*2)+10){
	     $half_length=sprintf("%.0f",$length/2);
	     print "\nTrimming and splitting fastq sequences from $length to $half_length nt\n";
	     system "$binPath/Trim-twiceOVER.pl $fq $half_length";
	     print "Trimming and splitting fastq sequences to $le nt sequences\n";
	     system "$binPath/Trim-twiceOVER.pl $root-$length-$half_length.fq $le";
	     system "rm $root-$length-$half_length.fq";
	     system "mv $root-$length-$half_length-$le.fq $root-$le.fq";
	 }
	 else {
	     print "\nTrimming and splitting fastq sequences to $le nt sequences\n";
	     system "$binPath/Trim-twiceOVER.pl $fq $le";
	     system "mv $root-$length-$le.fq $root-$le.fq";
	 }
     }
     elsif ($trim eq "once"){
	 print "\nTrimming fastq sequences to $le nt sequences\n";
	 system "$binPath/Trim-once.pl $fq $le";
	 system "mv $root-$length-$le.fq $root-$le.fq";
     }
 }
 print "Compressing raw fastq file\n" unless $length==50 || $length==36;
 #system "gzip $fq" unless $length==50 || $length==36;
####
 
#### Get effective reads (i.e. genome substraction).
 print "\nDoing genome substraction\n";
 system "$bowtie -p $cores -m 1 -v 2 --un $root-$le-e.fq --max /dev/null $dbDir/FILES/gDNA $root-$le.fq /dev/null";
 print "Compressing trimmed reads\n";
 system "$zip $root-$le.fq";
####
}

#### Map to the EEJ:
print "\nMapping reads to the \"splice site-based\" (aka \"a posteriori\") EEJ library\n";
system "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_COMBI-M-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - > $dbDir/OUTs/$species"."COMBI-M-$le-$root-e_s.out";
print "Mapping reads to the \"transcript-based\" (aka \"a priori\") SIMPLE EEJ library\n";
system "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/EXSK-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 - > $dbDir/OUTs/$species"."EXSK-$le-$root-e_s.out";  
print "Mapping reads to the \"transcript-based\" (aka \"a priori\") MULTI EEJ library\n";
system "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/MULTI-$le $root-$le-e.fq | cut -f 1-4,8 - | sort -u -k 1,1 > $dbDir/OUTs/$species"."MULTI-$le-$root-e.out";
print "Mapping reads to microexon EEJ library\n";
system "$bowtie -p $cores -m 1 -v 2 $dbDir/FILES/$species"."_MIC-$le $root-$le-e.fq $species"."MIC-$le-$root-e.bam";
print "Compressing genome-substracted reads\n";
system "$zip $root-$le-e.fq";
####

## Analyze MIC
print "\nStarting EEJ analyses:\n";
print "Analyzing microexons\n";
system "$binPath/Analyze_MIC.pl $dbDir/OUTs/$species"."MIC-$le-$root-e.out";

## Analyze MULTI and EXSK (A priori pipeline)
#print "Sorting a priori outputs\n";
#system "$binPath/sort_outs.pl $dbDir/OUTs/$species"."MULTI-$le-$root-e.out";
#system "$binPath/sort_outs.pl $dbDir/OUTs/$species"."EXSK-$le-$root-e.out";
print "Analyzing a priori outputs\n";
system "$binPath/Analyze_EXSK.pl $dbDir/OUTs/$species"."EXSK-$le-$root-e_s.out";
system "$binPath/Analyze_MULTI.pl $dbDir/OUTs/$species"."MULTI-$le-$root-e_s.out";

## Analyze a posteriori pipeline
#print "Sorting a posteriori output\n";
#system "$binPath/sort_outs.pl $dbDir/OUTs/$species"."COMBI-M-$le-$root-e.out";
print "Analyzing a posteriori output for exon skippings\n";
system "$binPath/Analyze_COMBI.pl $dbDir/OUTs/$species"."COMBI-M-$le-$root-e_s.out $dbDir/COMBI/$species/$species"."_COMBI-M-$le-gDNA.eff";
##

system "$zip $dbDir/OUTs/$species/*.out";

