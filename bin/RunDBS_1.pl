#!/usr/bin/perl
use Cwd;
$cwd = getcwd;
($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

if (!$ARGV[0] || $ARGV[0]=~/help/i){
    print "\nCommand: \nperl RunDBS_1.pl fastq_file_1 [fastq_file_2] -sp Mmu/Hsa [-expr/-exprONLY -trim Once/Twice -PE -c N]\n\n";
    
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

$cores=1;
foreach $i (1..$#ARGV){ #options set-up
    if ($ARGV[$i] eq "-sp"){
	if ($ARGV[$i+1]=~/Mm/ || $ARGV[$i+1]=~/Mouse/i){
	    $sp="Mmu";
	}
	elsif ($ARGV[$i+1]=~/Hs/ || $ARGV[$i+1]=~/Human/i){
	    $sp="Hsa";
	}
	else {
	    die "Can't recognize the species\n";
	}
    }
    elsif ($ARGV[$i] eq "-expr"){
	$expr=1;
    }
    elsif ($ARGV[$i] eq "-exprONLY"){
	$exprONLY=1;
    }
    elsif ($ARGV[$i] eq "-trim"){
	if ($ARGV[$i+1]=~/twice/i){
	    $trim="twice";
	}
	elsif ($ARGV[$i+1]=~/once/i){
	    $trim="once";
	}
	else {
	    die "Unrecognized trimming option\n";
	}
    }
    elsif ($ARGV[$i] eq "-PE"){ # paired-end option
	$PE=1;
	$file2=$ARGV[1];
    }
    elsif ($ARGV[$i] eq "-c"){ # number of threds
	$cores=$ARGV[$i+1];
	chomp($cores);
    }
}
die "Needs species\n" if !$sp;

### Getting sample name and length:
$file=$ARGV[0];
if ($file=~/\-e\.f/){
    $genome_sub=1;
    ($root,$length)=$file=~/(.+?)\-(.+?)\-e\.fq/;
    $fq=$&;
    die "Only for 50nt or 36nt if genome substracted\n" if $length!=36 && $length!=50;
}
else {
    if ($PE){
	($root,$length)=$file=~/(.+?)\_1\-(.+?)\.fq/;
	$fq1=$&;
	$fq2=$file2;
	$fq2=~s/\.gz//;
	$fq="$root-$length.fq";
    }
    else {
	($root,$length)=$file=~/(.+?)\-(.+?)\.fq/;
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
die "Reads <50nt not available for Human\n" if $le==36 && $sp eq "Hsa";
#####

system "gunzip $file" if $file=~/\.gz/;
system "gunzip $file2" if $file2=~/\.gz/ && $PE;

if (!$genome_sub){
#### Expression analysis (it maps only the first $le nucleotides of the read)
 if ($expr || $exprONLY){
     print "Mapping reads against mRNA sequences\n";
     system "$dir/bin/bowtie -p $cores -m 1 -v 2 -3 $difLE $dir/$sp/EXPRESSION/mRNA $fq $dir/$sp/EXPRESSION/OUTs/$sp"."mRNA-$le-$root.out" if !$PE;
     system "$dir/bin/bowtie -p $cores -m 1 -v 2 -3 $difLE $dir/$sp/EXPRESSION/mRNA $fq1 $dir/$sp/EXPRESSION/OUTs/$sp"."mRNA-$le-$root.out" if $PE;
     print "Calculating cRPKMs\n";
     system "$dir/bin/expr_RPKM.pl $dir/$sp/EXPRESSION/$sp"."mRNA-$le-$root.out $dir/$sp/EXPRESSION/OUTs/$sp"."_mRNA-$le.eff"; 
 }
 if ($exprONLY){
     print "Compressing raw fastq files\n";
     system "gzip $fq" if !$PE;
     system "gzip $fq1 $fq2" if $PE;
     die "Expression analysis done\n";
 }
###

#### Merge PE
 if ($PE){
     print "Concatenating paired end reads\n";
     system "cat $fq1 $fq2 > $fq";
     system "gzip $fq1 $fq2";
 }
 
#### Trimming
 if ($difLE>=10){
     if ($trim eq "twice" || !$trim){
	 if ($length > ($le*2)+10){
	     $half_length=sprintf("%.0f",$length/2);
	     print "\nTrimming and splitting fastq sequences from $length to $half_length nt\n";
	     system "$dir/bin/Trim-twiceOVER.pl $fq $half_length";
	     print "Trimming and splitting fastq sequences to $le nt sequences\n";
	     system "$dir/bin/Trim-twiceOVER.pl $root-$length-$half_length.fq $le";
	     system "rm $root-$length-$half_length.fq";
	     system "mv $root-$length-$half_length-$le.fq $root-$le.fq";
	 }
	 else {
	     print "\nTrimming and splitting fastq sequences to $le nt sequences\n";
	     system "$dir/bin/Trim-twiceOVER.pl $fq $le";
	     system "mv $root-$length-$le.fq $root-$le.fq";
	 }
     }
     elsif ($trim eq "once"){
	 print "\nTrimming fastq sequences to $le nt sequences\n";
	 system "$dir/bin/Trim-once.pl $fq $le";
	 system "mv $root-$length-$le.fq $root-$le.fq";
     }
 }
 print "Compressing raw fastq file\n" unless $length==50 || $length==36;
 system "gzip $fq" unless $length==50 || $length==36;
####
 
#### Get effective reads (i.e. genome substraction).
 print "\nDoing genome substraction\n";
 system "$dir/bin/bowtie -p $cores -m 1 -v 2 --un $root-$le-e.fq --max /dev/null $dir/$sp/FILES/gDNA $root-$le.fq /dev/null";
 print "Compressing trimmed reads\n";
 system "gzip $root-$le.fq";
####
}

#### Map to the EEJ:
print "\nMapping reads to the \"splice site-based\" (aka \"a posteriori\") EEJ library\n";
system "$dir/bin/bowtie -p $cores -m 1 -v 2 $dir/$sp/FILES/$sp"."_COMBI-M-$le $root-$le-e.fq $dir/$sp/OUTs/$sp"."COMBI-M-$le-$root-e.out";
print "Mapping reads to the \"transcript-based\" (aka \"a priori\") SIMPLE EEJ library\n";
system "$dir/bin/bowtie -p $cores -m 1 -v 2 $dir/$sp/FILES/EXSK-$le $root-$le-e.fq $dir/$sp/OUTs/$sp"."EXSK-$le-$root-e.out";  
print "Mapping reads to the \"transcript-based\" (aka \"a priori\") MULTI EEJ library\n";
system "$dir/bin/bowtie -p $cores -m 1 -v 2 $dir/$sp/FILES/MULTI-$le $root-$le-e.fq $dir/$sp/OUTs/$sp"."MULTI-$le-$root-e.out";
print "Mapping reads to microexon EEJ library\n";
system "$dir/bin/bowtie -p $cores -m 1 -v 2 $dir/$sp/FILES/$sp"."_MIC-$le $root-$le-e.fq $dir/$sp/OUTs/$sp"."MIC-$le-$root-e.out";
print "Compressing genome-substracted reads\n";
system "gzip $root-$le-e.fq";
####

## Analyze MIC
print "\nStarting EEJ analyses:\n";
print "Analyzing microexons\n";
system "$dir/bin/Analyze_MIC.pl $dir/$sp/OUTs/$sp"."MIC-$le-$root-e.out";

## Analyze MULTI and EXSK (A priori pipeline)
print "Sorting a priori outputs\n";
system "$dir/bin/sort_outs.pl $dir/$sp/OUTs/$sp"."MULTI-$le-$root-e.out";
system "$dir/bin/sort_outs.pl $dir/$sp/OUTs/$sp"."EXSK-$le-$root-e.out";
print "Analyzing a priori outputs\n";
system "$dir/bin/Analyze_EXSK.pl $dir/$sp/OUTs/$sp"."EXSK-$le-$root-e_s.out";
system "$dir/bin/Analyze_MULTI.pl $dir/$sp/OUTs/$sp"."MULTI-$le-$root-e_s.out";

## Analyze a posteriori pipeline
print "Sorting a posteriori output\n";
system "$dir/bin/sort_outs.pl $dir/$sp/OUTs/$sp"."COMBI-M-$le-$root-e.out";
print "Analyzing a posteriori output for exon skippings\n";
system "$dir/bin/Analyze_COMBI.pl $dir/$sp/OUTs/$sp"."COMBI-M-$le-$root-e_s.out $dir/COMBI/$sp/$sp"."_COMBI-M-$le-gDNA.eff";
##

