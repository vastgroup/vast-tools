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
use File::Path qw(make_path);

# INITIALIZE PATH AND FLAGS--TSW
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

my ($rc1,$rc2,$nrc1,$nrc2)=(undef,undef,undef,undef);  # if non of these specified by user, VTS tried to infer strandness of reads
my $helpFlag = 0;
my $bowtie = "bowtie"; # by default;
my $species = "Hsa"; # by default;
my $dbDir; #default
my $pairedEnd = 0; # no by default
my $notstrandaware=0;   # no by default
my $runExprFlag = 0; # no by default
my $onlyExprFlag = 0; # no by default
my $trim;
my $cores = 1; #default
my $readLength = ""; # default.
my $outdir;
my $noIRflag = 0;  # don't run intron retention (for speed..)
my $onlyIRflag = 0; # only run intron retention
my $stringentIRflag = 0; # Run extra genome/eej subtraction step
my $IR_version = 2; # IR version [09/Nov/2015] [new default 01/04/16]
my $minReadNum;
my $print_EEJs = 0; # print out a file with EEJ counts for MIC, MULTI and EXSK
my $trim5 = 0; # number of nt skipped at the 5 prime on R1 for strand specificity and expr

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

my $resume = 0;   # if this flag is set, vast-tools tries to resume a previous run 

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(		  "bowtieProg=s" => \$bowtie,
			  "sp=s" => \$species,
			  "dbDir=s" => \$dbDir,
			  "c=i" => \$cores, 
			  "cores=i" => \$cores,
			  "expr" => \$runExprFlag,
			  "ns" => \$notstrandaware,
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
			  "onlyIR" => \$onlyIRflag,
			  "stringentIR" => \$stringentIRflag,
			  "IR_version=i" => \$IR_version, 
			  "keep" => \$keepFlag,
			  "minReadDepth=i" => \$minReadNum, #to do
			  "EEJ_counts" => \$print_EEJs,
			  "ec" => \$print_EEJs,
			  "tmpDir=s" => \$tmpDir,
			  "stepSize=i" => \$trimStep,
			  "trimOnce" => \$trimOnceFlag,
			  "findSubtracted" => \$useGenSub,
                          "trimLen=i" => \$trimLen,
                          "mismatchNum=i" => \$bowtieV,
                          "preTrimmed" => \$trimmed,
                          "useFastq" => \$fastaOnly,
			  "riboFoot" => \$ribofoot,
			  "trim5=i" => \$trim5,
			  "resume" =>\$resume,
			  "rc1" => \$rc1,
			  "rc2" => \$rc2,
			  "nrc1" => \$nrc1,
			  "nrc2" => \$nrc2,
			  );

our $EXIT_STATUS = 0;

if(defined($rc1) && defined($nrc1)){errPrintDie("Only one of these two arguments --rc1 and --nrc1 can be given at a time.")}
if(defined($rc2) && defined($nrc2)){errPrintDie("Only one of these two arguments --rc2 and --nrc2 can be given at a time.")}
if(defined($rc1)){$rc1=1;}
if(defined($nrc1)){$rc1=0;}
if(defined($rc2)){$rc2=1;}
if(defined($nrc2)){$rc2=0;}
if(defined($rc1) && !defined($rc2)){$rc2=0;}
if(defined($rc2) && !defined($rc1)){$rc1=0;}

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
	if($fn =~ /(\S+)\.(fastq|fq)(\.gz)?/){  # FASTQ files
	    while(<$fh>){if($check){if(substr($_,0,1) ne "@"){errPrintDie("Sequence data must be provided in FASTQ format but file $fn does not look like FASTQ format (first line does not start with @).");}$check=0;}
			 chomp;
			 $c++;
			 if($c % 4==2){
			     unless($readL){
				 $readL=length($_);
				 $tally_readL{$readL}=1;
			     }
			     else{ 
				 $readL=length($_);
				 $tally_readL{$readL}++ if defined $tally_readL{$readL};
				 $tally_readL{$readL}=1 if !defined $tally_readL{$readL};
				 $total++;
			     }
			 }
			 if($c/4 > $maxN){last;}
	    }
	}elsif($fn =~ /(\S+)\.(fasta|fa)(\.gz)?/){  # FASTA files
	    my $read="";
	    my $str_tmp;
	    while(<$fh>){
		if($c>$maxN){last;}
		chomp;
		$str_tmp=substr($_,0,1);
		if(length($_)==0){next;}
		if($str_tmp eq "#"){next;}
		if($str_tmp eq ">"){
		    if($read ne ""){
			$c++;
			unless($readL){
			    $readL=length($read);
			    $tally_readL{$readL}=1;
			}else{ 
			    $readL=length($read);
			    $tally_readL{$readL}++ if defined $tally_readL{$readL};
			    $tally_readL{$readL}=1 if !defined $tally_readL{$readL};
			    $total++;
			}
		    }
		    $read="";next;
		}
		$read.=$_;
	    }
	}else{
	    errPrintDie("Type of sequence file $fn cannot be deduced from file name (must end with fastq|fq|fasta|fa (.gz).");
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
    if($resume){print STDERR "... --resumed!\n";return(1);}
#  print "\n".join(" ",@sysCommand)."\n";
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

sub checkResumeOption{
    my @files_to_be_checked=@_;   # if any of these files exist, we can skip the next computation 
    if($resume == 0){return(1);}
    $resume=0;   # assume we need to execute the next command
    foreach my $file_to_be_checked (@files_to_be_checked){
  	if(-e $file_to_be_checked){$resume=1;last;}  # except if one of these files exist
    }
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
	--sp Hsa/Mmu/etc	Three letter code for the database (default Hsa)
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
        --onlyIR                Only run intron retention pipeline (default off) 
	--stringentIR		Don't run first filtering step of IR 
				(this will increase speed a little) (default off)
        --IR_version 1/2        Version of the Intron Retention analysis (default 2)
	--keep			Don't remove trimmed and genome-subtracted reads 
				after use. (default off)
        -ec, --EEJ_counts       Print out files with EEJ counts for MIC and APR (default off)
	--findSubtracted	Set this flag to start alignment from genome-subtracted
				reads (default off). If enabled, must supply *-e.fq as input
	--trimOnce		Only use first 50bp of reads, if paired, only use 
				50 from fwd and 50 from rev (default off)
	--stepSize i		Trim 50bp every --stepSize (default is 25)
	--preTrimmed		If you are trying to use pre-trimmed fasta/q files 
				(only output from Trim.pl, default off)
        --trim5 i               Skips i nt of read R1 for expression and strand detection (default 0)
	--useFastq		This option is only necessary if you have pre-trimmed reads 
				in fastq not fasta format (default off)
	--resume		Resume a previous run using previous intermediate results

	By default, VAST-TOOLs will try to detect the strandness of the RNA-seq data automatically,
	and reverse-complement reads if necessary such that all map to the forward strand. If you
	specify any of the following arguments, the strandness will not be tested and reads will
	be treated as defined. This may be useful, if the strandness test fails. 	
	--rc1,  --rc2		Reverse-complement reads1 / reads2
	--nrc1, --nrc2		Do not reverse-complement reads1 / reads2
	--ns			RNA-seq reads will be treated as if they were strand-unspecific  

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

verbPrint "Using VASTDB -> $dbDir";
# change directories
make_path($outdir) unless (-e $outdir);
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
unless(-e "$tmpDir/tmp_read_files"){
	mkdir("$tmpDir/tmp_read_files") or die "$!";
}

# Quality control for trim5
if ($trimLen){
    if ($length-$trim5 < $trimLen){
	errPrint "Trim5 ($trim5) cannot be applied for reads of length $trimLen\n";
    }
}
else {
    if ($length-$trim5 < 50){
	errPrint "Trim5 ($trim5) cannot be applied for reads of length $length\n";
    }
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

my $fh_info;
my $bt_norc="";  # Bowtie option:will map only to fwd strand if set to --norc 
my $mapcorr_fileswitch="";  # change of file names with mappability correction (needs to change in strand-aware mode)
my $resumed=0;
my @files_to_be_deleted=();
# resume?
if($resume && -e "to_combine/$root".".info" && open($fh_info,"to_combine/$root".".info")){
    my $line=<$fh_info>; chomp($line); close($fh_info); my @fs=split("\t",$line);
    if(@fs[@fs-1] eq "done"){ # done in the end means read check and potential creation of reverse-complement has finished already successfully
	$bt_norc=$fs[@fs-3];
	$mapcorr_fileswitch=$fs[@fs-2];
	$fq1=$fs[@fs-4]; if($fs[0] eq "paired"){ ($fq1,$fq2)=($fs[@fs-5],$fs[@fs-4]); }
	$resumed=1;
    }
}
unless($resumed){
    # create info file for this RNAseq data set
    open($fh_info,">to_combine/$root".".info") or die "$!";
    if($fq2){print $fh_info "paired\t$fq1\t$fq2";}else{print $fh_info "single\t$fq1";}
    if($notstrandaware){
	print $fh_info "\t$species\t--ns (should be treated as strand-unspecific data)";
    }
    else{
	print $fh_info "\t$species\t--s (check if data is strand-specific)";
    }
    
    #### Check if paired-end reads are strand specific. If paired-end reads are strand-specific, all first/second reads get reverse-complemented if the majority of them maps to strand - of mRNA reference sequences.
    if($notstrandaware){
	print $fh_info "\tdue to given argument -ns, data are treated as being not strand-specific NA NA NA NA\t\t$bt_norc\t$mapcorr_fileswitch\tdone";
    }
    else{
	my $minNMappingReads=500;   # at least so many reads from all 10000 reads must get mapped
	my $minThresh=0.35;         # If fraction of reads mapping to strand - is larger than this threshold, we assume the data is indeed strand-specific.
	my $maxThresh=0.65; 
	my ($fh,$fh2);
	sub rvcmplt{ $_=$_[0]; tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy\[\]/TVGHCDKNYSAABWXRtvghcdknysaabwxr\]\[/; return(reverse($_));} 
	
       	my $N=400000; # check 100K fastq reads 
       	my $bowtie_fa_fq_flag="-q";  if($fq1 =~ /fasta$|fasta\.gz$|fa$|fa\.gz$/){$bowtie_fa_fq_flag="-f";$N=200000;}
       	my ($p1,$n1,$p2,$n2)=(0,0,0,0);   # number of reads 1 mapping to strand + and - , number of reads 2 mapping to strand + and -
       	my ($percR1p,$percR1n,$percR2p,$percR2n)=("NA","NA","NA","NA");

       	if(!defined($rc1)){
       		open($fh, "".getPrefixCmd($fq1)." | head -n $N | $bowtie $bowtie_fa_fq_flag -5 $trim5 -p $cores -m 1 -v $bowtieV $dbDir/EXPRESSION/mRNA - | cut -f 2 |") or errPrintDie "$!";  while(<$fh>){chomp;if($_ eq "-"){$n1++}else{$p1++}}; close($fh);
       	}else{ # user has defined if reads1 should be reversed or not
       		if($rc1==1){($p1,$n1)=(0,5000)}else{($p1,$n1)=(5000,0)}
       	}
       	if($p1+$n1<500){die "Too few first reads (<500) could be mapped to mRNA library for detecting strand-orientation of reads. Maybe wrong species selected?";}
       	($percR1p,$percR1n)=( ($p1/($p1+$n1)),($n1/($p1+$n1)) );
		$percR1p = sprintf("%.4f",$percR1p);
		$percR1n = sprintf("%.4f",$percR1n);
		verbPrint "   fraction of first reads mapping to fwd / rev strand : $percR1p / $percR1n";

       	if($pairedEnd){  # same for second reads of each pair
       		if(!defined($rc2)){	
       			open($fh, "".getPrefixCmd($fq2)." | head -n $N | $bowtie $bowtie_fa_fq_flag -5 $trim5 -p $cores -m 1 -v $bowtieV $dbDir/EXPRESSION/mRNA - | cut -f 2 |") or errPrintDie "$!";  while(<$fh>){chomp;if($_ eq "-"){$n2++}else{$p2++}}; close($fh);
       		}else{ # user has defined if reads1 should be reversed or not
       			if($rc2==1){($p2,$n2)=(0,5000)}else{($p2,$n2)=(5000,0)}
       		}
        	if($p2+$n2<500){die "Too few second reads (<500) could be mapped to mRNA library for detecting strand-orientation of reads. Maybe wrong species selected?";}
       		($percR2p,$percR2n)=( ($p2/($p2+$n2)),($n2/($p2+$n2)) );
		$percR2p = sprintf("%.4f",$percR2p);
		$percR2n = sprintf("%.4f",$percR2n);
       		verbPrint "   fraction of second reads mapping to fwd / rev strand : $percR2p / $percR2n";
        }

		if(($percR2n eq "NA" && $percR1n>=$minThresh && $percR1n<=$maxThresh) || ($percR1n>=$minThresh && $percR1n<=$maxThresh && $percR2n>=$minThresh && $percR2n<=$maxThresh)){
			print $fh_info "\tdata assumed to be not strand-specific $percR1p $percR1n $percR2p $percR2n";
			$notstrandaware=1;
		}else{
			my $fn;
			print $fh_info "\tdata assumed to be strand-specific $percR1p $percR1n $percR2p $percR2n";
			for(my $i=0;$i<2;$i++){
				if($i==0){
					unless($percR1n>=$maxThresh){print $fh_info "\t$fq1";next;}
					open($fh,"".getPrefixCmd($fq1)." |");
					$fn="$tmpDir/tmp_read_files/".pop(@{[split("/",$fq1)]});
					if(isZipped($fq1)){open($fh2,"| gzip -c > $fn" ) or die "$!";}else{open($fh2,">$fn") or die "$!";}
					verbPrint "   reverse-complementing reads from $fq1; writing into $fn";
					$fq1=$fn;
					push(@files_to_be_deleted,$fn);
				}
				if($i==1){
					if($percR2n eq "NA"){next;}  # single-end data
					unless($percR2n>=$maxThresh){print $fh_info "\t$fq2";next;}
					open($fh,"".getPrefixCmd($fq2)." |");
					$fn="$tmpDir/tmp_read_files/".pop(@{[split("/",$fq2)]});
					if(isZipped($fq2)){open($fh2,"| gzip -c > $fn" ) or die "$!";}else{open($fh2,">$fn") or die "$!";}
					verbPrint "   reverse-complementing reads from $fq2; writing into $fn";
					$fq2=$fn;
					push(@files_to_be_deleted,$fn);
				}
				print $fh_info "\t$fn";
				my $c=0; while(<$fh>){chomp;my $l=$_;$c++;
					if($c==1){print $fh2 "$l\n";}
					if($c==2){print $fh2 rvcmplt($l)."\n";  if($bowtie_fa_fq_flag eq "-f"){$c=0;}}
					if($c==3){print $fh2 "$l\n";}
					if($c==4){print $fh2 reverse($l)."\n"; $c=0;}
				}
				close($fh);close($fh2);
			} # for-i
			$bt_norc="--norc";          # set Bowtie argument --norc for strandaware mode
			$mapcorr_fileswitch="-SS";  # ending of files with mappability correction factors for strand-AWARE mode
		}
		print $fh_info "\t$bt_norc\t$mapcorr_fileswitch\tdone";
	}
	close($fh_info);
}


if (!$genome_sub and !$useGenSub){
 my $cmd;
#### Expression analysis (it maps only the first $le nucleotides of the read)
 if ($runExprFlag || $onlyExprFlag){
     verbPrint "Mapping RNAseq reads against mRNA sequences";

     $cmd = "$fq1";
     $cmd = getPrefixCmd($cmd);
     my $bowtie_fa_fq_flag="-q";
     if($fq1 =~ /fasta$|fasta\.gz$|fa$|fa\.gz$/){$bowtie_fa_fq_flag="-f";}

#    Different $cmd (24/12/16 --MI)
     if (defined($trimLen)){
	 $cmd .= " | $binPath/Trim.pl --once --targetLen $trimLen -v | $bowtie $bt_norc $bowtie_fa_fq_flag -5 $trim5 -p $cores -m 1 -v $bowtieV $dbDir/EXPRESSION/mRNA -"; 
     }
     else {     	
	 $cmd .= " | $binPath/Trim.pl --once --targetLen 50 -v | $bowtie $bt_norc $bowtie_fa_fq_flag -5 $trim5 -p $cores -m 1 -v $bowtieV $dbDir/EXPRESSION/mRNA -"; 
     }
     
     verbPrint "Calculating cRPKMs\n";
     checkResumeOption("$root-$le.fq.gz","$root-$le.fa.gz","to_combine/$root.eej2","to_combine/$root.IR.summary.txt","to_combine/$root.IR.summary_v2.txt");
     sysErrMsg "$cmd | $binPath/expr_RPKM.pl - $dbDir/EXPRESSION/$species"."_mRNA-${le}${mapcorr_fileswitch}.eff expr_out/$root > expr_out/$root\.cRPKM";
 }
 if ($onlyExprFlag){
     print STDERR "Expression analysis done\n";
     exit 0;
 }
}
###

#### Merge PE => happens in Trim
$fq = $fq1; # just assigning name

#### Trimming
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
    
    verbPrint "Trimming RNAseq reads to $le nt sequences";
    checkResumeOption("to_combine/$root.eej2","to_combine/$root.IR.summary.txt","to_combine/$root.IR.summary_v2.txt");
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
unless ($onlyIRflag){
    verbPrint "Doing genome subtraction\n";
    # Force bash shell to support process substitution
    $cmd = getPrefixCmd($fq);
    $cmd .= " | $bowtie -p $cores $inpType -m 1 -v 2 --un >(gzip > $subtractedFq) --max /dev/null $dbDir/FILES/gDNA - /dev/null";
    checkResumeOption("to_combine/$root.eej2","to_combine/$root.IR.summary.txt","to_combine/$root.IR.summary_v2.txt");
    sysErrMsg("bash", "-c", $cmd);
}

####

if ($EXIT_STATUS) {
    exit $EXIT_STATUS;
}

#### Map to the EEJ:
my $runArgs = "-dbDir=$dbDir -sp=$species -readLen=$le -root=$root"; 
     unless ($notstrandaware){$runArgs .=" -s";}
     if ($print_EEJs){$runArgs.=" -ec";}         
my $preCmd = getPrefixCmd($subtractedFq);
unless ($onlyIRflag){
    verbPrint "Mapping reads to the \"splice site-based\" (aka \"a posteriori\") EEJ library and Analyzing...\n";
    checkResumeOption("to_combine/$root.exskX");
    sysErrMsg "$preCmd | $bowtie $bt_norc $inpType -p $cores -m 1 -v $bowtieV " .
	"$dbDir/FILES/$species"."_COMBI-M-$le - | " .
	"cut -f 1-4,8 - | sort -T $tmpDir -k 1,1 | " .
	"$binPath/Analyze_COMBI.pl deprecated " .
	"$dbDir/COMBI/$species/$species"."_COMBI-M-$le-gDNA${mapcorr_fileswitch}.eff $runArgs";   # produces to_combine/$root.eej2
    
    verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") SIMPLE EEJ library and Analyzing...\n";
    checkResumeOption("to_combine/$root.MULTI3X");
    sysErrMsg "$preCmd | $bowtie $bt_norc $inpType -p $cores -m 1 -v $bowtieV " .
	"$dbDir/FILES/EXSK-$le - | " .
	"cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
	"$binPath/Analyze_EXSK.pl $runArgs";                                 # produces to_combine/$root.exskX
    
    verbPrint "Mapping reads to the \"transcript-based\" (aka \"a priori\") MULTI EEJ library and Analyzing...\n";
    checkResumeOption("to_combine/$root.micX");
    sysErrMsg "$preCmd | $bowtie $bt_norc $inpType -p $cores -m 1 -v $bowtieV " .
	"$dbDir/FILES/MULTI-$le - | " .
	"cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
	"$binPath/Analyze_MULTI.pl $runArgs";                                # produces to_combine/$root.MULTI3X

    verbPrint "Mapping reads to microexon EEJ library and Analyzing...\n";
    checkResumeOption("to_combine/$root.IR.summary.txt","to_combine/$root.IR.summary_v2.txt","$tmpDir/".pop(@{[split("/",$fq1)]}).".resume");
    sysErrMsg "$preCmd | $bowtie $bt_norc $inpType -p $cores -m 1 -v $bowtieV " .
	"$dbDir/FILES/$species"."_MIC-$le - | ".
	" cut -f 1-4,8 - | sort -T $tmpDir -k 1,1 | " .
	" $binPath/Analyze_MIC.pl $runArgs";                                 # produces to_combine/$root.micX
}

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
  checkResumeOption("to_combine/$root.IR","to_combine/$root.IR2");
  sysErrMsg "$preCmd | $bowtie $bt_norc $inpType -p $cores -m 1 -v $bowtieV " .
              "$dbDir/FILES/$species.IntronJunctions.$type.$le.8 - | " .
              "cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
              "$binPath/MakeSummarySAM.pl | " .
              "$binPath/RI_summarize$v.pl - $runArgs";                       # produces to_combine/$root.IR.summary.txt or to_combine/$root.IR.summary_v2.txt
  checkResumeOption("$tmpDir/".pop(@{[split("/",$fq1)]}).".resume");
  sysErrMsg "$preCmd | $bowtie $bt_norc $inpType -p $cores -m 1 -v $bowtieV " .
                  "$dbDir/FILES/$species.Introns.sample.200 - | " .
              "cut -f 1-4,8 | sort -T $tmpDir -k 1,1 | " .
              "$binPath/MakeSummarySAM.pl | " .
              "$binPath/RI_summarize_introns$v.pl - $runArgs";               # produces /to_combine/$root.IR or /to_combine/$root.IR2 
} else {
  verbPrint "Skipping intron retention step...\n";
}

unless($keepFlag or $keep_trimmed) {
  verbPrint "Cleaning $fq files!";
  sysErrMsg "rm -f $fq";
}

unless($keepFlag) {
    unless ($onlyIRflag){
	verbPrint "Cleaning up $subtractedFq!";
	sysErrMsg "rm -f $subtractedFq";
    }
}

unless($noIRflag || $IR_version == 2) {  # --UB
    my $juncAnnotationFile = "./to_combine/$root.IR.summary.txt";
    verbPrint "Cleaning up $juncAnnotationFile!";
    sysErrMsg "rm -f $juncAnnotationFile";
}

# Generate a control file for resume option. Allows us a complete resume if previous run did complete successfully.
open(my $fh,">$tmpDir/".pop(@{[split("/",$fq1)]}).".resume");print $fh "finished successfully";close($fh);

# remove files with reverse-complemented reads
verbPrint "Deleting temporary files with reverse-complemented reads.";
foreach my $file (@files_to_be_deleted){
	sysErrMsg "rm $file";
}

verbPrint "Completed " . localtime;
exit $EXIT_STATUS;
