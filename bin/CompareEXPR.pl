#!/usr/bin/perl
### General script to get differentially spliced events based on dPSI differences 

## LOG:
# 03/08/16: added fold and mean cRPKM to output table (local)
# 29/10/16: outputs a BG list with fold changes (for enrichment analysis)
# 31/12/16: added option to normalize cRPKMs using 'normalizebetweenarrays' from limma

use Getopt::Long;
use warnings;
use strict;
use Cwd qw(abs_path);
use List::Util qw(sum);

# INITIALIZE PATH AND FLAGS
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;


### Setting global variables:
my $helpFlag = 0;
my $verboseFlag = 1; 
my $dbDir; # directory of VASTDB
my $species; # to get the GO conversion file
my $Q = "O[KW]\,.+?\,.+?\,.+?\,.+?\@"; # quality search
my $input_file = $ARGV[0];
#my $output_file; => #changed for out_root
my $min_fold_av = 2; # min fold change
my $min_fold_r = 1.5; # min fold diff between max and min of each group
my $min_fold_log_av;
my $max_CV = 100; # so far too high
my $min_reads = 50; # min number of reads in 1 sample
my $min_cRPKM = 2; # min cRPKM in all samples from at least 1 group
my $min_cRPKM_loose;
my $samplesA;
my $samplesB;
my @samplesA;
my @samplesB;
my $repA; # number of replicates per type
my $repB; # number of replicates per type
my $get_GO;
my $paired;
my $use_names;
my $folder;
my $no_plot = 1;
my $plot_only_samples;
my $print_all;
my $normalize; # does quantile normalization (using limma)
my $install_limma;
my $out_root;
my @data_fold;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "min_fold_av=f" => \$min_fold_av,
			  "min_fold_r=f" => \$min_fold_r,
                          "max_CV=f" => \$max_CV,
                          "min_reads=i" => \$min_reads,
                          "min_cRPKM=f" => \$min_cRPKM,
			  "min_cRPKM_loose" => \$min_cRPKM_loose,
			  "a=s" => \$samplesA,
			  "samplesA=s" => \$samplesA,
			  "b=s" => \$samplesB,
			  "samplesB=s" => \$samplesB,
			  "outRoot=s" => \$out_root,
			  "output=s" => \$folder,
			  "o=s" => \$folder,
			  "help" => \$helpFlag,
			  "norm" => \$normalize,
			  "install_limma" => \$install_limma,
			  "GO" => \$get_GO,
			  "use_names" => \$use_names,
			  "print_all" => \$print_all,
			  "paired" => \$paired
    );

our $EXIT_STATUS = 0;

sub errPrint {
    my $errMsg = shift;
    print STDERR "[vast compare error]: $errMsg\n";
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
	print STDERR "[vast compare_expr]: $verbMsg\n";
    }
}

sub time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
    $year += 1900;
    $mon += 1;
    my $datetime = sprintf "%04d-%02d-%02d (%02d:%02d)", $year, $mday, $mon, $hour, $min;
    return $datetime;
}

### Gets the version
my $version;
open (VERSION, "$binPath/../VERSION");
$version=<VERSION>;
chomp($version);
$version="No version found" if !$version;

if (!defined($ARGV[0]) || $helpFlag){
    die "
VAST-TOOLS v$version

Usage: vast-tools compare_expr cRPKMS_AND_COUNTS-SpN.tab -a sample_a1,sample_a2 -b sample_b1,sample_b2 [options]

Compare two sample sets to find differentially expressed genes based on fold changes of cRPKM values

[General options] 
        --min_fold_av f          Minimum fold change of the averages (default 2)
        --min_fold_r f           Minimum fold change between max and min of the two groups (default 1.5)
        --min_cRPKM f            Minimum expression value (cRPKM) in all samples for at least one group (default 2)
        --min_cRPKM_loose        min_cRPKM applies to only one sample, not the whole group (default OFF [strict])
        --min_reads i            Minimum number of raw reads in at least one sample (default 50) 
        --outRoot root           Output file root (default based on option parameters)
        --norm                   Normalize cRPKMs using \'normalizeBetweenArrays\' from limma (default OFF, recommended)
        --install_limma          Installs limma package if needed (default OFF)
        -a/--samplesA sA1,sA2    Required, 1:n sample names or column_\# separated by , (mandatory)
        -b/--samplesB sB1,sB2    Required, 1:n sample names or column_\# separated by , (mandatory)
        --print_all              Print all samples (default OFF\; prints only the tested samples) 
        --paired                 Does a paired comparison (A1 vs B1, A2 vs B2, etc.)
                                   - It uses min_fold_av as the minimum fold change for the averages of individual fold changes
                                   - It uses min_fold_r as the minumum fold change for each paired comparison 
    
[GO options]
        --GO                     Generates gene lists for GO analysis (from default gene_ID)
        --use_names              Uses gene names (first column in INCLUSION table)


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

#### SANITY CHECKS
errPrintDie "Need to provide at least one replicate for each group\n" if (!defined $samplesA or !defined $samplesB);
# make the array with samples
@samplesA = split(/\,/,$samplesA);
@samplesB = split(/\,/,$samplesB);
$repA = $#samplesA+1;
$repB = $#samplesB+1;
# for paired (repA must be the same as repB)
errPrintDie "If paired comparison, the number of replicates must be the same\n" if (defined $paired) && ($repA != $repB);


### To get the folder in which the input file is
($folder) = $input_file =~/(.+)\//; # empty if no match (i.e. local folder)
$folder = "." unless (defined $folder);


open (GE, $input_file) or errPrintDie "Needs a cRPKM + COUNTs table $!\n";

# prints version (05/05/19)
verbPrint "VAST-TOOLS v$version";

### Creates the LOG
open (LOG, ">>$folder/VTS_LOG_commands.txt"); # || die "Cannot open the LOG file";
my $all_args="-a $samplesA -b $samplesB -min_fold_av $min_fold_av -min_fold_r $min_fold_r -min_reads $min_reads -min_cRPKM $min_cRPKM";
$all_args.=" -paired" if $paired;
$all_args.=" -min_cRPKM_loose" if $min_cRPKM_loose;
$all_args.=" -norm" if $normalize;
$all_args.=" -GO" if $get_GO;
$all_args.=" -use_names" if $use_names;
$all_args.=" -print_all" if $print_all;
$all_args.=" -outRoot" if $out_root;

print LOG "[VAST-TOOLS v$version, ".&time."] vast-tools compare_expr $all_args\n";

### Common for all numbers of replicates
# preparing the head
my $head_row=<GE>;
chomp($head_row);
$head_row =~ s/\-R\t/\t/g; # to remove the "-R" that says it's a cRPKM column
$head_row =~ s/\-cRPKM\t/\t/g; # to remove the "-R" that says it's a cRPKM column
my @head=split(/\t/,$head_row);
foreach my $i (2..$#head){
    if ($i%2==0){ # to match sample names with column number
	foreach my $j (0..$#samplesA){
	    $samplesA[$j] = $i if $samplesA[$j] eq $head[$i];
	}
	foreach my $j (0..$#samplesB){
	    $samplesB[$j] = $i if $samplesB[$j] eq $head[$i];
	}
    }
}

### Make normalized table
my %norm_cRPKMs;
if (defined $normalize){
    open (GE_2, $input_file) or errPrintDie "Needs a cRPKM + COUNTs table\n";
    my ($input_path,$root_input);
    if ($input_file=~/\//){
	($input_path,$root_input) = $input_file =~/(.+)\/(.+)\./;
    }
    else {
	$input_path=".";
	($root_input) = $input_file =~/(.+)\./;
    }
    $root_input=~s/cRPKM_AND_COUNTS/cRPKM/;
    
    open (TEMP, ">$input_path/temp_cRPKMs.tab");
    while (<GE_2>){
	chomp($_);
	my @t=split(/\t/,$_);
	print TEMP "$t[0]";
	foreach my $i (2..$#t){
	    if ($i%2==0){
		print TEMP "\t$t[$i]";
	    }
	}
	print TEMP "\n";
    }
    close TEMP;
    close GE_2;

    open (Temp_R, ">$input_path/temp.R");
    print Temp_R "source(\"https://bioconductor.org/biocLite.R\")
biocLite(\"limma\")\n" if (defined $install_limma);

    print Temp_R "
library(limma)
setwd(\"$input_path/\")
matrix=as.matrix(read.table(\"temp_cRPKMs.tab\", row.names=1, header=TRUE,sep=\"\\t\"))
Nmatrix=normalizeBetweenArrays(as.matrix(matrix))
NmatrixF=cbind(Names=row.names(matrix),Nmatrix)
write.table(NmatrixF,\"$root_input-NORM.tab\",
            sep=\"\\t\",col.names=T,row.names=F,quote=F)";
    close Temp_R;

    system "Rscript $input_path/temp.R";
    system "rm $input_path/temp*";

    open (GE_NORM, "$input_path/$root_input-NORM.tab") || die "Can't open the normalized cRPKMs";
    <GE_NORM>;
    while (<GE_NORM>){
	chomp($_);
	my @t = split(/\t/,$_);
	
	for my $i (1..$#t){
	    # keeps the normalized value as if there were raw counts too
	    $norm_cRPKMs{$t[0]}[$i*2]=sprintf("%.2f",$t[$i]) if $t[$i] ne "NA"; 
	    $norm_cRPKMs{$t[0]}[$i*2]=$t[$i] if $t[$i] eq "NA"; 
	}
    }
    close GE_NORM;
}

# check that columns provided are 0-based OR
# if names were provided, that all columns were properly matched
my $kill_0based;
my $kill_2lower;
foreach my $s (@samplesA){
    $kill_0based = 1 if $s%2 != 0;
    $kill_2lower = 1 if $s < 2;
}
foreach my $s (@samplesB){
    $kill_0based = 1 if $s%2 != 0;
    $kill_2lower = 1 if $s < 2;
}
errPrintDie "Column numbers do not seem 0-based or conversion did not work properly\n" if (defined $kill_0based);
errPrintDie "Column numbers do not seem to correspond to cRPKM samples\n" if (defined $kill_2lower);

my $short_head = "GENE\tNAME";
foreach my $j (0..$#samplesA){
    $short_head.= "\t$head[$samplesA[$j]]";
}
foreach my $j (0..$#samplesB){
    $short_head.= "\t$head[$samplesB[$j]]";
}
$short_head.= "\tCV_A\tCV_B\tAv_A\tAv_B\tLog2_Fold_Ch";

# representative names
my $name_A=$head[$samplesA[0]];
my $name_B=$head[$samplesB[0]];
$name_A=~s/(.+)\_.+/$1/ unless $repA == 1; # usually the rep number/id is encoded as "_a"
$name_B=~s/(.+)\_.+/$1/ unless $repA == 1;

# defining default output file name
my ($root)=$ARGV[0]=~/.+\-(.+?)\./;
my $tail = ""; # to be added to the output name
$tail.="_$name_A-vs-$name_B";
$tail.="-minReads$min_reads";
$tail.="-minRPKM$min_cRPKM";
$tail.="-loose" if (defined $min_cRPKM_loose);
$tail.="-strict" if (!defined $min_cRPKM_loose);
$tail.="-range$min_fold_r" if (defined $min_fold_r); 
$tail.="-paired" if (defined $paired);
$tail.="-NORM" if (defined $normalize);
$out_root="$root-fold$min_fold_av$tail" unless (defined $out_root);

my $output_file="DiffGE-$out_root.tab";

open (O, ">$folder/$output_file") or errPrintDie "Can't open the output file (do not provide a path)\n"; # output file

print O "$short_head\n" unless (defined $print_all); 
print O "$head_row\tLog2_Fold_Ch\n" if (defined $print_all);


# prepare to obtain gene IDs for GO analyses
my %ID_gene;
if (defined $get_GO){
    open (BG, ">$folder/GE_BG-$out_root.txt") or errPrintDie "Can't open GO output files";
    open (BG_FOLD, ">$folder/GE_BG_FOLD-$out_root.rnk.txt") or errPrintDie "Can't open GO output files";
    open (UP, ">$folder/GE_UP-$out_root.txt") or errPrintDie "Can't open GO output files";
    open (DOWN, ">$folder/GE_DOWN-$out_root.txt") or errPrintDie "Can't open GO output files";
    print BG_FOLD "Feature Name\tScore\n";
}

# Global variables for PSI analysis & GO
my %tally;
$tally{DOWN}=0; $tally{UP}=0;
$min_fold_log_av=log($min_fold_av)/log(2);

verbPrint "Doing comparisons of GE profiles ($name_A vs $name_B)\n";

while (<GE>){
    chomp($_);
    my @t=split(/\t/,$_);
    my @GE_A = ();
    my @GE_B = ();
    my $fold_BG;

    next if ($t[2] eq "NA" || $t[3] eq "NA"); # removes the lines with NA

    # coverage check (requires good coverage for ALL replicates)
    my $kill_coverage = 1;
#    my $kill_RPKM = 1; # no longer just one sample, but the whole group (30/03/16 --MI)
    foreach my $s (@samplesA){
	$t[$s+1] = 0 if ($t[$s+1] eq "");
	$t[$s+1] = 0 if ($t[$s+1] eq "NA");
	$kill_coverage = 0 if ($t[$s+1] >= $min_reads); # kill if NONE of the samples meets the coverage criteria
#	$kill_RPKM = 0 if ($t[$s] >= $min_cRPKM);

	### normalization added 31/12/16
	push(@GE_A,$t[$s]) if (!defined $normalize);
	push(@GE_A,$norm_cRPKMs{$t[0]}[$s]) if (defined $normalize);
    }
    foreach my $s (@samplesB){
	$t[$s+1] = 0 if ($t[$s+1] eq "");
	$t[$s+1] = 0 if ($t[$s+1] eq "NA");
	$kill_coverage = 0 if ($t[$s+1] >= $min_reads); # kill if NONE of the samples meets the coverage criteria
#	$kill_RPKM = 0 if ($t[$s] >= $min_cRPKM);

	### normalization added 31/12/16
	push(@GE_B,$t[$s]) if (!defined $normalize);
	push(@GE_B,$norm_cRPKMs{$t[0]}[$s]) if (defined $normalize);
    }
    next if ($kill_coverage == 1);
#    next if ($kill_RPKM == 1);

    ### Applicable to paired and non-paired (30/03/16) --MI
    # get RPKMs
    my $sum_A = sum(@GE_A);
    my $sum_B = sum(@GE_B);
    my $av_GE_A = sprintf("%.2f",$sum_A/$repA);
    my $av_GE_B = sprintf("%.2f",$sum_B/$repB);
    # get min and max (for ranges)
    my $min_A = (sort{$b<=>$a}@GE_A)[-1];
    my $max_A = (sort{$a<=>$b}@GE_A)[-1];
    my $min_B = (sort{$b<=>$a}@GE_B)[-1];
    my $max_B = (sort{$a<=>$b}@GE_B)[-1];

    next if ($max_A < $min_cRPKM && $max_B < $min_cRPKM) && defined ($min_cRPKM_loose); # to implement the new min expr cut-off (in all samples for 1 group)
    next if ($min_A < $min_cRPKM && $min_B < $min_cRPKM) && !defined ($min_cRPKM_loose);

    # NOT PAIRED: gets the average cRPKMs for A and B and the lowest (min) and highest (max) cRPKM for each replicate
    if (!defined $paired){
	# get fold change
	my $fold_change = ($av_GE_B+0.001)/($av_GE_A+0.001);
	my $log2_fold_change = sprintf("%.2f",log($fold_change)/log(2));
	$fold_BG = $log2_fold_change;

	# get coef variation
	my $std_dev_A = &std_dev(@GE_A);
	my $std_dev_B = &std_dev(@GE_B);
	my $CV_A = sprintf("%.2f",($std_dev_A+0.001)/($av_GE_A+0.001));
	my $CV_B = sprintf("%.2f",($std_dev_B+0.001)/($av_GE_B+0.001));

	# gets what to print
	# normalization added (31/12/16)
	my $short_line = "$t[0]\t$t[1]";
	foreach my $j (0..$#samplesA){
	    $short_line.= "\t$t[$samplesA[$j]]" if (!defined $normalize);
	    $short_line.= "\t$norm_cRPKMs{$t[0]}[$samplesA[$j]]" if (defined $normalize);
	}
	foreach my $j (0..$#samplesB){
	    $short_line.= "\t$t[$samplesB[$j]]" if (!defined $normalize);
	    $short_line.= "\t$norm_cRPKMs{$t[0]}[$samplesB[$j]]" if (defined $normalize);
	}
	$short_line.= "\t$CV_A\t$CV_B\t$av_GE_A\t$av_GE_B\t$log2_fold_change";
	
#	next if ($CV_A > 2**abs($log2_fold_change)/4 || $CV_B > 2**abs($log2_fold_change)/4); # at this point, 1/4 of the fold
	
	# does the diff GE test:
	if ($log2_fold_change > $min_fold_log_av){
	    # get fold change of ranges
	    my $fold_change_r = ($min_B+0.001)/($max_A+0.001);
	    next if $fold_change_r < $min_fold_r;
	    
	    $tally{UP}++;
	    print O "$short_line\n" unless (defined $print_all);
	    print O "$_\t$log2_fold_change\n" if (defined $print_all); 
	    
	    # print for GO
	    if (defined$get_GO){
		unless ($use_names){
		    print UP "$t[0]\n";
		}
		else {
		    print UP "$t[1]\n";
		}
	    }
	}
	if ($log2_fold_change < -1*$min_fold_log_av){
	    # get fold change of ranges
	    my $fold_change_r = ($min_A+0.001)/($max_B+0.001);
	    next if $fold_change_r < $min_fold_r;

	    $tally{DOWN}++;
	    print O "$short_line\n" unless (defined $print_all);
	    print O "$_\t$log2_fold_change\n" if (defined $print_all); 
	    
	    #print for GO
	    if (defined$get_GO){
		unless ($use_names){
		    print DOWN "$t[0]\n";
		}
		else {
		    print DOWN "$t[1]\n";
		}
	    }
	}
    }
    else { # if paired: calculates each pair's dPSI & the lower/higher of these
	# get each dPSI and the average
	my @GE_pairs = ();
	for my $i (0..$#samplesA){ # ~ $repA-1
	    $GE_pairs[$i] = sprintf("%.2f",log(($GE_B[$i]+0.001)/($GE_A[$i]+0.001))/log(2));
	}
	my $sum_paired_GE_fold = sum(@GE_pairs);
	my $av_paired_GE_fold = sprintf("%.2f",$sum_paired_GE_fold/$repA);
	$fold_BG = $av_paired_GE_fold;

	# get min and max (for min dPSI, pos and neg values)
	my $min_indiv_GE_fold = (sort{$b<=>$a}@GE_pairs)[-1];
	my $max_indiv_GE_fold = (sort{$a<=>$b}@GE_pairs)[-1];

	# get coef variation
	my $std_dev_A = &std_dev(@GE_A);
	my $std_dev_B = &std_dev(@GE_B);
	my $CV_A = sprintf("%.2f",($std_dev_A+0.001)/($av_GE_A+0.001));
	my $CV_B = sprintf("%.2f",($std_dev_B+0.001)/($av_GE_B+0.001));

	# gets what to print
	# normalization added (31/12/16)
	my $short_line = "$t[0]\t$t[1]";
	foreach my $j (0..$#samplesA){
	    $short_line.= "\t$t[$samplesA[$j]]" if (!defined $normalize);
	    $short_line.= "\t$norm_cRPKMs{$t[0]}[$samplesA[$j]]" if (defined $normalize);
	}
	foreach my $j (0..$#samplesB){
	    $short_line.= "\t$t[$samplesB[$j]]" if (!defined $normalize);
	    $short_line.= "\t$norm_cRPKMs{$t[0]}[$samplesB[$j]]" if (defined $normalize);
	}
	$short_line.= "\t$CV_A\t$CV_B\t$av_GE_A\t$av_GE_B\t$av_paired_GE_fold";
	
	if ($av_paired_GE_fold > $min_fold_log_av && $min_indiv_GE_fold > log($min_fold_r)/log(2)){ 
	    $tally{UP}++;
	    print O "$short_line\n" unless (defined $print_all);
	    print O "$_\t$av_paired_GE_fold\n" if (defined $print_all); 
	    
	    # print for GO
	    if (defined $get_GO){
		unless ($use_names){
		    print UP "$t[0]\n";
		}
		else {
		    print UP "$t[1]\n";
		}
	    }
	}
	if ($av_paired_GE_fold < -$min_fold_log_av && $max_indiv_GE_fold < -log($min_fold_r)/log(2)){ 
	    $tally{DOWN}++;
	    print O "$short_line\n" unless (defined $print_all);
	    print O "$_\t$av_paired_GE_fold\n" if (defined $print_all); 
	    
	    #print for GO
	    if (defined $get_GO){
		unless ($use_names){
		    print DOWN "$t[0]\n";
		}
		else {
		    print DOWN "$t[1]\n";
		}
	    }
	}
    }
    # prints out the genes for BG
    if (defined $get_GO){
	unless ($use_names){
	    print BG "$t[0]\n";
#	    print BG_FOLD "$t[0]\t$fold_BG\n" if (defined $fold_BG);
	    push(@data_fold, "$fold_BG=$t[0]") if (defined $fold_BG);
	}
	else {
	    print BG "$t[1]\n";
#	    print BG_FOLD "$t[1]\t$fold_BG\n" if (defined $fold_BG);
	    push(@data_fold, "$fold_BG=$t[1]") if (defined $fold_BG);
	}
    }
}
close GE;
close O;

if (defined $get_GO){
    verbPrint "Preparing files for GO analysis\n";

    no warnings;
    @data_fold=sort{$b<=>$a}(@data_fold);
    foreach my $temp (@data_fold){
	my ($fold,$g)=split(/\=/,$temp);
	print BG_FOLD "$g\t$fold\n";
    }
    use warnings;
    sleep(1);
    close BG;
    close BG_FOLD;
    close DOWN;
    close UP;
}

unless (defined $no_plot){
    verbPrint "Plotting differentially spliced AS events\n";
    my $config_file = $output_file;
    $config_file =~ s/\..+$//; # ~ getting a root
    $config_file.= ".config.txt";
    open (CONFIG, ">$folder/$config_file");
    
    print CONFIG "Order\tSampleName\tGroupName\tRColorCode\n";
    my $order = 0;
    my %column_seen;
    foreach my $i (@samplesA){
	$order++;
	$head[$i] = "X$head[$i]" if $head[$i]=~/^\d/;
	print CONFIG "$order\t$head[$i]\t$name_A\tred\n";
	$column_seen{$i} = 1;
    }
    foreach my $i (@samplesB){
	$order++;
	$head[$i] = "X$head[$i]" if $head[$i]=~/^\d/;
	print CONFIG "$order\t$head[$i]\t$name_B\tblue\n";
	$column_seen{$i} = 1;
    }
    unless (defined $plot_only_samples){
	for my $i (6..$#head){
	    if ($i%2 == 0){
		if (!defined $column_seen{$i}){
		    $order++;
		    $head[$i] = "X$head[$i]" if $head[$i]=~/^\d/;
		    print CONFIG "$order\t$head[$i]\tOthers\tblack\n";
		}
	    }
	}
    }
    close CONFIG;

    # defining default width (default ~3 OK for 6 events)
    my $mult = 1.3;
    my $width = sprintf("%.1f",$mult*(($order/2)+2));    
    my $heigth = sprintf("%.1f", $mult*(3.2));
    system "$binPath/../R/psiplotter.R $folder/$output_file -c $folder/$config_file -W $width -H $heigth -u TRUE";
}

verbPrint "Printing summary statistics\n";
my $extras = "";
$extras.=" (loose)" if (defined $min_cRPKM_loose);
$extras.=" (strict)" if (!defined $min_cRPKM_loose);
$extras.=", paired" if (defined $paired);
$extras.=", normalized" if (defined $normalize);

print "\n*** Options: Fold_change_average=$min_fold_av, Fold_change_range=$min_fold_r, Min_reads=$min_reads, Min_cRPKM=$min_cRPKM$extras\n";
print "*** Summary statistics:\n";
print "\tHigher_in_$name_A\tHigher_in_$name_B\n";
print "\t$tally{DOWN}\t$tally{UP}\n";
print "\n";


########################
sub average{
    my @data = @_;

    if ($#data==0) {
        die ("Empty array");
    }
    my $total = 0;
    foreach (@data) {
        $total += $_;
    }
    my $average = sprintf("%.2f",$total / ($#data+1));
    return $average;
}
sub std_dev{
    my @data = @_;
    if($#data == 0){
        return 0;
    }
    my $average = &average(@data);
    my $sqtotal = 0;
    foreach(@data) {
        $sqtotal += ($average-$_) ** 2;
    }
    my $std = sprintf("%.2f",($sqtotal / ($#data)) ** 0.5); 
    return $std;
}
