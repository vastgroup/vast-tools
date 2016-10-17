#!/usr/bin/perl
### General script to get differentially spliced events based on dPSI differences 

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
my $output_file;
my $min_dPSI = 15; # min dPSI difference (def=15)
my $samplesA;
my $samplesB;
my @samplesA;
my @samplesB;
my $repA; # number of replicates per type
my $repB; # number of replicates per type
my $min_range = 5; # min dPSI between ranges
my $noVLOW;
my $p_IR;
my $ID_file;
my $get_GO;
my $paired;
my $use_names;
my $folder;
my $no_plot;
my $plot_only_samples;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "min_dPSI=i" => \$min_dPSI,
                          "min_range=i" => \$min_range,
			  "a=s" => \$samplesA,
			  "samplesA=s" => \$samplesA,
			  "b=s" => \$samplesB,
			  "samplesB=s" => \$samplesB,
			  "outFile=s" => \$output_file,
			  "output=s" => \$folder,
			  "o=s" => \$folder,
			  "help" => \$helpFlag,
			  "p_IR" => \$p_IR,
			  "GO" => \$get_GO,
			  "sp=s" => \$species,
			  "species=s" => \$species,
			  "dbDir=s" => \$dbDir,
			  "GO_file=s" => \$ID_file,
			  "use_names" => \$use_names,
			  "paired" => \$paired,
			  "no_plot" => \$no_plot,
			  "only_samples" => \$plot_only_samples,
			  "noVLOW" => \$noVLOW
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
	print STDERR "[vast compare]: $verbMsg\n";
    }
}

# Check database directory and set up ID file
if (defined $get_GO){
    errPrintDie "Needs to provide species (\-\-sp) OR file with gene ID conversions (--GO_file) OR activate the --use_names flag\n" if (!defined $species && !defined $ID_file && !defined $use_names);
    
    if (!defined $ID_file && !defined $use_names){ # will use default file
	unless(defined($dbDir)){ # always the default VASTDB
	    $dbDir = "$binPath/../VASTDB";
	}
	$dbDir = abs_path($dbDir);
	$dbDir .= "/$species";
	errPrintDie "The database directory $dbDir does not exist" unless (-e $dbDir or $helpFlag);
	$ID_file = "$dbDir/FILES/$species.Event-Gene.IDs.txt";
	errPrintDie "The file with gene ID conversions does not exist in VASTDB" unless (-e $ID_file or $helpFlag);
    }
    elsif (defined $ID_file && !defined $use_names) { # will use user provided file
	errPrintDie "The file with gene ID conversions does not exist in VASTDB" unless (-e $ID_file or $helpFlag);
    }
}

if (!defined($ARGV[0]) || $helpFlag){
    die "\nUsage: vast-tools compare INCLUSION_LEVELS_FULL-root.tab -a sample_a1,sample_a2 -b sample_b1,sample_b2 [options]

Compare two sample sets to find differentially regulated AS events

[General options] 
        --min_dPSI i             Minimum delta PSI of the averages (default 15)
        --min_range i            Minimum distance between the ranges of both groups (default 5)
        --outFile file           Output file name (default based on option parameters)
        -a/--samplesA sA1,sA2    Required, 1:n sample names or column_\# separated by , (mandatory)
        -b/--samplesB sB1,sB2    Required, 1:n sample names or column_\# separated by , (mandatory)
        --noVLOW                 Does not use samples with VLOW coverage (default OFF)
        --p_IR                   Filter IR by the p-value of the binomial test (default OFF)
        --no_plot                Does NOT plot the DS events using \'plot\' (default OFF)
        --only_samples           Plots only the compared samples, otherwise the whole table (default OFF)
        --paired                 Does a paired comparison (A1 vs B1, A2 vs B2, etc.)
                                   - It uses min_dPSI as the minimum average of each paired dPSI
                                   - It uses min_range as the minumum dPSI for each paired comparison 
    
[GO options]
        --GO                     Generates gene lists for GO analysis (from default gene_ID)
        --GO_file file           To provide an alternative gene ID file for each event (default OFF)
        --use_names              Uses gene names (first column in INCLUSION table)
        --species/-sp Hsa/etc    Three letter code for the database 
        --dbDir db               Database directory (default VASTDB)


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


#### opens INCLUSION TABLE
open (PSI, $input_file) or errPrintDie "Needs a PSI INCLUSION table\n";

# Common for all numbers of replicates
# preparing the head
my $head_row=<PSI>;
chomp($head_row);
my @head=split(/\t/,$head_row);
foreach my $i (6..$#head){
    if ($i%2==0){ # to match sample names with column number
	foreach my $j (0..$#samplesA){
	    $samplesA[$j] = $i if $samplesA[$j] eq $head[$i];
	}
	foreach my $j (0..$#samplesB){
	    $samplesB[$j] = $i if $samplesB[$j] eq $head[$i];
	}
    }
}
# check that columns provided are 0-based OR
# if names were provided, that all columns were properly matched
my $kill_0based;
my $kill_6lower;
foreach my $s (@samplesA){
    $kill_0based = 1 if $s%2 != 0;
    $kill_6lower = 1 if $s < 6;
}
foreach my $s (@samplesB){
    $kill_0based = 1 if $s%2 != 0;
    $kill_6lower = 1 if $s < 6;
}
errPrintDie "Column numbers do not seem 0-based or conversion did not work properly\n" if (defined $kill_0based);
errPrintDie "Column numbers do not seem to correspond to INCLUSION samples\n" if (defined $kill_6lower);

# gets representative names
my $name_A=$head[$samplesA[0]];
my $name_B=$head[$samplesB[0]];
$name_A=~s/(.+)\_.+/$1/ unless $repA == 1; # usually the rep number/id is encoded as "_a" or "_1", but if it's only one, it's left as is
$name_B=~s/(.+)\_.+/$1/ unless $repB == 1;

####### Output file
# defining default output file name
my ($root)=$ARGV[0]=~/.+\-(.+?)\./;
my $tail = ""; # to be added to the output name
$tail.="-range$min_range" if (defined $min_range); 
$tail.="-noVLOW" if (defined $noVLOW);
$tail.="-p_IR" if (defined $p_IR);
$tail.="-paired" if (defined $paired);
$tail.="_$name_A-vs-$name_B";
my $out_root="$root-dPSI$min_dPSI$tail";
$output_file="DiffAS-$out_root.tab" unless (defined $output_file);
open (O, ">$folder/$output_file") or errPrintDie "Can't open the output file (do not provide a path)\n"; # output file

#### prepare to obtain gene IDs for GO analyses
my %ID_gene;
if (defined $get_GO){
    
    if (defined $ID_file){
	open (KEY, $ID_file) || errPrintDie "Can't open ID keys file ($ID_file)\n";
	while (<KEY>){
	    chomp($_);
	    my @temp=split(/\t/,$_);
	    $ID_gene{$temp[0]}=$temp[1]; # stores the geneID for each eventID
	}
	close KEY;
    }

    open (BG, ">$folder/BG-$out_root.txt") or errPrintDie "Can't open GO output files";
    open (IR_UP, ">$folder/IR_UP-$out_root.txt") or errPrintDie "Can't open GO output files";
    open (IR_DOWN, ">$folder/IR_DOWN-$out_root.txt") or errPrintDie "Can't open GO output files";
    open (EXSK, ">$folder/AltEx-$out_root.txt") or errPrintDie "Can't open GO output files";
}

#### Global variables for PSI analysis & GO
my %doneIR_UP;
my %doneIR_DOWN;
my %doneEXSK;
my %doneBG;
my %tally;
my %tally_total; # to count the total number of AS events with good coverage
my %tally_total_AS; # to count the total number of AS events within the compared samples

$tally{MIC}{DOWN}=0; $tally{MIC}{UP}=0; $tally_total_AS{MIC}=0; $tally_total{MIC}=0;
$tally{AltEx}{DOWN}=0; $tally{AltEx}{UP}=0; $tally_total_AS{AltEx}=0; $tally_total{AltEx}=0;
$tally{IR}{DOWN}=0; $tally{IR}{UP}=0; $tally_total_AS{IR}=0; $tally_total{IR}=0;
$tally{Alt3}{DOWN}=0; $tally{Alt3}{UP}=0; $tally_total_AS{Alt3}=0; $tally_total{Alt3}=0;
$tally{Alt5}{DOWN}=0; $tally{Alt5}{UP}=0; $tally_total_AS{Alt5}=0; $tally_total{Alt5}=0;

verbPrint "Doing comparisons of AS profiles ($name_A vs $name_B)\n";
print O "$head_row\n"; # it will print the original data (for plot later)

### starts the actual analysis
while (<PSI>){
    $_ =~ s/VLOW/N/g if (defined $noVLOW);
    chomp($_);
    my @t=split(/\t/,$_);
    my @PSI_A = ();
    my @PSI_B = ();

    next if $t[3] == 0; # removes the internal Alt3 and Alt5 splice sites to avoid double counting
    
    # defines AS type
    my $type="";
    $type="MIC" if ($t[1]=~/EX/ || ($t[5]!~/IR/ && $t[5]!~/Alt/)) && $t[3]<=27;
    $type="AltEx" if ($t[1]=~/EX/ || ($t[5]!~/IR/ && $t[5]!~/Alt/)) && $t[3]>27;
    $type="IR" if $t[1]=~/INT/ || $t[5]=~/IR/;
    $type="Alt3" if $t[1]=~/ALTA/ || $t[5]=~/Alt3/;
    $type="Alt5" if $t[1]=~/ALTD/ || $t[5]=~/Alt5/;
	
    # coverage check (requires good coverage for ALL replicates)
    my $kill_coverage = 0;
    foreach my $s (@samplesA){
	$kill_coverage = 1 if ($t[$s+1]!~/$Q/); # kill if ANY of the samples does not meet the coverage criteria
	push(@PSI_A,$t[$s]);
    }
    foreach my $s (@samplesB){
	$kill_coverage = 1 if ($t[$s+1]!~/$Q/); # kill if ANY of the samples does not meet the coverage criteria
	push(@PSI_B,$t[$s]);
    }
    next if ($kill_coverage == 1);

    # IR check (only checks the p if the p_IR is active)
    if (($type eq "IR") && (defined $p_IR)){ # only checks the p if the p_IR is active
	my $kill_pIR = 0;
	foreach my $s (@samplesA){
	    my ($temp_p)=$t[$s+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    $kill_pIR = 1 if $temp_p < 0.05;
	}
	foreach my $s (@samplesB){
	    my ($temp_p)=$t[$s+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    $kill_pIR = 1 if $temp_p < 0.05;
	}
	next if ($kill_pIR == 1);
    }

    # get PSIs
    my $sum_A = sum(@PSI_A);
    my $sum_B = sum(@PSI_B);
    my $av_PSI_A = sprintf("%.2f",$sum_A/$repA);
    my $av_PSI_B = sprintf("%.2f",$sum_B/$repB);
    # get min and max (for ranges)
    my $min_A = (sort{$b<=>$a}@PSI_A)[-1];
    my $max_A = (sort{$a<=>$b}@PSI_A)[-1];
    my $min_B = (sort{$b<=>$a}@PSI_B)[-1];
    my $max_B = (sort{$a<=>$b}@PSI_B)[-1];

    ### To count the total number of AS events considered 
    $tally_total_AS{$type}++ if ($av_PSI_A>10 && $av_PSI_A<90) || ($av_PSI_B>10 && $av_PSI_B<90) || abs($av_PSI_A-$av_PSI_B)>10;
    $tally_total{$type}++;
    
    # NOT PAIRED: gets the average PSI for A and B and the lowest (min) and highest (max) PSI for each replicate
    if (!defined $paired){
	# get dPSI
	my $dPSI = $av_PSI_B-$av_PSI_A;
	
	# does the diff AS test:
	if ($dPSI > $min_dPSI && $min_B > $max_A+$min_range){ # if rep1 it will always meet the criteria
	    $tally{$type}{UP}++;
	    print O "$_\n"; # dPSI is not printed so it can the be run with plot
	    
	    # print for GO
	    if (defined$get_GO){
		unless ($use_names){
		    if (defined $ID_gene{$t[1]}){
			print IR_UP "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_UP{$ID_gene{$t[1]}});
			$doneIR_UP{$ID_gene{$t[1]}}=1;
			print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			$doneEXSK{$ID_gene{$t[1]}}=1;
		    }
		}
		else {
		    print IR_UP "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_UP{$t[0]}) && (defined $t[0]);
		    $doneIR_UP{$t[0]}=1 if (defined $t[0]);
		    print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
		    $doneEXSK{$t[0]}=1 if (defined $t[0]);
		}
	    }
	}
	if ($dPSI < -1*$min_dPSI && $min_A > $max_B+$min_range){
	    $tally{$type}{DOWN}++;
	    print O "$_\n";
	    
	    #print for GO
	    if (defined$get_GO){
		unless ($use_names){
		    if (defined $ID_gene{$t[1]}){
			print IR_DOWN "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$ID_gene{$t[1]}});
			$doneIR_DOWN{$ID_gene{$t[1]}}=1;
			print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			$doneEXSK{$ID_gene{$t[1]}}=1;	
		    }
		}
		else {
		    print IR_DOWN "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$t[0]}) && (defined $t[0]);
		    $doneIR_DOWN{$t[0]}=1 if (defined $t[0]);
		    print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
		    $doneEXSK{$t[0]}=1 if (defined $t[0]);
		}
	    }
	}
    }
    else { # if paired: calculates each pair's dPSI & the lower/higher of these
	# get each dPSI and the average
	my @dPSI_pairs = ();
	for my $i (0..$#samplesA){ # ~ $repA-1
	    $dPSI_pairs[$i] = sprintf("%.2f",$PSI_B[$i]-$PSI_A[$i]);
	}
	my $sum_paired_dPSI = sum(@dPSI_pairs);
	my $av_paired_dPSI = sprintf("%.2f",$sum_paired_dPSI/$repA);

	# get min and max (for min dPSI, pos and neg values)
	my $min_indiv_dPSI = (sort{$b<=>$a}@dPSI_pairs)[-1];
	my $max_indiv_dPSI = (sort{$a<=>$b}@dPSI_pairs)[-1];
	
	if ($av_paired_dPSI > $min_dPSI && $min_indiv_dPSI > $min_range){ 
	    $tally{$type}{UP}++;
	    print O "$_\n";
	    
	    # print for GO
	    if (defined $get_GO){
		unless ($use_names){
		    if (defined $ID_gene{$t[1]}){
			print IR_UP "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_UP{$ID_gene{$t[1]}});
			$doneIR_UP{$ID_gene{$t[1]}}=1;
			print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			$doneEXSK{$ID_gene{$t[1]}}=1;
		    }
		}
		else {
		    print IR_UP "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_UP{$t[0]}) && (defined $t[0]);
		    $doneIR_UP{$t[0]}=1 if (defined $t[0]);
		    print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
		    $doneEXSK{$t[0]}=1 if (defined $t[0]);
		}
	    }
	}
	if ($av_paired_dPSI < -$min_dPSI && $max_indiv_dPSI < -$min_range){ 
	    $tally{$type}{DOWN}++;
	    print O "$_\n";
	    
	    #print for GO
	    if (defined $get_GO){
		unless ($use_names){
		    if (defined $ID_gene{$t[1]}){
			print IR_DOWN "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$ID_gene{$t[1]}});
			$doneIR_DOWN{$ID_gene{$t[1]}}=1;
			print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			$doneEXSK{$ID_gene{$t[1]}}=1;	
		    }	
		}
		else {
		    print IR_DOWN "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$t[0]}) && (defined $t[0]);
		    $doneIR_DOWN{$t[0]}=1 if (defined $t[0]);
		    print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
		    $doneEXSK{$t[0]}=1 if (defined $t[0]);
		}
	    }
	}
    }
    # prints out the genes for BG
    if (defined $get_GO){
	unless ($use_names){
	    if (defined $ID_gene{$t[1]}){
		print BG "$ID_gene{$t[1]}\n" if (!defined $doneBG{$ID_gene{$t[1]}});
		$doneBG{$ID_gene{$t[1]}}=1;
	    }
	}
	else {
	    print BG "$t[0]\n" if (!defined $doneBG{$t[0]}) && (defined $t[0]);
	    $doneBG{$t[0]}=1 if (defined $t[0]);
	}
    }
}
close PSI;
close O;

if (defined $get_GO){
    verbPrint "Preparing files for GO analysis\n";
    sleep(1);
    close BG;
    close IR_DOWN;
    close IR_UP;
    close EXSK;
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
$extras.=", noVLOW" if (defined $noVLOW);
$extras.=", p_IR" if (defined $p_IR);
$extras.=", paired" if (defined $paired);

print "\n*** Options: dPSI=$min_dPSI, range_dif=$min_range$extras\n";
print "*** Summary statistics:\n";
print "\tAS_TYPE\tHigher_in_$name_A\tHigher_in_$name_B\tTOTAL_EV\tTOTAL_AS(10<PSI<90)\n";
print "\tMicroexons\t$tally{MIC}{DOWN}\t$tally{MIC}{UP}\t$tally_total{MIC}\t$tally_total_AS{MIC}\n";
print "\tLong_AltEx\t$tally{AltEx}{DOWN}\t$tally{AltEx}{UP}\t$tally_total{AltEx}\t$tally_total_AS{AltEx}\n";
print "\tIntron_ret\t$tally{IR}{DOWN}\t$tally{IR}{UP}\t$tally_total{IR}\t$tally_total_AS{IR}\n";
print "\tAlt_3ss\t$tally{Alt3}{DOWN}\t$tally{Alt3}{UP}\t$tally_total{Alt3}\t$tally_total_AS{Alt3}\n";
print "\tAlt_5ss\t$tally{Alt5}{DOWN}\t$tally{Alt5}{UP}\t$tally_total{Alt5}\t$tally_total_AS{Alt5}\n";
print "\n";
