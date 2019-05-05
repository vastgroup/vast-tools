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
my $plot; # changed from no_plot to plot (01/04/18)
my $plot_only_samples;
my $print_dPSI;
my $print_sets;
my $max_dPSI;
my $print_all_ev;
my $print_AS_ev;
my $use_int_reads;
my $fr_int_reads = 0.4;
my $min_ALT_use = 25;

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
			  "use_int_reads" => \$use_int_reads,
			  "fr_int_reads=f" => \$fr_int_reads,
			  "GO" => \$get_GO,
			  "sp=s" => \$species,
			  "species=s" => \$species,
			  "dbDir=s" => \$dbDir,
			  "GO_file=s" => \$ID_file,
			  "use_names" => \$use_names,
			  "paired" => \$paired,
			  "print_dPSI" => \$print_dPSI,
			  "print_sets" => \$print_sets,
			  "print_all_ev" => \$print_all_ev,
			  "print_AS_ev" => \$print_AS_ev,
			  "max_dPSI=i"   => \$max_dPSI,
			  "plot_PSI" => \$plot,
			  "only_samples" => \$plot_only_samples,
			  "noVLOW" => \$noVLOW,
			  "min_ALT_use=i" => \$min_ALT_use
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

sub time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
    $year += 1900;
    $mon += 1;
    my $datetime = sprintf "%04d-%02d-%02d (%02d:%02d)", $year, $mday, $mon, $hour, $min;
    return $datetime;
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


### Gets the version
my $version;
open (VERSION, "$binPath/../VERSION");
$version=<VERSION>;
chomp($version);
$version="No version found" if !$version;

if (!defined($ARGV[0]) || $helpFlag){
    die "
VAST-TOOLS v$version

Usage: vast-tools compare /path/to/INCLUSION_LEVELS_FULL-root.tab -a sample_a1,sample_a2 -b sample_b1,sample_b2 [options]

Compare two sample sets to find differentially regulated AS events. 
INCLUSION_LEVELS_FULL-root.tab is final table produced by VAST-TOOLs command combine.

[General options] 
        --min_dPSI i             Minimum delta PSI of the averages (default >15)
        --min_range i            Minimum distance between the ranges of both groups (default >5)
        --outFile file           Output file name (default based on option parameters)
        -a/--samplesA sA1,sA2    Required, 1:n sample names or column_\# separated by , (mandatory)
        -b/--samplesB sB1,sB2    Required, 1:n sample names or column_\# separated by , (mandatory)
        --noVLOW                 Does not use samples with VLOW coverage (default OFF)
        --p_IR                   Filter IR by the p-value of the binomial test (default OFF)
        --use_int_reads          Requires a minimum fraction of intron body reads (--fr_int_reads) with respect to 
                                   those in the EIJ for IR (default OFF)(combine >= v2.1.3)
        --fr_int_reads fr        Minimum fraction of reads in the intron bodies respect to the 
                                   average of the EI/IE junctions (default 0.4). Only useful with --use_int_reads
        --min_ALT_use i          Minimum inclusion of the exon in which the Alt3/Alt5 is located across all 
                                   compared samples (default 25) (combine >= v2.2.1)
        --print_dPSI             Prints the mean dPSI (PSI_B-PSI_A) as last column (default OFF)
                                   - It does not allow ploting.
        --print_sets             Prints files with different sets for comparisons in Matt (http://matt.crg.eu):
                                   - CS: all events with coverage and constitutively spliced (PSI>95 for AltEx, PSI<5 for IR)
                                   - CR: all events with coverage and cryptically spliced (PSI<5 for AltEx, PSI>95 for IR)
                                   - AS_NC: all events with coverage, alternative (10 < av_PSI < 90 in a group)
                                            and that do not change between the two conditions (abs(dPSI)< max_dPSI)
        --print_all_ev           Prints a table with all event IDs that pass the coverage filters and dPSI (default OFF)
        --print_AS_ev            Prints a table with all AS events that pass the coverage filters and dPSI (default OFF)  
        --max_dPSI i             Maximum dPSI to consider an AS non-changing (default min_dPSI/5)
        --plot_PSI               Plots the DS events using \'plot\' (default OFF)
        --only_samples           Plots only the compared samples, otherwise the whole table (default OFF)
        --paired                 Does a paired comparison (A1 vs B1, A2 vs B2, etc.)
                                   - min_dPSI: minimum value for the average of each paired dPSI
                                   - min_range: minumum value the smallest paired dPSI can take
    
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

#### Check if plot and print_dPSI are active together
errPrintDie "print_dPSI cannot be used with plot\n" if (defined $print_dPSI) && (defined $plot);

#### opens INCLUSION TABLE
open (PSI, $input_file) or errPrintDie "Needs a PSI INCLUSION table\n";

# prints version (05/05/19)
verbPrint "VAST-TOOLS v$version";

### Creates the LOG
open (LOG, ">>$folder/VTS_LOG_commands.txt");
my $all_args="-o $folder -min_dPSI $min_dPSI -min_range $min_range -a $samplesA -b $samplesB -min_ALT_use $min_ALT_use";
$all_args.=" -paired" if $paired;
$all_args.=" -noVLOW" if $noVLOW;
$all_args.=" -p_IR" if $p_IR;
$all_args.=" -use_int_reads" if $use_int_reads;
$all_args.=" -fr_int_reads $fr_int_reads" if defined $fr_int_reads;
$all_args.=" -print_dPSI" if $print_dPSI;
$all_args.=" -print_sets" if $print_sets;
$all_args.=" -print_all_ev" if $print_all_ev;
$all_args.=" -print_AS_ev" if $print_AS_ev;
$all_args.=" -max_dPSI=i"   if defined $max_dPSI;
$all_args.=" -GO" if $get_GO;
$all_args.=" -sp $species" if $species;
$all_args.=" -outFile $output_file" if $output_file;
$all_args.=" -plot_PSI" if $plot;
$all_args.=" -only_samples" if $plot_only_samples;

print LOG "[VAST-TOOLS v$version, ".&time."] vast-tools compare $all_args\n";

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
my ($root)=$ARGV[0]=~/.+?\-([^\/]+?)\./;
my $tail = ""; # to be added to the output name
$tail.="-range$min_range" if (defined $min_range); 
$tail.="-noVLOW" if (defined $noVLOW);
$tail.="-p_IR" if (defined $p_IR);
$tail.="-IR_reads" if (defined $use_int_reads);
$tail.="-min_ALT_use$min_ALT_use";
$tail.="-paired" if (defined $paired);
$tail.="_$name_A-vs-$name_B";
$tail.="-with_dPSI" if (defined $print_dPSI);
my $out_root="$root-dPSI$min_dPSI$tail" unless (defined $output_file);
($out_root)=$output_file=~/([^\/]+)\./ if (defined $output_file && $output_file=~/[^\/]+\./);
($out_root)=$output_file=~/([^\/]+)/ if (defined $output_file && $output_file!~/[^\/]+\./);
$out_root=~s/DiffAS\-// if (defined $output_file);

$output_file="DiffAS-$out_root.tab" unless (defined $output_file);
open (O, ">$folder/$output_file") or errPrintDie "Can't open the output file (do not provide a path)\n"; # output file

my $all_ev_file;
if (defined $print_all_ev){
    $all_ev_file="AllEvents-$output_file";
    $all_ev_file=~s/DiffAS\-//g;
    open (O_ALL, ">$folder/$all_ev_file") or errPrintDie "Can't open the output file for all events (do not provide a path)\n"; # output file  
}
my $all_AS_file;
if (defined $print_AS_ev){
    $all_AS_file="ASEvents-$output_file";
    $all_AS_file=~s/DiffAS\-//g;
    open (O_AS, ">$folder/$all_AS_file") or errPrintDie "Can't open the output file for all AS events (do not provide a path)\n"; # output file  
}

####### Other sets file
if (defined $print_sets){
    $max_dPSI = $min_dPSI/5 if (!defined $max_dPSI);
    my $CS_file="CS-$out_root.tab";
    my $CR_file="CR-$out_root.tab";
    my $AS_file="AS_NC-$out_root-Max_dPSI$max_dPSI.tab";
    
    open (SET_CS, ">$folder/$CS_file") or errPrintDie "Can't open the CS file\n"; # output file 
    open (SET_CR, ">$folder/$CR_file") or errPrintDie "Can't open the CR file\n"; # output file 
    open (SET_AS, ">$folder/$AS_file") or errPrintDie "Can't open the AS_NC file\n"; # output file 
}

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
    open (ALL_EV, ">$folder/All_Ev-$out_root.txt") or errPrintDie "Can't open GO output files";
}

#### Global variables for PSI analysis & GO
my %doneIR_UP;
my %doneIR_DOWN;
my %doneEXSK;
my %doneALL;
my %doneBG;
my %tally;
my %tally_total; # to count the total number of AS events with good coverage
my %tally_total_AS; # to count the total number of AS events within the compared samples
my %tally_extra; # to count the number of other types

### Setting all to 0
$tally{MIC}{DOWN}=0; $tally{MIC}{UP}=0; $tally_total_AS{MIC}=0; $tally_total{MIC}=0;
$tally{AltEx}{DOWN}=0; $tally{AltEx}{UP}=0; $tally_total_AS{AltEx}=0; $tally_total{AltEx}=0;
$tally{IR}{DOWN}=0; $tally{IR}{UP}=0; $tally_total_AS{IR}=0; $tally_total{IR}=0;
$tally{Alt3}{DOWN}=0; $tally{Alt3}{UP}=0; $tally_total_AS{Alt3}=0; $tally_total{Alt3}=0;
$tally{Alt5}{DOWN}=0; $tally{Alt5}{UP}=0; $tally_total_AS{Alt5}=0; $tally_total{Alt5}=0;

$tally_extra{MIC}{CS}=0; $tally_extra{MIC}{CR}=0; $tally_extra{MIC}{AS_NC}=0;
$tally_extra{AltEx}{CS}=0; $tally_extra{AltEx}{CR}=0; $tally_extra{AltEx}{AS_NC}=0;
$tally_extra{IR}{CS}=0; $tally_extra{IR}{CR}=0; $tally_extra{IR}{AS_NC}=0;
$tally_extra{Alt3}{CS}=0; $tally_extra{Alt3}{CR}=0; $tally_extra{Alt3}{AS_NC}=0;
$tally_extra{Alt5}{CS}=0; $tally_extra{Alt5}{CR}=0; $tally_extra{Alt5}{AS_NC}=0;

verbPrint "Doing comparisons of AS profiles ($name_A vs $name_B)\n";
unless (defined $print_dPSI){
    print O "$head_row\n"; # it will print the original data (for plot later)
    if (defined $print_sets){
	print SET_AS "$head_row\n";
	print SET_CR "$head_row\n";
	print SET_CS "$head_row\n";
    }
}
else {
    print O "$head_row\tdPSI\n"; # it will print the original data + dPSI
    if (defined $print_sets){
	print SET_AS "$head_row\tdPSI\n";
	print SET_CR "$head_row\tdPSI\n";
	print SET_CS "$head_row\tdPSI\n";
    }
}
print O_ALL "EventID\tPSI_A\tPSI_B\tdPSI\tCATEGORY\n" if (defined $print_all_ev);
print O_AS "EventID\tPSI_A\tPSI_B\tdPSI\tCATEGORY\n" if (defined $print_AS_ev);

### starts the actual analysis
while (<PSI>){
    $_ =~ s/VLOW/N/g if (defined $noVLOW);
    chomp($_);
    my @t=split(/\t/,$_);
    my @PSI_A = ();
    my @PSI_B = ();
    my %int_body_reads = ();
    my %int_junct_reads = ();
    my %av_int_reads = ();
    my %av_junct_reads = ();
    my $event = $t[1];

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
	$kill_coverage = 1 if ($t[$s] eq "NA"); # for some MULTI cases without mappability (more strict in V2)
	push(@PSI_A,$t[$s]);
    }
    foreach my $s (@samplesB){
	$kill_coverage = 1 if ($t[$s+1]!~/$Q/); # kill if ANY of the samples does not meet the coverage criteria
	$kill_coverage = 1 if ($t[$s] eq "NA"); # for some MULTI cases without mappability (more strict in V2)
	push(@PSI_B,$t[$s]);
    }
    next if ($kill_coverage == 1);

    # min PSI-like usage for ALT3/5 (min_ALT_use) 04/05/19
    if ($type eq "Alt3" || $type eq "Alt5"){
	my $kill_ALT = 0;
        foreach my $s (@samplesA){
            my ($temp_ALT)=$t[$s+1]=~/O[KW]\,.+?\,(.+?)\,.+?\,.+?\@/;
	    if ($temp_ALT=~/\d/){
		$kill_ALT = 1 if $temp_ALT < $min_ALT_use;
	    }
	    else {
		$min_ALT_use = "NA (older version)";
	    }
        }
        foreach my $s (@samplesB){
            my ($temp_ALT)=$t[$s+1]=~/O[KW]\,.+?\,(.+?)\,.+?\,.+?\@/;
	    if ($temp_ALT=~/\d/){
		$kill_ALT = 1 if $temp_ALT < $min_ALT_use;
	    }
	    else {
		$min_ALT_use = "NA (older version)";
	    }
        }
        next if ($kill_ALT == 1);
    }

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
    # Gets the number of intron body reads (from v2.1.3)
    if (($type eq "IR") && (defined $use_int_reads)){
	foreach my $s (@samplesA){
            my ($temp_r_ib,$le_body,$temp_r_EI,$temp_r_IE)=$t[$s+1]=~/O[KW]\,.+?\,(.+?)\=(.+?)\,(.+?)\=(.+?)\=.+?\,.+?\@/;
	    if (defined $temp_r_ib && $temp_r_ib =~ /\d/ && $le_body >= 5){
		$int_body_reads{A}+= $temp_r_ib;
		$int_junct_reads{A}+= sprintf("%.2f",($temp_r_EI+$temp_r_IE)/2);
	    }
	    else {
		$int_body_reads{A} = "NA";
		$int_junct_reads{A} = "NA";
	    }
        }
        foreach my $s (@samplesB){
            my ($temp_r_ib,$le_body,$temp_r_EI,$temp_r_IE)=$t[$s+1]=~/O[KW]\,.+?\,(.+?)\=(.+?)\,(.+?)\=(.+?)\=.+?\,.+?\@/;
	    if (defined $temp_r_ib && $temp_r_ib =~ /\d/ && $le_body >= 5){
		$int_body_reads{B}+= $temp_r_ib;
		$int_junct_reads{B}+= sprintf("%.2f",($temp_r_EI+$temp_r_IE)/2);
	    }
	    else {
		$int_body_reads{B} = "NA";
		$int_junct_reads{B} = "NA";
	    }
        }
	### Does the averages
	if ($int_body_reads{A} ne "NA" && $int_body_reads{B} ne "NA"){
	    $av_int_reads{A}=sprintf("%.2f",$int_body_reads{A}/($#samplesA+1));
	    $av_int_reads{B}=sprintf("%.2f",$int_body_reads{B}/($#samplesB+1));
	    $av_junct_reads{A}=sprintf("%.2f",$int_junct_reads{A}/($#samplesA+1));
	    $av_junct_reads{B}=sprintf("%.2f",$int_junct_reads{B}/($#samplesB+1));
	}
	else {
	    $av_int_reads{A} = "NA";
	    $av_int_reads{B} = "NA";
	    $av_junct_reads{A} = "NA";
	    $av_junct_reads{B} = "NA";
	}
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

    # NOT PAIRED: gets the average PSI for A and B and the lowest (min) and highest (max) PSI for each replicate
    if (!defined $paired){
	# get dPSI
	my $dPSI = $av_PSI_B-$av_PSI_A;

	### To count the total number of AS events considered and print sets if required
	if (($av_PSI_A>10 && $av_PSI_A<90) || ($av_PSI_B>10 && $av_PSI_B<90) || abs($av_PSI_A-$av_PSI_B)>10){
	    $tally_total_AS{$type}++;
	    if (defined $print_AS_ev){
		print O_AS "$event\t$av_PSI_A\t$av_PSI_B\t$dPSI\tAS_EV\n";
	    }
	}
	$tally_total{$type}++;
	if (defined $print_all_ev){
	    print O_ALL "$event\t$av_PSI_A\t$av_PSI_B\t$dPSI\tBG\n";
	}
	
	
	# does the diff AS test:
	if ($dPSI > $min_dPSI && $min_B > $max_A+$min_range){ # if rep1 it will always meet the criteria
	    if (($type ne "IR") || (!defined $use_int_reads) || ($type eq "IR" && $av_int_reads{B} eq "NA")  || ($type eq "IR" && $av_int_reads{B}/$av_junct_reads{B} >= $fr_int_reads)){
		$tally{$type}{UP}++;
		unless (defined $print_dPSI){
		    print O "$_\n"; # dPSI is not printed so it can the be run with plot
		}
		else {
		    print O "$_\t$dPSI\n";
		}
		
		# print for GO
		if (defined $get_GO){
		    unless ($use_names){
			if (defined $ID_gene{$t[1]}){
			    print IR_UP "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_UP{$ID_gene{$t[1]}});
			    $doneIR_UP{$ID_gene{$t[1]}}=1;
			    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			    $doneEXSK{$ID_gene{$t[1]}}=1;
			    print ALL_EV "$ID_gene{$t[1]}\n" if !defined $doneALL{$ID_gene{$t[1]}};
			    $doneALL{$ID_gene{$t[1]}}=1;
			}
		    }
		    else {
			print IR_UP "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_UP{$t[0]}) && (defined $t[0]);
			$doneIR_UP{$t[0]}=1 if (defined $t[0]);
			print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
			$doneEXSK{$t[0]}=1 if (defined $t[0]);
			print ALL_EV "$t[0]\n" if (!defined $doneALL{$t[0]}) && (defined $t[0]);
			$doneALL{$t[0]}=1 if (defined $t[0]);
		    }
		}
	    }
	}
	if ($dPSI < -1*$min_dPSI && $min_A > $max_B+$min_range){
	    if (($type ne "IR") || (!defined $use_int_reads) || ($type eq "IR" && $av_int_reads{A} eq "NA")  || ($type eq "IR" && $av_int_reads{A}/$av_junct_reads{A} >= $fr_int_reads)){
		$tally{$type}{DOWN}++;
		unless (defined $print_dPSI){
		    print O "$_\n"; # dPSI is not printed so it can the be run with plot
		}
		else {
		    print O "$_\t$dPSI\n";
		}
		
		#print for GO
		if (defined$get_GO){
		    unless ($use_names){
			if (defined $ID_gene{$t[1]}){
			    print IR_DOWN "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$ID_gene{$t[1]}});
			    $doneIR_DOWN{$ID_gene{$t[1]}}=1;
			    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			    $doneEXSK{$ID_gene{$t[1]}}=1;	
			    print ALL_EV "$ID_gene{$t[1]}\n" if (!defined $doneALL{$ID_gene{$t[1]}});
			    $doneALL{$ID_gene{$t[1]}}=1;	
			}
		    }
		    else {
			print IR_DOWN "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$t[0]}) && (defined $t[0]);
			$doneIR_DOWN{$t[0]}=1 if (defined $t[0]);
			print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
			$doneEXSK{$t[0]}=1 if (defined $t[0]);
			print ALL_EV "$t[0]\n" if (!defined $doneALL{$t[0]}) && (defined $t[0]);
			$doneALL{$t[0]}=1 if (defined $t[0]);
		    }
		}
	    }
	}
	### Prints the extra sets
	if (defined $print_sets){
	    ### Set of AS with no change
	    if (abs($dPSI) < $max_dPSI && (($av_PSI_A>10 && $av_PSI_A<90) || ($av_PSI_B>10 && $av_PSI_B<90) || abs($av_PSI_A-$av_PSI_B)>10)){
		unless (defined $print_dPSI){
		    print SET_AS "$_\n"; # dPSI is not printed so it can the be run with plot
		}
		else {
		    print SET_AS "$_\t$dPSI\n";
		}
		$tally_extra{$type}{AS_NC}++;
	    }
	    ### Prints cryptic and constitutive
	    if ($type eq "IR"){
		unless (defined $print_dPSI){
		    print SET_CS "$_\n" if $max_A < 5 && $max_B < 5;
		    print SET_CR "$_\n" if $min_A > 95 && $min_B > 95;
		}
		else {
		    print SET_CS "$_\t$dPSI\n" if $max_A < 5 && $max_B < 5;
		    print SET_CR "$_\t$dPSI\n" if $min_A > 95 && $min_B > 95;
		}
		$tally_extra{$type}{CS}++ if $max_A < 5 && $max_B < 5;
		$tally_extra{$type}{CR}++ if $min_A > 95 && $min_B > 95;
	    }
	    else {
		unless (defined $print_dPSI){
		    print SET_CR "$_\n" if $max_A < 5 && $max_B < 5;
		    print SET_CS "$_\n" if $min_A > 95 && $min_B > 95;
		}
		else {
		    print SET_CR "$_\t$dPSI\n" if $max_A < 5 && $max_B < 5;
		    print SET_CS "$_\t$dPSI\n" if $min_A > 95 && $min_B > 95;
		}
		$tally_extra{$type}{CR}++ if $max_A < 5 && $max_B < 5;
		$tally_extra{$type}{CS}++ if $min_A > 95 && $min_B > 95;
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

	### To count the total number of AS events considered and print sets if required
	if (($av_PSI_A>10 && $av_PSI_A<90) || ($av_PSI_B>10 && $av_PSI_B<90) || abs($av_PSI_A-$av_PSI_B)>10){
	    $tally_total_AS{$type}++;
	    if (defined $print_AS_ev){
		print O_AS "$event\t$av_PSI_A\t$av_PSI_B\t$av_paired_dPSI\tAS_EV\n";
	    }
	}
	$tally_total{$type}++;
	if (defined $print_all_ev){
	    print O_ALL "$event\t$av_PSI_A\t$av_PSI_B\t$av_paired_dPSI\tBG\n";
	}
	
	### Does the diff tests
	if ($av_paired_dPSI > $min_dPSI && $min_indiv_dPSI > $min_range){ 
	    if (($type ne "IR") || (!defined $use_int_reads) || ($type eq "IR" && $av_int_reads{B} eq "NA") || ($type eq "IR" && $av_int_reads{B}/$av_junct_reads{B} >= $fr_int_reads)){
		$tally{$type}{UP}++;
		unless (defined $print_dPSI){
		    print O "$_\n";
		}
		else {
		    print O "$_\t$av_paired_dPSI\n";
		}
		
		# print for GO
		if (defined $get_GO){
		    unless ($use_names){
			if (defined $ID_gene{$t[1]}){
			    print IR_UP "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_UP{$ID_gene{$t[1]}});
			    $doneIR_UP{$ID_gene{$t[1]}}=1;
			    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			    $doneEXSK{$ID_gene{$t[1]}}=1;
			    print ALL_EV "$ID_gene{$t[1]}\n" if (!defined $doneALL{$ID_gene{$t[1]}});
			    $doneALL{$ID_gene{$t[1]}}=1;
			}
		    }
		    else {
			print IR_UP "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_UP{$t[0]}) && (defined $t[0]);
			$doneIR_UP{$t[0]}=1 if (defined $t[0]);
			print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
			$doneEXSK{$t[0]}=1 if (defined $t[0]);
			print ALL_EV "$t[0]\n" if (!defined $doneALL{$t[0]}) && (defined $t[0]);
			$doneALL{$t[0]}=1 if (defined $t[0]);
		    }
		}
	    }
	}
	if ($av_paired_dPSI < -$min_dPSI && $max_indiv_dPSI < -$min_range){ 
	    if (($type ne "IR") || (!defined $use_int_reads) || ($type eq "IR" && $av_int_reads{A} eq "NA") || ($type eq "IR" && $av_int_reads{A}/$av_junct_reads{A} >= $fr_int_reads)){
		$tally{$type}{DOWN}++;
		unless (defined $print_dPSI){
		    print O "$_\n";
		}
		else {
		    print O "$_\t$av_paired_dPSI\n";
		}
		
		#print for GO
		if (defined $get_GO){
		    unless ($use_names){
			if (defined $ID_gene{$t[1]}){
			    print IR_DOWN "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$ID_gene{$t[1]}});
			    $doneIR_DOWN{$ID_gene{$t[1]}}=1;
			    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
			    $doneEXSK{$ID_gene{$t[1]}}=1;	
			    print ALL_EV "$ID_gene{$t[1]}\n" if (!defined $doneALL{$ID_gene{$t[1]}});
			    $doneALL{$ID_gene{$t[1]}}=1;	
			}	
		    }
		    else {
			print IR_DOWN "$t[0]\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$t[0]}) && (defined $t[0]);
			$doneIR_DOWN{$t[0]}=1 if (defined $t[0]);
			print EXSK "$t[0]\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$t[0]}) && (defined $t[0]);
			$doneEXSK{$t[0]}=1 if (defined $t[0]);
			print ALL_EV "$t[0]\n" if (!defined $doneALL{$t[0]}) && (defined $t[0]);
			$doneALL{$t[0]}=1 if (defined $t[0]);
		    }
		}
	    }
	}
	### Prints the extra sets
	if (defined $print_sets){
	    ### Set of AS with no change
	    if (abs($av_paired_dPSI) < $max_dPSI && (($av_PSI_A>10 && $av_PSI_A<90) || ($av_PSI_B>10 && $av_PSI_B<90) || abs($av_PSI_A-$av_PSI_B)>10)){
		unless (defined $print_dPSI){
		    print SET_AS "$_\n"; # dPSI is not printed so it can the be run with plot
		}
		else {
		    print SET_AS "$_\t$av_paired_dPSI\n";
		}
		$tally_extra{$type}{AS_NC}++;
	    }
	    ### Prints cryptic and constitutive
	    if ($type eq "IR"){
		unless (defined $print_dPSI){
		    print SET_CS "$_\n" if $max_A < 5 && $max_B < 5;
		    print SET_CR "$_\n" if $min_A > 95 && $min_B > 95;
		}
		else {
		    print SET_CS "$_\t$av_paired_dPSI\n" if $max_A < 5 && $max_B < 5;
		    print SET_CR "$_\t$av_paired_dPSI\n" if $min_A > 95 && $min_B > 95;
		}
		$tally_extra{$type}{CS}++ if $max_A < 5 && $max_B < 5;
		$tally_extra{$type}{CR}++ if $min_A > 95 && $min_B > 95;
	    }
	    else {
		unless (defined $print_dPSI){
		    print SET_CR "$_\n" if $max_A < 5 && $max_B < 5;
		    print SET_CS "$_\n" if $min_A > 95 && $min_B > 95;
		}
		else {
		    print SET_CR "$_\t$av_paired_dPSI\n" if $max_A < 5 && $max_B < 5;
		    print SET_CS "$_\t$av_paired_dPSI\n" if $min_A > 95 && $min_B > 95;
		}
		$tally_extra{$type}{CR}++ if $max_A < 5 && $max_B < 5;
		$tally_extra{$type}{CS}++ if $min_A > 95 && $min_B > 95;
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
    close ALL_EV;
}

if (defined $plot){
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
$extras.=", use_int_reads" if (defined $use_int_reads);
$extras.=", paired" if (defined $paired);

print "\n*** Options: dPSI=$min_dPSI, range_dif=$min_range$extras, min_ALT_use=$min_ALT_use\n";
print "*** Summary statistics:\n";
print "\tAS_TYPE\tHigher_in_$name_A\tHigher_in_$name_B\tTOTAL_EV\tTOTAL_AS(10<PSI<90)\n";
print "\tMicroexons\t$tally{MIC}{DOWN}\t$tally{MIC}{UP}\t$tally_total{MIC}\t$tally_total_AS{MIC}\n";
print "\tLong_AltEx\t$tally{AltEx}{DOWN}\t$tally{AltEx}{UP}\t$tally_total{AltEx}\t$tally_total_AS{AltEx}\n";
print "\tIntron_ret\t$tally{IR}{DOWN}\t$tally{IR}{UP}\t$tally_total{IR}\t$tally_total_AS{IR}\n";
print "\tAlt_3ss\t$tally{Alt3}{DOWN}\t$tally{Alt3}{UP}\t$tally_total{Alt3}\t$tally_total_AS{Alt3}\n";
print "\tAlt_5ss\t$tally{Alt5}{DOWN}\t$tally{Alt5}{UP}\t$tally_total{Alt5}\t$tally_total_AS{Alt5}\n";
print "\n";

if (defined $print_sets){
    print 
"*** Summary statistics for extra event sets (Max_dPSI=$max_dPSI):
\tAS_TYPE\tConstitutive\tCryptic\tAS_non_change
\tMicroexons\t$tally_extra{MIC}{CS}\t$tally_extra{MIC}{CR}\t$tally_extra{MIC}{AS_NC}
\tLong_AltEx\t$tally_extra{AltEx}{CS}\t$tally_extra{AltEx}{CR}\t$tally_extra{AltEx}{AS_NC}
\tIntron_ret\t$tally_extra{IR}{CS}\t$tally_extra{IR}{CR}\t$tally_extra{IR}{AS_NC}
\tAlt_3ss\t$tally_extra{Alt3}{CS}\t$tally_extra{Alt3}{CR}\t$tally_extra{Alt3}{AS_NC}
\tAlt_5ss\t$tally_extra{Alt5}{CS}\t$tally_extra{Alt5}{CR}\t$tally_extra{Alt5}{AS_NC}

"
}
