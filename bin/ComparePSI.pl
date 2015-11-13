#!/usr/bin/perl
### General script to get differentially spliced events based on dPSI differences 

use Getopt::Long;
use warnings;
use strict;
use Cwd qw(abs_path);

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
my $a1; my $a2=0; my $a3=0; my $a4=0;
my $b1; my $b2=0; my $b3=0; my $b4=0;
my $rep1 = 2; # number of replicates per type (def=2);
my $rep2 = 2; # number of replicates per type (def=2);
my $min_range = 5; # min dPSI between ranges
my $noVLOW;
my $p_IR;
my $ID_file;
my $get_GO;
my $paired;
my $use_names;
my $folder;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "min_dPSI=i" => \$min_dPSI,
                          "repA=i" => \$rep1,
                          "repB=i" => \$rep2,
                          "min_range=i" => \$min_range,
                          "a1=s" => \$a1,
                          "a2=s" => \$a2,
                          "a3=s" => \$a3,
                          "a4=s" => \$a4,
                          "b1=s" => \$b1,
                          "b2=s" => \$b2,
                          "b3=s" => \$b3,
                          "b4=s" => \$b4,
			  "outFile=s" => \$output_file,
			  "output=s" => \$folder,
			  "o=s" => \$folder,
			  "help" => \$helpFlag,
			  "p_IR" => \$p_IR,
			  "GO" => \$get_GO,
			  "sp=s" => \$species,
			  "dbDir=s" => \$dbDir,
			  "GO_file=s" => \$ID_file,
			  "use_names" => \$use_names,
			  "paired" => \$paired,
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
    die "\nUsage: vast-tools compare INCLUSION_LEVELS_FULL-root.tab -a1 sample_a1 (-a2 sample_a2) -b1 sample_b1 (-b2 sample_b2) [options]

Compare two sample sets with 1,2, 3 or 4 replicates.

[General options] 
        --min_dPSI i             Minimum delta PSI of the averages (default 15)
        --min_range i            Minimum distance between the ranges of both groups (default 5)
        --outFile file           Output file name (default based on option parameters)
        --repA i                 Number of replicates for group A (default 2)
        --repB i                 Number of replicates for group B (default 2)
        -a1                      Sub-sample name for rep 1 of sample A OR column number (0-based in INCLUSION table)
        -a2                      Sub-sample name for rep 2 of sample A OR column number (0-based in INCLUSION table)
        -a3                      Sub-sample name for rep 3 of sample A OR column number (0-based in INCLUSION table)
        -a4                      Sub-sample name for rep 4 of sample A OR column number (0-based in INCLUSION table)
        -b1                      Sub-sample name for rep 1 of sample B OR column number (0-based in INCLUSION table)
        -b2                      Sub-sample name for rep 2 of sample B OR column number (0-based in INCLUSION table)
        -b3                      Sub-sample name for rep 3 of sample B OR column number (0-based in INCLUSION table)
        -b4                      Sub-sample name for rep 4 of sample B OR column number (0-based in INCLUSION table)
        --noVLOW                 Do not use samples with VLOW coverage (default OFF)
        --p_IR                   Filter IR b the p-value of the binomial test (default OFF)
        --paired                 Does a paired comparison (A1 vs B1, A2 vs B2, etc.)
                                   - It uses min_dPSI as the minimum average of each paired dPSI
                                   - It uses min_range as the minumum dPSI for each paired comparison 
    
[GO options]
        --GO                     Generates gene lists for GO analysis (from default gene_ID)
        --GO_file file           To provide an alternative gene ID file for each event (default OFF)
        --use_names              Uses gene names (first column in INCLUSION table)
        --species Hsa/etc        Three letter code for the database 
        --dbDir db               Database directory (default VASTDB)


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

#### SANITY CHECKS
# for paired (repA must be the same as repB)
errPrintDie "If paired comparison, the number of replicates must be the same\n" if (defined $paired) && ($rep1 != $rep2);
errPrintDie "Need to provide at least one replicate for each group\n" if (!defined $a1 or !defined $b1);

# defining default output file name
my ($root)=$ARGV[0]=~/.+\-(.+?)\./;
my $tail = ""; # to be added to the output name
$tail.="-range$min_range" if (defined $min_range); 
$tail.="-noVLOW" if (defined $noVLOW);
$tail.="-p_IR" if (defined $p_IR);
$tail.="-paired" if (defined $paired);
my $out_root="$root-reps$rep1"."_$rep2-dPSI$min_dPSI$tail";
$output_file="DiffAS-$out_root.tab" unless (defined $output_file);


### To get the folder in which the input file is
($folder) = $input_file =~/(.+)\//; # empty if no match (i.e. local folder)
$folder = "." if (!defined $folder);

# prepare to obtain gene IDs for GO analyses
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

open (PSI, $input_file) or errPrintDie "Needs a PSI INCLUSION table\n";
open (O, ">$folder/$output_file") or errPrintDie "Can't open the output file\n"; # output file

### Common for all numbers of replicates
# preparing the head
my $head_row=<PSI>;
chomp($head_row);
my @head=split(/\t/,$head_row);
foreach my $i (6..$#head){
    if ($i%2==0){ # to match sample names with column number
	$a1=$i if ($a1 eq $head[$i]);
	if ($rep1 >= 2){
	    $a2=$i if ($a2 eq $head[$i]);
	}
	if ($rep1 >= 3){
	    $a3=$i if ($a3 eq $head[$i]);
	}
	if ($rep1 >= 4){
	    $a4=$i if ($a4 eq $head[$i]);
	}
	$b1=$i if ($b1 eq $head[$i]);
	if ($rep2 >= 2){
	    $b2=$i if ($b2 eq $head[$i]);
	}
	if ($rep2 >= 3){
	    $b3=$i if ($b3 eq $head[$i]);
	}
	if ($rep2 >= 4){
	    $b4=$i if ($b4 eq $head[$i]);
	}
    }
}

# check that columns provided are 0-based OR
# if names were provided, that all columns were properly matched

errPrintDie "Column numbers do not seem 0-based or conversion did not work properly\n" if ($a1%2 != 0 || ($a2%2 != 0 && $a2 != 0) || ($a3%2 != 0 && $a3 != 0) || ($a4%2 != 0 && $a4 != 0) 
											       || $b1%2 != 0 || ($b2%2 != 0 && $b2 != 0) || ($b3%2 != 0 && $b3 != 0) || ($b4%2 != 0 && $b4 != 0));
errPrintDie "Column numbers do not seem to correspond to INCLUSION samples\n" if (($a1< 6 && $a1 != 0) || ($a2 < 6 && $a2 != 0) || ($a3 < 6 && $a3 != 0) || ($a4 < 6 && $a4 != 0) 
									|| ($b1 < 6 && $b1 != 0) || ($b2 < 6 && $b2 != 0) || ($b3 < 6 && $b3 != 0) || ($b4 < 6 && $b4 != 0));

#print O "$head_row\tdPSI\n"; # it will print the original data + the dPSI of the averages
print O "$head_row\n"; # it will print the original data (for plot later)
# representative names
my $name_A=$head[$a1];
my $name_B=$head[$b1];
$name_A=~s/(.+)\_.+/$1/; # usually the rep number/id is encoded as "_a"
$name_B=~s/(.+)\_.+/$1/;

# Global variables for PSI analysis & GO
my %doneIR_UP;
my %doneIR_DOWN;
my %doneEXSK;
my %doneBG;
my %tally;
$tally{MIC}{DOWN}=0; $tally{MIC}{UP}=0;
$tally{AltEx}{DOWN}=0; $tally{AltEx}{UP}=0;
$tally{IR}{DOWN}=0; $tally{IR}{UP}=0;
$tally{Alt3}{DOWN}=0; $tally{Alt3}{UP}=0;
$tally{Alt5}{DOWN}=0; $tally{Alt5}{UP}=0;

while (<PSI>){
    $_ =~ s/VLOW/N/g if (defined $noVLOW);
    chomp($_);
    my @t=split(/\t/,$_);
    
    next if $t[3] == 0; # removes the internal Alt3 and Alt5 splice sites to avoid double counting
    
    # defines AS type
    my $type="";
    $type="MIC" if ($t[1]=~/EX/ || ($t[5]!~/IR/ && $t[5]!~/Alt/)) && $t[3]<=27;
    $type="AltEx" if ($t[1]=~/EX/ || ($t[5]!~/IR/ && $t[5]!~/Alt/)) && $t[3]>27;
    $type="IR" if $t[1]=~/INT/ || $t[5]=~/IR/;
    $type="Alt3" if $t[1]=~/ALTA/ || $t[5]=~/Alt3/;
    $type="Alt5" if $t[1]=~/ALTD/ || $t[5]=~/Alt5/;
	
    # coverage check (requires good coverage for ALL replicates)
    my $OK_1=0;
    my $OK_2=0;
    if ($rep1 == 1){
	$OK_1=1 if ($t[$a1+1]=~/$Q/);
    }
    elsif ($rep1 == 2){
	$OK_1=1 if ($t[$a1+1]=~/$Q/ && $t[$a2+1]=~/$Q/);
    }
    elsif ($rep1 == 3){
	$OK_1=1 if ($t[$a1+1]=~/$Q/ && $t[$a2+1]=~/$Q/ && $t[$a3+1]=~/$Q/);
    }
    elsif ($rep1 == 4){
	$OK_1=1 if ($t[$a1+1]=~/$Q/ && $t[$a2+1]=~/$Q/ && $t[$a3+1]=~/$Q/ && $t[$a4+1]=~/$Q/);
    }
    if ($rep2 == 1){
	$OK_2=1 if ($t[$b1+1]=~/$Q/);
    }
    elsif ($rep2 == 2){
	$OK_2=1 if ($t[$b1+1]=~/$Q/ && $t[$b2+1]=~/$Q/);
    }
    elsif ($rep2 == 3){
	$OK_2=1 if ($t[$b1+1]=~/$Q/ && $t[$b2+1]=~/$Q/ && $t[$b3+1]=~/$Q/);
    }
    elsif ($rep2 == 4){
	$OK_2=1 if ($t[$b1+1]=~/$Q/ && $t[$b2+1]=~/$Q/ && $t[$b3+1]=~/$Q/ && $t[$b4+1]=~/$Q/);
    }

    next if ($OK_1 ==0 || $OK_2 == 0);
    
    # IR check (only checks the p if the p_IR is active)
    if (($type eq "IR") && (defined $p_IR)){ # only checks the p if the p_IR is active
	my ($p_a1)=$t[$a1+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	next if $p_a1 < 0.05;
	if ($rep1 >= 2){
	    my ($p_a2)=$t[$a2+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    next if $p_a2 < 0.05;
	}
	if ($rep1 >= 3){
	    my ($p_a3)=$t[$a3+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    next if $p_a3 < 0.05;
	}
	if ($rep1 >= 4){
	    my ($p_a4)=$t[$a4+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    next if $p_a4 < 0.05;
	}
	my ($p_b1)=$t[$b1+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	next if $p_b1 < 0.05;
	if ($rep2 >= 2){
	    my ($p_b2)=$t[$b2+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    next if $p_b2 < 0.05;
	}
	if ($rep2 >= 3){
	    my ($p_b3)=$t[$b3+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    next if $p_b3 < 0.05;	
	}
	if ($rep2 >= 4){ 
	    my ($p_b4)=$t[$b4+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	    next if $p_b4 < 0.05;
	}
    }
    
    # NOT PAIRED: gets the average PSI for A and B and the lowest (min) and highest (max) PSI for each replicate
    if (!defined $paired){
	my $PSI_A = "";
	my $min_A = "";
	my $max_A = "";
	my $PSI_B = "";
	my $min_B = "";
	my $max_B = "";

	if ($rep1 == 1) {
	    $PSI_A=sprintf("%.2f",$t[$a1]);
	    $min_A=$t[$a1];
	    $max_A=$t[$a1];
	}
	elsif ($rep1 == 2) {
	    $PSI_A=sprintf("%.2f",($t[$a1]+$t[$a2])/2);
	    $min_A=(sort{$b<=>$a} ($t[$a1],$t[$a2]))[-1];
	    $max_A=(sort{$a<=>$b} ($t[$a1],$t[$a2]))[-1];
	}
	elsif ($rep1 == 3) {
	    $PSI_A=sprintf("%.2f",($t[$a1]+$t[$a2]+$t[$a3])/3);
	    $min_A=(sort{$b<=>$a} ($t[$a1],$t[$a2],$t[$a3]))[-1];
	    $max_A=(sort{$a<=>$b} ($t[$a1],$t[$a2],$t[$a3]))[-1];
	}
	elsif ($rep1 == 4) {
	    $PSI_A=sprintf("%.2f",($t[$a1]+$t[$a2]+$t[$a3]+$t[$a4])/4);
	    $min_A=(sort{$b<=>$a} ($t[$a1],$t[$a2],$t[$a3],$t[$a4]))[-1];
	    $max_A=(sort{$a<=>$b} ($t[$a1],$t[$a2],$t[$a3],$t[$a4]))[-1];
	}
	if ($rep2 == 1) {
	    $PSI_B=sprintf("%.2f",$t[$b1]);
	    $min_B=$t[$b1];
	    $max_B=$t[$b1];
	}
	elsif ($rep2 == 2) {
	    $PSI_B=sprintf("%.2f",($t[$b1]+$t[$b2])/2);
	    $min_B=(sort{$b<=>$a} ($t[$b1],$t[$b2]))[-1];
	    $max_B=(sort{$a<=>$b} ($t[$b1],$t[$b2]))[-1];
	}
	elsif ($rep2 == 3) {
	    $PSI_B=sprintf("%.2f",($t[$b1]+$t[$b2]+$t[$b3])/3);
	    $min_B=(sort{$b<=>$a} ($t[$b1],$t[$b2],$t[$b3]))[-1];
	    $max_B=(sort{$a<=>$b} ($t[$b1],$t[$b2],$t[$b3]))[-1];
	}
	elsif ($rep2 == 4) {
	    $PSI_B=sprintf("%.2f",($t[$b1]+$t[$b2]+$t[$b3]+$t[$b4])/4);
	    $min_B=(sort{$b<=>$a} ($t[$b1],$t[$b2],$t[$b3],$t[$b4]))[-1];
	    $max_B=(sort{$a<=>$b} ($t[$b1],$t[$b2],$t[$b3],$t[$b4]))[-1];
	}
	
	# get dPSI
	my $dPSI = $PSI_B-$PSI_A;
	
	# does the diff AS test:
	if ($dPSI > $min_dPSI && $min_B > $max_A+$min_range){ # if rep1 it will always meet the criteria
	    $tally{$type}{UP}++;
#	    print O "$_\t$dPSI\n"; 
	    print O "$_\n"; # dPSI is not printed so it can the be run with plot
	    
	    # print for GO
	    if (defined$get_GO){
		unless ($use_names){
		    print IR_UP "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_UP{$ID_gene{$t[1]}});
		    $doneIR_UP{$ID_gene{$t[1]}}=1;
		    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
		    $doneEXSK{$ID_gene{$t[1]}}=1;
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
#	    print O "$_\t$dPSI\n";
	    print O "$_\n";
	    
	    #print for GO
	    if (defined$get_GO){
		unless ($use_names){
		    print IR_DOWN "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$ID_gene{$t[1]}});
		    $doneIR_DOWN{$ID_gene{$t[1]}}=1;
		    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
		    $doneEXSK{$ID_gene{$t[1]}}=1;		
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
    else { # if paired
	my $dPSI_pair1;
	my $dPSI_pair2;
	my $dPSI_pair3;
	my $dPSI_pair4;
	my $min_indiv_dPSI;
	my $max_indiv_dPSI;
	my $av_paired_dPSI;

	if ($rep1 == 2){
	    $dPSI_pair1 = $t[$b1]-$t[$a1];
	    $dPSI_pair2 = $t[$b2]-$t[$a2];
	    $min_indiv_dPSI = (sort{$b<=>$a} ($dPSI_pair1,$dPSI_pair2))[-1];
	    $max_indiv_dPSI = (sort{$a<=>$b} ($dPSI_pair1,$dPSI_pair2))[-1];
	    $av_paired_dPSI = sprintf("%.2f",($dPSI_pair1+$dPSI_pair2)/2);
	}
	elsif ($rep1 == 3){
	    $dPSI_pair1 = $t[$b1]-$t[$a1];
	    $dPSI_pair2 = $t[$b2]-$t[$a2];
	    $dPSI_pair3 = $t[$b3]-$t[$a3];
	    $min_indiv_dPSI = (sort{$b<=>$a} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3))[-1];
	    $max_indiv_dPSI = (sort{$a<=>$b} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3))[-1];
	    $av_paired_dPSI = sprintf("%.2f",($dPSI_pair1+$dPSI_pair2+$dPSI_pair3)/3);
	}
	elsif ($rep1 == 4){
	    $dPSI_pair1 = $t[$b1]-$t[$a1];
	    $dPSI_pair2 = $t[$b2]-$t[$a2];
	    $dPSI_pair3 = $t[$b3]-$t[$a3];
	    $dPSI_pair4 = $t[$b4]-$t[$a4];
	    $min_indiv_dPSI = (sort{$b<=>$a} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3,$dPSI_pair4))[-1];
	    $max_indiv_dPSI = (sort{$a<=>$b} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3,$dPSI_pair4))[-1];
	    $av_paired_dPSI = sprintf("%.2f",($dPSI_pair1+$dPSI_pair2+$dPSI_pair3+$dPSI_pair4)/4);
	}
	
	if ($av_paired_dPSI > $min_dPSI && $min_indiv_dPSI > $min_range){ 
	    $tally{$type}{UP}++;
#	    print O "$_\t$av_paired_dPSI\n";
	    print O "$_\n";
	    
	    # print for GO
	    if (defined $get_GO){
		unless ($use_names){
		    print IR_UP "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_UP{$ID_gene{$t[1]}});
		    $doneIR_UP{$ID_gene{$t[1]}}=1;
		    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
		    $doneEXSK{$ID_gene{$t[1]}}=1;
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
#	    print O "$_\t$av_paired_dPSI\n";
	    print O "$_\n";
	    
	    #print for GO
	    if (defined $get_GO){
		unless ($use_names){
		    print IR_DOWN "$ID_gene{$t[1]}\n" if ($type eq "IR") && (!defined $doneIR_DOWN{$ID_gene{$t[1]}});
		    $doneIR_DOWN{$ID_gene{$t[1]}}=1;
		    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && (!defined $doneEXSK{$ID_gene{$t[1]}});
		    $doneEXSK{$ID_gene{$t[1]}}=1;		
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
	    print BG "$ID_gene{$t[1]}\n" if (!defined $doneBG{$ID_gene{$t[1]}});
	    $doneBG{$ID_gene{$t[1]}}=1;
	}
	else {
	    print BG "$t[0]\n" if (!defined $doneBG{$t[0]}) && (defined $t[0]);
	    $doneBG{$t[0]}=1 if (defined $t[0]);
	}
    }
}

my $extras = "";
$extras.=", noVLOW" if (defined $noVLOW);
$extras.=", p_IR" if (defined $p_IR);
$extras.=", paired" if (defined $paired);

print "\n*** Options: dPSI=$min_dPSI, range_dif=$min_range$extras\n";
print "*** Summary statistics:\n";
print "\tAS_TYPE\tHigher_in_$name_A\tHigher_in_$name_B\n";
print "\tMicroexons\t$tally{MIC}{DOWN}\t$tally{MIC}{UP}\n";
print "\tLong_AltEx\t$tally{AltEx}{DOWN}\t$tally{AltEx}{UP}\n";
print "\tIntron_ret\t$tally{IR}{DOWN}\t$tally{IR}{UP}\n";
print "\tAlt_3ss\t$tally{Alt3}{DOWN}\t$tally{Alt3}{UP}\n";
print "\tAlt_5ss\t$tally{Alt5}{DOWN}\t$tally{Alt5}{UP}\n";
print "\n";
