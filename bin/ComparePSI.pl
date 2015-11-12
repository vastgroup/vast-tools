#!/usr/bin/perl
### General script to get event based on dPSI differences 

### v1: (04/09/15)
###### For N=1,2 or 3 replicates

### v2 (14/10/15): 
###### It allows different number of replicates in each group
###### It allows to provide a file to generate gene_ID sets for GO
###### It allows for 4 replicates

### v3 (07/11/15): 
###### It allows word definition of samples to be compared
###### It allows paired comparisons of dPSIs

use Getopt::Long;

### Setting global variables:
$Q="O[KW]\,.+?\,.+?\,.+?\,.+?\@"; # quality search
$input_file=$ARGV[0];
$min_dPSI=15; # min dPSI difference (def=15)
$rep1=2; # number of replicates per type (def=2);
$rep2=2; # number of replicates per type (def=2);
$min_range=5; # min dPSI between ranges
$noVLOW=0;
$p_IR=0;
$ID_file=0;
$paired=0;

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
			  "output=s"    => \$output_file,
			  "o=s"         => \$output_file,
			  "help" => \$helpFlag,
			  "p_IR" => \$p_IR,
			  "GO=s" => \$ID_file,
			  "paired" => \$paired,
			  "noVLOW" => \$noVLOW
    );

# defining default output file name
($root)=$ARGV[0]=~/.+\-(.+?)\./;
$tail.="-range$min_range" if $min_range; # to be added to the output
$tail.="-noVLOW" if $noVLOW;
$tail.="-p_IR" if $p_IR;
$tail.="-paired" if $paired;
$out_root="$root-reps$rep1"."_$rep2-dPSI$min_dPSI$tail";
$output_file="DiffAS-$out_root.tab" if !$output_file;

if (!defined($ARGV[0]) || $helpFlag){
    die "\nUsage: CompareAS_General_v2.pl INCLUSION_LEVELS_FULL-root.tab --min_dPSI min_dPSI --rep1 N_rep1 --rep2 N_rep2 -a1 col_a1  -b1 col_b2 [options]

Compare two sample sets with 1,2, 3 or 4 replicates.

OPTIONS: 
        --min_dPSI i             Minimum delta PSI of the averages (default 15)
        --repA i                 Number of replicates for group A (default 2)
        --repB i                 Number of replicates for group B (default 2)
        -a1                      Column # for rep 1 of sample A (0-based in INCLUSION table); also sub-sample name
        -a2                      Column # for rep 2 of sample A (0-based in INCLUSION table); also sub-sample name
        -a3                      Column # for rep 3 of sample A (0-based in INCLUSION table); also sub-sample name
        -a4                      Column # for rep 4 of sample A (0-based in INCLUSION table); also sub-sample name
        -b1                      Column # for rep 1 of sample B (0-based in INCLUSION table); also sub-sample name
        -b2                      Column # for rep 2 of sample B (0-based in INCLUSION table); also sub-sample name
        -b3                      Column # for rep 3 of sample B (0-based in INCLUSION table); also sub-sample name
        -b4                      Column # for rep 4 of sample B (0-based in INCLUSION table); also sub-sample name
        --min_range i            Minimum distance between the ranges of both groups (default 5)
        --noVLOW                 Do not use samples with VLOW coverage (default OFF)
        --p_IR                   Filter IR b the p-value of the binomial test (default OFF)
        --paired                 Does a paired comparison (A1 vs B1, A2 vs B2, etc.)
                                   - It uses min_dPSI as the minimum average of each paired dPSI
                                   - It uses min_range as the minumum dPSI for each paired comparison 
        --GO file                Generates gene lists for GO when provided a key file (e.g. Hsa.Event-Gene.IDs.txt)
        -o, --output             Output file name (default based on option parameters)


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

# prepare to obtain gene IDs for GO analyses
if ($ID_file){
    open (KEY, $ID_file) || die "Can't open ID keys file ($ID_file)\n";
    while (<KEY>){
	chomp;
	@t=split(/\t/);
	$ID_gene{$t[0]}=$t[1];
    }
    close KEY;
    open (BG, ">BG-$out_root.txt");
    open (IR_UP, ">IR_UP-$out_root.txt");
    open (IR_DOWN, ">IR_DOWN-$out_root.txt");
    open (EXSK, ">AltEx-$out_root.txt");
}

open (PSI, $input_file) || die "Needs a PSI INCLUSION table\n";
open (O, ">$output_file"); # output file

### Common for all numbers of replicates (v2 -- 14/10/15)
# preparing the head
$h=<PSI>;
chomp($h);
@head=split(/\t/,$h);
foreach $i (6..$#head){
    if ($i%2==0){
	$a1=$i if ($a1 eq $head[$i]);
	$a2=$i if ($a2 eq $head[$i]);
	$a3=$i if ($a3 eq $head[$i]);
	$a4=$i if ($a4 eq $head[$i]);
	$b1=$i if ($b1 eq $head[$i]);
	$b2=$i if ($b2 eq $head[$i]);
	$b3=$i if ($b3 eq $head[$i]);
	$b4=$i if ($b4 eq $head[$i]);
    }
}

#### SANITY CHECKS
# for paired (repA must be the same as repB)
die "*** If paired comparison, the number of replicates must be the same\n" if $paired && $rep1!=$rep2;

# check that columns provided are 0-based (added 23/09/15) --MI
# or, if names were provided, that all columns were properly matched
die "*** Column numbers do not seem 0-based or conversion did not work properly\n" if ($a1%2!=0 || $a2%2!=0 || $a3%2!=0 || $a4%2!=0 || $b1%2!=0 || $b2%2!=0 || $b3%2!=0 || $b4%2!=0);


print O "$h\tdPSI\n"; # it will print the original data + the dPSI of the averages
# representative names
$name_A=$head[$a1];
$name_B=$head[$b1];
$name_A=~s/(.+)\_.+/\1/; # usually the rep number/id is encoded as "_a"
$name_B=~s/(.+)\_.+/\1/;

while (<PSI>){
    s/VLOW/N/g if $noVLOW;
    chomp;
    @t=split(/\t/);
    
    next if $t[3]==0; # removes the internal Alt3 and Alt5 splice sites
    
    # defines AS type
    $type="";
    $type="MIC" if ($t[1]=~/EX/ || ($t[5]!~/IR/ && $t[5]!~/Alt/)) && $t[3]<=27;
    $type="AltEx" if ($t[1]=~/EX/ || ($t[5]!~/IR/ && $t[5]!~/Alt/)) && $t[3]>27;
    $type="IR" if $t[1]=~/INT/ || $t[5]=~/IR/;
    $type="Alt3" if $t[1]=~/ALTA/ || $t[5]=~/Alt3/;
    $type="Alt5" if $t[1]=~/ALTD/ || $t[5]=~/Alt5/;
	
    # coverage check (requires good coverage for ALL replicates)
    $OK_1=$OK_2="";
    $OK_1=1 if ($t[$a1+1]=~/$Q/ && $rep1 == 1);
    $OK_1=1 if ($t[$a1+1]=~/$Q/ && $t[$a2+1]=~/$Q/ && $rep1 == 2);
    $OK_1=1 if ($t[$a1+1]=~/$Q/ && $t[$a2+1]=~/$Q/ && $t[$a3+1]=~/$Q/ && $rep1 == 3);
    $OK_1=1 if ($t[$a1+1]=~/$Q/ && $t[$a2+1]=~/$Q/ && $t[$a3+1]=~/$Q/ && $t[$a4+1]=~/$Q/ && $rep1 == 4);
    $OK_2=1 if ($t[$b1+1]=~/$Q/ && $rep2 == 1);
    $OK_2=1 if ($t[$b1+1]=~/$Q/ && $t[$b2+1]=~/$Q/ && $rep2 == 2);
    $OK_2=1 if ($t[$b1+1]=~/$Q/ && $t[$b2+1]=~/$Q/ && $t[$b3+1]=~/$Q/ && $rep2 == 3);
    $OK_2=1 if ($t[$b1+1]=~/$Q/ && $t[$b2+1]=~/$Q/ && $t[$b3+1]=~/$Q/ && $t[$b4+1]=~/$Q/ && $rep2 == 4);
    next if !$OK_1 || !$OK_2;
    
    # IR check (only checks the p if the p_IR is active)
    if ($type eq "IR" && $p_IR){ # only checks the p if the p_IR is active
	($p_a1)=$t[$a1+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	($p_a2)=$t[$a2+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	($p_a3)=$t[$a3+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	($p_a4)=$t[$a4+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	($p_b1)=$t[$b1+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	($p_b2)=$t[$b2+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	($p_b3)=$t[$b3+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	($p_b4)=$t[$b4+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
	next if $p_a1<0.05;
	next if $p_a2<0.05 && $rep1>=2;
	next if $p_a3<0.05 && $rep1>=3;
	next if $p_a4<0.05 && $rep1>=4;
	next if $p_b1<0.05;
	next if $p_b2<0.05 && $rep2>=2;
	next if $p_b3<0.05 && $rep2>=3;
	next if $p_b4<0.05 && $rep2>=4;
    }
    
    # NOT PAIRED: gets the average PSI for A and B and the lowest (min) and highest (max) PSI for each replicate
    if (!$paired){
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
	$dPSI=$PSI_B-$PSI_A;
	
	# does the diff AS test:
	if ($dPSI > $min_dPSI && $min_B > $max_A+$min_range){ # if rep1 it will always meet the criteria
	    $tally{$type}{UP}++;
	    print O "$_\t$dPSI\n";
	    
	    # print for GO
	    print IR_UP "$ID_gene{$t[1]}\n" if $type eq "IR" && !$doneIR_UP{$ID_gene{$t[1]}};
	    $doneIR_UP{$ID_gene{$t[1]}}=1;
	    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && !$doneEXSK{$ID_gene{$t[1]}};
	    $doneEXSK{$ID_gene{$t[1]}}=1;
	}
	if ($dPSI < -$min_dPSI && $min_A > $max_B+$min_range){
	    $tally{$type}{DOWN}++;
	    print O "$_\t$dPSI\n";
	    
	    #print for GO
	    print IR_DOWN "$ID_gene{$t[1]}\n" if $type eq "IR" && !$doneIR_DOWN{$ID_gene{$t[1]}};
	    $doneIR_DOWN{$ID_gene{$t[1]}}=1;
	    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && !$doneEXSK{$ID_gene{$t[1]}};
	    $doneEXSK{$ID_gene{$t[1]}}=1;
	}
    }
    else { # if paired
	if ($rep1 == 2){
	    $dPSI_pair1=$t[$b1]-$t[$a1];
	    $dPSI_pair2=$t[$b2]-$t[$a2];
	    $min_indiv_dPSI=(sort{$b<=>$a} ($dPSI_pair1,$dPSI_pair2))[-1];
	    $max_indiv_dPSI=(sort{$a<=>$b} ($dPSI_pair1,$dPSI_pair2))[-1];
	    $av_paired_dPSI=sprintf("%.2f",($dPSI_pair1+$dPSI_pair2)/2);
	}
	elsif ($rep1 == 3){
	    $dPSI_pair1=$t[$b1]-$t[$a1];
	    $dPSI_pair2=$t[$b2]-$t[$a2];
	    $dPSI_pair3=$t[$b3]-$t[$a3];
	    $min_indiv_dPSI=(sort{$b<=>$a} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3))[-1];
	    $max_indiv_dPSI=(sort{$a<=>$b} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3))[-1];
	    $av_paired_dPSI=sprintf("%.2f",($dPSI_pair1+$dPSI_pair2+$dPSI_pair3)/3);
	}
	elsif ($rep1 == 4){
	    $dPSI_pair1=$t[$b1]-$t[$a1];
	    $dPSI_pair2=$t[$b2]-$t[$a2];
	    $dPSI_pair3=$t[$b3]-$t[$a3];
	    $dPSI_pair4=$t[$b4]-$t[$a4];
	    $min_indiv_dPSI=(sort{$b<=>$a} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3,$dPSI_pair4))[-1];
	    $max_indiv_dPSI=(sort{$a<=>$b} ($dPSI_pair1,$dPSI_pair2,$dPSI_pair3,$dPSI_pair4))[-1];
	    $av_paired_dPSI=sprintf("%.2f",($dPSI_pair1+$dPSI_pair2+$dPSI_pair3+$dPSI_pair4)/4);
	}
	
	if ($av_paired_dPSI > $min_dPSI && $min_indiv_dPSI > $min_range){ 
	    $tally{$type}{UP}++;
	    print O "$_\t$av_paired_dPSI\n";
	    
	    # print for GO
	    print IR_UP "$ID_gene{$t[1]}\n" if $type eq "IR" && !$doneIR_UP{$ID_gene{$t[1]}};
	    $doneIR_UP{$ID_gene{$t[1]}}=1;
	    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && !$doneEXSK{$ID_gene{$t[1]}};
	    $doneEXSK{$ID_gene{$t[1]}}=1;
	}
	if ($av_paired_dPSI < -$min_dPSI && $max_indiv_dPSI < -$min_range){ 
	    $tally{$type}{DOWN}++;
	    print O "$_\t$av_paired_dPSI\n";
	    
	    #print for GO
	    print IR_DOWN "$ID_gene{$t[1]}\n" if $type eq "IR" && !$doneIR_DOWN{$ID_gene{$t[1]}};
	    $doneIR_DOWN{$ID_gene{$t[1]}}=1;
	    print EXSK "$ID_gene{$t[1]}\n" if ($type eq "AltEx" || $type eq "MIC") && !$doneEXSK{$ID_gene{$t[1]}};
	    $doneEXSK{$ID_gene{$t[1]}}=1;
	}
    }
    # prints out the genes for BG
    print BG "$ID_gene{$t[1]}\n" unless $doneBG{$ID_gene{$t[1]}};
    $doneBG{$ID_gene{$t[1]}}=1;
}

print "AS_TYPE\tHigher_in_$name_A\tHigher_in_$name_B\n";
print "Microexons\t$tally{MIC}{DOWN}\t$tally{MIC}{UP}\n";
print "Long_AltEx\t$tally{AltEx}{DOWN}\t$tally{AltEx}{UP}\n";
print "Intron_ret\t$tally{IR}{DOWN}\t$tally{IR}{UP}\n";
print "Alt_3ss\t$tally{Alt3}{DOWN}\t$tally{Alt3}{UP}\n";
print "Alt_5ss\t$tally{Alt5}{DOWN}\t$tally{Alt5}{UP}\n";
