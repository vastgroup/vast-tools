#!/usr/bin/env perl 
# This script takes an annotation template and quantifies the complex PSI for all exons that can match with coverage
# It does NOT do a new call for AS events

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Getopt::Long;

my $dbDir;
my $sp;
my $samLen;
my $verboseFlag;
my $ALL_EXC_EEJ; # non-active variable
my $min_eff_complex = 2;  # cut-off for the minimum number of mappable position a "complex" eej can have (before 1)

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "len=i" => \$samLen,
	   "verbose=i" => \$verboseFlag);

sub verbPrint {
    my $verbMsg = shift;
    if($verboseFlag) {
	chomp($verbMsg);
	print STDERR "[vast combine annot]: $verbMsg\n";
    }
}

die "Needs Species key\n" if !defined($sp);

$COMB="M"; # Only available version

@EEJ=glob("to_combine/*.ee*"); # uses the same input files as COMBI
@EFF=glob("$dbDir/FILES/$sp*-$COMB-*-gDNA.ef*");
die "[vast combine annot error] Needs effective from database!\n" if !@EFF;

###
verbPrint "Loading Mappability for each EEJ and length:\n";
foreach $file (@EFF){
    ($length)=$file=~/COMBI\-[A-Z]\-(\d+?)\-/;
    verbPrint "Loading: $file\tLength: $length\n";
    open (MAPPABILITY, $file);
    while (<MAPPABILITY>){
	chomp;
	@t=split(/\t/);
	($gene,$donor,$acceptor,$donor_coord,$acceptor_coord)=$t[0]=~/(.+?)\-(\d+?)\_(\d+?)\-(\d+?)\_(\d+)/;
	$eej="$gene-$donor-$acceptor";
	$eff{$length}{$eej}=$t[1];

	# keeps the coordinate for each donor/acceptor
	$D_CO{$gene}{$donor}=$donor_coord;
        $A_CO{$gene}{$acceptor}=$acceptor_coord;
    }
    close MAPPABILITY;
}

###
verbPrint "Parsing Template file\n";
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.ANNOT.Template.txt") || die "Can't find the template file for ANNOT\n";
$head=<TEMPLATE>;
chomp($head);
$head_PSIs=$head; # no file with read counts
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $event=$t[1];
    $ALL{$event}=$_; # also acts as holder for all event ids
}
close TEMPLATE;

###
verbPrint "Loading EEJ read counts data\n";
foreach $file (@EEJ){
    my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
    $head_PSIs.="\t$sample\t$sample-Q";
    
    open (EEJ, $file);
    while (<EEJ>){
	chomp;
	@t=split(/\t/);
	$gene=$t[0];
	$eej=$t[1];
	$event="$gene-$eej";

        $reads{$sample}{$event}=$t[2];
        ($donor,$acceptor)=$eej=~/(\d+?)\-(\d+)/;
        $last_acceptor{$gene}=$acceptor if $last_acceptor{$gene}<$acceptor || !$last_acceptor{$gene}; # keeps track of the last used acceptor
        $last_donor{$gene}=$donor if $last_donor{$gene}<$donor || !$last_donor{$gene}; # keeps track of the last used acceptora
    }
    close EEJ;
}

### Setting output files
my $NUM=$#EEJ+1; # number of samples
open (PSIs, ">raw_incl/INCLUSION_LEVELS_ANNOT-$sp$NUM-n.tab");
print PSIs "$head_PSIs\n";


verbPrint "Quantifying PSIs\n";
foreach $event (sort (keys %ALL)){
    print PSIs "$ALL{$event}" if $event;
    
    ($gene,$e,$i1,$i2)=$event=~/(.+?)\-\'*(\d+?\-\d+?)\,\'*(\d+?\-\d+?)\,(\d+?\-\d+)/; # ANNOT have a \= right after
    $eej_exc="$gene-$e"; # EEJ for exclusion
    $eej_inc1="$gene-$i1"; # EEJ for inclusion 1 (upstream)
    $eej_inc2="$gene-$i2"; # EEJ for inclusion 2 (downstream)
    
    @DATA=split(/\t/,$ALL{$event});
    
    # Infer strand (encoded in C1 and C2 coordinates)
    ($c1_from_fullCO,$c2_from_fullCO)=$DATA[4]=~/(\d+)\,.+?\,(\d+)/;
    $strand="+" if $c2_from_fullCO >= $c1_from_fullCO; 
    $strand="-" if $c2_from_fullCO < $c1_from_fullCO; 
    
    # Obtain coordinates for Donor (d2) and Acceptor (a1)
    ($acceptor_coord,$donor_coord)=$DATA[2]=~/\:(\d+?)\-(\d+)/ if $strand eq "+";
    ($donor_coord,$acceptor_coord)=$DATA[2]=~/\:(\d+?)\-(\d+)/ if $strand eq "-";

    foreach $file (@EEJ){
	my $fname = $file;
	$fname =~ s/^.*\///;
	($sample)=$fname=~/^(.*)\..*$/;
	$length = $samLen;
	$exc=$inc1=$inc2=$Rexc=$Rinc1=$Rinc2=0; # empty temporary variables for read counts
	
	# data from the reference EEJ (C1A, AC2, C1C2)
	# corrected read counts
	$exc=$reads{$sample}{$eej_exc}/$eff{$length}{$eej_exc} if $eff{$length}{$eej_exc};
	$inc1=$reads{$sample}{$eej_inc1}/$eff{$length}{$eej_inc1} if $eff{$length}{$eej_inc1};
	$inc2=$reads{$sample}{$eej_inc2}/$eff{$length}{$eej_inc2} if $eff{$length}{$eej_inc2};
	# raw read counts
	$Rexc=$reads{$sample}{$eej_exc} if $eff{$length}{$eej_exc};
	$Rinc1=$reads{$sample}{$eej_inc1} if $eff{$length}{$eej_inc1};
	$Rinc2=$reads{$sample}{$eej_inc2} if $eff{$length}{$eej_inc2};
	
	# Donor1 (d1), Donor2 (d2), Acceptor1 (a1) and Acceptor2 (a2)
	($d1,$a2)=$eej_exc=~/\-(\d+?)\-(\d+)/;
	($a1)=$eej_inc1=~/\-\d+?\-(\d+)/;
	($d2)=$eej_inc2=~/\-(\d+?)\-\d+/;

	#### To calculate "complex" PSIs
	$Rexc1C=$Rexc2C=$exc1C=$exc2C=$inc1C=$inc2C=$excC=$Rinc1C=$Rinc2C=$RexcC=0; # empty temporary variables for read counts for each sample
	
	#### Quantifying COMPLEX READS: doing it "very" locally (+/-5 acc/donors)
	### Inclusion reads
	for $i ($d1-5..$d2-1){
	    if ((($D_CO{$gene}{$i} < $acceptor_coord && $strand eq "+") || ($D_CO{$gene}{$i} > $acceptor_coord && $strand eq "-")) && $D_CO{$gene}{$i} && $i != $d1 && $i>=0){
		$temp_eej="$gene-$i-$a1";
		if ($eff{$length}{$temp_eej} >= $min_eff_complex){
		    $inc1C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
		    $Rinc1C+=$reads{$sample}{$temp_eej};
		} 
	    }
	}
	for $i ($a1+1..$a2+5){
	    if ((($A_CO{$gene}{$i} > $donor_coord && $strand eq "+") || ($A_CO{$gene}{$i} < $donor_coord && $strand eq "-")) && $A_CO{$gene}{$i} && $i != $a2 && $i<=$last_acceptor{$gene}){
		$temp_eej="$gene-$d2-$i";
		if ($eff{$length}{$temp_eej} >= $min_eff_complex){
		    $inc2C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
		    $Rinc2C+=$reads{$sample}{$temp_eej};
		}
	    }
	}
	### Exclusion reads
	### COMBI-like:
	### Exclusion reads (It does NOT take all EEJs around the alternative exon, but only those including C1 or C2.)
	if (!$ALL_EXC_EEJ){
	    for $i ($d1-5..$d1-1){
		if ($i>=0){
		    $temp_eej="$gene-$i-$a2";
		    if ($eff{$length}{$temp_eej} >= $min_eff_complex){
			$exc1C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
			$Rexc1C+=$reads{$sample}{$temp_eej};
		    }
		}
	    }
	    for $i ($d1+1..$d1+5){
		if (($D_CO{$gene}{$i} < $A_CO{$gene}{$a1} && $strand eq "+") || ($D_CO{$gene}{$i} > $A_CO{$gene}{$a1} && $strand eq "-") && $i <= $last_donor{$gene}){
		    $temp_eej="$gene-$i-$a2";
		    if ($eff{$length}{$temp_eej} >= $min_eff_complex){
			$exc1C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
			$Rexc1C+=$reads{$sample}{$temp_eej};
		    }
		}
	    }
	    for $i ($a2+1..$a2+5){
		if ($i <= $last_acceptor{$gene}){
		    $temp_eej="$gene-$d1-$i";
		    if ($eff{$length}{$temp_eej} >= $min_eff_complex){
			$exc2C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
			$Rexc2C+=$reads{$sample}{$temp_eej};
		    }
		}
	    }
	    for $i ($a2-5..$a2-1){
		if (($A_CO{$gene}{$i} > $D_CO{$gene}{$d2} && $strand eq "+") || ($A_CO{$gene}{$i} < $D_CO{$gene}{$d2} && $strand eq "-") && $i >= 0){
		    $temp_eej="$gene-$d1-$i";
		    if ($eff{$length}{$temp_eej} >= $min_eff_complex){
			$exc2C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
			$Rexc2C+=$reads{$sample}{$temp_eej};
		    }
		}
	    }
            $excC=$exc1C+$exc2C;
	    $RexcC=$Rexc1C+$Rexc2C;
	}
	### ALL exclusion (prone to more false positives)
	elsif ($ALL_EXC_EEJ){
	    for $i ($d1-5..$d2-1){ # The only true ANNOT-specific thing
		for $j ($a1+1..$a2+5){
		    if ((($D_CO{$gene}{$i} < $acceptor_coord && $A_CO{$gene}{$j} > $donor_coord && $strand eq "+") || ($D_CO{$gene}{$i} > $acceptor_coord && $A_CO{$gene}{$j} < $donor_coord && $strand eq "-")) && $D_CO{$gene}{$i} && $A_CO{$gene}{$j} && ($i != $d1 || $j != $a2) && $i >= 0 && $j <= $last_acceptor{$gene}){ # either of the two or both are not the cannonical
			$temp_eej="$gene-$i-$j";
			if ($eff{$length}{$temp_eej} >= $min_eff_complex){
			    $excC+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
			    $RexcC+=$reads{$sample}{$temp_eej};
			}
		    }
		}
	    }
	}
	############# END ANNOT QUANT


#### Sum of corrected reads
	$inc1F=$inc1C+$inc1;
	$inc2F=$inc2C+$inc2;
	$excF=$excC+$exc;
	# Correction by maximum mappability
	$Creads_exc=$excF*($length-15);
	$Creads_inc1=$inc1F*($length-15);
	$Creads_inc2=$inc2F*($length-15);
	$total_Creads=$Creads_exc+$Creads_inc1+$Creads_inc2;
	
# Sum of actual reads
	$reads_inc1=$Rinc1C+$Rinc1;
	$reads_inc2=$Rinc2C+$Rinc2;
	$reads_exc=$RexcC+$Rexc;
	$total_reads=$reads_exc+$reads_inc1+$reads_inc2;

#### Coverage scores. Score 1: Using all raw reads
	if (($reads_exc >= 20 || (($reads_inc1 >= 20 && $reads_inc2 >= 15) || ($reads_inc1 >= 15 && $reads_inc2 >= 20))) && $total_reads>=100){
	    $Q="SOK";
	}
	elsif (($reads_exc >= 20 || (($reads_inc1 >= 20 && $reads_inc2 >= 15) || ($reads_inc1 >= 15 && $reads_inc2 >= 20))) && $total_reads<100){
	    $Q="OK";
	}
	elsif ($reads_exc >= 15 || (($reads_inc1 >= 15 && $reads_inc2 >= 10) || ($reads_inc1 >= 10 && $reads_inc2 >= 15))){
	    $Q="LOW";
	}
	elsif ($reads_exc >= 10 || (($reads_inc1 >= 10 && $reads_inc2 >= 5) || ($reads_inc1 >= 5 && $reads_inc2 >= 10))){
	    $Q="VLOW";
	}
	else {
	    $Q="N";
	}
	
#### Score 2: Using all corrected reads
	if (($Creads_exc >= 20 || (($Creads_inc1 >= 20 && $Creads_inc2 >= 15) || ($Creads_inc1 >= 15 && $Creads_inc2 >= 20))) && $total_Creads>=100){
	    $Q.=",SOK";
	}
	elsif (($Creads_exc >= 20 || (($Creads_inc1 >= 20 && $Creads_inc2 >= 15) || ($Creads_inc1 >= 15 && $Creads_inc2 >= 20))) && $total_Creads<100){
	    $Q.=",OK";
	}
	elsif ($Creads_exc >= 15 || (($Creads_inc1 >= 15 && $Creads_inc2 >= 10) || ($Creads_inc1 >= 10 && $Creads_inc2 >= 15))){
	    $Q.=",LOW";
	}
	elsif ($Creads_exc >= 10 || (($Creads_inc1 >= 10 && $Creads_inc2 >= 5) || ($Creads_inc1 >= 5 && $Creads_inc2 >= 10))){
	    $Q.=",VLOW";
	}
	else {
	    $Q.=",N";
	}
#### Score 3: Using simple (=reference, C1A, AC2, C1C2) raw reads
	$total_simple_reads=$Rexc+$Rinc1+$Rinc2;
	if (($Rexc >= 20 || (($Rinc1 >= 20 && $Rinc2 >= 15) || ($Rinc1 >= 15 && $Rinc2 >= 20))) && $total_simple_reads >= 100){
	    $Qs="SOK";
	    $Q.=",SOK";
	}
	elsif (($Rexc >= 20 || (($Rinc1 >= 20 && $Rinc2 >= 15) || ($Rinc1 >= 15 && $Rinc2 >= 20))) && $total_simple_reads < 100){
	    $Qs="OK";
	    $Q.=",OK";
	}
	elsif ($Rexc >= 15 || (($Rinc1 >= 15 && $Rinc2 >= 10) || ($Rinc1 >= 10 && $Rinc2 >= 15))){
	    $Qs="LOW";
	    $Q.=",LOW";
	}
	elsif ($Rexc >= 10 || (($Rinc1 >= 10 && $Rinc2 >= 5) || ($Rinc1 >= 5 && $Rinc2 >= 10))){
	    $Qs="VLOW";
	    $Q.=",VLOW";
	}
	else {
	    $Qs="N";
	    $Q.=",N";
	}
###

### Score 4: Calculate imbalance between inclusion EEJs (OK<B1<B2; Bl=not enough inclusion reads):
	if ($Creads_inc1 && $Creads_inc2){
	    if (($Creads_inc1/$Creads_inc2 > 2 && $Creads_inc1/$Creads_inc2 <= 5) || ($Creads_inc2/$Creads_inc1 > 2 && $Creads_inc2/$Creads_inc1 <= 5)){
		$Q.=",B1";
		$Qs.=",B1";
	    }
	    elsif (($Creads_inc1/$Creads_inc2 > 5) || ($Creads_inc2/$Creads_inc1 > 5)){
		$Q.=",B2";
		$Qs.=",B2";
	    }
	    else {
		$Q.=",OK";
		$Qs.=",OK";
	    }
	}
	else {
	    if (!$Creads_inc1){
		$Q.=",B2" if $Creads_inc2>=20;
		$Q.=",Bl" if ($Creads_inc2<20 && $Creads_inc2>=15);
		$Q.=",Bn" if $Creads_inc2<15 && $Creads_inc2>0;
		$Qs.=",B2" if $Creads_inc2>=20;
		$Qs.=",Bl" if ($Creads_inc2<20 && $Creads_inc2>=15);
		$Qs.=",Bn" if $Creads_inc2<15 && $Creads_inc2>0;
	    }
	    if (!$Creads_inc2){
		$Q.=",B2" if $Creads_inc1>=20;
		$Q.=",Bl" if ($Creads_inc1<20 && $Creads_inc1>=15);
		$Q.=",Bn" if $Creads_inc1<15;
		$Qs.=",B2" if $Creads_inc1>=20;
		$Qs.=",Bl" if ($Creads_inc1<20 && $Creads_inc1>=15);
		$Qs.=",Bn" if $Creads_inc1<15;
	    }
	}
### 
	
### Score 5: Complexity score (S<C1<C2<C3)
	$from_C=$excC+$inc1C+$inc2C;
	$from_S=$exc+$inc1+$inc2;

        if ($from_C > ($from_C+$from_S)/2) {$Q.=",C3";}
        elsif ($from_C > ($from_C+$from_S)/5 && $from_C <= ($from_C+$from_S)/2){$Q.=",C2";}
        elsif ($from_C > ($from_C+$from_S)/20 && $from_C <= ($from_C+$from_S)/5){$Q.=",C1";}
        else {$Q.=",S";}
###

### Final PSI value  
	$PSI_complex=sprintf ("%.2f", (100*($inc1F+$inc2F))/($inc1F+$inc2F+(2*$excF))) if ($inc1F+$inc2F+(2*$excF))>0;
	$PSI_complex="NA" if ($inc1F+$inc2F+(2*$excF))==0;

	### DIFF OUTPUT ADDITION TO QUAL SCORE!  --TSW
	### Essentially adding the expected number of reads re-distributed to INC or EXC after normalization..
	### These values are added to the qual score and used to infer the posterior distribution
	my $totalN = $total_reads; # no point in re-assigning really.
	my($pPSI, $exValOfInc, $exValOfExc) = (0, 0, 0);
	unless($PSI_complex eq "NA") {
	    $pPSI = $PSI_complex / 100;
	    $exValOfInc = sprintf("%.2f", $pPSI * $totalN);
	    $exValOfExc = sprintf("%.2f", (1-$pPSI) * $totalN);
	}
	$Q.= "\@$exValOfInc,$exValOfExc";

	print PSIs "\t$PSI_complex\t$Q" if $event;
    }
    print PSIs "\n" if $event;
}
close PSIs;

