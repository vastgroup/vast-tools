#!/usr/bin/env perl 
# This script takes an a posteriori template and uses it to get PSIs for exactly those events. 
# It does NOT do a new call for AS events

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Getopt::Long;

my $dbDir;
my $sp;
my $samLen;
my $verboseFlag;
my $legacyFlag;
my $min_eff_complex=2; # cut-off for the minimum number of mappable position a "complex" eej can have (before 1)
my $use_all_excl_eej;
my $extra_eej; # only used if $ALL_EXC_EEJ is active

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "len=i" => \$samLen, "extra_eej=i" => \$extra_eej, "use_all_excl_eej" => \$use_all_excl_eej,
			  "verbose=i" => \$verboseFlag, "legacy" => \$legacyFlag);

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine combi]: $verbMsg\n";
  }
}

die "Needs Species key\n" if !defined($sp);

$COMB="M"; # Only available version

@EEJ=glob("to_combine/*.eej2"); # from v2, only works with eej2
@EFF_ns=glob("$dbDir/FILES/$sp*-$COMB-*-gDNA.ef*"); die "[vast combine combi error] Needs effective (strand-unspecific) from database!\n" if !@EFF_ns;
@EFF_ss=glob("$dbDir/FILES/$sp*-$COMB-*-gDNA-SS.ef*"); die "[vast combine combi error] Needs effective (strand-specific) from database!\n" if !@EFF_ss;

###
verbPrint "Loading strand-unspecific mappability for each EEJ and length:\n";
foreach $file (@EFF_ns){
    ($length)=$file=~/COMBI\-[A-Z]\-(\d+?)\-/;
    verbPrint "Loading: $file\tLength: $length\n";
    open (MAPPABILITY, $file);
    while (<MAPPABILITY>){
	chomp;
	@t=split(/\t/);
	($gene,$donor,$acceptor,$donor_coord,$acceptor_coord)=$t[0]=~/(.+?)\-(\d+?)\_(\d+?)\-(\d+?)\_(\d+)/;
	$eej="$gene-$donor-$acceptor";
	$eff_ns{$length}{$eej}=$t[1];

	# keeps the coordinate for each donor/acceptor
	$D_CO_ns{$gene}{$donor}=$donor_coord;
        $A_CO_ns{$gene}{$acceptor}=$acceptor_coord;
    }
    close MAPPABILITY;
}

verbPrint "Loading strand-specific mappability for each EEJ and length:\n";
foreach $file (@EFF_ss){
    ($length)=$file=~/COMBI\-[A-Z]\-(\d+?)\-/;
    verbPrint "Loading: $file\tLength: $length\n";
    open (MAPPABILITY, $file);
    while (<MAPPABILITY>){
	chomp;
	@t=split(/\t/);
	($gene,$donor,$acceptor,$donor_coord,$acceptor_coord)=$t[0]=~/(.+?)\-(\d+?)\_(\d+?)\-(\d+?)\_(\d+)/;
	$eej="$gene-$donor-$acceptor";
	$eff_ss{$length}{$eej}=$t[1];

	# keeps the coordinate for each donor/acceptor
	$D_CO_ss{$gene}{$donor}=$donor_coord;
	$A_CO_ss{$gene}{$acceptor}=$acceptor_coord;
    }
    close MAPPABILITY;
}

###
verbPrint "Parsing Template file\n";
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.COMBI.Template.txt") || die "Can't find the template file for COMBI\n";
$head=<TEMPLATE>;
chomp($head);
$head_PSIs=$head_ReadCounts=$head; # setting headings for both output files
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $event=$t[1];
    $ALL{$event}=$_; # also acts as holder for all event ids
}
close TEMPLATE;

###
my %is_ss=();  # get strand-specific samples / groups
verbPrint "Loading EEJ read counts data\n";
foreach $file (@EEJ){
    my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
    $head_PSIs.="\t$sample\t$sample-Q";
    $head_ReadCounts.="\t$sample-Re\t$sample-Ri1\t$sample-Ri2\t$sample-ReC\t$sample-Ri1C\t$sample-Ri2C\t$sample-Q";
    
    unless(-e "to_combine/${sample}.info"){ verbPrint "   $sample: do not find to_combine/${sample}.info. Sample will be treated as being not strand-specific.";
    }else{
    	open(my $fh_info,"to_combine/${sample}.info") or die "$!"; my $line=<$fh_info>; close($fh_info);
    	my @fs=split("\t",$line);
    	if($fs[@fs-2] eq "-SS"){
	    $is_ss{$sample}=1;
	    verbPrint "   $sample: found to_combine/${sample}.info. Sample will be treated as being strand-specific."
    	}
	else{
	    verbPrint "   $sample: found to_combine/${sample}.info. Sample will be treated as being not strand-specific."
	}
    }
    
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
	$last_donor{$gene}=$donor if $last_donor{$gene}<$donor || !$last_donor{$gene}; # keeps track of the last used donor
    }
    close EEJ;
}

### Setting output files
$NUM=$#EEJ+1; # number of samples
open (PSIs, ">raw_incl/INCLUSION_LEVELS_COMBI-$sp$NUM-n.tab");
open (COUNTs, ">raw_reads/RAW_READS_COMBI-$sp$NUM-n.tab");
print PSIs "$head_PSIs\n";
print COUNTs "$head_ReadCounts\n";

my ($eff_href,$D_CO_href,$A_CO_href);
verbPrint "Quantifying PSIs\n";
foreach $event (sort (keys %ALL)){
    print PSIs "$ALL{$event}" if $event;
    print COUNTs "$ALL{$event}" if $event;
    
    ($gene,$e,$i1,$i2)=$event=~/(.+?)\-\'*(\d+?\-\d+?)\,\'*(\d+?\-\d+?)\,(\d+?\-\d+)/;
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
	
	# set mappability correction accordingly
	if($is_ss{$sample}){$eff_href=\%eff_ss;$D_CO_href=\%D_CO_ss;$A_CO_href=\%A_CO_ss;}else{$eff_href=\%eff_ns;$D_CO_href=\%D_CO_ns;$A_CO_href=\%A_CO_ns;}
	
	# data from the reference EEJ (C1A, AC2, C1C2)
	# corrected read counts
	$exc=$reads{$sample}{$eej_exc}/$eff_href->{$length}{$eej_exc} if $eff_href->{$length}{$eej_exc};
	$inc1=$reads{$sample}{$eej_inc1}/$eff_href->{$length}{$eej_inc1} if $eff_href->{$length}{$eej_inc1};
	$inc2=$reads{$sample}{$eej_inc2}/$eff_href->{$length}{$eej_inc2} if $eff_href->{$length}{$eej_inc2};
	# raw read counts
	$Rexc=$reads{$sample}{$eej_exc} if $eff_href->{$length}{$eej_exc};
	$Rinc1=$reads{$sample}{$eej_inc1} if $eff_href->{$length}{$eej_inc1};
	$Rinc2=$reads{$sample}{$eej_inc2} if $eff_href->{$length}{$eej_inc2};
	# Simple PSI (only C1A, AC2 and C1C2 EEJs); Not used
	$PSI_simple=sprintf("%.2f",100*($inc1+$inc2)/(($exc*2)+($inc1+$inc2))) if (($exc*2)+($inc1+$inc2))>0;
	$PSI_simple="NA" if (($exc*2)+($inc1+$inc2))==0;
	
	# Donor1 (d1), Donor2 (d2), Acceptor1 (a1) and Acceptor2 (a2)
	($d1,$a2)=$eej_exc=~/\-(\d+?)\-(\d+)/;
	($a1)=$eej_inc1=~/\-\d+?\-(\d+)/;
	($d2)=$eej_inc2=~/\-(\d+?)\-\d+/;
	
	#### To calculate "complex" PSIs
	$Rexc1C=$Rexc2C=$exc1C=$exc2C=$inc1C=$inc2C=$excC=$Rinc1C=$Rinc2C=$RexcC=0; # empty temporary variables for read counts for each sample

	### Inclusion reads
	for $i (0..$d2-1){ # loops through all donors in the gene from the first to d2-1
	    if ((($D_CO_href->{$gene}{$i} < $acceptor_coord && $strand eq "+") || ($D_CO_href->{$gene}{$i} > $acceptor_coord && $strand eq "-")) && $i != $d1){
		$temp_eej="$gene-$i-$a1";
		if ($eff_href->{$length}{$temp_eej} >= $min_eff_complex){
		    $inc1C+=$reads{$sample}{$temp_eej}/$eff_href->{$length}{$temp_eej};
		    $Rinc1C+=$reads{$sample}{$temp_eej};
		}
	    }
	}
	for $i ($a1+1..$last_acceptor{$gene}){ # loops through all acceptors in the gene from a1+1 to the last
	    if ((($A_CO_href->{$gene}{$i} > $donor_coord && $strand eq "+") || ($A_CO_href->{$gene}{$i} < $donor_coord && $strand eq "-")) && $i != $a2){
		$temp_eej="$gene-$d2-$i";
		if ($eff_href->{$length}{$temp_eej} >= $min_eff_complex){
		    $inc2C+=$reads{$sample}{$temp_eej}/$eff_href->{$length}{$temp_eej};
		    $Rinc2C+=$reads{$sample}{$temp_eej};
		}
	    }
	}
	### Exclusion reads (It does NOT take all EEJs around the alternative exon, but only those including C1 or C2.)
	if (!$use_all_excl_eej){
	    for $i (0..$d1-1){
		$temp_eej="$gene-$i-$a2";
		if ($eff_href->{$length}{$temp_eej} >= $min_eff_complex){
		    $exc1C+=$reads{$sample}{$temp_eej}/$eff_href->{$length}{$temp_eej};
		    $Rexc1C+=$reads{$sample}{$temp_eej};
		}
	    }
	    for $i ($d1+1..$last_donor{$gene}){
		if (($D_CO_href->{$gene}{$i} < $A_CO_href->{$gene}{$a1} && $strand eq "+") || ($D_CO_href->{$gene}{$i} > $A_CO_href->{$gene}{$a1} && $strand eq "-")){
		    $temp_eej="$gene-$i-$a2";
		    if ($eff_href->{$length}{$temp_eej} >= $min_eff_complex){
			$exc1C+=$reads{$sample}{$temp_eej}/$eff_href->{$length}{$temp_eej};
			$Rexc1C+=$reads{$sample}{$temp_eej};
		    }
		}
	    }
	    for $i ($a2+1..$last_acceptor{$gene}){
		$temp_eej="$gene-$d1-$i";
		if ($eff_href->{$length}{$temp_eej} >= $min_eff_complex){
		    $exc2C+=$reads{$sample}{$temp_eej}/$eff_href->{$length}{$temp_eej};
		    $Rexc2C+=$reads{$sample}{$temp_eej};
		}
	    }
	    for $i (0..$a2-1){
		if (($A_CO_href->{$gene}{$i} > $D_CO_href->{$gene}{$d2} && $strand eq "+") || ($A_CO_href->{$gene}{$i} < $D_CO_href->{$gene}{$d2} && $strand eq "-")){
		    $temp_eej="$gene-$d1-$i";
		    if ($eff_href->{$length}{$temp_eej} >= $min_eff_complex){
			$exc2C+=$reads{$sample}{$temp_eej}/$eff_href->{$length}{$temp_eej};
			$Rexc2C+=$reads{$sample}{$temp_eej};
		    }
		}
	    }
            $excC=$exc1C+$exc2C;
	    $RexcC=$Rexc1C+$Rexc2C;
	}
	### Taking all EEJs around the alternative exon 
	elsif ($use_all_excl_eej){ ### NOT USED, NOT TESTED, added as user's option
#	    for $i (0..$d2-1){
#		for $j ($a1+1..$last_acceptor{$gene}){
	    for $i ($d2-$extra_eej..$d2-1){
		for $j ($a1+1..$a1+$extra_eej){
		    if (($D_CO_href{$gene}{$i} < $acceptor_coord && $A_CO_href->{$gene}{$j} > $donor_coord && $str eq "+") || ($D_CO_href{$gene}{$i} > $acceptor_coord && $A_CO_href->{$gene}{$j} < $donor_coord && $str eq "-")){
			$temp_eej="$gene-$i-$j";
			if ($eff_href->{$length}{$temp_eej} >= $min_eff_complex){
			    $excC+=$reads{$sample}{$temp_eej}/$eff_href->{$length}{$temp_eej};
			    $RexcC+=$reads{$sample}{$temp_eej};
			}
		    }
		}
	    }
	}
	
# Sum of corrected reads
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
	unless($legacyFlag) {
	    my $totalN = $total_reads; # no point in re-assigning really.
	    my($pPSI, $exValOfInc, $exValOfExc) = (0, 0, 0);
	    unless($PSI_complex eq "NA" or $totalN == 0) {
		$pPSI = $PSI_complex / 100;
		$exValOfInc = sprintf("%.2f", $pPSI * $totalN);
		$exValOfExc = sprintf("%.2f", (1-$pPSI) * $totalN);
	    }
	    # ALTER QUAL OUTPUT HERE>>
	    $Q .= "\@$exValOfInc,$exValOfExc";
	} 
	
	print PSIs "\t$PSI_complex\t$Q" if $event;
	print COUNTs "\t$Rexc\t$Rinc1\t$Rinc2\t$RexcC\t$Rinc1C\t$Rinc2C\t$PSI_complex=$Q" if $event;
    }
    
    print PSIs "\n" if $event;
    print COUNTs "\n" if $event;
}
close PSIs;
close COUNTs;
