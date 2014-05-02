#!/usr/bin/perl 
# This script takes an a posteriori template and uses it to get PSIs for exactly those events. 
# It does NOT do a new call for AS events

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

#use Cwd;
#$cwd = getcwd;
#($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

use Getopt::Long;

my $dbDir;
my $sp;
my $samLen;
my $verboseFlag;
my $legacyFlag;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "len=i" => \$samLen,
			  "verbose=i" => \$verboseFlag, "legacy" => \$legacyFlag);

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine combi]: $verbMsg\n";
  }
}

#$sp=$ARGV[0];
die "Needs Species key\n" if !defined($sp);

$COMB="M"; # Only available version

@EEJ=glob("spli_out/*.ee*"); # is this right? --TSW
@EFF=glob("$dbDir/FILES/$sp*-$COMB-*-gDNA.ef*");
die "[vast combine combi error] Needs effective from database!\n" if !@EFF;

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
verbPrint "Loading EEJ read counts data\n";
foreach $file (@EEJ){
	  my $fname = $file;
     $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
#    ($sample)=$file=~/COMBI\-[A-Z]\-\d+?\-(.+)\./;
    $head_PSIs.="\t$sample\t$sample-Q";
    $head_ReadCounts.="\t$sample-Re\t$sample-Ri1\t$sample-Ri2\t$sample-ReC\t$sample-Ri1C\t$sample-Ri2C\t$sample-Q";

    open (EEJ, $file);
    while (<EEJ>){
		chomp;
		@t=split(/\t/);
		$gene=$t[0];
		$eej=$t[1];
		$event="$gene-$eej";

        $reads{$sample}{$event}=$t[2];
        ($donor,$acceptor)=$eej=~/(\d+?)\-(\d+)/;
        $last_acceptor{$gene}=$acceptor if $last_acceptor{$gene}<$acceptor; # keeps track of the last used acceptor
    }
    close EEJ;
}

### Setting output files
$NUM=$#EEJ+1; # number of samples
open (PSIs, ">raw_incl/INCLUSION_LEVELS_COMBI-$sp$NUM-n.tab");
open (COUNTs, ">raw_reads/RAW_READS_COMBI-$sp$NUM-n.tab");
print PSIs "$head_PSIs\n";
print COUNTs "$head_ReadCounts\n";

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
#	($length,$sample)=$file=~/COMBI\-[A-Z]\-(\d+?)\-(.+)\./;
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
	    if ((($D_CO{$gene}{$i} < $acceptor_coord && $strand eq "+") || ($D_CO{$gene}{$i} > $acceptor_coord && $strand eq "-")) && $i != $d1){
		$temp_eej="$gene-$i-$a1";
		$inc1C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		$Rinc1C+=$reads{$sample}{$temp_eej} if $eff{$length}{$temp_eej}>0;
	    }
	}
	for $i ($a1+1..$last_acceptor{$gene}){ # loops through all acceptors in the gene from a1+1 to the last
	    if ((($A_CO{$gene}{$i} > $donor_coord && $strand eq "+") || ($A_CO{$gene}{$i} < $donor_coord && $strand eq "-")) && $i != $a2){
		$temp_eej="$gene-$d2-$i";
		$inc2C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		$Rinc2C+=$reads{$sample}{$temp_eej} if $eff{$length}{$temp_eej}>0;
	    }
	}
	### Exclusion reads (It does NOT take all EEJs around the alternative exon, but only those including C1 or C2.)
	if (!$ALL_EXC_EEJ){
	    for $i (0..$d1-1){
		$temp_eej="$gene-$i-$a2";
		$exc1C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		$Rexc1C+=$reads{$sample}{$temp_eej} if $eff{$length}{$temp_eej}>0;
	    }
	    for $i ($d1+1..$last_acceptor{$gene}){
		if (($D_CO{$gene}{$i} < $A_CO{$gene}{$a1} && $strand eq "+") || ($D_CO{$gene}{$i} > $A_CO{$gene}{$a1} && $strand eq "-")){
		    $temp_eej="$gene-$i-$a2";
		    $exc1C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		    $Rexc1C+=$reads{$sample}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		}
	    }
	    for $i ($a2+1..$last_acceptor{$gene}){
		$temp_eej="$gene-$d1-$i";
		$exc2C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		$Rexc2C+=$reads{$sample}{$temp_eej} if $eff{$length}{$temp_eej}>0;
	    }
	    for $i (0..$a2-1){
		if (($A_CO{$gene}{$i} > $D_CO{$gene}{$d2} && $strand eq "+") || ($A_CO{$gene}{$i} < $D_CO{$gene}{$d2} && $strand eq "-")){
		    $temp_eej="$gene-$d1-$i";
		    $exc2C+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		    $Rexc2C+=$reads{$sample}{$temp_eej} if $eff{$length}{$temp_eej}>0;
		}
	    }
            $excC=$exc1C+$exc2C;
	    $RexcC=$Rexc1C+$Rexc2C;
	}
	### Taking all EEJs around the alternative exon 
	elsif ($ALL_EXC_EEJ){ ### NOT USED, NOT TESTED
	    for $i (0..$d2-1){
		for $j ($a1+1..$last_acceptor{$gene}){
		    if (($D_CO{$gene}{$i} < $acceptor_coord && $A_CO{$gene}{$j} > $donor_coord && $str eq "+") || ($D_CO{$gene}{$i} > $acceptor_coord && $A_CO{$gene}{$j} < $donor_coord && $str eq "-")){
			$temp_eej="$gene-$i-$j";
			if ($eff{$length}{$temp_eej}>0){
			    $excC+=$reads{$sample}{$temp_eej}/$eff{$length}{$temp_eej};
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
     unless($PSI_complex eq "NA" or $totalN < 2) {
       $pPSI = $PSI_complex / 100;
       $exValOfInc = $pPSI * $totalN;
       $exValOfExc = (1-$pPSI) * $totalN;
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
