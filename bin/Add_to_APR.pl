#!/usr/bin/perl
# This script gets the PSIs for the APR events, plus qualities, etc.

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Getopt::Long;

use Cwd;
#$cwd = getcwd;
#($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

#$sp=$ARGV[0]; # DEPRECATED --TSW
#$type=$ARGV[1];

my $sp;
my $type;
my $dbDir;
my $samLen;

GetOptions("sp=s" => \$sp, "type=s" => \$type,
			  "dbDir=s" => \$dbDir, "len=i" => \$samLen);


#die "Needs Species (Hsa/Mmu) and type (exskX/MULTI3X)\n" if ($#ARGV<1);

$type_of_template="EXSK" if $type eq "exskX";
$type_of_template="MULTI" if $type eq "MULTI3X";

#@EXSK=glob("spli_out/$sp*$type");
my(@EXSK) = glob("spli_out/*$type");  # not sure what to do with this.  TEST PLZ --TSW

#die "@EXSK";

open (TEMPLATE, "$dbDir/TEMPLATES/$sp.$type_of_template.Template.2.txt") || die "Can't find $type_of_template template file for $sp\n";
$head=<TEMPLATE>;
chomp($head);
$head_reads=$head;
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $event=$t[1];
    $ALL{$event}=$_;
}
close TEMPLATE;

print "Loading and parsing data for each sample for $type_of_template\n";
foreach my $file (@EXSK){
	 my $fname = $file;
	 $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
    $head.="\t$sample\t$sample-Q";
    $head_reads.="\t$sample-Re\t$sample-Ri1\t$sample-Ri2\t$sample-e\t$sample-i1\t$sample-i2\t$sample-Q";

    open (I, $file);
    while (<I>){
	chomp;
	@t=split(/\t/);
	$event="$t[3]";
	$event="$t[3]=$t[17]-$t[18]" if $type eq "MULTI3X";
	
	$PSI{$event}{$sample}=$t[12];
	
	### raw reads # same for both types
	$Rexc{$event}{$sample}=$t[13];
	$Rinc1{$event}{$sample}=$t[14];
	$Rinc2{$event}{$sample}=$t[15];
	    
	if ($type eq "MULTI3X"){
	    ($cE,$Es,$cEs)=$t[19]=~/(.+?)\=(.+?)\=(.+)/; #corrected all (cE), raw reference (Es), corrected reference (cEs)
	    ($cI1,$I1s,$cI1s)=$t[20]=~/(.+?)\=(.+?)\=(.+)/;
	    ($cI2,$I2s,$cI2s)=$t[21]=~/(.+?)\=(.+?)\=(.+)/;

	    # All corrected reads	    
	    $exc{$event}{$sample}=$cE;
	    $inc1{$event}{$sample}=$cI1;
	    $inc2{$event}{$sample}=$cI2;
	    # Raw reads from reference
	    $RexcS{$event}{$sample}=$Es;
	    $Rinc1S{$event}{$sample}=$I1s;
	    $Rinc2S{$event}{$sample}=$I2s;
	    # Corrected reads from reference
	    $excS{$event}{$sample}=$cEs;
	    $inc1S{$event}{$sample}=$cI1s;
	    $inc2S{$event}{$sample}=$cI2s;
	}
	elsif ($type eq "exskX"){
	    # All corrected reads	
	    $exc{$event}{$sample}=$t[19];
	    $inc1{$event}{$sample}=$t[20];
	    $inc2{$event}{$sample}=$t[21];
	}
	$complexity{$event}{$sample}=$t[22]; # only in MULTI
    }
    close I;
}	

$NUM=$#EXSK+1;

open (PSIs, ">raw_incl/INCLUSION_LEVELS_$type_of_template-$sp$NUM-n.tab"); # change output directories --TSW
open (COUNTs, ">raw_reads/RAW_READS_$type_of_template-$sp$NUM-n.tab");

print PSIs "$head\n";
print COUNTs "$head_reads\n";

print STDERR "Parsing info and getting Quality scores (Q) for $type_of_template\n";
foreach $event (sort keys %ALL){
    print PSIs "$ALL{$event}";
    print COUNTs "$ALL{$event}";
    
    foreach $file (@EXSK){
	#($sample)=$file=~/$sp.+?\-\d+?\-(.+?)\./;
    my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;

	$PSI=sprintf("%.2f",$PSI{$event}{$sample});
	$total_raw_reads=$Rexc{$event}{$sample}+$Rinc1{$event}{$sample}+$Rinc2{$event}{$sample};
	$total_corr_reads=$exc{$event}{$sample}+$inc1{$event}{$sample}+$inc2{$event}{$sample};
	$total_ref_reads=$RexcS{$event}{$sample}+$Rinc1S{$event}{$sample}+$Rinc2S{$event}{$sample};

	if ($type eq "MULTI3X"){
#### Coverage scores. Score 1: Using all raw reads	
	    if (($Rexc{$event}{$sample}>=20 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=20) || ($Rinc1{$event}{$sample} >=20 && $Rinc2{$event}{$sample}>=15)) && $total_raw_reads >= 100){
		$Q="SOK";
	    }
	    elsif (($Rexc{$event}{$sample}>=20 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=20) || ($Rinc1{$event}{$sample} >=20 && $Rinc2{$event}{$sample}>=15)) && $total_raw_reads < 100){
		$Q="OK";
	    }
	    elsif ($Rexc{$event}{$sample}>=15 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=10) || ($Rinc1{$event}{$sample} >= 10 && $Rinc2{$event}{$sample}>=15)){	
		$Q="LOW";
	    }
	    elsif ($Rexc{$event}{$sample} >= 10 || ($Rinc1{$event}{$sample} >= 10 && $Rinc2{$event}{$sample} >= 5) || ($Rinc1{$event}{$sample} >= 5 && $Rinc2{$event}{$sample} >= 10)){	
		$Q="VLOW";
	    }
	    else {
		$Q="N";
	    }
#### Score 2: Using all corrected reads
	    if (($exc{$event}{$sample}>=20 || ($inc1{$event}{$sample} >=15 && $inc2{$event}{$sample}>=20) || ($inc1{$event}{$sample} >=20 && $inc2{$event}{$sample}>=15)) && $total_corr_reads>=100){
		$Q.=",SOK";
	    }
	    elsif (($exc{$event}{$sample}>=20 || ($inc1{$event}{$sample} >=15 && $inc2{$event}{$sample}>=20) || ($inc1{$event}{$sample} >=20 && $inc2{$event}{$sample}>=15)) && $total_corr_reads<100){
		$Q.=",OK";
	    }
	    elsif ($exc{$event}{$sample}>=15 || ($inc1{$event}{$sample} >=15 && $inc2{$event}{$sample}>=10) || ($inc1{$event}{$sample} >= 10 && $inc2{$event}{$sample}>=15)){	
		$Q.=",LOW";
	    }
	    elsif ($exc{$event}{$sample} >= 10 || ($inc1{$event}{$sample} >= 10 && $inc2{$event}{$sample} >= 5) || ($inc1{$event}{$sample} >= 5 && $inc2{$event}{$sample} >= 10)){	
		$Q.=",VLOW";
	    }
	    else {
		$Q.=",N";
	    }
#### Score 3: Using simple (=reference, C1A, AC2, C1C2) raw reads
	    if (($RexcS{$event}{$sample}>=20 || ($Rinc1S{$event}{$sample} >=15 && $Rinc2S{$event}{$sample}>=20) || ($Rinc1S{$event}{$sample} >=20 && $Rinc2S{$event}{$sample}>=15)) && $total_ref_reads>=100){
		$Q.=",SOK";
	    }
	    elsif (($RexcS{$event}{$sample}>=20 || ($Rinc1S{$event}{$sample} >=15 && $Rinc2S{$event}{$sample}>=20) || ($Rinc1S{$event}{$sample} >=20 && $Rinc2S{$event}{$sample}>=15)) && $total_ref_reads<100){
		$Q.=",OK";
	    }
	    elsif ($RexcS{$event}{$sample}>=15 || ($Rinc1S{$event}{$sample} >=15 && $Rinc2S{$event}{$sample}>=10) || ($Rinc1S{$event}{$sample} >= 10 && $Rinc2S{$event}{$sample}>=15)){	
		$Q.=",LOW";
	    }
	    elsif ($RexcS{$event}{$sample} >= 10 || ($Rinc1S{$event}{$sample} >= 10 && $Rinc2S{$event}{$sample} >= 5) || ($Rinc1S{$event}{$sample} >= 5 && $Rinc2S{$event}{$sample} >= 10)){	
		$Q.=",VLOW";
	    }
	    else {
		$Q.=",N";
	    }
	}
	
	if ($type eq "exskX"){
#### Coverage scores. Score 1: Using all raw reads	
	    if (($Rexc{$event}{$sample}>=20 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=20) || ($Rinc1{$event}{$sample} >=20 && $Rinc2{$event}{$sample}>=15)) && $total_raw_reads>=100){
		$Q="SOK";
	    }
	    elsif (($Rexc{$event}{$sample}>=20 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=20) || ($Rinc1{$event}{$sample} >=20 && $Rinc2{$event}{$sample}>=15)) && $total_raw_reads<100){
		$Q="OK";
	    }
	    elsif ($Rexc{$event}{$sample}>=15 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=10) || ($Rinc1{$event}{$sample} >= 10 && $Rinc2{$event}{$sample}>=15)){	
		$Q="LOW";
	    }
	    elsif ($Rexc{$event}{$sample} >= 10 || ($Rinc1{$event}{$sample} >= 10 && $Rinc2{$event}{$sample} >= 5) || ($Rinc1{$event}{$sample} >= 5 && $Rinc2{$event}{$sample} >= 10)){	
		$Q="VLOW";
	    }
	    else {
		$Q="N";
	    }
#### Score 2: Using all corrected reads
	    if (($exc{$event}{$sample}>=20 || ($inc1{$event}{$sample} >=15 && $inc2{$event}{$sample}>=20) || ($inc1{$event}{$sample} >=20 && $inc2{$event}{$sample}>=15)) && $total_corr_reads>=100){
		$Q.=",SOK";
	    }
	    elsif (($exc{$event}{$sample}>=20 || ($inc1{$event}{$sample} >=15 && $inc2{$event}{$sample}>=20) || ($inc1{$event}{$sample} >=20 && $inc2{$event}{$sample}>=15)) && $total_corr_reads<100){
		$Q.=",OK";
	    }
	    elsif ($exc{$event}{$sample}>=15 || ($inc1{$event}{$sample} >=15 && $inc2{$event}{$sample}>=10) || ($inc1{$event}{$sample} >= 10 && $inc2{$event}{$sample}>=15)){	
		$Q.=",LOW";
	    }
	    elsif ($exc{$event}{$sample} >= 10 || ($inc1{$event}{$sample} >= 10 && $inc2{$event}{$sample} >= 5) || ($inc1{$event}{$sample} >= 5 && $inc2{$event}{$sample} >= 10)){	
		$Q.=",VLOW";
	    }
	    else {
		$Q.=",N";
	    }
#### Score 3: Using simple (=reference, C1A, AC2, C1C2) raw reads	
	    if (($Rexc{$event}{$sample}>=20 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=20) || ($Rinc1{$event}{$sample} >=20 && $Rinc2{$event}{$sample}>=15)) && $total_raw_reads>=100){
		$Q.=",SOK";
	    }
	    elsif (($Rexc{$event}{$sample}>=20 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=20) || ($Rinc1{$event}{$sample} >=20 && $Rinc2{$event}{$sample}>=15)) && $total_raw_reads<100){
		$Q.=",OK";
	    }
	    elsif ($Rexc{$event}{$sample}>=15 || ($Rinc1{$event}{$sample} >=15 && $Rinc2{$event}{$sample}>=10) || ($Rinc1{$event}{$sample} >= 10 && $Rinc2{$event}{$sample}>=15)){	
		$Q.=",LOW";
	    }
	    elsif ($Rexc{$event}{$sample} >= 10 || ($Rinc1{$event}{$sample} >= 10 && $Rinc2{$event}{$sample} >= 5) || ($Rinc1{$event}{$sample} >= 5 && $Rinc2{$event}{$sample} >= 10)){	
		$Q.=",VLOW";
	    }
	    else {
		$Q.=",N";
	    }
	}

### Score 4: Calculate imbalance between inclusion EEJs (OK<B1<B2; Bl=not enough inclusion reads):
	$inc1=$inc1S{$event}{$sample} if $type eq "MULTI3X"; # only corrected reads from reference
	$inc2=$inc2S{$event}{$sample} if $type eq "MULTI3X"; # only corrected reads from reference
	$inc1=$inc1{$event}{$sample} if $type eq "exskX";
	$inc2=$inc2{$event}{$sample} if $type eq "exskX";

	if ($inc1 && $inc2){
	    if (($inc1/$inc2 > 2 && $inc1/$inc2 <=5) || ($inc2/$inc1> 2 && $inc2/$inc1 <=5)){
		$Q.=",B1";
	    }
	    elsif (($inc1/$inc2 > 5) || ($inc2/$inc1>5)){
		$Q.=",B2";
	    }
	    else {
		$Q.=",OK";
	    }
	}
	else {
	    if (!$inc1){
		$Q.=",B2" if $inc2>=20;
		$Q.=",Bl" if ($inc2<20 && $inc2>=15);
		$Q.=",Bn" if $inc2<15 && $inc2>0;
	    }
	    if (!$inc2){
		$Q.=",B2" if $inc1>=20;
		$Q.=",Bl" if ($inc1<20 && $inc1>=15);
		$Q.=",Bn" if $inc1<15;
	    }
	}
#### Score 5: Complexity score (S<C1<C2<C3)	
	$Q.=",S" if $type eq "exskX";
	$Q.=",$complexity{$event}{$sample}" if $type eq "MULTI3X";
	
	# Print out data
	print PSIs "\t$PSI\t$Q";
	print COUNTs "\t$Rexc{$event}{$sample}\t$Rinc1{$event}{$sample}\t$Rinc2{$event}{$sample}\t\t\t\t$PSI=$Q";
    }
    close I;
    print PSIs "\n";
    print COUNTs "\n";
}
close PSIs;
close COUNTs;

