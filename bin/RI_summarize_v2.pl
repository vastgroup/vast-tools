#!/usr/bin/env perl
# This script gets, for each event, the normalized read counts for the three associated junctions in the sample

# Modification of RI_summary.pl to run IR_version == 2
# By: Manuel Irimia [09/Nov/2015]
### IT CANNOT BE USED TO RUN v1

use warnings; 
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Getopt::Long;

my $dbDir;
my $sp;
my $rle;
my $root;
my $strandaware=0;
my $silent=0; # not used

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "ec" => \$silent,
           "readLen=i" => \$rle, "root=s" => \$root,  "s" => \$strandaware);

my $mapcorr_fileswitch=""; if($strandaware){$mapcorr_fileswitch="-SS"}
my $maxcount = $rle - 15;
my $type = "ALL"; # added to incorporate all COMBI-B EEJs.

# Getting mappability inexit bash scriptformation (from output of uniquecount.IJ.pl)
my %ucount;
my %transcripts;

# just to get the transcript per gene in COMBI
my $ucountFile = "$dbDir/FILES/$sp.IntronJunctions.$type${mapcorr_fileswitch}.$rle.8.uniquecount.txt"; 
open (EFF, $ucountFile);
while (<EFF>){
    chomp($_);
    my($junction0,$count) = split(/\t/,$_);
    
    if ($junction0 =~ /.+?\:.+?\:/){ # i.e. it's a EI junction
	my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
        my ($N)=$en=~/(\d+)\-/; # adapted to the special case of Hsa/Mmu [11/11/15]
        $transcripts{$gene}{$N}=$trans; # first should go the EI junctions # CHANGED 11/11/15 
    }
}
close EFF;

my $UC = openFileHandle($ucountFile);
while(<$UC>){
    chomp($_);
    my($junction0,$count) = split(/\t/,$_);

    if ($junction0 =~ /.+?\:.+?\:/){ # i.e. it's a EI junction
	my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
	my $junction = $gene.":".$trans.":".$en;
	$ucount{$junction} = $count; 
    }
    else {
#	my($gene,$en,$co) = split(/\-/,$junction0);
	my($gene,$en,$co)=$junction0=~/(.+)\-(.+?)\-(.+)/; #CHANGED 11/11/15 (genes with "-", e.g. D17H6S56E-3)
	my $junction = $gene."-".$en;

	my($e1,$e2)=$en=~/(\d+)\_(\d+)/;
	$e1++; # they were 0-based
	$e2++;
	my $new_en="$e1-EE";
	if (defined $transcripts{$gene}{$e1}){ #CHANGED 11/11/15 (some Hsa/Mmu "introns" do not exist if multitr)
	    my $final_junction = $gene.":".$transcripts{$gene}{$e1}.":".$new_en; # requires all IEs going first #CHANGED 11/11/15
	    $ucount{$final_junction} = $count if $e1 eq $e2; # stores the mappability for the simple EE
	}
	$ucount{$junction} = $count; 
    }
}
close $UC;

# Getting raw and corrected read counts
my %rcount;
my %corrected_count;
my $RC = openFileHandle($ARGV[0]);
while(<$RC>){
    chomp($_);
    my($junction0,$count,$pos) = split(/\t/,$_);
 
    if ($junction0 =~ /.+?\:.+?\:/){ # i.e. it's a EI junction
	my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
	my $junction = $gene.":".$trans.":".$en;
	$rcount{$junction} = $count if $ucount{$junction}; # changed 07/09/15
	$corrected_count{$junction} = ($count / $ucount{$junction}) * $maxcount if $ucount{$junction};
    }
    else {
#	my($gene,$en,$co) = split(/\-/,$junction0);
	my($gene,$en,$co)=$junction0=~/(.+)\-(.+?)\-(.+)/; #CHANGED 11/11/15 (genes with "-", e.g. D17H6S56E-3)
	my $junction = $gene."-".$en;
	my($e1,$e2)=$en=~/(\d+)\_(\d+)/;
	$e1++; # they were 0-based
	$e2++;

	for my $i ($e1..$e2){ # it counts as splice out for all exons in between
	    my $new_en="$i-EE"; # copies the 2-EE structures (EE for intron 2)

	    if (defined $transcripts{$gene}{$e1}){ #CHANGED 11/11/15 (some Hsa/Mmu "introns" do not exist if multitr)
		my $final_junction = $gene.":".$transcripts{$gene}{$e1}.":".$new_en;  #CHANGED 11/11/15
		if (!defined ($rcount{$final_junction}) && defined ($ucount{$junction})){
		    $rcount{$final_junction}=0;
		}
		if (!defined ($corrected_count{$final_junction}) && defined ($ucount{$junction})){
		    $corrected_count{$final_junction} = 0;
		}
		
		$rcount{$final_junction} = $rcount{$final_junction} + $count if $ucount{$junction};
		$corrected_count{$final_junction} = $corrected_count{$final_junction} + (($count / $ucount{$junction}) * $maxcount) if $ucount{$junction}; 
	    }
	}
    }
# Examples of junctions: 
    #geneID:transcriptID:3-EI1:6:42
    #geneID-0_1-872_443
}
close $RC;

# Getting junction annotation (i.e. the IDs of the 3 junctions associated with each event) and generating output file
my %eventseen;
my $juncAnnotationFile = "$dbDir/FILES/$sp.IntronJunctions.new.annotation.txt"; # still maintained
my $ANOT = openFileHandle($juncAnnotationFile);
# summary_v2: contains corrected and raw read counts. It should not be deleted (used to build coverage key)
my $outfile = "./to_combine/$root.IR.summary_v2.txt";
open(OUT,">$outfile") or die "Failed to open $outfile: $!\n";
my $head = <$ANOT>;
print OUT "Event\tcEIJ1\tcEIJ2\tcEEJ\trEIJ1\trEIJ2\trEEJ\n";
while(<$ANOT>){
    chomp($_);
    my ($event,$EI1j,$EI2j,$EEj,@rest) = split(/\t/,$_); #names modified from v1 for IR consitency [11/11/15]

    if(! defined $eventseen{$event}){	
	my $EI1_cor="";
	my $EI1_raw="";
	$ucount{$EI1j}=0 if (!defined $ucount{$EI1j}); #CHANGED 11/11/15 (for missing junctions in conversion)  
	if ($ucount{$EI1j}>0){
	    if (!defined $corrected_count{$EI1j}){
		$corrected_count{$EI1j} = 0;
	    }
	    if (!defined $rcount{$EI1j}){
		$rcount{$EI1j} = 0;
	    }
	    $EI1_cor=$corrected_count{$EI1j};
	    $EI1_raw=$rcount{$EI1j};
	}
	else {
	    $EI1_cor="ne";
	    $EI1_raw="ne";
	}

	my $EI2_cor="";
	my $EI2_raw="";
	$ucount{$EI2j}=0 if (!defined $ucount{$EI2j}); #CHANGED 11/11/15 (for missing junctions in conversion) 
	if ($ucount{$EI2j}>0){
	    if (!defined $corrected_count{$EI2j}){
		$corrected_count{$EI2j} = 0;
	    }
	    if (!defined $rcount{$EI2j}){
		$rcount{$EI2j} = 0;
	    }
	    $EI2_cor=$corrected_count{$EI2j};
	    $EI2_raw=$rcount{$EI2j};
	}
	else {
	    $EI2_cor="ne";
	    $EI2_raw="ne";
	}

	my $EE_cor="";
	my $EE_raw="";
	$ucount{$EEj}=0 if (!defined $ucount{$EEj}); #CHANGED 11/11/15 (for missing junctions in conversion) 
	if ($ucount{$EEj}>0){
	    if (!defined $corrected_count{$EEj}){
		$corrected_count{$EEj} = 0;
	    }
	    if (!defined $rcount{$EEj}){
		$rcount{$EEj} = 0;
	    }
	    $EE_cor=$corrected_count{$EEj};
	    $EE_raw=$rcount{$EEj};
	}
	else {
	    $EE_cor="ne";
	    $EE_raw="ne";
	}

	print OUT "$event\t$EI1_cor\t$EI2_cor\t$EE_cor\t$EI1_raw\t$EI2_raw\t$EE_raw\n";
	$eventseen{$event} = "Y";
    }
}
close $ANOT;
close OUT;
