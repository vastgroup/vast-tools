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

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp,
           "readLen=i" => \$rle, "root=s" => \$root);

my $maxcount = $rle - 15;
my $type = "ALL"; # added to incorporate all COMBI-B EEJs.

# Getting mappability inexit bash scriptformation (from output of uniquecount.IJ.pl)
my %ucount;
my %transcripts;

# just to get the transcript per gene in COMBI
my $ucountFile = "$dbDir/IR/$sp/$sp.IntronJunctions.$type.$rle.8.uniquecount.txt"; 
open (EFF, $ucountFile);
while (<EFF>){
    chomp($_);
    my($junction0,$count) = split(/\t/,$_);
    
    if ($junction0 =~ /.+?\:.+?\:/){ # i.e. it's a EI junction
	my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
	$transcripts{$gene}=$trans; # first should go the EI junctions
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
	my($gene,$en,$co) = split(/\-/,$junction0);
	my $junction = $gene."-".$en;

	my($e1,$e2)=$en=~/(\d+)\_(\d+)/;
	$e1++; # they were 0-based
	$e2++;
	my $new_en="$e1-EE";
	my $final_junction = $gene.":".$transcripts{$gene}.":".$new_en; # requires all IEs going first
	die "$gene\n" if !$transcripts{$gene};
	
	$ucount{$final_junction} = $count if $e1 eq $e2; # stores the mappability for the simple EE
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
	my($gene,$en,$co) = split(/\-/,$junction0);
	my $junction = $gene."-".$en;
	my($e1,$e2)=$en=~/(\d+)\_(\d+)/;
	$e1++; # they were 0-based
	$e2++;

	for my $i ($e1..$e2){ # it counts as splice out for all exons in between
	    my $new_en="$i-EE"; # copies the 2-EE structures (EE for intron 2)
	    my $final_junction = $gene.":".$transcripts{$gene}.":".$new_en;

	    if (!defined ($rcount{$final_junction}) && defined ($ucount{$junction})){
		$rcount{$final_junction}=0;
	    }
	    if (!defined ($corrected_count{$final_junction}) && defined ($ucount{$junction})){
		$corrected_count{$final_junction} = 0;
	    }

	    $rcount{$final_junction} = $rcount{$final_junction} + $count if $ucount{$junction};
	    $corrected_count{$final_junction} = $corrected_count{$final_junction} + (($count / $ucount{$junction}) * $maxcount) if $ucount{$junction}; 
	    # $ucount{$junction} and not final_junction!
	}
    }
# Examples of junctions: 
    #geneID:transcriptID:3-EI1:6:42
    #geneID-0_1-872_443
}
close $RC;

# Getting junction annotation (i.e. the IDs of the 3 junctions associated with each event) and generating output file
my %eventseen;
my $juncAnnotationFile = "$dbDir/IR/$sp/$sp.IntronJunctions.new.annotation.txt"; # still maintained
my $ANOT = openFileHandle($juncAnnotationFile);
# summary_v2: contains corrected and raw read counts. It should not be deleted (used to build coverage key)
my $outfile = "$dbDir/SAM/$sp/$root.IR.summary_v2.txt";
open(OUT,">$outfile") or die "Failed to open $outfile: $!\n";
my $head = <$ANOT>;
print OUT "Event\tcEIJ1\tcEIJ2\tcEEJ\trEIJ1\trEIJ2\trEEJ\n";
while(<$ANOT>){
    chomp($_);
    my ($event,$C1Aj,$AC2j,$C1C2j,@rest) = split(/\t/,$_);

    if(! defined $eventseen{$event}){	
	my $C1A_cor="";
	my $C1A_raw="";
	if ($ucount{$C1Aj}>0){
	    if (!defined $corrected_count{$C1Aj}){
		$corrected_count{$C1Aj} = 0;
	    }
	    if (!defined $rcount{$C1Aj}){
		$rcount{$C1Aj} = 0;
	    }
	    $C1A_cor=$corrected_count{$C1Aj};
	    $C1A_raw=$rcount{$C1Aj};
	}
	else {
	    $C1A_cor="ne";
	    $C1A_raw="ne";
	}
	my $AC2_cor="";
	my $AC2_raw="";
	if ($ucount{$AC2j}>0){
	    if (!defined $corrected_count{$AC2j}){
		$corrected_count{$AC2j} = 0;
	    }
	    if (!defined $rcount{$AC2j}){
		$rcount{$AC2j} = 0;
	    }
	    $AC2_cor=$corrected_count{$AC2j};
	    $AC2_raw=$rcount{$AC2j};
	}
	else {
	    $AC2_cor="ne";
	    $AC2_raw="ne";
	}
	my $C1C2_cor="";
	my $C1C2_raw="";
	if ($ucount{$C1C2j}>0){
	    if (!defined $corrected_count{$C1C2j}){
		$corrected_count{$C1C2j} = 0;
	    }
	    if (!defined $rcount{$C1C2j}){
		$rcount{$C1C2j} = 0;
	    }
	    $C1C2_cor=$corrected_count{$C1C2j};
	    $C1C2_raw=$rcount{$C1C2j};
	}
	else {
	    $C1C2_cor="ne";
	    $C1C2_raw="ne";
	}

	print OUT "$event\t$C1A_cor\t$AC2_cor\t$C1C2_cor\t$C1A_raw\t$AC2_raw\t$C1C2_raw\n";
	$eventseen{$event} = "Y";
    }
}
close $ANOT;
close OUT;
