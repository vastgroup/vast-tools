#!/usr/bin/env perl
# This script gets, for each event, the normalized read counts for the three associated junctions in the sample

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

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp,
           "readLen=i" => \$rle, "root=s" => \$root, "s" => \$strandaware);

my $mapcorr_fileswitch=""; if($strandaware){$mapcorr_fileswitch="-SS"}
my $maxcount = $rle - 15;

# Getting mappability inexit bash scriptformation (from output of uniquecount.IJ.pl)
my %ucount;

my $ucountFile = "$dbDir/FILES/$sp.IntronJunctions.new${mapcorr_fileswitch}.$rle.8.uniquecount.txt";
my $UC = openFileHandle($ucountFile);
while(<$UC>){
    chomp($_);
    my($junction0,$count) = split(/\t/,$_);
    my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
    my $junction = $gene.":".$trans.":".$en;
    $ucount{$junction} = $count;
}
close $UC;

# Getting raw read counts
my %rcount;
my $RC = openFileHandle($ARGV[0]);
while(<$RC>){
    chomp($_);
    my($junction0,$count,$pos) = split(/\t/,$_);
    my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
    my $junction = $gene.":".$trans.":".$en;
    $rcount{$junction} = $count;
}
close $RC;

# Getting junction annotation (i.e. the IDs of the 3 junctions associated with each event) and generating output file
my %eventseen;
my $juncAnnotationFile = "$dbDir/FILES/$sp.IntronJunctions.new.annotation.txt";
my $ANOT = openFileHandle($juncAnnotationFile);
my $outfile = "./to_combine/$root.IR.summary.txt";
open(OUT,">$outfile") or die "Failed to open $outfile: $!\n";
my $head = <$ANOT>;
print OUT "Event\tEIJ1\tEIJ2\tEEJ\n";
while(<$ANOT>){
    chomp($_);
    my ($event,$C1Aj,$AC2j,$C1C2j,@rest) = split(/\t/,$_);
    if(! defined $eventseen{$event} && 
        (defined $rcount{$C1Aj} || defined $rcount{$AC2j} || 
            defined $rcount{$C1C2j})){
	my $C1A = 0;
	if(defined $rcount{$C1Aj} && defined $ucount{$C1Aj} && $ucount{$C1Aj} > 0){
	    $C1A = $rcount{$C1Aj} / $ucount{$C1Aj} * $maxcount;
	}
	my $AC2 = 0;
	if(defined $rcount{$AC2j} && defined $ucount{$AC2j} && $ucount{$AC2j} > 0){
	    $AC2 = $rcount{$AC2j} / $ucount{$AC2j} * $maxcount;
	}
	my $C1C2 = 0;
	if(defined $rcount{$C1C2j} && defined $ucount{$C1C2j} && $ucount{$C1C2j} > 0){
	    $C1C2 = $rcount{$C1C2j} / $ucount{$C1C2j} * $maxcount;
	}
	print OUT "$event\t$C1A\t$AC2\t$C1C2\n";
	$eventseen{$event} = "Y";
    }
}
close $ANOT;
close OUT;
