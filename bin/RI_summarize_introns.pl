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

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp,
           "readLen=i" => \$rle, "root=s" => \$root);

my $maxcount = $rle - 15;

# Getting mappability information (from output of uniquecount.IJ.pl)
my %ucount;

my $ucountFile = "$dbDir/FILES/$sp.Introns.sample.200.$rle.uniquecount.txt";
my $UC = openFileHandle($ucountFile);
while(<$UC>){
    chomp($_);
    my($junction0,$count,$pos) = split(/\t/,$_);
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
    my($junction0,$count) = split(/\t/,$_);
    my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
    my $junction = $gene.":".$trans.":".$en;
    $rcount{$junction} = $count;
}
close $RC;

# Getting junction annotation (i.e. the IDs of the 3 junctions associated with each event) and generating output file
my %eventseen;
my $juncAnnotationFile = "./to_combine/$root.IR.summary.txt";
my $ANOT = openFileHandle($juncAnnotationFile);
my $outfile = "./to_combine/$root.IR";
open(OUT,">$outfile") or die "Failed to open $outfile: $!\n";
my $head = <$ANOT>;
print OUT "Event\tEIJ1\tEIJ2\tEEJ\tI\n";
while(<$ANOT>){
    chomp($_);
    my ($event,$C1A,$AC2,$C1C2) = split(/\t/,$_);
    my $I = 0;

    if(defined $rcount{$event} && defined $ucount{$event} && $ucount{$event} > 0){
	    $I = $rcount{$event} / $ucount{$event} * $maxcount;
    }
    print OUT "$event\t$C1A\t$AC2\t$C1C2\t$I\n";
}
close $ANOT;
close OUT;
