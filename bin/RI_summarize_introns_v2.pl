#!/usr/bin/env perl
# This script gets, for each event, the normalized read counts for the three associated junctions in the sample

# Version of RI_summarize_introns.pl to run IR on IR_version == 2
# By Manuel Irimia [09/Nov/2015]

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
           "readLen=i" => \$rle, "root=s" => \$root, "s" => \$strandaware);

my $mapcorr_fileswitch=""; if($strandaware){$mapcorr_fileswitch="-SS"}

my $maxcount = $rle - 15;
my $type = "ALL"; # ALL or new

# Getting mappability information (from output of uniquecount.IJ.pl)
my %ucount;

my $ucountFile = "$dbDir/FILES/$sp.Introns.sample${mapcorr_fileswitch}.200.$rle.uniquecount.txt";
my $UC = openFileHandle($ucountFile);
while(<$UC>){
    chomp($_);
    my($junction0,$count,$pos) = split(/\t/,$_);
    my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
    my $junction = $gene.":".$trans.":".$en;
    $ucount{$junction} = $count;
}
close $UC;

# Getting raw read counts for intron bodies
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
# summary_v2 has raw and corrected read counts, and should not be deleted
my $juncAnnotationFile = "./to_combine/$root.IR.summary_v2.txt"; # new input
my $ANOT = openFileHandle($juncAnnotationFile);
my $outfile = "./to_combine/$root.IR2"; 
open(OUT,">$outfile") or die "Failed to open $outfile: $!\n";
my $head = <$ANOT>;
print OUT "Event\tEIJ1\tEIJ2\tEEJ\tI\n"; # File to be used in RI_MakeTablePIR.R
while(<$ANOT>){
    chomp($_);
    my ($event,$C1A,$AC2,$C1C2) = split(/\t/,$_); # 7 elements; event and next three are the corrected counts
    my $I = 0;
    
    if(defined $rcount{$event} && defined $ucount{$event} && $ucount{$event} > 0){
	$I = sprintf("%.2f",($rcount{$event} / $ucount{$event}) * $maxcount);
    }    
    # added on 02/10/15
    $C1A = 0 if !$C1A;
    $AC2 = 0 if !$AC2;
    $C1C2 = 0 if !$C1C2;
    
    print OUT "$event\t$C1A\t$AC2\t$C1C2\t$I\n" if $C1A ne "ne" && $AC2 ne "ne" && $C1C2 ne "ne"; # to avoid ne's in subsequent R processing
}
close $ANOT;
close OUT;
