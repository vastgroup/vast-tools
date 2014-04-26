#!/usr/bin/perl
# This scripts splits reads into two X-nt reads. 
# If the length of the original reads is < 2X, it splits in an overlapping manner.

BEGIN {push @INC, '../lib'}
use FuncBasics qw(:all);

($root)=$ARGV[0]=~/(.+?)\.f/;
$length=$ARGV[1];
die "You need to provide length as ARGV[1]\n" if !$ARGV[1];

$file=$ARGV[0];

### Obtains the begining of the read to set \$/
open (TEMP,$file);
$head=<TEMP>;
close TEMP;
($/)=$head=~/(\@.{3})/;
$del=$/;

### Parses the original reads
open (OUTPUT, ">$root-$length.fq");
$INPUT = openFileHandle ($ARGV[0]);
<$INPUT>; #invalid bit
while (<$INPUT>){
    /\n(.+?)\n(.+?)\n(.+?)\n/;
    $name=$`;
    $seq=$1;
    $name2=$2;
    $rest=$3;
    $total_reads++;

    $rest=~s/$del$//; # removes the heading bit at the end
    
    if (length($seq)>=$length && length($rest)>=$length) { #make sure the length of the sequence AND quality scores are >= length
        ($S1)=$seq=~/^(.{$length})/;
	($S2)=$seq=~/(.{$length})$/;
        ($R1)=$rest=~/^(.{$length})/;
        ($R2)=$rest=~/(.{$length})$/;
        print OUTPUT "$del$name-1\n$S1\n$name2-1\n$R1\n";
        print OUTPUT "$del$name-2\n$S2\n$name2-2\n$R2\n";
        $total_reads_accepted++;
    }
}
close $INPUT;

print "Total processed reads: $total_reads\n";
print "Total valid reads: $total_reads_accepted\n";

