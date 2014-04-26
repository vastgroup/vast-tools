#!/usr/bin/perl
# This scripts trims reads into one X-nt read. 
# Discards shorter reads, if any

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
open (INPUT, $ARGV[0]);
<INPUT>; #invalid line
while (<INPUT>){
    /\n(.+?)\n(.+?)\n(.+?)\n/;
    $name=$`;
    $seq=$1;
    $name2=$2;
    $rest=$3;
    $total_reads++;

    $rest=~s/$del$//; # removes the heading bit at the end
    
    if (length($seq)>=$length && length($rest)>=$length) { #make sure the length of the sequence AND quality scores are >= length
        ($S1)=$seq=~/^(.{$length})/;
        ($R1)=$rest=~/^(.{$length})/;
        print OUTPUT "$del$name-1\n$S1\n$name2-1\n$R1\n";
        $total_reads_accepted++;
    }
}

print "Total processed reads: $total_reads\n";
print "Total valid reads: $total_reads_accepted\n";

