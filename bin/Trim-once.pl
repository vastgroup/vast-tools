#!/usr/bin/perl
# This scripts trims reads into one X-nt read. 
# Discards shorter reads, if any

BEGIN {push @INC, '../lib'}
use FuncBasics qw(:all);

($root)=$ARGV[0]=~/(\S+?)\.f/; #fixed regex -TSW
$length=$ARGV[1];
die "You need to provide length as ARGV[1]\n" if !$ARGV[1];

$file=$ARGV[0];

### Obtains the begining of the read to set \$/
my $TMP = openFileHandle($ARGV[0]);
$head=<$TMP>;
close $TMP;
($/)=$head=~/(\@.{3})/;
$del=$/;

### Parses the original reads
#open (STDOUT, ">$root-$length.fq");
my $INPUT = openFileHandle($ARGV[0]);
<$INPUT>; #invalid line
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
        ($R1)=$rest=~/^(.{$length})/;
        print STDOUT "$del$name-1\n$S1\n$name2-1\n$R1\n";
        $total_reads_accepted++;
    }
}
close $INPUT;

print STDERR "[vast align trim]: Total processed reads: $total_reads\n";
print STDERR "[vast align trim]: Total valid reads: $total_reads_accepted\n";

if($total_reads <= 1 or $total_reads_accepted <= 1) { exit 1; }

