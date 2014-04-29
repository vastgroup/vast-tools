#!/usr/bin/perl
# This scripts splits reads into two X-nt reads. 
# If the length of the original reads is < 2X, it splits in an overlapping manner.

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

my ($root)=$ARGV[0]=~/(.+?)\.f/;
my $length=$ARGV[1];
die "You need to provide length as ARGV[1]\n" if !$ARGV[1];

my $stdPrefix = '[vast align trim]:';

my $file=$ARGV[0];

### Obtains the begining of the read to set \$/
print STDERR "$stdPrefix Checking first read... \n";
print STDERR "$stdPrefix Ignore gzip broken pipe warning, if present\n";
my $TMP = openFileHandle($file);
my $head=<$TMP>;
close $TMP;
($/)=$head=~/(\@.{3})/;
my $del=$/;

### Initialize variables
my $total_reads = 0;
my $total_reads_accepted = 0;
my $name;
my $name2;
my $seq;
my $rest;
my ($S1, $S2, $R1, $R2);

### Parses the original reads
my $INPUT = openFileHandle ($ARGV[0]);
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
        print STDOUT "$del$name-1\n$S1\n$name2-1\n$R1\n";
        print STDOUT "$del$name-2\n$S2\n$name2-2\n$R2\n";
        $total_reads_accepted++;
    }
}
close $INPUT;

print STDERR "$stdPrefix Total processed reads: $total_reads\n";
print STDERR "$stdPrefix Total valid reads: $total_reads_accepted\n";

if($total_reads <= 1 or $total_reads_accepted <= 1) { exit 1; }
