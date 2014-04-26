#!/usr/bin/perl
# Script to sort bowtie outputs by reads name (to then not double-count the same reads)
# It also discards part of the output, saving around half of the space. For more saving, use MakeSummarySAM.pl
# It is based on pure brute force, and can take a lot of memmory for large outputs (~20GB for 1 lane).

BEGIN {push @INC, '../lib'}
use FuncBasics qw(:all);

$file=$ARGV[0];
$IN = open ($ARGV[0]) || die "Can't open the input file\n";
($root)=$ARGV[0]=~/(.+)\./;
open (OUT,">$root"."_s.out") || die "Can't open the output file\n";

while (<$IN>){
    @t=split(/\t/);
    $data{$t[0]}="$t[1]\t$t[2]\t$t[3]\t$t[-1]"; # for the default bowtie1 output
}
close $IN;

foreach $read (sort (keys %data)){ # massive sorting of the keys (~read name)
    print OUT "$read\t$data{$read}";
}
close OUT;

system "rm $file"; # Removes the original file
