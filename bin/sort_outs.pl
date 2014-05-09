#!/usr/bin/perl
#script to sort bowtie outputs by reads name (to then not double-count the same reads)
#it also discards part of the output, saving around half of the space. For more saving, use MakeSummarySAM.pl

#use English; # supposedly slow on old perl
use English qw( -no_match_vars );

$f=$ARGV[0];
open (I,$ARGV[0]) || die "Can't open the input file\n";
($root)=$ARGV[0]=~/(.+)\./;
open (O,">$root"."_s.out") || die "Can't open the output file\n";

while (<I>){
    @t=split(/\t/);
    $h{$t[0]}="$t[1]\t$t[2]\t$t[3]\t$t[-1]";
}
close I;

foreach $r (sort (keys %h)){
    print O "$r\t$h{$r}";
}
close O;

system "rm $f";
