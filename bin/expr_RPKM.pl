#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Cwd;
#$cwd = getcwd;
#($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

use Getopt::Long;

my $sp;
my $dbDir;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp);

if ($ARGV[0] eq "--help" || !$ARGV[0]){
    print "expr_RPKM.pl mRNA.out mRNA.eff\n\n";
    die "Script to calculate cRPKMs from read counts per gene, provided a file with mappability\n";
}

### Loads effective lengths (i.e. mappable positions per transcript).
open (MAPPABILITY, $ARGV[1]) || die "Needs file with effective lengths (mappability file) as \$ARGV[1]\n";
while (<MAPPABILITY>){
    chomp;
    @t=split(/\t/);
    $effective_length{$t[0]}=$t[1]; # gets the mappability for each gene
}
close MAPPABILITY;

### Read count
$file=$ARGV[0]; #Output from bowtie
system "gunzip $file" if $file=~/\.gz/;
$file=~s/\.gz//;
open (INPUT, $file) || die "Can't find the input file\n"; 

while (<INPUT>){ #analyzes the bowtie output
    chomp;
    @d=split(/\t/);
    $gene=$d[2];
    $tally{$gene}++;
    $total_mapped++;
}
close INPUT;

### Calculating cRPKMs and finalizing
($root)=$file=~/(.+?)\.out/;
open (OUTPUT, ">$root"."_exprRPKM.txt"); #output file (cRPKM per gene)

foreach $gene (sort (keys %effective_length)){
    $cRPKM="ne";
    $cRPKM=sprintf("%.2f",1000000*(1000*$tally{$gene}/$effective_length{$gene})/$total_mapped) if $effective_length{$gene};
    print OUTPUT "$gene\t$cRPKM\t$tally{$gene}\n" if $gene;
}
close OUTPUT;

system "gzip $file";
