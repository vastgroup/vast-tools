#!/usr/bin/env perl
# This script puts all cRPKMs together in two tables (one alone, the other one with the raw read counts for each sample).

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Cwd;
use Getopt::Long;

my $dbDir;
my $sp;
my $cRPKMCounts = 0;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "C" => \$cRPKMCounts);

die "[vast combine cRPKM error] Needs Species\n" if !$sp;
my @files=glob("expr_out/*.cRPKM");
my $index=$#files+1;

# creates output files where run
my $OUTPUT;
open ($OUTPUT, ">cRPKM_AND_COUNTS-$sp$index.tab") if $cRPKMCounts;
open (RPKM, ">cRPKM-$sp$index.tab");

my %names;
open (NAMES, "$dbDir/FILES/$sp.ID.names.txt");
while (<NAMES>){
    chomp;
    my @t=split(/\t/);
    $names{$t[0]}=$t[1];
}
close NAMES;

my %data;
my %RPKM;
my $head="ID\tNAME";
my $headRPKM=$head;
foreach my $f (@files){
    my ($root)=$f=~/([^\/]+).cRPKM/;
    $head.="\t$root-cRPKM\t$root-Counts";
    $headRPKM.="\t$root";
    open (INPUT, $f);
    while (<INPUT>){
        chomp;
        my @t=split(/\t/);
        my $cRPKM = sprintf("%.2f", 0);
        my $raw_count = 0;

        if ($t[1] eq 'NA') {
            $cRPKM = 'NA';
            $raw_count = 'NA';
        } elsif ($t[1] != 0) {
            $cRPKM=sprintf("%.2f",$t[1]);
            $raw_count=$t[2];
        }
        
        $data{$t[0]}.="\t$cRPKM\t$raw_count";
        $RPKM{$t[0]}.="\t$cRPKM";
    }
    close INPUT;
}

print $OUTPUT "$head\n" if $cRPKMCounts;
print RPKM "$headRPKM\n";

foreach my $gene (sort (keys %data)){
    print $OUTPUT "$gene\t$names{$gene}$data{$gene}\n" if $cRPKMCounts;
    print RPKM "$gene\t$names{$gene}$RPKM{$gene}\n";
}

close $OUTPUT if $cRPKMCounts;
close RPKM;
