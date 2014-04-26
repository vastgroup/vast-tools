#!/usr/bin/perl
# This script puts all cRPKMs together in two tables (one alone, the other one with the raw read counts for each sample).

BEGIN {push @INC, '../lib'}
use FuncBasics qw(:all);

use Cwd;
$cwd = getcwd;
($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

$sp=$ARGV[0];
die "Needs Species\n" if !$sp;
@files=glob("$sp*_exprRPKM.txt");
$index=$#files+1;

# creates output files where run
open (OUTPUT, ">$sp"."_cRPKMs-$index.tab");
open (RPKM, ">$sp"."_cRPKMs_only-$index.tab");

open (NAMES, "$dir/$sp/FILES/$sp.ID.names.txt");
while (<NAMES>){
    chomp;
    @t=split(/\t/);
    $names{$t[0]}=$t[1];
}
close NAMES;

$headRPKM=$head="ID\tNAME";
foreach $f (@files){
    ($root)=$f=~/\-\d+?\-(.+?)\_exprRP/;
    $head.="\t$root-R\t$root-#";
    $headRPKM.="\t$root";
    open (INPUT, $f);
    while (<INPUT>){
	chomp;
	@t=split(/\t/);
	$cRPKM=sprintf("%.2f",$t[1]);
	$raw_count=$t[2];
	
	$data{$t[0]}.="\t$cRPKM\t$raw_count";
	$RPKM{$t[0]}.="\t$cRPKM";
    }
    close INPUT;
}

print OUTPUT "$head\n";
print RPKM "$headRPKM\n";

foreach $gene (sort (keys %data)){
    print OUTPUT "$gene\t$names{$gene}$data{$gene}\n";
    print RPKM "$gene\t$names{$gene}$RPKM{$gene}\n";
}
