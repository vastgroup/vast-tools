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
my $normalize = 0;
my $install_limma = 0;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "C" => \$cRPKMCounts, "norm" => \$normalize, "install_limma" => \$install_limma);

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
my $first_sample_count=0;
foreach my $f (@files){
    my ($root)=$f=~/([^\/]+).cRPKM/;
    $head.="\t$root-cRPKM\t$root-Counts";
    $headRPKM.="\t$root";
    my $sample_count=0;
    
    open (INPUT, $f);
    while (<INPUT>){
        chomp;
        my @t=split(/\t/);
        my $cRPKM = sprintf("%.2f", 0);
        my $raw_count = 0;
	$sample_count++;

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

    if ($sample_count >0 && $first_sample_count==0){
	$first_sample_count=$sample_count;
    }
    elsif ($sample_count > 0 && $sample_count != $first_sample_count){
	die "[vast combine cRPKM error] Count files do not have the same number of genes (incorrect versions?)\n";
    }
}

print $OUTPUT "$head\n" if $cRPKMCounts;
print RPKM "$headRPKM\n";

foreach my $gene (sort (keys %data)){
    $names{$gene}="NA" if (!defined $names{$gene});
    print $OUTPUT "$gene\t$names{$gene}$data{$gene}\n" if $cRPKMCounts;
    print RPKM "$gene\t$names{$gene}$RPKM{$gene}\n";
}
close $OUTPUT if $cRPKMCounts;
close RPKM;


### Makes normalized table:
my %norm_cRPKMs;
if ($normalize){
    my $input_file = "cRPKM-$sp$index.tab";
    open (GE_2, $input_file) || die "[vast combine cRPKM error] Needs a cRPKM table\n";
    my ($input_path,$root_input);
    if ($input_file=~/\//){
	($input_path,$root_input) = $input_file =~/(.+)\/(.+)\./;
    }
    else {
	$input_path=".";
	($root_input) = $input_file =~/(.+)\./;
    }
    open (TEMP, ">$input_path/temp_cRPKMs.tab");
    while (<GE_2>){
	chomp($_);
	my @t=split(/\t/,$_);
	print TEMP "$t[0]";
	foreach my $i (2..$#t){ # not count table
	    print TEMP "\t$t[$i]";
	}
	print TEMP "\n";
    }
    close TEMP;
    close GE_2;
    
    open (Temp_R, ">$input_path/temp.R");
    print Temp_R "source(\"https://bioconductor.org/biocLite.R\")
biocLite(\"limma\")\n" if ($install_limma);
    
    print Temp_R "
library(limma)
setwd(\"$input_path/\")
matrix=as.matrix(read.table(\"temp_cRPKMs.tab\", row.names=1, header=TRUE,sep=\"\\t\"))
Nmatrix=normalizeBetweenArrays(as.matrix(matrix))
NmatrixF=cbind(Names=row.names(matrix),Nmatrix)
write.table(NmatrixF,\"$root_input-NORM.tab\",
            sep=\"\\t\",col.names=T,row.names=F,quote=F)";
    close Temp_R;
    
    system "Rscript $input_path/temp.R";
    system "rm $input_path/temp*";
}
