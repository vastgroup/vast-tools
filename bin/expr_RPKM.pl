#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Cwd;
use Getopt::Long;

my $sp;
my $dbDir;
my $root;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp);

if ($ARGV[0] eq "--help" || !$ARGV[0]){
    print "expr_RPKM.pl mRNA.out mRNA.eff out_dir/root\n\n";
    die "Script to calculate cRPKMs from read counts per gene, provided a file with mappability\n";
}
$root=$ARGV[2];

### Loads effective lengths (i.e. mappable positions per transcript).
open (MAPPABILITY, $ARGV[1]) || die "Needs file with effective lengths (mappability file) as \$ARGV[1]\n";
while (<MAPPABILITY>){
    chomp;
    @t=split(/\t/);
    $effective_length{$t[0]}=$t[1]; # gets the mappability for each gene
}
close MAPPABILITY;

### Read count
open (TEMP, ">$root"); # to avoid error from issue #40
my $outDir = getcwd;
while (<>){ #analyzes the bowtie output
    chomp;
    @d=split(/\t/);
    $gene=$d[2];
    if ($gene){ #ARGV[1] does not have d[2]
	$tally{$gene}++;
	$total_mapped++;
	$positions{$gene}[$d[3]]++;
    }
}
system "rm $outDir/$root"; # to avoid error from issue #40

### Calculating 3'bias of the sample
($sub_root)=$root=~/.+\/(.+)/;
$sub_root=~s/\_R1.+//;
$sub_root=~s/\-.+//;

open (BIAS3, ">$root.3bias") || die "can't open the output for 3'bias\n";

foreach $gene (sort keys %positions){
    $last_pos=$#{$positions{$gene}};
    if ($tally{$gene}>200 && $#{$positions{$gene}}>2500){ # more than 200 reads and >2500 in covered length
	### A: 500nt-500nt-500nt-500nt-500nt
	$A_tot=$A_1=$A_2=$A_3=$A_4=$A_5=0;
	for $i ($last_pos-500..$last_pos){
	    $A_1+=$positions{$gene}[$i];
	    $A_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-1000..$last_pos-501){
	    $A_2+=$positions{$gene}[$i];
	    $A_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-1500..$last_pos-1001){
	    $A_3+=$positions{$gene}[$i];
	    $A_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-2000..$last_pos-1501){
	    $A_4+=$positions{$gene}[$i];
	    $A_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-2500..$last_pos-2001){
	    $A_5+=$positions{$gene}[$i];
	    $A_tot+=$positions{$gene}[$i];
	}
	$A_perc[0]+=100*$A_1/$A_tot;
	$A_perc[1]+=100*$A_2/$A_tot;
	$A_perc[2]+=100*$A_3/$A_tot;
	$A_perc[3]+=100*$A_4/$A_tot;
	$A_perc[4]+=100*$A_5/$A_tot;
	$A_tot_genes++;
    }
    if ($tally{$gene}>500 && $#{$positions{$gene}}>5000){ # more than 500 reads and >5000 in covered length
	### B: 1000nt-1000nt-1000nt-1000nt-1000nt
	$B_tot=$B_1=$B_2=$B_3=$B_4=$B_5=0;
	for $i ($last_pos-1000..$last_pos){
	    $B_1+=$positions{$gene}[$i];
	    $B_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-2000..$last_pos-1001){
	    $B_2+=$positions{$gene}[$i];
	    $B_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-3000..$last_pos-2001){
	    $B_3+=$positions{$gene}[$i];
	    $B_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-4000..$last_pos-3001){
	    $B_4+=$positions{$gene}[$i];
	    $B_tot+=$positions{$gene}[$i];
	}
	for $i ($last_pos-5000..$last_pos-4001){
	    $B_5+=$positions{$gene}[$i];
	    $B_tot+=$positions{$gene}[$i];
	}
	$B_perc[0]+=100*$B_1/$B_tot;
	$B_perc[1]+=100*$B_2/$B_tot;
	$B_perc[2]+=100*$B_3/$B_tot;
	$B_perc[3]+=100*$B_4/$B_tot;
	$B_perc[4]+=100*$B_5/$B_tot;
	$B_tot_genes++;
    }
}
$A_av_1=sprintf("%.2f",$A_perc[0]/($A_tot_genes+0.001));
$A_av_2=sprintf("%.2f",$A_perc[1]/($A_tot_genes+0.001));
$A_av_3=sprintf("%.2f",$A_perc[2]/($A_tot_genes+0.001));
$A_av_4=sprintf("%.2f",$A_perc[3]/($A_tot_genes+0.001));
$A_av_5=sprintf("%.2f",$A_perc[4]/($A_tot_genes+0.001));
$B_av_1=sprintf("%.2f",$B_perc[0]/($B_tot_genes+0.001));
$B_av_2=sprintf("%.2f",$B_perc[1]/($B_tot_genes+0.001));
$B_av_3=sprintf("%.2f",$B_perc[2]/($B_tot_genes+0.001));
$B_av_4=sprintf("%.2f",$B_perc[3]/($B_tot_genes+0.001));
$B_av_5=sprintf("%.2f",$B_perc[4]/($B_tot_genes+0.001));

print BIAS3 "$sub_root\t500nt\t$A_av_1\t$A_av_2\t$A_av_3\t$A_av_4\t$A_av_5\t$A_tot_genes\n";
print BIAS3 "$sub_root\t1000nt\t$B_av_1\t$B_av_2\t$B_av_3\t$B_av_4\t$B_av_5\t$B_tot_genes\n";

### Calculating cRPKMs and finalizing
foreach $gene (sort (keys %effective_length)){
    $cRPKM="NA";
    $tally{$gene}=0 if !$tally{$gene};
    $tally{$gene}="NA" if $effective_length{$gene}==0;
    $cRPKM=sprintf("%.2f",1000000*(1000*$tally{$gene}/$effective_length{$gene})/$total_mapped) if $effective_length{$gene}>0;
    print STDOUT "$gene\t$cRPKM\t$tally{$gene}\n" if $gene;
}
#close OUTPUT;

#system "gzip $file";
