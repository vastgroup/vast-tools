#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Getopt::Long;

# In this simplified version, it only aims to get an exon-exon junction (eej) read count. It does not identify cassette events.
#use Cwd;
#$cwd = getcwd;
#($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

my $dbDir;
my $sp;
my $length;
my $root;
my $strandaware=0;  # dummy; to prevent GetOptions getting confused between -sp and -s 
my $silent=0; # non used option

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "readLen=i" => \$length, 
	   "root=s" => \$root, , "s" => \$strandaware, "ec" => \$silent);


### Parsesread counts
while (<STDIN>) {
    $read="";
    @t=split(/\t/);
    
    ($read)=$t[0]=~/(.+)\-/;
    $read=$t[0] if !$read;
    
    $hit=$t[2];
    ($gene, $donor, $acceptor, $donor_coord, $acceptor_coord)=$hit=~/(.+?)\-(\d+?)\_(\d+?)\-(\d+?)\_(\d+)/;
    $eej="$donor-$acceptor";
    $gene_eej="$gene-$eej";

    if ($read ne $previous_read){ 
	$EEJ{$gene}{$eej}++;
	$POS{$gene_eej}{$t[3]}++;
    }
    $previous_read=$read;
}

### Prints the read counts
open (OUTPUT, ">to_combine/$root.eej2"); # read counts for all EEJ
foreach $gene (keys %EEJ) {  
    foreach $eej (sort {$a<=>$b}(keys %{$EEJ{$gene}})){
	$p_p="";
	$gene_eej="$gene-$eej";
	foreach $POS (sort {$a<=>$b}(keys %{$POS{$gene_eej}})){
	    $p_p.="$POS:$POS{$gene_eej}{$POS},"; # read count at each EEJ position (POS)
	}
	chop($p_p);
	print OUTPUT "$gene\t$eej\t$EEJ{$gene}{$eej}\tNA\t$p_p\n";
    }
}
close OUTPUT;
