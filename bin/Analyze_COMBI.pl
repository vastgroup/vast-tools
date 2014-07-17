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

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp,
           "readLen=i" => \$length, "root=s" => \$root);

# ARGV[0] is now a dummy variable.
#system "gunzip $ARGV[0]" if $ARGV[0]=~/\.gz/;

### it accepts both bowtie outputs and summarySAMs
#($input)=$ARGV[0]=~/(.+\.out)/;

#$INPUT = openFileHandle ($input) || die "Needs a bowtie COMBI output\n";
#($root)=$ARGV[0]=~/(.+?COMBI\-.+?\-\d+?\-.+?)\-e/; # analysis of a combination of boundaries (COMBI)   DEPRECATED --TSW

### Parsesread counts
while (<STDIN>) {
    $read="";
    @t=split(/\t/);
    
    # General ranked search pattern for gene name. It may vary in specific cases.
#    ($read)=$t[0]=~/(.+) /;
#    ($read)=$t[0]=~/(.+)\#/ if !$read;
#    ($read)=$t[0]=~/(.+)\:/ if !$read;
#    ($read)=$t[0]=~/(.+)\// if !$read;
    ($read)=$t[0]=~/(.+)\-\d+/;
    $read=$t[0] if !$read;
    
    $hit=$t[2];
    ($gene, $donor, $acceptor, $donor_coord, $acceptor_coord)=$hit=~/(.+?)\-(\d+?)\_(\d+?)\-(\d+?)\_(\d+)/;
    $eej="$donor-$acceptor";
    
    if ($read ne $previous_read){ 
	$EEJ{$gene}{$eej}++;
	$POS{$gene}{$eej}{$t[3]}++;
    }
    $previous_read=$read;
}
#close $INPUT;

### Prints the read counts
open (OUTPUT, ">to_combine/$root.eej2"); # read counts for all EEJ
foreach $gene (keys %EEJ) {  
    foreach $eej (sort {$a<=>$b}(keys %{$EEJ{$gene}})){
	$p_p="";
	foreach $POS (sort {$a<=>$b}(keys %{$POS{$gene}{$eej}})){
	    $p_p.="$POS:$POS{$gene}{$eej}{$POS},"; # read count at each EEJ position (POS)
	}
	chop($p_p);	
	print OUTPUT "$gene\t$eej\t$EEJ{$gene}{$eej}\tNA\t$p_p\n";
    }
}
close OUTPUT;
