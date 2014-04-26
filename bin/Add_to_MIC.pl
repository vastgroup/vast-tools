#!/usr/bin/perl

BEGIN {push @INC, '../lib'}
use FuncBasics qw(:all);

# This script is to calculate PSIs for MIC and produce a table with them for all samples.
use Cwd;
$cwd = getcwd;
($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

$sp=$ARGV[0];
die "Needs the 3-letter species key\n" if !$sp;

@files=glob("$dir/$sp/SAMPLES/$sp*micX");

$head_counts=$head_PSI="GENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX";
foreach $file (@files){
    ($sample)=$file=~/MIC\-\d+?\-(.+)\./;
    $head_PSI.="\t$sample\t$sample-Q";
    $head_counts.="\t$sample-Rexc\t$sample-Rinc\t$sample-exc\t$sample-inc\tPSI=Q";
    open (INPUT, $file);
    while (<INPUT>){
	chomp;
	@t=split(/\t/);
	$event=$t[1];
	$pre_data{$event}=join("\t",@t[0..5]) if !$pre_data{$event};

	# Quality scores (only coverage)
	$Q="";
#### Score 1: Using all raw reads
	if (($t[7]+$t[8])>=100){$Q="SOK";}
	elsif ($t[7]>=20 || $t[8]>=20){$Q="OK";}
	elsif ($t[7]>=15 || $t[8]>=15){$Q="LOW";}
	elsif ($t[7]>=10 || $t[8]>=10){$Q="VLOW";}
	else {$Q="N";}
#### Score 2: Using all corrected reads
	if (($t[9]+$t[10])>=100){$Q.=",SOK";}
	elsif ($t[9]>=20 || $t[10]>=20){$Q.=",OK";}
	elsif ($t[9]>=15 || $t[10]>=15){$Q.=",LOW";}
	elsif ($t[9]>=10 || $t[10]>=10){$Q.=",VLOW";}
	else {$Q.=",N";}
#### Score 3: Using simple (=reference, C1A, AC2, C1C2) raw reads
	if (($t[7]+$t[8])>=100){$Q.=",SOK";}
	elsif ($t[7]>=20 || $t[8]>=20){$Q.=",OK";}
	elsif ($t[7]>=15 || $t[8]>=15){$Q.=",LOW";}
	elsif ($t[7]>=10 || $t[8]>=10){$Q.=",VLOW";}
	else {$Q.=",N";}
#### No scores 4 and 5 for microexon pipeline.	
	$Q.=",na,na";

	$PSI{$event}.="\t$t[6]\t$Q";
	$read_counts{$event}.="\t".join("\t",@t[7..10])."\t$t[6]=$Q";
    }
    close INPUT;
}
$N=$#files+1;
open (PSIs, ">$dir/$sp/SAMPLES/INCLUSION_LEVELS_MIC-$sp$N-n.tab");
open (COUNTs, ">$dir/$sp/RAW_READS/RAW_READS_MIC-$sp$N-n.tab");

print PSIs "$head_PSI\n";
print COUNTs "$head_counts\n";

foreach $event (sort keys %pre_data){
    print PSIs "$pre_data{$event}"."$PSI{$event}\n";
    print COUNTs "$pre_data{$event}"."$read_counts{$event}\n";
}

close PSIs;
close COUNTs;

