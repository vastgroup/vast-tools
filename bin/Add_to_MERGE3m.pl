#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Cwd;
$cwd = getcwd;
($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

$EXSK=$ARGV[0]; # Cassette
$MULTI=$ARGV[1]; # MULTI
$COMBI=$ARGV[2]; # COMBI
$MIC=$ARGV[3]; # MIC

# checks if files exist and opens them
open (I1, $EXSK) || die "Can't find the EXSK file\n";
open (I2, $MULTI) || die "Can't find the MULTI file\n";
open (I3, $COMBI) || die "Can't find the exsk5 file\n";
open (I4, $MIC) || die "Can't find the micX file\n";

print "Loading EXSK data\n";
$head1=<I1>;
@head1=split(/\t/,$head1);
$head1=join("\t",@head1[6..$#head1]);
while (<I1>){
    chomp;
    @t=split(/\t/);
    $data{$t[1]}=join("\t",@t[6..$#t]);
}
close I1;

print "Loading MULTI data\n";
$head2=<I2>;
@head2=split(/\t/,$head2);
$head2=join("\t",@head2[6..$#head2]);
while (<I2>){
    chomp;
    @t=split(/\t/);
    $data{$t[1]}=join("\t",@t[6..$#t]);
}
close I2;

print "Loading COMBI data\n";
$head3=<I3>;
@head3=split(/\t/,$head3);
$head3=join("\t",@head3[6..$#head3]);
while (<I3>){
    chomp;
    @t=split(/\t/);
    $data{$t[1]}=join("\t",@t[6..$#t]);
}
close I3;

print "Loading MIC data\n";
$head4=<I4>;
@head4=split(/\t/,$head4);
$head4=join("\t",@head4[6..$#head4]);
while (<I4>){
    chomp;
    @t=split(/\t/);
    $data{$t[1]}=join("\t",@t[6..$#t]);
}
close I4;

# Checks that all four files have the same samples in the same order
die "Some of the headings are not identical\n" if ($head1 ne $head2 || $head1 ne $head3 || $head1 ne $head4);

### Output file
($sp,$N_samples)=$EXSK=~/\-(.{3})(\d+)\-n/;
$output="INCLUSION_LEVELS_MERGE3m-$sp$N_samples-n.tab";
open (OUTPUT, ">$output");

# Loads the template 
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.MERGE3m.Template.txt") || die "Can't find the MERGE3m Template";
$head=<TEMPLATE>;
chomp($head);
@t=split(/\t/,$head);
$pre_data=join("\t",@t[0..5]);
print OUTPUT "$pre_data\t$head1";
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $NAMES=join("\t",@t[0..5]);
    $event=$t[1];
    print OUTPUT "$NAMES\t$data{$event}\n" unless $done{$t[2]};
    $done{$t[2]}=1;
}
close TEMPLATE;
close OUTPUT;
