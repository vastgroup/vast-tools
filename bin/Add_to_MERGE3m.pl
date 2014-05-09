#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Cwd;
#$cwd = getcwd;
#($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

use Getopt::Long;

my $dbDir;
my $sp;
my $samLen;
my $verboseFlag;

GetOptions("sp=s" => \$sp, "dbDir=s" => \$dbDir, "len=i" => \$samLen,
			  "verbose=i" => \$verboseFlag);


sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine merge]: $verbMsg\n";
  }
}


verbPrint "Loading data for merge..\n";
my $head1 = <STDIN>;
chomp($head1);
my(@header);
push(@header, $head1);

my(@head1) = split(/\t/,$head1);
$head1 = join("\t",@head1[6..$#head1]);

while (my $l = <STDIN>){
    chomp($l);
    if($l =~ /^GENE/) { push(@header, $l); next; }
    @t=split(/\t/);
    $data{$t[1]}=join("\t",@t[6..$#t]);
}

# Checks that all four files have the same samples in the same order
die "[vast combine merge error]: Some of the headings are not identical!\n" if ($header[0] ne $header[1] || $header[0] ne $header[2] || $header[0] ne $header[3]);

### Output file
($sp,$N_samples)=$EXSK=~/\-(.{3})(\d+)\-n/;
$output="raw_incl/INCLUSION_LEVELS_MERGE3m-$sp$N_samples-n.tab";
open (OUTPUT, ">$output");

# Loads the template 
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.MERGE3m.Template.txt") || die "Can't find the MERGE3m Template";
$head=<TEMPLATE>;
chomp($head);
@t=split(/\t/,$head);
$pre_data=join("\t",@t[0..5]);
print OUTPUT "$pre_data\t$head1\n";
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
