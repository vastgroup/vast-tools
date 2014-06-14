#!/usr/bin/perl -w
# Script to trim reads
#
# Trim settings:
# Once:
# Trims reads into one X-nt read
#
# Twice:
# Splits reads into two X-nt reads. 
# If the length of the original reads is < 2X, it splits in an overlapping manner.

use strict;
use Getopt::Long;

my $verboseFlag = 1;
my $stepSize = 25;

my $targetLength = 50;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions("verbose" => \$verboseFlag,
		    "v" => \$verboseFlag,
			 "step=i" => \$stepSize
);

sub sysErrMsg {
  my $sysCommand = shift;
  not system($sysCommand) or die "[vast trim]: $sysCommand Failed in $0!";
}

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast trim]: $verbMsg\n";
  }
}


### Initialize variables
my $total_reads = 0;
my $total_reads_accepted = 0;
my $name;
my $name2;
my $seq;
my $rest;

### Parses the original reads
my $lineCounter = 1;
while (<>){
  chomp;

  my $mod = $lineCounter % 4;
  if ($mod == 1) {
      $name = $_;
  } elsif ($mod == 2) {
      $seq = $_;
  } elsif ($mod == 3) {
      $name2 = $_;
  } elsif ($mod == 0) {
      $rest = $_;
  
      $total_reads++;
    
      my $trimNum = 1; 
      for(my $off=0; length($seq) >= $off + $targetLength; $off += $stepSize) { 
          # make sure the length of the sequence AND quality scores are >= length
 
          my ($S1) = substr($seq, $off, $targetLength);
          my ($R1) = substr($rest, $off, $targetLength);

          print STDOUT "$name-$trimNum\n$S1\n$name2-$trimNum\n$R1\n";

			 $trimNum++;
      }
      $total_reads_accepted++;
  }
  $lineCounter++;
}

verbPrint "Total processed reads: $total_reads\n";
verbPrint "Total valid reads: $total_reads_accepted\n";

if($total_reads <= 1 or $total_reads_accepted <= 1) { exit 1; }

exit $total_reads_accepted;

