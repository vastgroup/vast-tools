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
use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Getopt::Long;

my $trim = "once";
my $verboseFlag = 1;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions("trim=s" => \$trim,
		    "verbose" => \$verboseFlag,
		    "v" => \$verboseFlag
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

if (! ($trim eq "once" || $trim eq "twice")) {
  sysErrMsg "Invalid trim option\n";    
}

my $length=$ARGV[1];
sysErrMsg "Not a valid length\n" if ($length !~ /^[0-9]+$/);

### Initialize variables
my $total_reads = 0;
my $total_reads_accepted = 0;
my $name;
my $name2;
my $seq;
my $rest;

### Parses the original reads
my $INPUT = openFileHandle ($ARGV[0]);
my $lineCounter = 1;
while (<$INPUT>){
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

      if (length($seq)>=$length && length($rest)>=$length) { 
          # make sure the length of the sequence AND quality scores are >= length

          my ($S1)=$seq=~/^(.{$length})/;
          my ($R1)=$rest=~/^(.{$length})/;
          print STDOUT "$name-1\n$S1\n$name2-1\n$R1\n";

          if ($trim eq "twice") {
              my ($S2)=$seq=~/(.{$length})$/;
              my ($R2)=$rest=~/(.{$length})$/;
              print STDOUT "$name-2\n$S2\n$name2-2\n$R2\n";
          }
          $total_reads_accepted++;
      }
  }
  $lineCounter++;
}

verbPrint "Total processed reads: $total_reads\n";
verbPrint "Total valid reads: $total_reads_accepted\n";

if($total_reads <= 1 or $total_reads_accepted <= 1) { exit 1; }
