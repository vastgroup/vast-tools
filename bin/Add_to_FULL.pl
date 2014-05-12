#!/usr/bin/perl
#
# Convert old event IDs with new IDs using the conversation table found in the
# database directory.
#
# e.g. Hsa/FILES/New_ID-Hsa.txt.gz or Mmu/FILES/New_ID-Mmu.txt.gz

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
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
    print STDERR "[vast combine convert]: $verbMsg\n";
  }
}

sub loadKeyVal {
  my $fileName = shift;
  my $hndl = openFileHandle($fileName);
  my(%retHash);
  while(my $l = <$hndl>) {
    chomp($l);
    my(@a) = split(/\t/, $l);
    if(defined($retHash{$a[1]})) { die "non-unique key value pair!\n"; }
    $retHash{$a[1]} = $a[0];
  }
  return(\%retHash);
}

# Load conversation table file to memory
my $newID = "$dbDir/FILES/New_ID-$sp.txt.gz";
my %newIDs = %{loadKeyVal($newID)};


my $header = undef;

while (<STDIN>) {
  chomp;
  
  # Check headers
  if (/^GENE\tEVENT/) {
    if (! defined $header) {
      $header = $_;
      print STDOUT $header . "\n";
    }   
    next;
  }

  my @l = split(/\t/);

  # Replace ID in column 2 in INCLUSION_LEVELS
  if (defined $newIDs{$l[1]}) {
    $l[1] = $newIDs{$l[1]};
  } else {
    verbPrint "Could not finding matching ID for " . $l[1] . "\n";
  }

  print STDOUT join("\t", @l)."\n";
}
