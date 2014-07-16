#!/usr/bin/env perl
#
# Final "Add_to_*" script that combines all previous AS event type-specific PSI
# tables into one final ("FULL") table.
#
# Event IDs are also converted to the "new" IDs, which are specified in the
# library file "New_ID-*.txt.gz" in VASTDB/FILES/ directory

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

# Load conversation table file to memory
my $NEWID = openFileHandle("$dbDir/FILES/New_ID-$sp.txt.gz");
my %newIDs;
while (<$NEWID>) {
  chomp;
  my @l = split("\t");
  if (defined $newIDs{$l[1]}) {
      die "Non-unique key value pair!\n";
  }
  $newIDs{$l[1]} = $l[0];
}
close $NEWID;

# Loads the template 
my $TEMPLATE = openFileHandle("$dbDir/TEMPLATES/$sp.FULL.Template.txt.gz");
my $h = <$TEMPLATE>;
chomp $h;
my @header = split(/\t/, $h);
@header = @header[0..5];

my %template;
while (<$TEMPLATE>){
  chomp;
  my @l = split(/\t/);
  if (defined $template{$l[1]}) {
    die "Non-unique key value pair!\n";
  }
  $template{$l[1]} = \@l;
}
close $TEMPLATE;

# Load input data
my %done;
my $sawHeader = 0;

while (<STDIN>) {
  chomp;
  
  my @l = split(/\t/);

  # Check headers
  if (/^GENE\tEVENT/) {
    if (!$sawHeader) {
      push @header, @l[6..$#l];
      $sawHeader = @l;  # store number of expected columns
      print STDOUT join("\t", @header) . "\n";
    } elsif ($sawHeader != @l) {
      die "Number of columns in subsequent header does not match. Terminating!!\n";
    }  
    next;
  }

  # Check if input is found in template
  if ($template{$l[1]}) {
    
    my @prefix = @{$template{$l[1]}};

    if ($newIDs{$prefix[1]}) {

      $prefix[1] = $newIDs{$prefix[1]};

      print STDOUT join("\t", (@prefix, @l[6..$#l])) . "\n" 
          unless $done{$l[2]};

      $done{$l[2]} = 1;
    }
  }
}

