#!/usr/bin/env perl
#
# Final "Add_to_*" script that combines all previous AS event type-specific PSI
# tables into one final ("FULL") table.
#
# Event IDs are also converted to the "new" IDs, which are specified in the
# library file "New_ID-*.txt.gz" in VASTDB/FILES/ directory

use strict;
use warnings;
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

###############################################################################

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine convert]: $verbMsg\n";
  }
}

sub simplifyComplex {
  # Ad hoc routine to simplify COMPLEX types
  # (should eventually be simplified in the template source files)
  my $type = shift;
  $type =~ s/\*//;
  if ($type =~ /^ME\(.*\)$/) {
      $type = "C3";
  } elsif ($type =~ /MIC/) {
      $type = "MIC";
  }
  return $type;
}

sub reorderColumns {
  # Re-order columns if input files don't have the same sample ordering
  my $columns = shift;
  my $refOrder = shift;
  my @newOrder = (0) x keys %{$refOrder};
  for (my $i = 0; $i < @{$columns}; $i++) {
   # Iterate through each column, find out it's actual position based on the
   # original sample ordering from previous input file
   my $pos = $refOrder->{$columns->[$i]};
   $newOrder[$pos] = $i; 
  }
  return @newOrder;
}

###############################################################################

# Load conversation table file to memory
my $NEWID = openFileHandle("$dbDir/FILES/New_ID-$sp.txt.gz");
my %newIDs;
while (<$NEWID>) {
  chomp;
  my @l = split("\t");
  if (defined $newIDs{$l[1]}) {
      die "Non-unique key value pair in $NEWID!\n";
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
    die "Non-unique key value pair in $TEMPLATE!\n";
  }

  $template{$l[1]} = \@l;
}
close $TEMPLATE;

# Load input data
my %done;
my $sawHeader = 0;
my @prevSampleCols;     # remember last header
my $headerCount = 0;    # count number of columns
my %headerOrder;        # store order of samples
my @newOrder;           # used for fixing out-of-order headers

while (<STDIN>) {
  chomp;
  
  my @l = split(/\t/);

  my @sampleCols = @l[6..$#l];

  # Check headers
  if (/^GENE\tEVENT/) {
    $headerCount++;
    if (!$sawHeader) {
      push @header, @sampleCols;
      $sawHeader = @l;  # store number of expected columns
      print STDOUT join("\t", @header) . "\n";

      for (my $i=0; $i < @sampleCols; $i++) {
        $headerOrder{$sampleCols[$i]} = $i;
        push @newOrder, $i;
      }
    } elsif ($sawHeader != @l) {
      die "Number of columns in subsequent header of input file $headerCount" .
      " does not match. Terminating!\n";
    } elsif (!(@sampleCols ~~ @prevSampleCols)) {
      print STDERR "Inconsistent ordering of samples in input file $headerCount!" .
      " Re-ordering columns.\n";
      @newOrder = reorderColumns(\@sampleCols, \%headerOrder);
    }
    @prevSampleCols = @sampleCols;
    next;
  }

  # Check if input is found in template
  if ($template{$l[1]}) {
    
    my @prefix = @{$template{$l[1]}};

    if ($newIDs{$prefix[1]}) {

      $prefix[1] = $newIDs{$prefix[1]};
      $prefix[5] = simplifyComplex($prefix[5]);     # simplify complex codes

      print STDOUT join("\t", (@prefix, @sampleCols[@newOrder])) . "\n" 
          unless $done{$l[2]};

      $done{$l[2]} = 1;
    }
  }
}

