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

our $EXIT_STATUS = 0;

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine convert]: $verbMsg\n";
  }
}

sub errPrint {
    my $errMsg = shift;
    print STDERR "[vast combine error]: $errMsg\n";
    $EXIT_STATUS++; 
}

sub errPrintDie {
    my $errMsg = shift;
    errPrint $errMsg;
    exit $EXIT_STATUS if ($EXIT_STATUS != 0);
}

sub simplifyComplex {
    # Ad hoc routine to simplify COMPLEX types
    # (should eventually be simplified in the template source files)
    my $type = shift;
    $type =~ s/\*//;
    if ($type =~ /^ME\(.*\)$/) {
	$type = "C3";
    } 
    elsif ($type =~ /MIC/) {
	$type = "MIC";
    }
    elsif ($type =~ /A\_/){
	$type = "ANN";
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
    
# Removed in V2 (18/01/18): if repeated (due to assembly conversion) are re-rewritten
#  if (defined $newIDs{$l[1]}) {
#      die "Non-unique key value pair in $NEWID!\n";
#   }
#  else {
#      $newIDs{$l[1]} = $l[0]; # old_ID => new_ID
#  }
    
    $newIDs{$l[1]} = $l[0]; # old_ID => new_ID
    
    # to correct a small discordance in old human/mouse IDs --MI [23/12/15]
    if ($l[0] =~ /INT/ && $l[1] =~ /^\-/){
	my $temp_ID = "NA".$l[1];
	$newIDs{$temp_ID} = $l[0];
    }
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

# Requirement removed in V2 (if repeated, overwrites)
#  if (defined $template{$l[1]}) {
#    die "Non-unique key value pair in $TEMPLATE!\n";
#  }

  ### Added on v2.4.1 to sort coordinates consistently with Alt5 (14/04/20)
  if ($l[5] eq "Alt3"){
      my ($t_A,$t_B,$t_C,$t_D) = $l[4] =~ /(.+?)\:(.*)\,(.*)\-(.*)/;
      my $t_strand;
      $t_strand = "+" if $t_C =~ /\+/;
      $t_strand = "-" if $t_D =~ /\+/;
      
      if ($t_strand eq "-"){
	  my $new_coord = join("+", (sort{$a<=>$b}(split(/\+/,$t_D))));
	  $l[4]="$t_A:$t_B,$t_C-$new_coord";
      }
  }

  @l = @l[0..5]; # to make sure unexpected extras are not included --MI [30/12/15] (old: \@l)
  $template{$l[1]} = \@l; 
  
  # to correct a small discordance in old human/mouse IDs --MI [23/12/15]
  if ($l[5] =~ /IR/ && $l[1] =~ /^\-/){
      my $temp_ID = "NA".$l[1];
      $template{$temp_ID} = \@l; 
  }
}
close $TEMPLATE;

# Load input data
my %done;
my $sawHeader = 0;
my @samples;     
my $headerCount = 0;    # count number of columns
my %headerOrder;        # store order of samples
my @newOrder;           # used for fixing out-of-order headers
my %seen_sample;       # keeps the samples that were seen in the first file; if not matched in the rest, die 

while (<STDIN>) {
  chomp;
  
  my @l = split(/\t/);

  my @sampleCols = @l[6..$#l];

  # Check headers
  if (/^GENE\tEVENT/) {
    $headerCount++;
    
    #### Check that samples are the exact ones, irrespective of the order (31/01/18 --MI)
    foreach my $temp_sample (@sampleCols){
	if ($headerCount == 1){ # first table
	    $seen_sample{$temp_sample}=1;
	}
	else {
	    errPrintDie("Sample $temp_sample is not in all tables\n") if (!defined $seen_sample{$temp_sample});
	}
    }
    
    if (!$sawHeader) {
      push @header, @sampleCols;
      $sawHeader = @l;  # store number of expected columns
      print STDOUT join("\t", @header) . "\n";

      for (my $i=0; $i < @sampleCols; $i++) {
        $headerOrder{$sampleCols[$i]} = $i;
        push @newOrder, $i;
      }

      @samples = @sampleCols;
    } elsif ($sawHeader != @l) {
      die "Number of columns in subsequent header of input file $headerCount" .
      " does not match. Terminating!\n";
    } elsif (!(@samples ~~ @sampleCols)) {
      print STDERR "Inconsistent ordering of samples in input file $headerCount!" .
      " Re-ordering columns.\n";
      @newOrder = reorderColumns(\@sampleCols, \%headerOrder);
    } else {
      @newOrder = sort {$a <=> $b} values %headerOrder; 
    }
    next;
  }

  # Check if input is found in template
  if ($template{$l[1]}) {
    
    my @prefix = @{$template{$l[1]}};
    my $eventType = $prefix[5];

    if ($newIDs{$prefix[1]}) {

      $prefix[1] = $newIDs{$prefix[1]};
      $prefix[5] = simplifyComplex($prefix[5]);     # simplify complex codes

      print STDOUT join("\t", (@prefix, @sampleCols[@newOrder])) . "\n"
          unless $done{$eventType}{$l[2]};

      $done{$eventType}{$l[2]} = 1;
    }
  }
}
