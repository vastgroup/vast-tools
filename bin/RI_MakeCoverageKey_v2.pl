#!/usr/bin/env perl
# New version to account for complex EEJ (v2) 
# Fairly different from v1
# By: Manuel Irimia [09/Nov/2015]
# mirimia@gmail.com

use warnings;
use strict;
use Getopt::Long;

my $helpFlag = 0;
my $dbDir=".";
my $species = $ARGV[0]; 

GetOptions("help"     => \$helpFlag,
            "dbDir=s" => \$dbDir,
            "sp=s"    => \$species,
);

if ($helpFlag) {
  print STDERR "Usage: $0 [options] combine_folder

OPTIONS:
  -dbDir DBDIR        : Database directory
  -sp Hsa/Mmu/etc     : $species
  -h, --help          : Print this message
";
  exit 1;
}

my $combineFolder = $ARGV[0];

# One intermediary file that contains the raw and the corrected read counts
my @files_counts=glob($combineFolder . "/*.summary_v2.txt"); # This contains raw and corrected counts
my $N=$#files_counts+1;

die "Error in $0: No samples found in folder $combineFolder\n" if ($N == 0);

### Mappability is no longer needed in v2 (02/10/15)

### For each file with counts
my %samples_counts;
my %corrected_reads;
my %raw_reads;
my %ne_events;

foreach my $file (@files_counts){
    my ($sample)=$file=~/([^\/]+)\.IR\.summary_v2\.txt$/;
    
    $samples_counts{$sample}=1;
    open (IN, $file);
    <IN>;
    while (<IN>){ # Format:  Event_ID cEIJ1 cEIJ2 cEEJ rEIJ1 rEIJ2 rEEJ
        chomp;
        my @temp=split(/\t/);

        $corrected_reads{$temp[0]}{'EI1'}{$sample}=$temp[1];
        $corrected_reads{$temp[0]}{'EI2'}{$sample}=$temp[2];
        $corrected_reads{$temp[0]}{'EE'}{$sample}=$temp[3];
	
        $raw_reads{$temp[0]}{'EI1'}{$sample}=$temp[4];
        $raw_reads{$temp[0]}{'EI2'}{$sample}=$temp[5];
        $raw_reads{$temp[0]}{'EE'}{$sample}=$temp[6];

	$ne_events{$temp[0]}="ne" if $temp[1] eq "ne" || $temp[2] eq "ne" || $temp[3] eq "ne"; 
    }
    close IN;
}

my $outputFile = "$combineFolder/Coverage_key_v2-$species$N.IRQ";
open (OUT, ">$outputFile") or die "Failed to open $outputFile";

print OUT "EVENT";
foreach my $sample (sort keys %samples_counts){
    print OUT "\t$sample";
}
print OUT "\n";

foreach my $event (sort keys %corrected_reads){
    print OUT "$event";
    foreach my $sample (sort keys %samples_counts){
        my $eEI=$corrected_reads{$event}{EI1}{$sample};
        my $eIE=$corrected_reads{$event}{EI2}{$sample};
        my $eEE=$corrected_reads{$event}{EE}{$sample};

        my $rEI=$raw_reads{$event}{EI1}{$sample};
        my $rIE=$raw_reads{$event}{EI2}{$sample};
        my $rEE=$raw_reads{$event}{EE}{$sample};

        $eEI=0 if !$corrected_reads{$event}{EI1}{$sample};
        $eIE=0 if !$corrected_reads{$event}{EI2}{$sample};
        $eEE=0 if !$corrected_reads{$event}{EE}{$sample};

        $rEI=0 if !$raw_reads{$event}{EI1}{$sample};
        $rIE=0 if !$raw_reads{$event}{EI2}{$sample};
        $rEE=0 if !$raw_reads{$event}{EE}{$sample};

        my $reads="$rEI=$rIE=$rEE"; # string with the raw read counts for the three junctions

        my $Q;

        ### If any of the three junctions is "ne" (=no mappability) the whole KEY is turned into N's
        if ($eIE eq "ne" || $eEI eq "ne" || $eEE eq "ne" || $ne_events{$event}) { # added last bit (02/10/15)
            $Q="N,N,NA,ne";  # score 3 set to NA for consistency --UB
        } 
	else {
            ### Corresponds to Q1 in other events
            $Q="N";
            $Q="VLOW" if (($rEI>=5 && $rIE>=10) || ($rEI>=10 && $rIE>=5) || $rEE>=10);
            $Q="LOW" if (($rEI>=10 && $rIE>=15) || ($rEI>=15 && $rIE>=10) || $rEE>=15);
            $Q="OK" if (($rEI>=15 && $rIE>=20) || ($rEI>=20 && $rIE>=15) || $rEE>=20);
            $Q="SOK" if ((($rEI>=15 && $rIE>=20) || ($rEI>=20 && $rIE>=15) || $rEE>=20) && ($rEI+$rIE+$rEE)>=100);

            ### Includes keys 2,3 and 4 of other events, but only 2 is real. 3 is always
            ### "NA" and "Reads" is just all the reads.
            ### The 5th key would be Ulrich's imbalanced test p-value and will be added later
            if (! (defined $eEI && defined $eIE && defined $eEE)) {
                $Q.=",N,NA,$reads";
            }
            elsif ((($eEI>=15 && $eIE>=20) || ($eEI>=20 && $eIE>=15) || $eEE>=20) && ($eEI+$eIE+$eEE)>=100){
                $Q.=",SOK,NA,$reads";
            }
            elsif (($eEI>=15 && $eIE>=20) || ($eEI>=20 && $eIE>=15) || $eEE>=20){
                $Q.=",OK,NA,$reads";
            }
            elsif (($eEI>=10 && $eIE>=15) || ($eEI>=15 && $eIE>=10) || $eEE>=15){
                $Q.=",LOW,NA,$reads";
            }
            elsif (($eEI>=5 && $eIE>=10) || ($eEI>=10 && $eIE>=5) || $eEE>=10){
                $Q.=",VLOW,NA,$reads";
            }
            else {
                $Q.=",N,NA,$reads";
            }
        }
        print OUT "\t$Q";
    }
    print OUT "\n";
}
