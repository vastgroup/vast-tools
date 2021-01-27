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
my $verboseFlag=1;

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


sub verbPrint {
    my $verbMsg = shift;
    if($verboseFlag) {
	chomp($verbMsg);
	print STDERR "[vast combine IR]: $verbMsg\n";
    }
}

# One intermediary file that contains the raw and the corrected read counts
my @files_counts1=glob($combineFolder . "/*.summary_v2.txt"); # This contains raw and corrected counts
my @files_counts2=glob($combineFolder . "/*.summary_v2.txt.gz"); # This contains raw and corrected counts
my @files_counts=(@files_counts1,@files_counts2); # This contains raw and corrected counts
my $N=$#files_counts+1;

die "Error in $0: No samples found in folder $combineFolder\n" if ($N == 0);

### Mappability is no longer needed in v2 (02/10/15)
### Addition on 15/10/18: to incorportate the intron body read number as third score.
my %intron_sample_mappability_noSS;
my %intron_sample_mappability_SS;
my %intron_sample_cor_counts;

open (MAP_I_NOSS, "$dbDir/FILES/$species.Introns.sample.200.50.uniquecount.txt") || die "    Error: Needs the intron body sample mappability for non-SS\n"; 
while (<MAP_I_NOSS>){
    chomp;
    my @temp=split(/\t/,$_);
    $intron_sample_mappability_noSS{$temp[0]}=$temp[1]; # intron => eff positions (max = 151)
}
close MAP_I_NOSS;

open (MAP_I_SS, "$dbDir/FILES/$species.Introns.sample-SS.200.50.uniquecount.txt") || die "    Error: Needs the intron body sample mappability for SS\n"; 
while (<MAP_I_SS>){
    chomp;
    my @temp=split(/\t/,$_);
    $intron_sample_mappability_SS{$temp[0]}=$temp[1]; # intron => eff positions (max = 151)
}
close MAP_I_SS;

my %is_ss;
my @files_IR2_1=glob($combineFolder . "/*.IR2"); # This contains the corrected counts, including the I (colum 5)
my @files_IR2_2=glob($combineFolder . "/*.IR2.gz"); # This contains the corrected counts, including the I (colum 5)
my @files_IR2=(@files_IR2_1,@files_IR2_2); # This contains the corrected counts, including the I (colum 5)
my $M = $#files_IR2+1;
die "Error in $0: Different number of IR2 and IR.summary_v2.txt files\n" if ($M != $N);
foreach my $file (@files_IR2){
    my ($sample)=$file=~/([^\/]+)\.IR2/;
    
    unless(-e "to_combine/${sample}.info" || -e "to_combine/${sample}.info.gz"){ verbPrint "   $sample: do not find to_combine/${sample}.info. Sample will be treated as being not strand-specific.";
    } else{
	my $fh_info;
        if (-e "to_combine/${sample}.info.gz"){
            open($fh_info, "gunzip -c to_combine/${sample}.info.gz | ") or die "$!";
        } else {
            open($fh_info, "to_combine/${sample}.info") or die "$!";
        }
	my $line=<$fh_info>; close($fh_info);
	my @fs=split("\t",$line);
	if($fs[@fs-2] eq "-SS"){
	    $is_ss{$sample}=1;
	    verbPrint "   $sample: found to_combine/${sample}.info. Sample will be treated as being strand-specific."
	}
	else{
	    verbPrint "   $sample: found to_combine/${sample}.info. Sample will be treated as being not strand-specific."
	}
    }
    
    if ($file=~/\.gz$/){
	open (IN, "gunzip -c $file | ") || die "It cannot open the $file\n";
    } else {
        open (IN, $file);
    }
    <IN>;
    while (<IN>){ #Format:  Event_ID EIJ1 EIJ2 EEJ I (all corrected)
        chomp;
        my @temp=split(/\t/,$_);
	$temp[4]=sprintf("%.1f",$temp[4]);
	$intron_sample_mappability_noSS{$temp[0]}=0 if (!defined $intron_sample_mappability_noSS{$temp[0]});
	$intron_sample_mappability_SS{$temp[0]}=0 if (!defined $intron_sample_mappability_SS{$temp[0]});
	
	unless ($is_ss{$sample}){
	    if ($intron_sample_mappability_noSS{$temp[0]}>0){
		$intron_sample_cor_counts{$temp[0]}{$sample}=$temp[4]."=".$intron_sample_mappability_noSS{$temp[0]};
	    }
	    else {
		$intron_sample_cor_counts{$temp[0]}{$sample}="NA=0";
	    }
	}
	else {
	    if ($intron_sample_mappability_SS{$temp[0]}>0){
		$intron_sample_cor_counts{$temp[0]}{$sample}=$temp[4]."=".$intron_sample_mappability_SS{$temp[0]};
	    }
	    else {
		$intron_sample_cor_counts{$temp[0]}{$sample}="NA=0";
	    }
	}
    }
    close IN;
}
#### Update finished

### For each file with counts
my %samples_counts;
my %corrected_reads;
my %raw_reads;
my %ne_events;

foreach my $file (@files_counts){
    my ($sample)=$file=~/([^\/]+)\.IR\.summary_v2\.txt/;
    
    $samples_counts{$sample}=1;

    if ($file=~/\.gz$/){
	open (IN, "gunzip -c $file | ") || die "It cannot open the $file\n";
    } else {
        open (IN, $file);
    }
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
	
        $eEI=sprintf("%.1f",$eEI) if $eEI=~/\d/;
        $eIE=sprintf("%.1f",$eIE) if $eIE=~/\d/;
        $eEE=sprintf("%.1f",$eEE) if $eEE=~/\d/;

        my $rEI=$raw_reads{$event}{EI1}{$sample};
        my $rIE=$raw_reads{$event}{EI2}{$sample};
        my $rEE=$raw_reads{$event}{EE}{$sample};

        $eEI=0 if !$corrected_reads{$event}{EI1}{$sample};
        $eIE=0 if !$corrected_reads{$event}{EI2}{$sample};
        $eEE=0 if !$corrected_reads{$event}{EE}{$sample};
	
        $rEI=0 if !$raw_reads{$event}{EI1}{$sample};
        $rIE=0 if !$raw_reads{$event}{EI2}{$sample};
        $rEE=0 if !$raw_reads{$event}{EE}{$sample};

        my $reads="$rEI=$rIE=$rEE"; # string with the corrected read counts for the three junctions (from v2.1.3)

        my $Q;

        ### If any of the three junctions is "ne" (=no mappability) the whole KEY is turned into N's
        if ($eIE eq "ne" || $eEI eq "ne" || $eEE eq "ne" || $ne_events{$event}) { # added last bit (02/10/15)
            $Q="N,N,NA=0,ne";  # score 3 set to NA for consistency --UB
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
	    # to solve a warning for some cases without mappability (likely older versions)
	    $intron_sample_cor_counts{$event}{$sample}="NA=0" if (!defined $intron_sample_cor_counts{$event}{$sample});
            if (! (defined $eEI && defined $eIE && defined $eEE)) {
                $Q.=",N,$intron_sample_cor_counts{$event}{$sample},$reads";
            }
            elsif ((($eEI>=15 && $eIE>=20) || ($eEI>=20 && $eIE>=15) || $eEE>=20) && ($eEI+$eIE+$eEE)>=100){
                $Q.=",SOK,$intron_sample_cor_counts{$event}{$sample},$reads";
            }
            elsif (($eEI>=15 && $eIE>=20) || ($eEI>=20 && $eIE>=15) || $eEE>=20){
                $Q.=",OK,$intron_sample_cor_counts{$event}{$sample},$reads";
            }
            elsif (($eEI>=10 && $eIE>=15) || ($eEI>=15 && $eIE>=10) || $eEE>=15){
                $Q.=",LOW,$intron_sample_cor_counts{$event}{$sample},$reads";
            }
            elsif (($eEI>=5 && $eIE>=10) || ($eEI>=10 && $eIE>=5) || $eEE>=10){
                $Q.=",VLOW,$intron_sample_cor_counts{$event}{$sample},$reads";
            }
            else {
                $Q.=",N,$intron_sample_cor_counts{$event}{$sample},$reads";
            }
        }
        print OUT "\t$Q";
    }
    print OUT "\n";
}
