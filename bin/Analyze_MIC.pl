#!/usr/bin/env perl
# this script analyzes the bowtie output for the microexon (MIC) junctions; EEJs and EEEJs
# comment

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Getopt::Long;

use Cwd qw(abs_path);
$cwd = abs_path($0);
($dir)=$cwd=~/(.+)\/bin/;

my $dbDir;
my $sp;
my $read_length;
my $root;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp,
			  "readLen=i" => \$read_length, "root=s" => \$root);

### DEPRECATED --TSW
#$file=$ARGV[0];

#system "gunzip $file" if $file=~/\.gz/;
#$file=~s/\.gz//;

#($sp,$read_length)=$file=~/.+\/(.+?)MIC\-(\d+)\-[^\-]+/;
#$root=$&;


### Loads mappability/effective lengths
open (EFF, "$dbDir/FILES/$sp"."_MIC-$read_length-gDNA.eff2") || die "Needs file with effective length for each EEEJ\n";
while (<EFF>){ #loads the effective length in the hash \%eff
    chomp;
    @t=split(/\t/);
    ($event,$coord,$inc_exc,$n)=$t[0]=~/(.+)\.(.+?)\.(.+?)\.(\d+)/;
    $eej="$coord=$n";
    $eff{$event}{$inc_exc}{$eej}=$t[1];
}
close EFF;

### Loads bowtie output file and does the read count
#$INPUT = openFileHandle ($file) || die "Needs the bowtie output for MIC file\n";
while (<STDIN>){
    chomp;
    @t=split(/\t/);
    ($event,$coord,$inc_exc,$n)=$t[2]=~/(.+)\.(.+?)\.(.+?)\.(\d+)/;
    $eej="$coord=$n";
    
    ($read)=$t[0]=~/(.+)\-/;
    $read=$t[0] if !$read;

    if ($event ne $previous_event || $read ne $previous_read){	#avoid multiple counting
    	$reads{$event}{$inc_exc}{$eej}++;
    }
    #keeps a 1-line memmory in the loop
    $previous_read=$read;
    $previous_event=$event;
}
#close STDIN;

### Loads the template information (first 6 columns)
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.MIC.Template.txt") || die "Can't find the Template for MIC\n";
$head=<TEMPLATE>;
chomp($head);
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $pre_data{$t[1]}=$_;
}
close TEMPLATE;

open (O, ">to_combine/$root.micX"); # summary info (output)
print O "GENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tPSI\tRaw_reads_exc\tRaw_reads_inc\tCorr_reads_exc\tCorr_reads_inc\n";

### Calculates the PSIs
foreach $event (sort (keys %eff)){
    ### Single microexons
    if ($event!~/MICROEX\d+\-/){
	$OKinc=$OKexc=$r_inc=$r_exc=$inc=$exc=$PSI=""; # empty temporary variables for each event
	foreach $n (sort (keys %{$eff{$event}{INC}})){ # loops through each inclusion EEEJ (when alternative C1, C2, splice sites, etc)
	    ($coord,$n2)=$n=~/(.+?)\=(.+)/;
	    $full="$event.$coord.INC.$n2";
	    $OKinc=1 if $eff{$event}{INC}{$n}>0;
	    $inc+=sprintf("%.2f",($read_length-15)*$reads{$event}{INC}{$n}/$eff{$event}{INC}{$n}) if $eff{$event}{INC}{$n}>0;
	    $r_inc+=$reads{$event}{INC}{$n};
	}
	foreach $n (sort (keys %{$eff{$event}{EXC}})){ # loops through each exclusion EEJ (when alternative C1, C2, splice sites, etc)
	    ($coord,$n2)=$n=~/(.+?)\=(.+)/;
	    $full="$event.$coord.EXC.$n2";
	    $OKexc=1 if $eff{$event}{EXC}{$n}>0;
	    $exc+=sprintf("%.2f",($read_length-15)*$reads{$event}{EXC}{$n}/$eff{$event}{EXC}{$n}) if $eff{$event}{EXC}{$n}>0;
	    $r_exc+=$reads{$event}{EXC}{$n};
	}
	
	if ($OKinc && $OKexc){ # i.e. if there is any mappable position for at least one inclusion (INC) and one exclusion (EXC) junction
	    $PSI=sprintf("%.2f",100*$inc/($inc+$exc)) if ($inc+$exc)>0;
	    $PSI="NA" if ($inc+$exc)==0;
	    print O "$pre_data{$event}\t$PSI\t$r_exc\t$r_inc\t$exc\t$inc\n" if (defined $pre_data{$event});
	}
	else {
	    print O "$pre_data{$event}\tNA\tNA\tNA\tNA\tNA\n" if (defined $pre_data{$event});
	}
    }
    
    ### Events with multiple microexons
    elsif ($event=~/MICROEX\d+\-/){
	%OKinc=%OKexc=%r_inc=%r_exc=%inc=%exc=(); # empty temporary hashes for each event
	($total_exones)=$event=~/.+\-(\d+)/; # Total number of exons in the event 
	foreach $ie (sort (keys %{$eff{$event}})){
	    foreach $n (sort (keys %{$eff{$event}{$ie}})){
		($coord,$n2)=$n=~/(.+?)\=(.+)/;
		$full="$event.$coord.$ie.$n2";

		# Modified in the new version (18/01/16) --MI
		@exons=split(/\_/,$ie); # array with all (micro)exons in the junctions
		for $z (1..$total_exones){ # looping through all MICs in the group
		    $INC_read="";
		    foreach $exon (@exons){
			$INC_read=1 if $z eq $exon; # i.e. $z is one of the exons in the junction
		    }
		    if ($INC_read){
			$OKinc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			$inc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			$r_inc{$z}+=$reads{$event}{$ie}{$n};
		    }
		    else {
			if (($exons[0] eq "C1" || $exons[0]<$z) && ($exons[$#exons] eq "C2" || $exons[$#exons] > $z)){
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
		    }
		} #### End of update
	    }
	}
	for $z (1..$total_exones){ # Generates the actual event
	    ($eventN)=$event=~/(.+)\-/;
	    $eventN="$eventN-$z/$total_exones";

	    if ($OKinc{$z} && $OKexc{$z}){
		$PSI=sprintf("%.2f",100*$inc{$z}/($inc{$z}+$exc{$z})) if ($inc{$z}+$exc{$z})>0;
		$PSI="NA" if ($inc{$z}+$exc{$z})==0;
		print O "$pre_data{$eventN}\t$PSI\t$r_exc{$z}\t$r_inc{$z}\t$exc{$z}\t$inc{$z}\n" if (defined $pre_data{$eventN});
	    }
	    else {
		print O "$pre_data{$eventN}\tNA\tNA\tNA\tNA\tNA\n" if (defined $pre_data{$eventN});
	    }
	}
    }
}
close O;
