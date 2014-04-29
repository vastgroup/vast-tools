#!/usr/bin/perl
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
    $reads{$event}{$inc_exc}{$eej}++;
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

open (O, ">spli_out/$root.micX"); # summary info (output)
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
	    print O "$pre_data{$event}\t$PSI\t$r_exc\t$r_inc\t$exc\t$inc\n";
	}
	else {
	    print O "$pre_data{$event}\tNA\tNA\tNA\tNA\tNA\n";
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
		
		if (($i,$j,$k,$l)=$ie=~/(.+?)\_(.+?)\_(.+?)\_(.+)/){ # for EEEEJ (i.e. junction of 4 exons)
		    for $z (1..$total_exones){ # then it loops through each exon and compares to those in the EEEEJ
			if ($z eq $i || $z eq $j || $z eq $k || $z eq $l){
			    $OKinc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $inc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_inc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($i eq "C1" && $l eq "C2") {
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($i eq "C1" && $z<$l) {
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($z>$i && $l eq "C2") {
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($z>$i && $z<$l){
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
		    }
		}
		elsif (($i,$j,$k)=$ie=~/(.+?)\_(.+?)\_(.+)/){ # for EEEEJ (i.e. junction of 3 exons)
		    for $z (1..$total_exones){
			if ($z eq $i || $z eq $j || $z eq $k){
			    $OKinc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $inc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_inc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($i eq "C1" && $k eq "C2") {
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($i eq "C1" && $z<$k) {
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($z>$i && $k eq "C2") {
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
			elsif ($z>$i && $z<$k){
			    $OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			    $exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			    $r_exc{$z}+=$reads{$event}{$ie}{$n};
			}
		    }
		}
		elsif ($ie=~/C1\_C2/){
		    for $z (1..$total_exones){
			$OKexc{$z}=1 if $eff{$event}{$ie}{$n}>0;
			$exc{$z}+=sprintf("%.2f",($read_length-15)*$reads{$event}{$ie}{$n}/$eff{$event}{$ie}{$n}) if $eff{$event}{$ie}{$n}>0;
			$r_exc{$z}+=$reads{$event}{$ie}{$n};
		    }
		}
	    }
	}
	for $z (1..$total_exones){ # Generates the actual event
	    ($eventN)=$event=~/(.+)\-/;
	    $eventN="$eventN-$z/$total_exones";

	    if ($OKinc{$z} && $OKexc{$z}){
		$PSI=sprintf("%.2f",100*$inc{$z}/($inc{$z}+$exc{$z})) if ($inc{$z}+$exc{$z})>0;
		$PSI="NA" if ($inc{$z}+$exc{$z})==0;
		print O "$pre_data{$eventN}\t$PSI\t$r_exc{$z}\t$r_inc{$z}\t$exc{$z}\t$inc{$z}\n";
	    }
	    else {
		print O "$pre_data{$eventN}\tne\tne\tne\tne\tne\n";
	    }
	}
    }
}
