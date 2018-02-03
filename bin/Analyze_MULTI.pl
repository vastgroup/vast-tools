#!/usr/bin/env perl
# This script parses the bowtie output against the EEJ for MULTIEX events (a priori complex events).

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Getopt::Long;

use Cwd qw(abs_path);
$cwd = abs_path($0);
($dir)=$cwd=~/(.+)\/bin/;

my $dbDir;
my $sp;
my $length;
my $root;
my $strandaware=0;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp,
           "readLen=i" => \$length, "root=s" => \$root,  "s" => \$strandaware);

my $mapcorr_fileswitch=""; if($strandaware){$mapcorr_fileswitch="-SS"}

#($sp,$length)=$ARGV[0]=~/(\S{3})MULTI\-(\d+?)\-/; # uniform format
#($file)=$ARGV[0]=~/(.+?\.out)/; # input file: bowtie output    # DEPRECATED --TSW
#system "gunzip $ARGV[0]" if $ARGV[0]=~/\.gz/;

### Loads template information
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.MULTI.Template.1.txt") || die "Can't find MULTI template 1\n";
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $event=$t[3];
    $pre_data{$event}=join("\t",@t[0..11]);
    $ref_C1{$event}=$t[12]; # Pre-defined Reference C1 for this exon
    $ref_C2{$event}=$t[13]; # Pre-defined Reference C2 for this exon
    $post_data{$event}=join("\t",@t[14..16]);

    # Keeps track of the number of exons in the event
    ($event_root,$N)=$event=~/(.+)\-\d+?\/(\d+)/;
    $Nex{$event_root}=$N;
}
close TEMPLATE;

### Loads mappability for each EEJ
open (MAPPABILITY, "$dbDir/FILES/MULTI-$length-gDNA${mapcorr_fileswitch}.eff") || die "Can't find mappability file for MULTI\n";
while (<MAPPABILITY>){
    chomp;
    @t=split(/\t/);
    ($event,$eej)=$t[0]=~/(.+)_.+?_.+?_(.+?\-.+?\..+)/;
    $eff{$event}{$eej}=$t[1];
}
close MAPPABILITY;

### Loads and counts raw reads
#$INPUT = openFileHandle($ARGV[0]); # DEPRECATED --TSW

while (<STDIN>){
    chomp;
    @t=split(/\t/);

    # obtains original read source
    ($read)=$t[0]=~/(.+)\-/;
    $read=$t[0] if !$read;

    ($event,$eej)=$t[2]=~/(.+)_.+?_.+?_(.+?\-.+?\..+)/;
    
    if ($event ne $previous_event || $read ne $previous_read){	
	$read_count{$event}{$eej}++;
    }

    #keeps a 1-line memmory in the loop                                                                             
    $previous_read=$read;
    $previous_event=$event;
}
#close $INPUT;

#($root)=$file=~/(.+\-$length\-.+?)\-e/;
open (O, ">to_combine/$root.MULTI3X");
print O "Ensembl_ID\tA_coord\tStrand\tEvent_ID\tFullCoord\tType\tLength\tC1_coord\tC2_coord\tLength_Int1\tLength_Int2\t3n\t";
print O "PSI\tReads_exc\tReads_inc1\tReads_inc2\tSum_all_Reads\tRef_C1\tRef_C2\tExcl_reads(corr_all=raw_ref=corr_ref)\tInc1_reads(corr_all=raw_ref=corr_ref)\tInc2_reads(corr_all=raw_ref=corr_ref)\tComplexity\tPre_C1\tPre_C2\tGene_name\n";

### Calculates PSIs
$cle=$length-15; #corrected length to correct reads.
foreach $event_root (sort (keys %eff)){
    ### Empty all temporary variables and hashes
    %I1=%rI1=%I2=%rI2=%EXCL=%rEXCL=(); # will keep read counts (raw and corrected) for each exon in event
    %refI1=%refI2=%refE=%RrefI1=%RrefI2=%RrefE=(); # will keep read counts (raw and corrected) for each exon in event FROM refence EEJs only

    foreach $eej (sort (keys %{$eff{$event_root}})){
	($u, $d, $n)=$eej=~/(.+?)\-(.+?)\.(.+)/; # structure of EEJ: upstream_exon-downstream_exon.variant => $u-$d.$n
	
	### C1-C2 are named so only in C1C2 EEJs. Otherwise, C1=0 and C2=$Nex{$event_root_root}+1
      	if ($eff{$event_root}{$eej}>0){ # Mappability of the EEJ > 0
	    # Corrected read counts
	    $I1{$d}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
	    $I2{$u}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
	    # Raw read counts
	    $rI1{$d}+=$read_count{$event_root}{$eej};
	    $rI2{$u}+=$read_count{$event_root}{$eej};

	    # Calculates inclusion from reference only
	    $u2=$u;
	    $d2=$d;
	    $u2="C1" if $u==0;
	    $d2="C2" if $d==$Nex{$event_root}+1;
	    $event1d="$event_root-$d/$Nex{$event_root}";
	    $event1u="$event_root-$u/$Nex{$event_root}";

	    if ($u2 eq $ref_C1{$event1d}){
		$refI1{$d}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
		$RrefI1{$d}+=$read_count{$event_root}{$eej}; 
	    }
	    if ($d2 eq $ref_C2{$event1u}){
		$refI2{$u}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
		$RrefI2{$u}+=$read_count{$event_root}{$eej};
	    }

	    if ($u eq "C1"){ # all skipped (only in C1C2)
		for $ii (1..$Nex{$event_root}){ # loops through all exons in the event
		    $EXCL{$ii}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
		    $rEXCL{$ii}+=$read_count{$event_root}{$eej};
		    
		    # Checks if these are the reference C1C2 for exclusion
		    $event1ii="$event_root-$ii/$Nex{$event_root}";
		    if ($ref_C2{$event1ii} eq "C2" && $ref_C1{$event1ii} eq "C1"){
			$refE{$ii}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
			$RrefE{$ii}+=$read_count{$event_root}{$eej};
		    }
		}
	    }
	    else {
		for $inter ($u+1..$d-1){ # exclusion for any exon between the two exons in the EEJ being analyzed
		    $EXCL{$inter}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
		    $rEXCL{$inter}+=$read_count{$event_root}{$eej};

		    # Checks if these are the reference C1C2 for exclusion
		    $event1inter="$event_root-$inter/$Nex{$event_root}";
		    if ($ref_C2{$event1inter} eq $d2 && $ref_C1{$event1inter} eq $u2){
			$refE{$inter}+=sprintf("%.2f",$cle*$read_count{$event_root}{$eej}/$eff{$event_root}{$eej});
			$RrefE{$inter}+=$read_count{$event_root}{$eej};
		    }
		}
	    }
	}
    }

    for $i (1..$Nex{$event_root}){
	$event_F="$event_root-$i/$Nex{$event_root}"; # rebuilds the final event ID
	next if !$pre_data{$event_F}; # Prefiltered events only
	
	# Complexity definition scores are explained in scripts from Step 2
	$from_S=$refI1{$i}+$refI2{$i}+$refE{$i}; # reads coming only from the reference EEJs (refI1, refI2 and refE)
	$from_C=($I1{$i}+$I2{$i}+$EXCL{$i})-$from_S; # all other reads
	if ($from_C > ($from_C+$from_S)/2) {$Q="C3";}
	elsif ($from_C > ($from_C+$from_S)/5 && $from_C <= ($from_C+$from_S)/2){$Q="C2";}
	elsif ($from_C > ($from_C+$from_S)/20 && $from_C <= ($from_C+$from_S)/5){$Q="C1";}
        else {$Q="S";}

	# calculates PSIs
	$PSI=sprintf("%.2f",100*($I1{$i}+$I2{$i})/(($I1{$i}+$I2{$i})+2*$EXCL{$i})) if (($I1{$i}+$I2{$i})+2*$EXCL{$i})>0;
	$PSI="NA" if (($I1{$i}+$I2{$i})+2*$EXCL{$i})==0;
	$sum_all_EEJ=$rI1{$i}+$rI2{$i}+$rEXCL{$i};

	### In case some of the EEJ don't have mappability
	### col 13-16 (0-based): $rEXCL{$i}\t$rI1{$i}\t$rI2{$i}\t$sum_all_EEJ
	$rEXCL{$i}="NA" if $rEXCL{$i}!~/\d/; $rI1{$i}="NA" if $rI1{$i}!~/\d/;
	$rI2{$i}="NA" if $rI2{$i}!~/\d/; $sum_all_EEJ="NA" if $sum_all_EEJ!~/\d/;
	if ($rEXCL{$i} eq "NA" || $rI1{$i} eq "NA" || $rI2{$i} eq "NA"){
	    $PSI="NA";
	}
	
	print O "$pre_data{$event_F}\t$PSI\t$rEXCL{$i}\t$rI1{$i}\t$rI2{$i}\t$sum_all_EEJ\t$ref_C1{$event_F}\t$ref_C2{$event_F}\t";
	print O "$EXCL{$i}=$RrefE{$i}=$refE{$i}\t$I1{$i}=$RrefI1{$i}=$refI1{$i}\t$I2{$i}=$RrefI2{$i}=$refI2{$i}\t$Q\t$post_data{$event_F}\n";	    
    }
}
