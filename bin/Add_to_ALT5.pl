#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Cwd;
#$cwd = getcwd;
#($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

use Getopt::Long;

my $dbDir;
my $sp;
my $samLen;
my $verboseFlag;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "len=i" => \$samLen,
			  "verbose=i" => \$verboseFlag);

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine alt5]: $verbMsg\n";
  }
}

#$sp=$ARGV[0];
die "[vast combine alt5]: Needs 3-letter species key\n" if !defined($sp);
$COMB="M"; # only version implemented.

verbPrint "Parsing Template file\n";
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.ALT5.Template.txt") || die "Can't find the ALT5 template for $sp\n";
$head=<TEMPLATE>;
chomp($head);
$head_PSIs=$head_ReadCounts=$head; # setting headings for both output files
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $event=$t[1];
    $pre_data{$event}=$_; # First 6 information columns
    ($event_root,$N_ss)=$event=~/(.+)\-\d+?\/(\d+)/;
    $ALL{$event_root}=$N_ss; # keeps the total number of alternative splice sites; placeholder for event_root
}
close TEMPLATE;

@EEJ=glob("to_combine/*.ee*");
@EFF=glob("$dbDir/FILES/$sp"."_COMBI-$COMB-*gDNA.ef*");
die "Needs effective\n" if !@EFF;

verbPrint "Loading Effective files:\n";
foreach $file (@EFF){
    ($length)=$file=~/COMBI\-[A-Z]\-(\d+?)\-/;
    verbPrint "Loading: $file\tLength: $length\n";
    open (MAPPABILITY, $file);
    while (<MAPPABILITY>){
	chomp;
	@t=split(/\t/);
	($gene,$donor,$acceptor)=$t[0]=~/(.+?)\-(\d+?)\_(\d+?)\-\d+?\_\d+/;
	$eej="$gene-$donor-$acceptor";
	$eff{$length}{$eej}=$t[1];
    }
    close MAPPABILITY;
}

verbPrint "Loading EEJ data for ALT5\n";
foreach $file (@EEJ){
   # ($sample)=$file=~/COMBI\-$COMB\-\d+?\-(.+)\./; DEPRECATED --TSW
   #($sample)=$file=~/^(.*)\..*$/;
    my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
     $length = $samLen;
    # generates headings
    $head_PSIs.="\t$sample\t$sample-Q";
    $head_ReadCounts.="\t$sample-Ri\t$sample-Rtot\t$sample-Q";

    # Loads EEJ files
    open (EEJ, $file);
    while (<EEJ>){
        chomp;
        @t=split(/\t/);
	     $gene=$t[0];
	     $eej=$t[1];
        $gene_eej="$gene-$t[1]";
        $reads{$sample}{$gene_eej}=$t[2];
        ($donor,$acceptor)=$eej=~/(\d+?)\-(\d+)/;
        $last_acceptor{$gene}=$acceptor if $last_acceptor{$gene}<$acceptor; # keeps track of the last acceptor used
    }
    close EEJ;
}

# Output files
$NUM=$#EEJ+1;
open (PSIs, ">raw_incl/INCLUSION_LEVELS_ALT5-$sp$NUM-n.tab");
open (COUNTs, ">raw_reads/RAW_READS_ALT5-$sp$NUM-n.tab");
print PSIs "$head_PSIs\n";
print COUNTs "$head_ReadCounts\n";

verbPrint "Starting PSI quantifications\n";
foreach $event_root (sort keys %ALL){
    ($gene,$junctions)=$event_root=~/(.+?)\-(.+)/;
    @junctions=split(/\,/,$junctions);
    @TPC1=@TPC2=(); # empty temporary array with data to print

    for $i (0..$#junctions){ # loops through each alt donor (minimun 2)
	$ind=$i+1;
	$event="$event_root-$ind/$ALL{$event_root}"; # actual event_id
	
	$TPC1[$i]="$pre_data{$event}\t"; # temporary array for PSIs
	$TPC2[$i]="$pre_data{$event}\t"; # temporary array for read counts
    }
    
    foreach $file (@EEJ){
	#($length,$sample)=$file=~/COMBI\-[A-Z]\-(\d+?)\-(.+)\./;   
	 my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
     $length = $samLen;  #replacement --TSW
	# Emptying variables and arrays with read counts per sample
	$total_raw_reads_S=$total_corr_reads_S=0; # total simple reads
	$total_raw_reads_ALL=$total_corr_reads_ALL=0; # total complex and simple reads
	@corr_inc_reads_S=@raw_inc_reads_S=(); # Simple reads for each donor
	@corr_inc_reads_ALL=@raw_inc_reads_ALL=(); # Complex reads for each donor
	@PSI=(); # Percent splice site usage for each donor
	
	$max_mappability=$length-15;
	for $i (0..$#junctions){ #does Simple read counts
	    $eej="$gene-$junctions[$i]";
	    $corr_inc_reads_S[$i]=$max_mappability*($reads{$sample}{$eej}/$eff{$length}{$eej}) if $eff{$length}{$eej};
	    $raw_inc_reads_S[$i]=$reads{$sample}{$eej} if $eff{$length}{$eej};
	    $total_corr_reads_S+=$corr_inc_reads_S[$i];
	    $total_raw_reads_S+=$raw_inc_reads_S[$i];
	}

	for $i (0..$#junctions){ # does Simple and Complex read counts
	    ($donor)=$junctions[$i]=~/(\d+?)\-\d+/;
	    for $j (0..$last_acceptor{$gene}){  # It includes the "reference" acceptor
		$eej="$gene-$donor-$j"; 
		$corr_inc_reads_ALL[$i]+=$max_mappability*($reads{$sample}{$eej}/$eff{$length}{$eej}) if $eff{$length}{$eej};
		$raw_inc_reads_ALL[$i]+=$reads{$sample}{$eej} if $eff{$length}{$eej};
	    }
	    $total_corr_reads_ALL+=$corr_inc_reads_ALL[$i];
	    $total_raw_reads_ALL+=$raw_inc_reads_ALL[$i];
	}
	 
#### QUALITY SCORES
	 $Q="";
	 ### Score 1
	 $Q="SOK" if $total_raw_reads_ALL >= 100;
	 $Q="OK" if $total_raw_reads_ALL >= 40 && $total_raw_reads_ALL < 100;
	 $Q="LOW" if $total_raw_reads_ALL >= 25 && $total_raw_reads_ALL < 40;
	 $Q="VLOW" if $total_raw_reads_ALL >= 15 && $total_raw_reads_ALL < 25;
	 $Q="N" if $total_raw_reads_ALL < 15;
	 ### Score 2
        $Q.=",SOK" if $total_corr_reads_ALL >= 100;
        $Q.=",OK" if $total_corr_reads_ALL >= 40 && $total_corr_reads_ALL < 100;
	 $Q.=",LOW" if $total_corr_reads_ALL >= 25 && $total_corr_reads_ALL < 40;
        $Q.=",VLOW" if $total_corr_reads_ALL >= 15 && $total_corr_reads_ALL < 25;
        $Q.=",N" if $total_corr_reads_ALL < 15;
        ### Score 3 (instead of Score 4, read count)
	 $Q.=",SOK,$total_raw_reads_ALL=$total_raw_reads_S" if $total_raw_reads_S >= 100;
	 $Q.=",OK,$total_raw_reads_ALL=$total_raw_reads_S" if $total_raw_reads_S >= 40 && $total_raw_reads_S < 100;
	 $Q.=",LOW,$total_raw_reads_ALL=$total_raw_reads_S" if $total_raw_reads_S >= 25 && $total_raw_reads_S < 40;
	 $Q.=",VLOW,$total_raw_reads_ALL=$total_raw_reads_S" if $total_raw_reads_S >= 15 && $total_raw_reads_S < 25;
	 $Q.=",N,$total_raw_reads_ALL=$total_raw_reads_S" if $total_raw_reads_S < 15;

        #### Score 5: COMPLEXITY
        $from_C=$total_corr_reads_ALL-$total_corr_reads_S; # All reads minus simple reads
        $from_S=$total_corr_reads_S;

		  if ($from_C > ($from_C+$from_S)/2) {$Q.=",C3"; $Qs.=",C3";}
        elsif ($from_C > ($from_C+$from_S)/5 && $from_C <= ($from_C+$from_S)/2){$Q.=",C2";$Qs.=",C2";}
        elsif ($from_C > ($from_C+$from_S)/20 && $from_C <= ($from_C+$from_S)/5){$Q.=",C1";$Qs.=",C1";}
		  else {$Q.=",S"; $Qs.=",S";}
        ####### 

	for $i (0..$#junctions){
	    $PSI[$i]=sprintf("%.2f",100*$corr_inc_reads_ALL[$i]/$total_corr_reads_ALL) if $total_corr_reads_ALL>0;
	    $PSI[$i]="NA" if $total_corr_reads_ALL==0;

	    $Q[$i] = $Q;
             ### DIFF OUTPUT ADDITION TO QUAL SCORE!  --TSW
             ### Essentially adding the expected number of reads re-distributed to INC or EXC after normalization..
             ### These values are added to the qual score and used to infer the posterior distribution
             unless($legacyFlag) {
               my $totalN = $total_raw_reads_ALL;
               my($pPSI, $exValOfInc, $exValOfExc) = (0, 0, 0);
               unless($PSI[$i] eq "NA" or $totalN == 0) {
                 $pPSI = $PSI[$i] / 100;
                 #$exValOfInc = $pPSI * $totalN;
                 #$exValOfExc = (1-$pPSI) * $totalN;
                 $exValOfInc = sprintf("%.2f", $pPSI * $totalN);
                 $exValOfExc = sprintf("%.2f", (1-$pPSI) * $totalN);
                }
                # ALTER QUAL OUTPUT HERE>>
	        $Q[$i] .= "\@$exValOfInc,$exValOfExc";
	    }

	    $ind=$i+1;
	    $event="$event_root-$ind/$ALL{$event_root}";
	    $TPC1[$i].="$PSI[$i]\t$Q[$i]\t";
	    $TPC2[$i].="\t$raw_inc_reads_ALL[$i]\t$total_raw_reads_ALL\t$PSI[$i]=$Q[$i]";
	}
    }
    
    for $i (0..$#junctions){
		$ind=$i+1;
      $event="$event_root-$ind/$ALL{$event_root}";
		if ($pre_data{$event}){
	    print PSIs "$TPC1[$i]\n";
	    print COUNTs "$TPC2[$i]\n";
		}
    }
}
close PSIs;
close COUNTs;
