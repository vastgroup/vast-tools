#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Cwd;
use Getopt::Long;

my $dbDir;
my $sp;
my $verboseFlag;
my $samLen;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "verbose=i" => \$verboseFlag,
			  "len=i" => \$samLen);

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine alt3]: $verbMsg\n";
  }
}

#$sp=$ARGV[0];
die "[vast combine alt3]: Needs 3-letter species key\n" if !defined($sp);
$COMB="M"; # only version implemented.

verbPrint "Parsing Template file\n";
open (TEMPLATE, "$dbDir/TEMPLATES/$sp.ALT3.Template.txt") || die "Can't find the ALT3 template for $sp\n";
$head=<TEMPLATE>;
chomp($head);
$head_PSIs=$head_ReadCounts=$head; # setting headings for both output files
while (<TEMPLATE>){
    chomp;
    @t=split(/\t/);
    $event=$t[1];
    $pre_data{$event}=$_; # also acts as holder for all event ids
    ($event_root,$N_ss)=$event=~/(.+)\-\d+?\/(\d+)/;
    $ALL{$event_root}=$N_ss; # keeps the total number of alternative splice sites

    ### gets the strand of gene
    ($gene)=$event=~/(.+?)\-/;
    ($C1)=$t[4]=~/\:(\d+)/;
    ($core)=$t[4]=~/\:.+?\,[^\d]*?(\d+)/;
    die "Can't identify a C1 ($C1) or core ($core) coordinate for $event\n" if !$C1 || !$core;
    $strand{$gene}="+" if $C1<$core;
    $strand{$gene}="-" if $C1>=$core;
}
close TEMPLATE;

@EEJ=glob("to_combine/*.eej2");

@EFF=glob("$dbDir/FILES/$sp"."_COMBI-$COMB-*gDNA.ef*");
die "[vast combine alt3]: Needs strand-unspecific effective from database!\n" if !@EFF;
verbPrint "Loading mappability files (strand-unspecific and strand-specific):\n";
foreach $file (@EFF){
    ($length)=$file=~/COMBI\-[A-Z]\-(\d+?)\-/;
    verbPrint "Loading: $file\tLength: $length\n";
    open (MAPPABILITY, $file);
    while (<MAPPABILITY>){
	chomp;
	@t=split(/\t/);
	($gene,$donor,$acceptor,$co_donor,$co_acceptor)=$t[0]=~/(.+?)\-(\d+?)\_(\d+?)\-(\d+?)\_(\d+)/;
	$eej="$gene-$donor-$acceptor";
	$eff_ns{$length}{$eej}=$t[1];
	
	### keeps the coordinates for global PSI-like (04/05/19)
	$co_donor{$gene}{$donor}=$co_donor;
	$co_acceptor{$gene}{$acceptor}=$co_acceptor;
    }
    close MAPPABILITY;
}

@EFF=glob("$dbDir/FILES/$sp"."_COMBI-$COMB-*gDNA-SS.ef*");
die "[vast combine alt3]: Needs strand-specific effective from database!\n" if !@EFF;
foreach $file (@EFF){
    ($length)=$file=~/COMBI\-[A-Z]\-(\d+?)\-/;
    verbPrint "Loading: $file\tLength: $length\n";
    open (MAPPABILITY, $file);
    while (<MAPPABILITY>){
	chomp;
	@t=split(/\t/);
	($gene,$donor,$acceptor)=$t[0]=~/(.+?)\-(\d+?)\_(\d+?)\-\d+?\_\d+/;
	$eej="$gene-$donor-$acceptor";
	$eff_ss{$length}{$eej}=$t[1];
    }
    close MAPPABILITY;
}

my %is_ss;
verbPrint "Loading EEJ data for ALT3\n";
foreach $file (@EEJ){
    my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
    $length = $samLen; # replacement --TSW
    
    # generates headings
    $head_PSIs.="\t$sample\t$sample-Q";
    $head_ReadCounts.="\t$sample-Ri\t$sample-Rtot\t$sample-Q";

    unless(-e "to_combine/${sample}.info"){ verbPrint "   $sample: do not find to_combine/${sample}.info. Sample will be treated as being not strand-specific.";
    }else{
    	open(my $fh_info,"to_combine/${sample}.info") or die "$!"; my $line=<$fh_info>; close($fh_info);
    	my @fs=split("\t",$line);
    	if($fs[@fs-2] eq "-SS"){
    		$is_ss{$sample}=1;
    		verbPrint "   $sample: found to_combine/${sample}.info. Sample will be treated as being strand-specific."
    	}else{
    		verbPrint "   $sample: found to_combine/${sample}.info. Sample will be treated as being not strand-specific."
    	}
    }
    ### Loads EEJ files
    open (EEJ, $file);
    while (<EEJ>){
        chomp;
        @t=split(/\t/);
 	$gene=$t[0];
	$eej=$t[1];
        $gene_eej="$gene-$t[1]";
        $reads{$sample}{$gene_eej}=$t[2];
        ($donor,$acceptor)=$eej=~/(\d+?)\-(\d+)/;
	$last_donor{$gene}=$donor if $last_donor{$gene}<$donor; # keeps track of the last donor used
	$last_acceptor{$gene}=$acceptor if $last_acceptor{$gene}<$acceptor; # keeps track of the last acceptor used 
    }
    close EEJ;
}

# Output files
$NUM=$#EEJ+1;
open (PSIs, ">raw_incl/INCLUSION_LEVELS_ALT3-$sp$NUM-n.tab"); #change output directories --TSW
open (COUNTs, ">raw_reads/RAW_READS_ALT3-$sp$NUM-n.tab");
print PSIs "$head_PSIs\n";
print COUNTs "$head_ReadCounts\n";

my $eff_href;
verbPrint "Starting PSI quantification\n";
foreach $event_root (sort (keys %ALL)){
    ($gene,$junctions)=$event_root=~/(.+?)\-(.+)/;
    @junctions=split(/\,/,$junctions);
    @TPC1=@TPC2=(); # empty temporary array with data to print
    
    for $i (0..$#junctions){ # loops through each alt acceptor (minimun 2)
	$ind=$i+1;
	$event="$event_root-$ind/$ALL{$event_root}"; # actual event_id

	$TPC1[$i]="$pre_data{$event}\t"; # temporary array for PSIs
        $TPC2[$i]="$pre_data{$event}\t"; # temporary array for read counts
    }
    
    foreach $file (@EEJ){
	my $fname = $file;
	$fname =~ s/^.*\///;
	($sample)=$fname=~/^(.*)\..*$/;
	$length = $samLen;  #replacement --TSW
	
	# Emptying variables and arrays with read counts per sample
	$total_raw_reads_S=$total_corr_reads_S=0; # total simple reads
	$total_raw_reads_ALL=$total_corr_reads_ALL=0; # total complex and simple reads
	@corr_inc_reads_S=@raw_inc_reads_S=(); # Simple reads for each acceptor
	@corr_inc_reads_ALL=@raw_inc_reads_ALL=(); # Complex reads for each acceptor
	@PSI=(); # Percent splice site usage for each acceptor
	
	# set mappability correction accordingly
	if($is_ss{$sample}){$eff_href=\%eff_ss}else{$eff_href=\%eff_ns}
	
	$max_mappability=$length-15;
	for $i (0..$#junctions){ #does Simple read counts
	    $eej="$gene-$junctions[$i]";
	    $corr_inc_reads_S[$i]=$max_mappability*($reads{$sample}{$eej}/$eff_href->{$length}{$eej}) if $eff_href->{$length}{$eej};
	    $raw_inc_reads_S[$i]=$reads{$sample}{$eej} if $eff_href->{$length}{$eej};
	    $total_corr_reads_S+=$corr_inc_reads_S[$i];
	    $total_raw_reads_S+=$raw_inc_reads_S[$i];
	}

	for $i (0..$#junctions){ # does Simple and Complex read counts
	    ($acceptor)=$junctions[$i]=~/\d+?\-(\d+)/;	    
	    for $j (0..$last_donor{$gene}){ # any donor to the tested acceptor
		$eej="$gene-$j-$acceptor";
		$corr_inc_reads_ALL[$i]+=$max_mappability*($reads{$sample}{$eej}/$eff_href->{$length}{$eej}) if $eff_href->{$length}{$eej};
		$raw_inc_reads_ALL[$i]+=$reads{$sample}{$eej} if $eff_href->{$length}{$eej};
	    }
	    # "internal usage" reads
	    $total_corr_reads_ALL+=$corr_inc_reads_ALL[$i];
	    $total_raw_reads_ALL+=$raw_inc_reads_ALL[$i];
	}
	# reads that jump over the event (any upstream donor to any downstream acceptor)
	$skipping_corr_reads=0; $skipping_raw_reads=0;

	($ext_acceptor)=$junctions[$#junctions]=~/\d+?\-(\d+)/;
	($int_acceptor)=$junctions[0]=~/\d+?\-(\d+)/;

	for $t_don (0..$last_donor{$gene}){
	    for $t_acc (0..$last_acceptor{$gene}){ # redundant call for ext/int, just in case
		if ($strand{$gene} eq "+" && 
		    $co_donor{$gene}{$t_don} < $co_acceptor{$gene}{$ext_acceptor} && $co_donor{$gene}{$t_don} < $co_acceptor{$gene}{$int_acceptor} &&
		    $co_acceptor{$gene}{$t_acc} > $co_acceptor{$gene}{$ext_acceptor} && $co_acceptor{$gene}{$t_acc} > $co_acceptor{$gene}{$int_acceptor}){

		    $eej="$gene-$t_don-$t_acc";
		    $skipping_corr_reads+=$max_mappability*($reads{$sample}{$eej}/$eff_href->{$length}{$eej}) if $eff_href->{$length}{$eej};
		    $skipping_raw_reads+=$reads{$sample}{$eej} if $eff_href->{$length}{$eej};		    
		}
		elsif ($strand{$gene} eq "-" && 
		    $co_donor{$gene}{$t_don} > $co_acceptor{$gene}{$ext_acceptor} && $co_donor{$gene}{$t_don} > $co_acceptor{$gene}{$int_acceptor} &&
		    $co_acceptor{$gene}{$t_acc} < $co_acceptor{$gene}{$ext_acceptor} && $co_acceptor{$gene}{$t_acc} < $co_acceptor{$gene}{$int_acceptor}){
		    
		    $eej="$gene-$t_don-$t_acc";
		    $skipping_corr_reads+=$max_mappability*($reads{$sample}{$eej}/$eff_href->{$length}{$eej}) if $eff_href->{$length}{$eej};
		    $skipping_raw_reads+=$reads{$sample}{$eej} if $eff_href->{$length}{$eej};		    
		}
	    }
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
        ### Score 3: PSI-like score for the event as a whole (instead of coverage score by raw reads) ## main diff from 04/05/19
	$PSI_like=sprintf("%.2f",100*$total_corr_reads_ALL/($total_corr_reads_ALL+$skipping_corr_reads)) if ($total_corr_reads_ALL+$skipping_corr_reads)>0;
	$PSI_like="NA" if ($total_corr_reads_ALL+$skipping_corr_reads) == 0;
	$Q.=",$PSI_like";

	### Scores 4 and 5 moved to splice site loop (v2.2.2, 11/05/19)

        ####### 

	for $i (0..$#junctions){
	    $PSI[$i]=sprintf("%.2f",100*$corr_inc_reads_ALL[$i]/$total_corr_reads_ALL) if $total_corr_reads_ALL>0;
	    $PSI[$i]="NA" if $total_corr_reads_ALL==0;

	    ### Completes the Q scores:
	    $Q[$i] = $Q; # only 3 scores so far	    

	    ### Score 4: now allows to have
	    $Q[$i].=",$raw_inc_reads_ALL[$i]=$total_raw_reads_ALL=$skipping_raw_reads";
	    
	    #### Score 5: COMPLEXITY
	    $from_C=$total_corr_reads_ALL-$total_corr_reads_S; # All reads minus simple reads
	    $from_S=$total_corr_reads_S;
	    
	    if ($from_C > ($from_C+$from_S)/2) {$Q[$i].=",C3";}
	    elsif ($from_C > ($from_C+$from_S)/5 && $from_C <= ($from_C+$from_S)/2){$Q[$i].=",C2";}
	    elsif ($from_C > ($from_C+$from_S)/20 && $from_C <= ($from_C+$from_S)/5){$Q[$i].=",C1";}
	    else {$Q[$i].=",S";}
	    
	    ### DIFF OUTPUT ADDITION TO QUAL SCORE!  --TSW
           ### Essentially adding the expected number of reads re-distributed to INC or EXC after normalization..
           ### These values are added to the qual score and used to infer the posterior distribution
           unless($legacyFlag) {
             my $totalN = $total_raw_reads_ALL;
             my($pPSI, $exValOfInc, $exValOfExc) = (0, 0, 0);
             unless($PSI[$i] eq "NA" or $totalN == 0) {
		 $pPSI = $PSI[$i] / 100;
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
