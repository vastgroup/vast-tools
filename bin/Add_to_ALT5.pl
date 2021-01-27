#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Cwd;
use Getopt::Long;

my $dbDir;
my $sp;
my $samLen;
my $verboseFlag;
my $extra_eej = 10; # hardcoded for now

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "len=i" => \$samLen,
			  "verbose=i" => \$verboseFlag);

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine alt5]: $verbMsg\n";
  }
}

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

    next if $N_ss > 15; # change 1 for speed

    $ALL{$event_root}=$N_ss; # keeps the total number of alternative splice sites; placeholder for event_root

    ### gets the strand of gene
    ($gene)=$event=~/(.+?)\-/;
    ($C2)=$t[4]=~/.+\,(\d+)/;
    ($core)=$t[4]=~/\:[^\d]*?(\d+)/;
    die "Can't identify a C2 ($C2) or core ($core) coordinate for $event\n" if !$C2 || !$core;
    $strand{$gene}="+" if $C2 > $core;
    $strand{$gene}="-" if $C2 <= $core;
}
close TEMPLATE;

@EEJ1=glob("to_combine/*.eej2");
@EEJ2=glob("to_combine/*.eej2.gz");
@EEJ=(@EEJ1,@EEJ2);

@EFF=glob("$dbDir/FILES/$sp"."_COMBI-$COMB-*gDNA.eff");
die "Needs strand-unspecific effective\n" if !@EFF;
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

@EFF=glob("$dbDir/FILES/$sp"."_COMBI-$COMB-*gDNA-SS.eff");
die "Needs strand-specific effective (strand-specific)\n" if !@EFF;
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
verbPrint "Loading EEJ data for ALT5\n";
foreach $file (@EEJ){
    my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\.eej/;
    $length = $samLen;
    # generates headings
    $head_PSIs.="\t$sample\t$sample-Q";
    $head_ReadCounts.="\t$sample-Ri\t$sample-Rtot\t$sample-Q";
    
    unless(-e "to_combine/${sample}.info" || -e "to_combine/${sample}.info.gz"){ verbPrint "   $sample: do not find to_combine/${sample}.info. Sample will be treated as being not strand-specific.";
    }else{
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

    # Loads EEJ files
    if ($file=~/\.gz$/){
	open (EEJ, "gunzip -c $file | ") || die"It cannot open the $file\n";
    } else {
        open (EEJ, $file);
    }
    while (<EEJ>){
        chomp;
        @t=split(/\t/);
	$gene=$t[0];
	$eej=$t[1];
        $gene_eej="$gene-$t[1]";

	# to correct for stack reads per sample
	if ($is_ss{$sample}){$eff_href=\%eff_ss}else{$eff_href=\%eff_ns}
	$new_count=0;
	$risky_pos=0;
	@temp_vals=();
	@temp_POS=split(/\,/,$t[4]);
	foreach $t_var (@temp_POS){
	    ($pos,$pos_count)=$t_var=~/(.+?)\:(.+)/;
	    push(@temp_vals,$pos_count);
	    $risky_pos++ if ($pos == 0 || ($pos == 34 || $pos == $eff_href->{$length}{$event}-1));
	}
	$median=median(@temp_vals);
	$positive_pos=$#temp_vals+1;

	if ($eff_href->{$length}{$gene_eej} < 4){ # 3 or fewer positions
	    $new_count = $t[2];
	}
        elsif ($t[2] >= 2 && $positive_pos == 1 && $risky_pos == $positive_pos){ # i.e. 2 or more reads stack into the same first or last position
	    $new_count = 1 if $t[2]<=3;
            $new_count = 0 if $t[2]>3;
        }
	elsif ($t[2] >= 3 && $positive_pos == 1){ # i.e. 3 or more reads stack into 1 position
	    $new_count = 1 if $t[2]<=4;
            $new_count = 0 if $t[2]>4;
	}
	elsif ($t[2] >= 6 && $positive_pos == 2){ # i.e. 6 or more reads stack into 2 position
	    $new_count = 2 if $t[2]<=10;
            $new_count = 0 if $t[2]>10;
	}
	else {
	    foreach $temp_val (@temp_vals){
		if ($temp_val > $median*4){ # this median can never be 0 by definition.
		    $new_count+=$median*4;
		}
		else {
		    $new_count+=$temp_val;
		}
	    }
	}
	
        $reads{$sample}{$gene_eej}=$new_count;
        ($donor,$acceptor)=$eej=~/(\d+?)\-(\d+)/;
        $last_acceptor{$gene}=$acceptor if $last_acceptor{$gene}<$acceptor; # keeps track of the last acceptor used
	$last_donor{$gene}=$donor if $last_donor{$gene}<$donor; # keeps track of the last donor used
    }
    close EEJ;
}

# Output files
$NUM=$#EEJ+1;
open (PSIs, ">raw_incl/INCLUSION_LEVELS_ALT5-$sp$NUM-n.tab");
open (COUNTs, ">raw_reads/RAW_READS_ALT5-$sp$NUM-n.tab");
print PSIs "$head_PSIs\n";
print COUNTs "$head_ReadCounts\n";

my $eff_href;
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
	my $fname = $file;
	$fname =~ s/^.*\///;
	($sample)=$fname=~/^(.*)\.eej/;
	$length = $samLen;  #replacement --TSW
	# Emptying variables and arrays with read counts per sample
	$total_raw_reads_S=$total_corr_reads_S=0; # total simple reads
	$total_raw_reads_ALL=$total_corr_reads_ALL=0; # total complex and simple reads
	@corr_inc_reads_S=@raw_inc_reads_S=(); # Simple reads for each donor
	@corr_inc_reads_ALL=@raw_inc_reads_ALL=(); # Complex reads for each donor
	@PSI=(); # Percent splice site usage for each donor
	
	# set mappability correction accordingly
	if($is_ss{$sample}){$eff_href=\%eff_ss}else{$eff_href=\%eff_ns}
	
	$max_mappability=$length-15;
	for $i (0..$#junctions){ # does Simple and Complex read counts
	    $eejS="$gene-$junctions[$i]";
	    $corr_inc_reads_S[$i]=$max_mappability*($reads{$sample}{$eejS}/$eff_href->{$length}{$eejS}) if $eff_href->{$length}{$eejS};
	    $raw_inc_reads_S[$i]=$reads{$sample}{$eejS} if $eff_href->{$length}{$eejS};
	    $total_corr_reads_S+=$corr_inc_reads_S[$i];
	    $total_raw_reads_S+=$raw_inc_reads_S[$i];

	    ### Multi Q
	    ($donor)=$junctions[$i]=~/(\d+?)\-\d+/;
	    for $j (0..$last_acceptor{$gene}){  # It includes the "reference" acceptor
		$eej="$gene-$donor-$j"; 
		$corr_inc_reads_ALL[$i]+=$max_mappability*($reads{$sample}{$eej}/$eff_href->{$length}{$eej}) if $eff_href->{$length}{$eej};
		$raw_inc_reads_ALL[$i]+=$reads{$sample}{$eej} if $eff_href->{$length}{$eej};
	    }
	    $total_corr_reads_ALL+=$corr_inc_reads_ALL[$i];
	    $total_raw_reads_ALL+=$raw_inc_reads_ALL[$i];
	}
	# reads that jump over the event (any upstream donor to any downstream acceptor)
	$skipping_corr_reads=0; $skipping_raw_reads=0;

	($ref_acceptor)=$junctions[0]=~/\d+?\-(\d+)/;
	($ext_donor)=$junctions[$#junctions]=~/(\d+?)\-\d+/;
	($int_donor)=$junctions[0]=~/(\d+?)\-\d+/;

	$min_donor=$int_donor-($extra_eej*2);
	$min_donor=0 if $min_donor<0;
	$max_donor=$int_donor-1;
	$max_donor=0 if $max_donor<0;

	$min_acceptor=$ref_acceptor-$extra_eej;
	$min_acceptor=0 if $min_acceptor<0;
	$max_acceptor=$ref_acceptor+$extra_eej;
	$max_acceptor=$last_acceptor{$gene} if $max_acceptor > $last_acceptor{$gene};
	
#	for $t_don (0..$last_donor{$gene}){
#	    for $t_acc (0..$last_acceptor{$gene}){ # redundant call for ext/int, just in case
	for $t_don ($min_donor..$max_donor){ # change 3 to improve speed 
	    for $t_acc ($min_acceptor..$max_acceptor){ # change 3 to improve speed
		if ($strand{$gene} eq "+" && 
		    $co_donor{$gene}{$t_don} < $co_donor{$gene}{$ext_donor} && $co_donor{$gene}{$t_don} < $co_donor{$gene}{$int_donor} &&
		    $co_acceptor{$gene}{$t_acc} > $co_donor{$gene}{$ext_donor} && $co_acceptor{$gene}{$t_acc} > $co_donor{$gene}{$int_donor}){
		    
		    $eej="$gene-$t_don-$t_acc";
		    $skipping_corr_reads+=$max_mappability*($reads{$sample}{$eej}/$eff_href->{$length}{$eej}) if $eff_href->{$length}{$eej};
		    $skipping_raw_reads+=$reads{$sample}{$eej} if $eff_href->{$length}{$eej};    
		}
		elsif ($strand{$gene} eq "-" && 
		       $co_donor{$gene}{$t_don} > $co_donor{$gene}{$ext_donor} && $co_donor{$gene}{$t_don} > $co_donor{$gene}{$int_donor} &&
		       $co_acceptor{$gene}{$t_acc} < $co_donor{$gene}{$ext_donor} && $co_acceptor{$gene}{$t_acc} < $co_donor{$gene}{$int_donor}){
		    
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
	    $Q[$i] = $Q; # only 3 scores here
	    
	    ### Score 4: from v2.2.2 is specific for each splice site
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



sub median {
    my @temp=@_;
    
    @temp=sort{$a<=>$b}(@temp);
    $ind=$#temp+1; # number of positions with counts
    if ($ind%2 != 0){
        $m = int($ind/2);
        $median=$temp[$m];
    }
    else {
        $m1 = int($ind/2)-1;
        $m2 = int($ind/2);
        $median = ($temp[$m1]+$temp[$m2])/2;
    }
    return $median;
}
