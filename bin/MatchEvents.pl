#!/usr/bin/perl
### General script to get differentially spliced events based on dPSI differences 

use Getopt::Long;
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);
use Cwd qw(abs_path);
use List::Util qw(sum);

# INITIALIZE PATH AND FLAGS
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

### Setting global variables:
my $helpFlag = 0;
my $verboseFlag = 1; 
my $dbDir; # directory of VASTDB
my $species;
my $input_file = $ARGV[0];
my $output_file;
my $outdir;
my $input_format = "vast";
my $use_names;
my @samples;

GetOptions(
    "dbDir=s" => \$dbDir,
    "o=s" => \$outdir,
    "format=s" => \$input_format,
    "f=s" => \$input_format,
    "sp=s" => \$species, 
    "verbose=i" => \$verboseFlag,
    "use_names" => \$use_names
    );

our $EXIT_STATUS = 0;

sub errPrint {
    my $errMsg = shift;
    print STDERR "[vast match]: $errMsg\n";
    $EXIT_STATUS++; 
}

sub errPrintDie {
    my $errMsg = shift;
    errPrint $errMsg;
    exit $EXIT_STATUS if ($EXIT_STATUS != 0);
}

sub verbPrint {
    my $verbMsg = shift;
    if($verboseFlag) {
	chomp($verbMsg);
	print STDERR "[vast match]: $verbMsg\n";
    }
}

### Option/help
if (!defined $ARGV[0] || !defined $species || !defined $outdir || $helpFlag){
    die "\nUsage: vast-tools match list_of_events -sp Hsa/Mmu/etc -o vast_folder [options]

Matches events from other tools to vast-tools internal EEJs and quantifies simple and complex PSIs (when possible)

OPTIONS: 
        -o, --outDir                Path to output folder of vast-tools align (default vast_out)
        --sp Hsa/Mmu/etc            Three letter code for the database
        --dbDir db                  Database directory (default VASTDB)
        --format/-f vast/suppa/SJ   Format of the events (default 'vast')
        --use_names                 Information based on gene symbol, not EnsemblIDs (default OFF)
        --help                      Prints this help message

FORMATS:
        - vast:   GENE_ID\\tchr:C1donor,Ai-Af,C2acceptor\\tAltEx
                  GENE_ID\\tchr:C1i-C1f=C2i-C2f:strand\\tIR        (not yet available)
                  GENE_ID\\tchr:C1d,A1+A2-X\\tAlt3                 (not yet available)
                  GENE_ID\\tchr:X-A1+A2,C2a\\tAlt5                 (not yet available)

        - suppa:  GENE_ID:SE:chr:C1d-Aa:Ad-C2a:+
                  GENE_ID:SE:chr:C2a-Ad:Aa-C1d:-
                  GENE_ID:A3:chr:C1d-A1:C1d-A2:+                 (not yet available)
                  GENE_ID:A3:chr:A1-C1d:A2-C1d:-                 (not yet available)
                  GENE_ID:A5:chr:A2-C2a:A1-C2a:+                 (not yet available)
                  GENE_ID:A5:chr:C2a-A2:C2a-A1:-                 (not yet available)

          - SJ:   chr_C1d_C2a_+\\t*\\tAi\\tAf\\t*\\tGENE_NAME
                  chr_C2a_C1d_-\\t*\\tAi\\tAf\\t*\\tGENE_NAME


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

# change directories
my $input_file_fullpath=abs_path($input_file);

errPrintDie "The output directory \"$outdir/to_combine\" does not exist" unless (-e "$outdir/to_combine");
chdir($outdir) or errPrint "Unable to change directories into output" and die;
verbPrint "Setting output directory to $outdir";

unless(defined($dbDir)) {
    $dbDir = "$binPath/../VASTDB/$species";
}
$dbDir = abs_path($dbDir);


###### Data collection
my %names_ID;
my %ID_names;
open (NAMES, "$dbDir/FILES/$species.ID.names.txt") || errPrintDie "Cannot find ID names";
while (<NAMES>){
    chomp ($_);
    my @t=split(/\t/);
    $ID_names{$t[0]}=$t[1];
    $names_ID{$t[1]}=$t[0];
}
close NAMES;

my @EEJ=glob("to_combine/*.eej*");
my @EFF=glob("$dbDir/FILES/$species*-M-*-gDNA.ef*");
errPrintDie "Needs effective from database!" if !@EFF;

###
verbPrint "Loading Mappability for each EEJ and length:\n";
my %eff;
my %D_CO;
my %A_CO;
my %CO_A;
my %CO_D;
my %last_acceptor;
foreach my $file (@EFF){
    my ($length)=$file=~/COMBI\-[A-Z]\-(\d+?)\-/;
    verbPrint "Loading: $file\tLength: $length\n";
    open (MAPPABILITY, $file);
    while (<MAPPABILITY>){
	chomp;
	my @t=split(/\t/);
	my ($gene,$donor,$acceptor,$donor_coord,$acceptor_coord)=$t[0]=~/(.+?)\-(\d+?)\_(\d+?)\-(\d+?)\_(\d+)/;
	my $eej="$gene-$donor-$acceptor";
	$eff{$length}{$eej}=$t[1];

	# keeps the coordinate for each donor/acceptor
	$D_CO{$gene}{$donor}=$donor_coord;
        $A_CO{$gene}{$acceptor}=$acceptor_coord;
	$CO_D{$gene}{$donor_coord}=$donor;
        $CO_A{$gene}{$acceptor_coord}=$acceptor;
	$last_acceptor{$gene} = 0 if !defined $last_acceptor{$gene};
        $last_acceptor{$gene}=$acceptor if $last_acceptor{$gene} < $acceptor; # keeps track of the last used acceptor
    }
    close MAPPABILITY;
}

###
verbPrint "Loading EEJ read counts data\n";
my $head_PSIs;
my %reads;
foreach my $file (@EEJ){
    my $fname = $file;
    $fname =~ s/^.*\///;
    my ($sample)=$fname=~/^(.*)\..*$/;
    verbPrint "  $sample\n";
    push(@samples,$sample);

    $head_PSIs.="\t$sample\t$sample-Q";

    open (EEJ, $file);
    while (<EEJ>){
	chomp;
	my @t=split(/\t/);
	my $gene=$t[0];
	my $eej=$t[1];
	my $event="$gene-$eej";

        $reads{$sample}{$event}=$t[2];
        my ($donor,$acceptor)=$eej=~/(\d+?)\-(\d+)/;
	$last_acceptor{$gene} = 0 if !defined $last_acceptor{$gene};
        $last_acceptor{$gene}=$acceptor if $last_acceptor{$gene} < $acceptor; # keeps track of the last used acceptor
    }
    close EEJ;
}

### to match vast IDs
# Load conversation table file to memory
my $NEWID = openFileHandle("$dbDir/FILES/New_ID-$species.txt.gz");
my %newIDs;
while (<$NEWID>) {
    chomp;
    my @l = split("\t");
    if (defined $newIDs{$l[1]}) {
	die "Non-unique key value pair in $NEWID!\n";
    }
    $newIDs{$l[1]} = $l[0];
  
  # to correct a small discordance in old human/mouse IDs --MI [23/12/15]
    if ($l[0] =~ /INT/ && $l[1] =~ /^\-/){
	my $temp_ID = "NA".$l[1];
	$newIDs{$temp_ID} = $l[0];
    }
}
close $NEWID;

# Loads ID => GeneID
my %ID_gene;
open (G_ID, "$dbDir/FILES/$species.Event-Gene.IDs.txt") || errPrintDie "Cannot open Event to GeneID key";
while (<G_ID>){
    chomp($_);
    my @t=split(/\t/,$_);
    $ID_gene{$t[0]} = $t[1];
}
close G_ID;

# Loads the template 
my $TEMPLATE = openFileHandle("$dbDir/TEMPLATES/$species.FULL.Template.txt.gz");
my %vast_events_A;
my %vast_events_D;
<$TEMPLATE>;
while (<$TEMPLATE>){
    chomp;
    my @t = split(/\t/);
    my $event;

    my $type;
    if ($t[5]=~/IR/){
	$type = "IR";
	if ($t[1] =~ /^\-/){
	    my $temp_ID = "NA".$t[1];
	}
    }
    elsif ($t[5]=~/Alt/){
	$type = $t[5];
    }
    else {
	$type = "AltEx";
	my ($C1)=$t[4]=~/\:(\d+)/;
	my ($C2)=$t[4]=~/\,.+?\,(\d+)/;
	my ($I,$F)=$t[4]=~/\,(.+?)\-(.+?)\,/;
	my @I = split(/\+/,$I);
	my @F = split(/\+/,$F);
	my $strand = "+" if $C2>=$C1;
	$strand = "-" if $C2<$C1;
	$event = $newIDs{$t[1]};
	my $gene = $ID_gene{$event};
	
	foreach my $i (@I){
	    my $t_co = "$gene-$i";
	    $vast_events_A{$type}{$t_co} = $event if $strand eq "+";
	    $vast_events_D{$type}{$t_co} = $event if $strand eq "-";
	}
	foreach my $f (@F){
	    my $t_co = "$gene-$f";
	    $vast_events_A{$type}{$t_co} = $event if $strand eq "-";
	    $vast_events_D{$type}{$t_co} = $event if $strand eq "+";
	}
    }
}
close $TEMPLATE;


#### Starts the event matching
verbPrint "Matching list of events and quantifying PSIs\n";

my ($root_file) = $input_file_fullpath =~ /.+\/(.+)/;
open (OUTs, ">$root_file-PSI_S.tab") or die "Cannot create file $root_file-PSI_S.tab for writing\n";
open (OUTc, ">$root_file-PSI_C.tab") or die "Cannot create file $root_file-PSI_C.tab for writing\n";

my $line = 0;
open (EVENTS, $input_file_fullpath) || errPrintDie "Can't open event list";

my $head = <EVENTS>;
chomp($head);
print OUTs "$head$head_PSIs\tVastID\n";
print OUTc "$head$head_PSIs\tVastID\n";

# get a string with as many NAs as head_PSIs has elements; this string will be output for all not-mappable exons
my $NAs="";
for(my $i=0;$i<@{[split("\t",$head_PSIs)]};$i++){$NAs.="\tNA";}


my %tally_matches;
while (<EVENTS>){
    s/\r//g;
    s/\"//g;
    chomp;
    my @t=split(/\t/);    
    my $gene;
    my $type_event;
    my $strand;
    my $event_coord;
    my $vastID;

    my $C1;
    my $C2;
    my $Aa;
    my $Ad;
    $line++;

    ### deals with different formats
    if (!defined $input_format || $input_format eq "vast"){
	$gene = $t[0];
	$event_coord = $t[1];
	$type_event = $t[2];
	if (defined $use_names){
	    $gene = $names_ID{$t[0]};
	}
    }
    elsif ($input_format eq "SUPPA" || $input_format eq "suppa"){
	my ($t1,$t2,$t3,$t4) = $t[0]=~/(.+?)\;(.+?)\:(.+)\:([\+\-])/;
	$gene = $t1;
	$event_coord = $t3;
	$strand = $t4;

	if (defined $use_names){
	    $gene = $names_ID{$t1};
	}
	$input_format = "suppa";

	# define event type
	$type_event = "AltEx" if $t2 eq "SE";
	$type_event = "IR" if $t2 eq "RI";
	$type_event = "Alt3" if $t2 eq "A3";
	$type_event = "Alt5" if $t2 eq "A5";
	$type_event = "Other" if !defined $type_event;
    }
    elsif ($input_format eq "SJ" || $input_format eq "SANJUAN" || $input_format eq "sanjuan"){
	my ($t1,$t2,$t3,$t4) = $t[0]=~/(.+?)\_(.+?)\_(.+?)\_([\+\-])/;
	$t3++;
	$t[3]++;
	$event_coord = "$t1:$t2-$t[2]:$t[3]-$t3";
	$strand = $t4;

	### sometime it contains multiple names
	my @temp_gene = split(/\,/,$t[5]);
	$gene = $names_ID{$temp_gene[0]};
    
	$input_format = "SJ";
	$type_event = "AltEx";
    }

    ### does the match by event type
    if ($type_event eq "AltEx"){
	my $chr;
	my $C1_co; my $C2_co;
	my $i_co; my $f_co;
	if (!defined $input_format || $input_format eq "vast"){
	    ($chr,$C1_co,$i_co,$f_co,$C2_co) = $event_coord =~ /(.+?)\:(.+?)\,(.+?)\-(.+?)\,(.+)/;
	    $strand = "+" if $C2_co >= $C1_co;
	    $strand = "-" if $C2_co < $C1_co;
	}
	elsif ($input_format eq "suppa" || $input_format eq "SJ"){
	    if ($strand eq "+"){
		($chr,$C1_co,$i_co,$f_co,$C2_co) = $event_coord=~/(.+?)\:(.+?)\-(.+?)\:(.+?)\-(.+)/;
	    }
	    else {
		($chr,$C2_co,$i_co,$f_co,$C1_co) = $event_coord=~/(.+?)\:(.+?)\-(.+?)\:(.+?)\-(.+)/;
	    }
	}

	# prints the root
	print OUTs $_;
	print OUTc $_;

	### match to vastID
	$gene = "NA" if !defined $gene;
	my $t_coA = "$gene-$i_co";
	my $t_coB = "$gene-$f_co";
	
	if ($strand eq "+"){
	    $vastID = $vast_events_A{AltEx}{$t_coA} if defined $vast_events_A{AltEx}{$t_coA};
	    $vastID = $vast_events_D{AltEx}{$t_coB} if defined $vast_events_D{AltEx}{$t_coB} && !defined $vast_events_A{AltEx}{$t_coA};
	}
	else {
	    $vastID = $vast_events_D{AltEx}{$t_coA} if defined $vast_events_D{AltEx}{$t_coA};
	    $vastID = $vast_events_A{AltEx}{$t_coB} if defined $vast_events_A{AltEx}{$t_coB} && !defined $vast_events_D{AltEx}{$t_coA};
	}
	$vastID = "NA" if !defined $vastID;

	### matching coordinates to junctions
	$C1 = $CO_D{$gene}{$C1_co};
	$C2 = $CO_A{$gene}{$C2_co};
	if ($strand eq "+"){
	    $Aa = $CO_A{$gene}{$i_co};
	    $Ad = $CO_D{$gene}{$f_co};
	}
	else {
	    $Aa = $CO_A{$gene}{$f_co};
	    $Ad = $CO_D{$gene}{$i_co};
	}
	
	if (!defined $C1 || !defined $C2 || !defined $Aa || !defined $Ad) {
#	    errPrint "Cannot map all coordinates for event in gene $gene \[line $line\]";
	    print OUTs "$NAs\n";
	    print OUTc "$NAs\n";
	    $tally_matches{AltEx}{NO}++;
	}
	else {
	    ### making three ref junctions
	    my $inc1 = "$gene-$C1-$Aa";
	    my $inc2 = "$gene-$Ad-$C2";
	    my $exc = "$gene-$C1-$C2";
	    
	    #i.e. any of the basic have no mappability
	    if (!$eff{50}{$inc1} || !$eff{50}{$inc2} || !$eff{50}{$exc}){
		print OUTs "$NAs\n";
		print OUTc "$NAs\n";
		$tally_matches{AltEx}{NO}++;
	    }
	    next if !$eff{50}{$inc1} || !$eff{50}{$inc2} || !$eff{50}{$exc};

	    $tally_matches{AltEx}{YES}++; # no mappablity already discounted
	    
	    foreach my $sample (@samples){
		my $reads_inc1 = $reads{$sample}{$inc1};
		my $reads_inc2 = $reads{$sample}{$inc2};
		my $reads_exc = $reads{$sample}{$exc};
		$reads_inc1 = 0 if !defined $reads_inc1;
		$reads_inc2 = 0 if !defined $reads_inc2;
		$reads_exc = 0 if !defined $reads_exc;
		my $total_reads = $reads_exc + $reads_inc1 + $reads_inc2;
		
		my $Creads_inc1 = sprintf("%.2f", 35 * $reads_inc1/$eff{50}{$inc1}) if (defined $eff{50}{$inc1} && $eff{50}{$inc1}>0);
		my $Creads_inc2 = sprintf("%.2f", 35 * $reads_inc2/$eff{50}{$inc2}) if (defined $eff{50}{$inc2} && $eff{50}{$inc2}>0);
		my $Creads_exc = sprintf("%.2f", 35 * $reads_exc/$eff{50}{$exc}) if (defined $eff{50}{$exc} && $eff{50}{$exc}>0);
		$Creads_inc1 = 0 if !defined $Creads_inc1;
		$Creads_inc2 = 0 if !defined $Creads_inc2;
		$Creads_exc = 0 if !defined $Creads_exc;
		my $total_Creads = $Creads_exc + $Creads_inc1 + $Creads_inc2;
		
		my $PSI_simple;
		if (($Creads_inc1 + $Creads_inc2 + 2*$Creads_exc)>0){
		    $PSI_simple = sprintf ("%.2f", 100 * ($Creads_inc1 + $Creads_inc2)/ ($Creads_inc1 + $Creads_inc2 + 2*$Creads_exc));
		}
		else {
		    $PSI_simple = "NA";
		}
		
		my $Q_simple;
		my $inc_average_simple;

		my $PSI_complex;
		my $Q_complex;
		my $inc_average_complex;

		### Complex PSI
		my $inc1C = 0; my $Rinc1C = 0;
		my $inc2C = 0; my $Rinc2C = 0;
		my $excC = 0; my $RexcC = 0;

		#### Inclusion
                for my $i (0..$Ad-1){ # loops through all donors in the gene from the first to d2-1
                    if ((($D_CO{$gene}{$i} < $i_co && $strand eq "+") || ($D_CO{$gene}{$i} > $f_co && $strand eq "-")) && $i != $C1){
                        my $temp_eej="$gene-$i-$Aa";
                        $reads{$sample}{$temp_eej} = 0 if !defined $reads{$sample}{$temp_eej};
			
                        if (defined $eff{50}{$temp_eej} && $eff{50}{$temp_eej}>0){
                            $inc1C+=$reads{$sample}{$temp_eej}/$eff{50}{$temp_eej};
                            $Rinc1C+=$reads{$sample}{$temp_eej};
                        }
                    }
                }
                for my $i ($Aa+1..$last_acceptor{$gene}){ # loops through all acceptors in the gene from a1+1 to the last
                    if ((($A_CO{$gene}{$i} > $f_co && $strand eq "+") || ($A_CO{$gene}{$i} < $i_co && $strand eq "-")) && $i != $C2){
                        my $temp_eej="$gene-$Ad-$i";
                        $reads{$sample}{$temp_eej} = 0 if !defined $reads{$sample}{$temp_eej};
                        if (defined $eff{50}{$temp_eej} && $eff{50}{$temp_eej}>0){
                            $inc2C+=$reads{$sample}{$temp_eej}/$eff{50}{$temp_eej};
                            $Rinc2C+=$reads{$sample}{$temp_eej};
                        }
                    }
                }

		### Exclusion reads (It does NOT take all EEJs around the alternative exon, but only those including C1 or C2.)
		for my $i (0..$C1-1){ # $d1 => $C1, $a1 => $Aa, $d2 => $Ad, $a2 => $C2 
		    my $temp_eej="$gene-$i-$C2";
		    if ($eff{50}{$temp_eej} >= 1){
			$reads{$sample}{$temp_eej} = 0 if (!defined $reads{$sample}{$temp_eej});
			$excC+=$reads{$sample}{$temp_eej}/$eff{50}{$temp_eej};
			$RexcC+=$reads{$sample}{$temp_eej};
		    }
		}
		for my $i ($C1+1..$Ad-1){
		    if (($D_CO{$gene}{$i} < $A_CO{$gene}{$Aa} && $strand eq "+") || ($D_CO{$gene}{$i} > $A_CO{$gene}{$Aa} && $strand eq "-")){
			my $temp_eej="$gene-$i-$C2";
			if ($eff{50}{$temp_eej} >= 1){
			    $reads{$sample}{$temp_eej} = 0 if (!defined $reads{$sample}{$temp_eej});
			    $excC+=$reads{$sample}{$temp_eej}/$eff{50}{$temp_eej};
			    $RexcC+=$reads{$sample}{$temp_eej};
			}
		    }
		}
		for my $i ($C2+1..$last_acceptor{$gene}){
		    my $temp_eej="$gene-$C1-$i";
		    if ($eff{50}{$temp_eej} >= 1){
			$reads{$sample}{$temp_eej} = 0 if (!defined $reads{$sample}{$temp_eej});
			$excC+=$reads{$sample}{$temp_eej}/$eff{50}{$temp_eej};
			$RexcC+=$reads{$sample}{$temp_eej};
		    }
		}
		for my $i ($Ad+1..$C2-1){
		    if (($A_CO{$gene}{$i} > $D_CO{$gene}{$Ad} && $strand eq "+") || ($A_CO{$gene}{$i} < $D_CO{$gene}{$Ad} && $strand eq "-")){
			my $temp_eej="$gene-$C1-$i";
			if ($eff{50}{$temp_eej} >= 1){
			    $reads{$sample}{$temp_eej} = 0 if (!defined $reads{$sample}{$temp_eej});
			    $excC+=$reads{$sample}{$temp_eej}/$eff{50}{$temp_eej};
			    $RexcC+=$reads{$sample}{$temp_eej};
			}
		    }
		}
#####
	    ### exclusion (21/06/17) --MI
#	    for my $i (0..$Ad-1){
#		for my $j ($Aa+1..$last_acceptor{$gene}){
#		    if (($D_CO{$gene}{$i} < $i_co && $A_CO{$gene}{$j} > $f_co && $strand eq "+") || 
#			($D_CO{$gene}{$i} > $f_co && $A_CO{$gene}{$j} < $i_co && $strand eq "-")) {
#			my $temp_eej="$gene-$i-$j";
#			if (defined $eff{50}{$temp_eej} && $eff{50}{$temp_eej} > 0  && ($i != $C1 || $j != $C2)){
#				$reads{$sample}{$temp_eej} = 0 if (!defined $reads{$sample}{$temp_eej});
#				$excC+=$reads{$sample}{$temp_eej}/$eff{50}{$temp_eej};
#				$RexcC+=$reads{$sample}{$temp_eej};
#			}
#		    }
#		}
#	    }
		
		### sum of reads
		my $reads_inc1_all = $reads_inc1 + $Rinc1C;
		my $reads_inc2_all = $reads_inc2 + $Rinc2C;
		my $reads_exc_all = $reads_exc + $RexcC;
		my $total_reads_all = $reads_inc1_all + $reads_inc2_all + $reads_exc_all;
		my $Creads_inc1_all = sprintf ("%.2f",$Creads_inc1 + $inc1C);
		my $Creads_inc2_all = sprintf ("%.2f",$Creads_inc2 + $inc2C);
		my $Creads_exc_all = sprintf ("%.2f",$Creads_exc + $excC);
		my $total_Creads_all = $Creads_inc1_all + $Creads_inc2_all + $Creads_exc_all;
		
		
#### Coverage scores. Score 1: Using all raw reads
		if (($reads_exc >= 20 || (($reads_inc1 >= 20 && $reads_inc2 >= 15) || ($reads_inc1 >= 15 && $reads_inc2 >= 20))) && $total_reads>=100){
		    $Q_simple="SOK";
		}
		elsif (($reads_exc >= 20 || (($reads_inc1 >= 20 && $reads_inc2 >= 15) || ($reads_inc1 >= 15 && $reads_inc2 >= 20))) && $total_reads<100){
		    $Q_simple="OK";
		}
		elsif ($reads_exc >= 15 || (($reads_inc1 >= 15 && $reads_inc2 >= 10) || ($reads_inc1 >= 10 && $reads_inc2 >= 15))){
		    $Q_simple="LOW";
		}
		elsif ($reads_exc >= 10 || (($reads_inc1 >= 10 && $reads_inc2 >= 5) || ($reads_inc1 >= 5 && $reads_inc2 >= 10))){
		    $Q_simple="VLOW";
		}
		else {
		    $Q_simple="N";
		}
		### COMPLEX
		if (($reads_exc_all >= 20 || (($reads_inc1_all >= 20 && $reads_inc2_all >= 15) || ($reads_inc1_all >= 15 && $reads_inc2_all >= 20))) && $total_reads_all>=100){
		    $Q_complex="SOK";
		}
		elsif (($reads_exc_all >= 20 || (($reads_inc1_all >= 20 && $reads_inc2_all >= 15) || ($reads_inc1_all >= 15 && $reads_inc2_all >= 20))) && $total_reads_all<100){
		    $Q_complex="OK";
		}
		elsif ($reads_exc_all >= 15 || (($reads_inc1_all >= 15 && $reads_inc2_all >= 10) || ($reads_inc1_all >= 10 && $reads_inc2_all >= 15))){
		    $Q_complex="LOW";
		}
		elsif ($reads_exc_all >= 10 || (($reads_inc1_all >= 10 && $reads_inc2_all >= 5) || ($reads_inc1_all >= 5 && $reads_inc2_all >= 10))){
		    $Q_complex="VLOW";
		}
		else {
		    $Q_complex="N";
		}
		
#### Score 2: Using all corrected reads
		if (($Creads_exc >= 20 || (($Creads_inc1 >= 20 && $Creads_inc2 >= 15) || ($Creads_inc1 >= 15 && $Creads_inc2 >= 20))) && $total_Creads>=100){
		    $Q_simple.=",SOK";
		}
		elsif (($Creads_exc >= 20 || (($Creads_inc1 >= 20 && $Creads_inc2 >= 15) || ($Creads_inc1 >= 15 && $Creads_inc2 >= 20))) && $total_Creads<100){
		    $Q_simple.=",OK";
		}
		elsif ($Creads_exc >= 15 || (($Creads_inc1 >= 15 && $Creads_inc2 >= 10) || ($Creads_inc1 >= 10 && $Creads_inc2 >= 15))){
		    $Q_simple.=",LOW";
		}
		elsif ($Creads_exc >= 10 || (($Creads_inc1 >= 10 && $Creads_inc2 >= 5) || ($Creads_inc1 >= 5 && $Creads_inc2 >= 10))){
		    $Q_simple.=",VLOW";
		}
		else {
		    $Q_simple.=",N";
		}
		### COMPLEX
		if (($Creads_exc_all >= 20 || (($Creads_inc1_all >= 20 && $Creads_inc2_all >= 15) || ($Creads_inc1_all >= 15 && $Creads_inc2_all >= 20))) && $total_Creads_all>=100){
		    $Q_complex.=",SOK";
		}
		elsif (($Creads_exc_all >= 20 || (($Creads_inc1_all >= 20 && $Creads_inc2_all >= 15) || ($Creads_inc1_all >= 15 && $Creads_inc2_all >= 20))) && $total_Creads_all<100){
		    $Q_complex.=",OK";
		}
		elsif ($Creads_exc_all >= 15 || (($Creads_inc1_all >= 15 && $Creads_inc2_all >= 10) || ($Creads_inc1_all >= 10 && $Creads_inc2_all >= 15))){
		    $Q_complex.=",LOW";
		}
		elsif ($Creads_exc_all >= 10 || (($Creads_inc1_all >= 10 && $Creads_inc2_all >= 5) || ($Creads_inc1_all >= 5 && $Creads_inc2_all >= 10))){
		    $Q_complex.=",VLOW";
		}
		else {
		    $Q_complex.=",N";
		}

##### score 3 is the same as 1 in SIMPLE
		if (($reads_exc >= 20 || (($reads_inc1 >= 20 && $reads_inc2 >= 15) || ($reads_inc1 >= 15 && $reads_inc2 >= 20))) && $total_reads>=100){
		    $Q_simple.=",SOK";
		    $Q_complex.=",SOK";
		}
		elsif (($reads_exc >= 20 || (($reads_inc1 >= 20 && $reads_inc2 >= 15) || ($reads_inc1 >= 15 && $reads_inc2 >= 20))) && $total_reads<100){
		    $Q_simple.=",OK";
		    $Q_complex.=",OK";
		}
		elsif ($reads_exc >= 15 || (($reads_inc1 >= 15 && $reads_inc2 >= 10) || ($reads_inc1 >= 10 && $reads_inc2 >= 15))){
		    $Q_simple.=",LOW";
		    $Q_complex.=",LOW";
		}
		elsif ($reads_exc >= 10 || (($reads_inc1 >= 10 && $reads_inc2 >= 5) || ($reads_inc1 >= 5 && $reads_inc2 >= 10))){
		    $Q_simple.=",VLOW";
		    $Q_complex.=",VLOW";
		}
		else {
		    $Q_simple.=",N";
		    $Q_complex.=",N";
		}
		
### Score 4: Calculate imbalance between inclusion EEJs (OK<B1<B2; Bl=not enough inclusion reads):
		if ((defined $Creads_inc1) && $Creads_inc1 > 0 && (defined $Creads_inc2) && $Creads_inc2 > 0){
		    if (($Creads_inc1/$Creads_inc2 > 2 && $Creads_inc1/$Creads_inc2 <= 5) || ($Creads_inc2/$Creads_inc1 > 2 && $Creads_inc2/$Creads_inc1 <= 5)){
			$Q_simple.=",B1";
			$Q_complex.=",B1";
		    }
		    elsif (($Creads_inc1/$Creads_inc2 > 5) || ($Creads_inc2/$Creads_inc1 > 5)){
			$Q_simple.=",B2";
			$Q_complex.=",B2";
		    }
		    else {
			$Q_simple.=",OK";
			$Q_complex.=",OK";
		    }
		}
		else {
		    if (!defined $Creads_inc1){
			$Q_simple.=",B2" if $Creads_inc2>=20;
			$Q_simple.=",Bl" if ($Creads_inc2<20 && $Creads_inc2>=15);
			$Q_simple.=",Bn" if $Creads_inc2<15 && $Creads_inc2>0;
			$Q_complex.=",B2" if $Creads_inc2>=20;
			$Q_complex.=",Bl" if ($Creads_inc2<20 && $Creads_inc2>=15);
v			$Q_complex.=",Bn" if $Creads_inc2<15 && $Creads_inc2>0;
		    }
		    if (!defined $Creads_inc2){
			$Q_simple.=",B2" if $Creads_inc1>=20;
			$Q_simple.=",Bl" if ($Creads_inc1<20 && $Creads_inc1>=15);
			$Q_simple.=",Bn" if $Creads_inc1<15;
			$Q_complex.=",B2" if $Creads_inc1>=20;
			$Q_complex.=",Bl" if ($Creads_inc1<20 && $Creads_inc1>=15);
			$Q_complex.=",Bn" if $Creads_inc1<15;
		    }
		}


### Score 5: Complexity score (S<C1<C2<C3) [only for complex]
		my $from_C=$excC+$inc1C+$inc2C;
		my $from_S=$Creads_exc+$Creads_inc1+$Creads_inc2;
		
		if ($from_C > ($from_C+$from_S)/2) {$Q_complex.=",C3";}
		elsif ($from_C > ($from_C+$from_S)/5 && $from_C <= ($from_C+$from_S)/2){$Q_complex.=",C2";}
		elsif ($from_C > ($from_C+$from_S)/20 && $from_C <= ($from_C+$from_S)/5){$Q_complex.=",C1";}
		else {$Q_complex.=",S";}


		### PSI and N of reads
		my($pPSI_simple, $exValOfInc_simple, $exValOfExc_simple) = (0, 0, 0);
		unless($PSI_simple eq "NA"){
		    $pPSI_simple = $PSI_simple / 100;
		    $exValOfInc_simple = sprintf("%.2f", $pPSI_simple * $total_reads);
		    $exValOfExc_simple = sprintf("%.2f", (1-$pPSI_simple) * $total_reads);
		}
#		$inc_average_simple=sprintf("%.2f",($Creads_inc1+$Creads_inc2)/2);
#		$Q_simple.=",NA\@$inc_average_simple,$Creads_exc";
		$Q_simple.=",NA\@$exValOfInc_simple,$exValOfExc_simple";

		### for complex
		if (($Creads_inc1_all+$Creads_inc2_all+2*$Creads_exc_all)>0){
		    $PSI_complex=sprintf("%.2f",100*($Creads_inc1_all+$Creads_inc2_all)/($Creads_inc1_all+$Creads_inc2_all+2*$Creads_exc_all));
		}
		else {
		    $PSI_complex = "NA";
		}

		my($pPSI_complex, $exValOfInc_complex, $exValOfExc_complex) = (0, 0, 0);
		unless($PSI_complex eq "NA"){
		    $pPSI_complex = $PSI_complex / 100;
		    $exValOfInc_complex = sprintf("%.2f", $pPSI_complex * $total_reads_all);
		    $exValOfExc_complex = sprintf("%.2f", (1-$pPSI_complex) * $total_reads_all);
		}
		$Q_complex.=",NA\@$exValOfInc_complex,$exValOfExc_complex";
#		$inc_average_complex=sprintf("%.2f",($Creads_inc1_all+$Creads_inc2_all)/2);
#		$Q_complex.="\@$inc_average_complex,$Creads_exc_all";
		
		print OUTs "\t$PSI_simple\t$Q_simple";		
		print OUTc "\t$PSI_complex\t$Q_complex";
	    }
	    print OUTs "\t$vastID\n";
	    print OUTc "\t$vastID\n";
	}
    }
    else {
#	print OUTs "\tNot_yet_available\n";
#	print OUTc "\tNot_yet_available\n";
    }
}
verbPrint "Summary output:\n";
print "\nTYPE\tTOTAL\tMAPPED\tNON_MAPPED\t\%MAPPED\n";
foreach my $type (sort keys %tally_matches){
    $tally_matches{$type}{YES} = 0 if !defined $tally_matches{$type}{YES};
    $tally_matches{$type}{NO} = 0 if !defined $tally_matches{$type}{NO};
    my $sum = $tally_matches{$type}{YES} + $tally_matches{$type}{NO};
    my $perc = sprintf("%.2f",100*$tally_matches{$type}{YES}/$sum);
    print "$type\t$sum\t$tally_matches{$type}{YES}\t$tally_matches{$type}{NO}\t$perc\n";
}
print "\n";
