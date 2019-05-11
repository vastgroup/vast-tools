#!/usr/bin/perl
# Script to prepare and filter vast-tools PSI tables for other analyses

use Getopt::Long;
use Cwd qw(abs_path);

### Setting global variables:
$Q="O[KW]\,.+?\,.+?\,.+?\,.+?\@"; # NEW quality search

my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

$input_file=$ARGV[0];
$min_SD=5; # min standard deviation of the event (def=5)
#$min_Fraction=0.8; # min fraction of samples with good coverage (def=0.8)
#$min_N=10; # min number of samples with good coverage (def=10)
$noVLOW=0;
$p_IR=0;
$print_all = "";
$samples = "";
$group_file;
$noB3;
$min_ALT_use = 25;
$verboseFlag = 1;  # on for debugging 
$command_line=join(" ", @ARGV);
@TYPES=("MIC","AltEx","IR","Alt3","Alt5");

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "min_SD=i" => \$min_SD,
                          "min_N=i" => \$min_N,
                          "min_Fr=f" => \$min_Fraction,
                          "outFile=s"    => \$output_file,
                          "help" => \$helpFlag,
                          "p_IR" => \$p_IR,
                          "noVLOW" => \$noVLOW,
			  "noB3" => \$noB3,
			  "min_ALT_use=i" => \$min_ALT_use,
			  "samples=s" => \$samples,
			  "groups=s" => \$group_file,
			  "log" => \$log,
			  "onlyEX" => \$onlyEXSK,
			  "verbose" => \$verboseFlag,
			  "add_names" => \$AddName
    );

sub errPrint {
    my $errMsg = shift;
    print STDERR "[vast tidy error]: $errMsg\n";
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
	print STDERR "[vast tidy]: $verbMsg\n";
    }
}

sub time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
    $year += 1900;
    $mon += 1;
    my $datetime = sprintf "%04d-%02d-%02d (%02d:%02d)", $year, $mday, $mon, $hour, $min;
    return $datetime;
}

# gets the array of samples
if ($samples){
    @samples=split(/\,/,$samples);
    $N_samples=$#samples+1;
}

# sanity checks for min_Fraction
errPrintDie "min_Fr has to be between 0 and 1\n" if $min_Fraction > 1;
errPrintDie "min_N ($min_N) cannot be higher than N of samples ($N_samples)\n" if $min_N > $N_samples && $samples;

# defining default output file name
($root)=$input_file=~/(.+)\./;
$root_out=$root;
$root_out.="-minN_$min_N" if $min_N;
$root_out.="-minFr_$min_Fraction" if $min_Fraction;
$root_out.="-minSD_$min_SD";
$root_out.="-noVLOW" if $noVLOW;
$root_out.="-noB3" if $noB3;
$root_out.="-p_IR" if $p_IR;
$root_out.="-min_ALT_use$min_ALT_use";
$root_out.="-onlyEX" if $onlyEXSK;
$root_out.="-samples$N_samples" if $samples;
$root_out.="-groups" if $groups;

$output_file="$root_out-Tidy.tab" unless $output_file;
$log_file="$root_out-Tidy.log" if $log;

### Gets the version
my $version;
open (VERSION, "$binPath/../VERSION");
$version=<VERSION>;
chomp($version);
$version="No version found" if !$version;

if (!$ARGV[0] || $helpFlag){
    die "
VAST-TOOLS v$version

Usage: vast-tools tidy INCLUSION_LEVELS_FULL-SpN.tab \(--min_N min_N_samples OR --min_Fr min_fraction\) --min_SD min_SD [options]

Prepares and filters a vast-tools output for general analyses.

[General options] 
        -min_N i                Minimum number of samples with good coverage. Alternatively, you can define it by fraction
        -min_Fr f               Minimum fraction of samples with good coverage. 
        -min_SD i               Minimum standard deviation of the event (def=5)
        -samples S1,S2,...      Samples to be considered (default all)
        -groups FILE            Provide a config file to set two groups. (default OFF)
                                   The number/fraction of minimal samples will be applied to EACH group.
                                   The output table will have samples ordered by group.
                                   Format: SAMPLE_NAME\\tGROUP_IDENTIFIER
        -outFile                Output file name (default based on option parameters)
        --noVLOW                Do not use samples with VLOW coverage (default OFF)
        --noB3                   Does not use AltEx events with B3 imbalance (default OFF)
        --p_IR                  Filter IR by the p-value of the binomial test (default OFF)
        --min_ALT_use i          Minimum inclusion of the exon in which the Alt3/Alt5 is located across all 
                                   compared samples (default 25) (combine >= v2.2.1)
        --onlyEX                Outputs only EXSK events (default OFF)
        --add_names             Adds gene name to the event_ID. E.g. Mta1\=MmuEX0029874 (default OFF)
        --log                   Print the summary stats into a file (default OFF)


*** Questions \& Bug Reports: Manuel Irimia \(mirimia\@gmail.com\)

";
}

errPrintDie "*** You can only define a minimum fraction or absolute number of samples with good coverage\n" if $min_N && $min_Fraction;
errPrintDie "*** You need to define either a minimum fraction or absolute number of samples with good coverage\n" if !$min_N && !$min_Fraction;

# prints version (05/05/19)
verbPrint "VAST-TOOLS v$version";

### Creates the LOG
($folder)=$root=~/(.+)\//;
open (LOG, ">>$folder/VTS_LOG_commands.txt");
my $all_args="-min_SD $min_SD -min_ALT_use $min_ALT_use";
$all_args.=" -min_N $min_N" if defined $min_N;
$all_args.=" -min_Fr $min_Fraction" if defined $min_Fraction;
$all_args.=" -p_IR" if $p_IR;
$all_args.=" -noVLOW" if $noVLOW;
$all_args.=" -noB3" if $noB3;
$all_args.=" -samples $samples" if $samples;
$all_args.=" -groups $group_file" if $group_file;
$all_args.=" -log" if $log;
$all_args.=" -onlyEX" if $onlyEXSK;
$all_args.=" -add_names" if $AddName;
$all_args.=" -outFile $output_file" if $output_file;

print LOG "[VAST-TOOLS v$version, ".&time."] vast-tools tidy $all_args\n";


open (I, $input_file) || die "Can't open input file\n";
open (O, ">$output_file") || die "Can't open output file ($output_file)\n";

### saying hi:
verbPrint "Cleaning and filtering $ARGV[0]\n";

### Open group_file
if ($group_file){
    open (GROUPS, $group_file) || die "Cannot open group config file\n";
    while (<GROUPS>){
	chomp;
	@t=split(/\t/);
	$sample_group{$t[0]}=$t[1];
	push(@{$group_samples{$t[1]}},$t[0]);
    }
    close GROUPS;
}

### Printing a new clean heading (only Event_ID and PSI)
print O "EVENT";
$head=<I>;
chomp($head);
@H=split(/\t/,$head);
for $i (6..$#H){
    if ($i%2==0){
	if ($samples){
	    foreach $j (0..$#samples){
		if ($samples[$j] eq $H[$i]){
		    $valid_sample[$i] = 1;
		    print O "\t$H[$i]";
		    $TISSUES{$H[$i]}=1;
		    $max_N++;
		}
	    }
	}
	elsif ($group_file){
	    if ($sample_group{$H[$i]}){
		$sample_index{$H[$i]}=$i;
		$TISSUES{$H[$i]}=1;
	    }
	}
	else {
	    $valid_sample[$i] = 1;
	    print O "\t$H[$i]";
	    $TISSUES{$H[$i]}=1;
	    $max_N++;
	}
    }
}
if ($group_file){ # prints them in order by groups
    foreach $group (sort keys %group_samples){
	foreach $sample (@{$group_samples{$group}}){
	    print O "\t$sample";
	    die "Sample ($sample) from $group not found\n" if !$sample_index{$sample};
	    $max_N_groups{$group}++;
	}
	verbPrint "   Group $group: $max_N_groups{$group} samples\n";
    }
}
print O "\n";

### Processing events
while (<I>){
    s/VLOW/N/g if $noVLOW;
    chomp;
    @t=split(/\t/);

    $gene_name=$t[0];
    $event=$t[1];
    $length=$t[3];
    if ($t[5]=~/Alt/){
	$type=$t[5];
    }
    elsif ($t[5]=~/IR/){
	$type="IR";
    }
    else {
	if ($t[3]<=27){
	    $type="MIC";
	}
	else {
	    $type="AltEx";
	}
    }

    
    next if $length==0; # to remove the internal ss in Alt3 and Alt5
    next if $onlyEXSK && ($type=~/Alt[35]/ || $type=~/IR/);
    
    $event_ID="$gene_name=$event" if $AddName;
    $event_ID=$event if !$AddName;

    %PRINT=();
    next if $done{$event}; # to avoid possible repeated IDs (in preliminary tables)
    
    #### STANDARD TIDY
    if (!$group_file){
	@PSIs=();
	foreach $i (6..$#t){
	    if ($i%2==0){
		next if !$valid_sample[$i];
		if ($t[$i+1]=~/$Q/){
		    ### For IR
		    if ($type=~/IR/ && $p_IR){ # only checks the p if the p_IR is active
			($p_i)=$t[$i+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
			if ($p_i<0.05){
			    $PRINT{$event}.="\tNA"; # storages the PSIs
			    $tallyNA{$event}{$H[$i]}=1; # for summary stats 
			}
			next if $p_i<0.05;
		    }
		    ### For Alt3 and Alt5
		    if ($type eq "Alt3" || $type eq "Alt5"){
			my $kill_ALT = 0;
			my ($temp_ALT)=$t[$i+1]=~/O[KW]\,.+?\,(.+?)\,.+?\,.+?\@/;
			if ($temp_ALT=~/\d/){
			    $kill_ALT = 1 if $temp_ALT < $min_ALT_use;
			}
			else {
			    $min_ALT_use = "NA (older version)";
			}
			next if $kill_ALT;
		    }
		    # B3 check for AltEx events
		    if (($type eq "AltEx" || $type eq "MIC") && $noB3){ 
			my $kill_B3 = 0;
			my ($score3,$temp_B3)=$t[$i+1]=~/O[KW]\,.+?\,(.+?)\,(.+?)\,.+?\@/;
			if ($score3 =~ /\=/){ # i.e. from v2.2.2 onwards
			    my ($t_i1,$t_i2)=$score3=~/(\d+?)\=(\d+?)\=/;
			    $kill_B3 = 1 if $temp_B3 eq "B3" && $t_i1+$t_i2 > 15;
			}
			else {
			    $noB3="NA (older version)";
			}
			next if $kill_B3;
		    }
			
		    $total_N{$event}++;
		    push(@PSIs,$t[$i]);
		    
		    $PRINT{$event}.="\t$t[$i]"; # storages the PSIs
		    $tallyNA{$event}{$H[$i]}=0; # for summary stats
		}
		else {
		    $PRINT{$event}.="\tNA"; # storages the PSIs
		    $tallyNA{$event}{$H[$i]}=1; # for summary stats
		}
	    }
	}
	next if $min_N && $total_N{$event} < $min_N; # check for absolute number
	next if $min_Fraction && $total_N{$event}/$max_N < $min_Fraction; # check for fraction
	
	$SD{$event}=&std_dev(@PSIs);
	next if $SD{$event} < $min_SD;
	
	print O "$event_ID"."$PRINT{$event}\n";
	$tally_type{$type}++;
	$OK{$event}=1;
    }
    else {
	%PSIs=();
	$OK_group=0;
	foreach $group (sort keys %group_samples){
	    foreach $sample (@{$group_samples{$group}}){
		$i=$sample_index{$sample};
		
		if ($t[$i+1]=~/$Q/){
		    if ($type=~/IR/ && $p_IR){ # only checks the p if the p_IR is active
			($p_i)=$t[$i+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
			if ($p_i<0.05){
			    $PRINT{$event}.="\tNA"; # storages the PSIs
			    $tallyNA{$event}{$H[$i]}=1; # for summary stats 
			}
			next if $p_i<0.05;
		    }
		    ### For Alt3 and Alt5
		    if ($type eq "Alt3" || $type eq "Alt5"){
			my $kill_ALT = 0;
			my ($temp_ALT)=$t[$i+1]=~/O[KW]\,.+?\,(.+?)\,.+?\,.+?\@/;
			if ($temp_ALT=~/\d/){
			    $kill_ALT = 1 if $temp_ALT < $min_ALT_use;
			}
			else {
			    $min_ALT_use = "NA (older version)";
			}
			next if $kill_ALT;
		    }
		    # B3 check for AltEx events
		    if (($type eq "AltEx" || $type eq "MIC") && $noB3){ 
			my $kill_B3 = 0;
			my ($score3,$temp_B3)=$t[$i+1]=~/O[KW]\,.+?\,(.+?)\,(.+?)\,.+?\@/;
			if ($score3 =~ /\=/){ # i.e. from v2.2.2 onwards
			    my ($t_i1,$t_i2)=$score3=~/(\d+?)\=(\d+?)\=/;
			    $kill_B3 = 1 if $temp_B3 eq "B3" && $t_i1+$t_i2 > 15;
			}
			else {
			    $noB3="NA (older version)";
			}
			next if $kill_B3;
		    }		    

		    $total_N{$event}{$group}++;
		    push(@{$PSIs{$group}},$t[$i]);
		    
		    $PRINT{$event}.="\t$t[$i]"; # storages the PSIs
		    $tallyNA{$event}{$H[$i]}=0; # for summary stats
		}
		else {
		    $PRINT{$event}.="\tNA"; # storages the PSIs
		    $tallyNA{$event}{$H[$i]}=1; # for summary stats
		}
	    }
	    ### At the level of group
	    next if $min_N && $total_N{$event}{$group} < $min_N; # check for absolute number in each group
	    next if $min_Fraction && $total_N{$event}{$group}/$max_N_groups{$group} < $min_Fraction; # check for fraction
	    
	    $SD{$event}{$group}=&std_dev(@{$PSIs{$group}});
	    next if $SD{$event}{$group} < $min_SD;
	    $OK_group++;
	}
	if ($OK_group==2){
	    print O "$event_ID"."$PRINT{$event}\n";
	    $tally_type{$type}++;
	    $OK{$event}=1;
	}
    }
    $done{$event}=1;
}

### this scores the number of events missing in each sample
foreach $ev (sort keys %OK){
    $total_events++;
    foreach $tis (sort (keys %TISSUES)){
	$CUENTA{$tis}++ if $tallyNA{$ev}{$tis};
    }
}

$min_N="NA" if !$min_N;
$min_Fraction="NA" if !$min_Fraction;

### Print summary and LOG
open (LOG, ">$log_file") if $log;
print LOG "OPTIONS: $command_line\n\n";

$extras="";
$extras.= " -noVOW" if $noVLOW;
$extras.= " -noB3" if $noB3 && $noB3 ne "NA (older version)";
$extras.= " -noB3=NA (older version)" if $noB3 && $noB3 eq "NA (older version)";
$extras.= " -min_ALT_use $min_ALT_use";
$extras.= " -p_IR" if $p_IR;
$extras.= " GROUPS" if $group_file;

print LOG "\nSettings: -Min_N $min_N -Min_Fr $min_Fraction -Min_SD $min_SD$extras\n";
verbPrint "Settings: -Min_N $min_N -Min_Fr $min_Fraction -Min_SD $min_SD$extras\n";
print LOG "TOTAL # of Events: $total_events\n";
print "\t\tTOTAL # of Events: $total_events\n";
foreach $type (@TYPES){
    print LOG "$type\t$tally_type{$type}\n";
    print "\t\t$type\t$tally_type{$type}\n";
}
print LOG "\nTISSUE\tMISSING\t\%\n";
print "\n\t\tTISSUE\tMISSING\t\%\n";
foreach $tis (sort (keys %TISSUES)){
    $CUENTA{$tis}=0 if !$CUENTA{$tis};
    $perc=sprintf("%.2f",100*$CUENTA{$tis}/$total_events);
    print LOG "$tis\t$CUENTA{$tis}\t$perc\n";
    print "\t\t$tis\t$CUENTA{$tis}\t$perc\n";
}
print LOG "\n";
print "\n";



########################
sub average{
    my @data = @_;

    if ($#data==0) {
	die("Empty array");
    }
    my $total = 0;
    foreach (@data) {
	$total += $_;
    }
    my $average = sprintf("%.2f",$total / ($#data+1));
    return $average;
}
sub std_dev{
    my @data = @_;
    if($#data == 0){
	return 0;
    }
    my $average = &average(@data);
    my $sqtotal = 0;
    foreach(@data) {
	$sqtotal += ($average-$_) ** 2;
    }
    my $std = sprintf("%.2f",($sqtotal / ($#data)) ** 0.5); ## total or total - 1
    return $std;
}
