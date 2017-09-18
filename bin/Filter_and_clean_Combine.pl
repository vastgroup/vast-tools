#!/usr/bin/perl
# Script to prepare and filter vast-tools PSI tables for other analyses

use Getopt::Long;

### Setting global variables:
$Q="O[KW]\,.+?\,.+?\,.+?\,.+?\@"; # NEW quality search

$input_file=$ARGV[0];
$min_SD=5; # min standard deviation of the event (def=5)
#$min_Fraction=0.8; # min fraction of samples with good coverage (def=0.8)
#$min_N=10; # min number of samples with good coverage (def=10)
$noVLOW=0;
$p_IR=0;
$command_line=join(" ", @ARGV);

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "min_SD=i" => \$min_SD,
                          "min_N=i" => \$min_N,
                          "min_Fr=f" => \$min_Fraction,
                          "outFile=s"    => \$output_file,
                          "help" => \$helpFlag,
                          "p_IR" => \$p_IR,
                          "noVLOW" => \$noVLOW,
			  "log" => \$log,
			  "onlyEXSK" => \$onlyEXSK,
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



# defining default output file name
($root)=$input_file=~/(.+)\./;
#$root=~s/.+\///;
$root_out=$root;
$root_out.="-minN_$min_N" if $min_N;
$root_out.="-minFr_$min_Fr" if $min_Fr;
$root_out.="-minSD_$min_SD";
$root_out.="-noVLOW" if $noVLOW;
$root_out.="-p_IR" if $p_IR;
$root_out.="-onlyEXSK" if $onlyEXSK;

$output_file="$root_out-Tidy.tab" unless $output_file;
$log_file="$root_out-Tidy.log" if $log;

if (!$ARGV[0] || $helpFlag){
    die "\nUsage: vast-tools tidy INCLUSION_LEVELS_FULL-SpN.tab (--min_N min_N_samples OR --min_Fr min_fraction) --min_SD min_SD [options]

Prepares and filters a vast-tools output for general analyses.

[General options] 
        -min_N i                Minimum number of samples with good coverage. Alternatively, you can define it by fraction
        -min_Fr f               Minimum fraction of samples with good coverage. 
        -min_SD i               Minimum standard deviation of the event (def=5)
        -outFile                Output file name (default based on option parameters)
        --noVLOW                Do not use samples with VLOW coverage (default OFF)
        --p_IR                  Filter IR by the p-value of the binomial test (default OFF)
        --onlyEXSK              Outputs only EXSK events (default OFF)
        --add_names             Adds gene name to the event_ID. E.g. Mta1=MmuEX0029874 (default OFF)
        --log                   Print the summary stats into a file (default OFF)


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

errPrintDie "*** You can only define a minimum fraction or absolute number of samples with good coverage\n" if $min_N && $min_Fraction;
errPrintDie "*** You need to define either a minimum fraction or absolute number of samples with good coverage\n" if !$min_N && !$min_Fraction;

open (I, $input_file) || die "Can't open input file\n";
open (O, ">$output_file") || die "Can't open output file ($output_file)\n";

### saying hi:
verbPrint "Cleaning and filtering $ARGV[0]\n";


### Printing a new clean heading (only Event_ID and PSI)
print O "EVENT";
$head=<I>;
chomp($head);
@H=split(/\t/,$head);
for $i (6..$#H){
    if ($i%2==0){
	print O "\t$H[$i]";
        $TISSUES{$H[$i]}=1;
	$max_N++;
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
    $type=$t[5];

    next if $length==0; # to remove the internal ss in Alt3 and Alt5
    next if $onlyEXSK && ($type=~/Alt/ || $type=~/IR/);
    
    $event_ID="$gene_name=$event" if $AddName;
    $event_ID=$event if !$AddName;

    %PRINT=();
    @PSIs=();
    next if $done{$event}; # to avoid possible repeated IDs (in preliminary tables)
    
    foreach $i (6..$#t){
	if ($i%2==0){
	    if ($t[$i+1]=~/$Q/){
		if ($type=~/IR/ && $p_IR){ # only checks the p if the p_IR is active
		    ($p_i)=$t[$i+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
		    if ($p_i<0.05){
			$PRINT{$event}.="\tNA"; # storages the PSIs
			$tallyNA{$event}{$H[$i]}=1; # for summary stats 
		    }
		    next if $p_i<0.05;
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
    $done{$event}=1;
    $OK{$event}=1;
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
$extras.= " -p_IR" if $p_IR;

print LOG "\nSettings: -Min_N $min_N -Min_Fr $min_Fraction -Min_SD $min_SD$extras\n";
print "\nSettings: -Min_N $min_N -Min_Fr $min_Fraction -Min_SD $min_SD$extras\n";
print LOG "TOTAL # of Events: $total_events\n";
print "TOTAL # of Events: $total_events\n";
foreach $type (sort keys %tally_type){
    print LOG "$type\t$tally_type{$type}\n";
    print "$type\t$tally_type{$type}\n";
}
print LOG "\nTISSUE\tMISSING\t\%\n";
print "\nTISSUE\tMISSING\t\%\n";
foreach $tis (sort (keys %TISSUES)){
    $CUENTA{$tis}=0 if !$CUENTA{$tis};
    $perc=sprintf("%.2f",100*$CUENTA{$tis}/$total_events);
    print LOG "$tis\t$CUENTA{$tis}\t$perc\n";
    print "$tis\t$CUENTA{$tis}\t$perc\n";
}
print LOG "\n";
print "\n";



########################
sub average{
    my @data = @PSIs;

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
