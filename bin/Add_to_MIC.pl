#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

# This script is to calculate PSIs for MIC and produce a table with them for all samples.
use Cwd;
#$cwd = getcwd;
#($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;
use Getopt::Long;

my $dbDir;
my $sp;
my $samLen;
my $legacyFlag;

GetOptions("dbDir=s" =>\$dbDir, "sp=s" => \$sp, "len=i" => \$samLen,
			  "verbose=i" => \$verboseFlag, "legacy" => \$legacyFlag);

sub verbPrint {
  my $verbMsg = shift;
  if($verboseFlag) {
    chomp($verbMsg);
    print STDERR "[vast combine micro]: $verbMsg\n";
  }
}

#$sp=$ARGV[0];
die "[vast combine micro]: Needs the 3-letter species key\n" if !defined($sp);

@files=glob("spli_out/*micX");

$head_counts=$head_PSI="GENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX";
foreach $file (@files){
 #   ($sample)=$file=~/MIC\-\d+?\-(.+)\./;
    my $fname = $file;
    $fname =~ s/^.*\///;
    ($sample)=$fname=~/^(.*)\..*$/;
    $head_PSI.="\t$sample\t$sample-Q";
    $head_counts.="\t$sample-Rexc\t$sample-Rinc\t$sample-exc\t$sample-inc\tPSI=Q";
    open (INPUT, $file);
    while (<INPUT>){
	chomp;
	@t=split(/\t/);
	$event=$t[1];
	$pre_data{$event}=join("\t",@t[0..5]) if !$pre_data{$event};

	# Quality scores (only coverage)
	$Q="";i
   my $raw_reads_exc=$t[7];
   my $raw_reads_inc=$t[8];
   my $corr_reads_exc=$t[9];
   my $corr_reads_inc=$t[10];
#### Score 1: Using all raw reads                                                                                                       
   if (($raw_reads_exc+$raw_reads_inc)>=100){$Q="SOK";}
   elsif ($raw_reads_exc>=20 || $raw_reads_inc >=20){$Q="OK";}
   elsif ($raw_reads_exc>=15 || $raw_reads_inc >=15){$Q="LOW";}
   elsif ($raw_reads_exc>=10 || $raw_reads_inc >=10){$Q="VLOW";}
   else {$Q="N";}
#### Score 2: Using all corrected reads                                                                                                 
   if (($corr_reads_exc+$corr_reads_inc)>=100){$Q.=",SOK";}
   elsif ($corr_reads_exc>=20 || $corr_reads_inc>=20){$Q.=",OK";}
   elsif ($corr_reads_exc>=15 || $corr_reads_inc>=15){$Q.=",LOW";}
   elsif ($corr_reads_exc>=10 || $corr_reads_inc>=10){$Q.=",VLOW";}
   else {$Q.=",N";}
#### Score 3: Using simple (=reference, C1A, AC2, C1C2) raw reads                                                                       
   if (($raw_reads_exc+$raw_reads_inc)>=100){$Q.=",SOK";}
   elsif ($raw_reads_exc>=20 || $raw_reads_inc>=20){$Q.=",OK";}
   elsif ($raw_reads_exc>=15 || $raw_reads_inc>=15){$Q.=",LOW";}
   elsif ($raw_reads_exc>=10 || $raw_reads_inc>=10){$Q.=",VLOW";}
   else {$Q.=",N";}
#### No scores 4 and 5 for microexon pipeline.	
	$Q.=",na,na";

   ### DIFF OUTPUT ADDITION TO QUAL SCORE!
   ### Essentially adding the expected number of reads re-distributed to INC or EXC after normalization..
   ### These values are added to the qual score and used to infer the posterior distribution
   unless($legacyFlag) {
     my $totalN = $raw_reads_inc + $raw_reads_exc;
     my($pPSI, $exValOfInc, $exValOfExc) = (0, 0, 0);
     unless($PSI eq "NA" or $totalN < 2) {
       $pPSI = $PSI / 100;
       $exValOfInc = $pPSI * $totalN;
       $exValOfExc = (1-$pPSI) * $totalN;
     }
     # ALTER QUAL OUTPUT HERE>>
     $Q .= "\@$exValOfInc,$exValOfExc";
   }

	$PSI{$event}.="\t$t[6]\t$Q";
	$read_counts{$event}.="\t".join("\t",@t[7..10])."\t$t[6]=$Q";
    }
    close INPUT;
}
$N=$#files+1;
open (PSIs, ">raw_incl/INCLUSION_LEVELS_MIC-$sp$N-n.tab");
open (COUNTs, ">raw_reads/RAW_READS_MIC-$sp$N-n.tab");

print PSIs "$head_PSI\n";
print COUNTs "$head_counts\n";

foreach $event (sort keys %pre_data){
    print PSIs "$pre_data{$event}"."$PSI{$event}\n";
    print COUNTs "$pre_data{$event}"."$read_counts{$event}\n";
}

close PSIs;
close COUNTs;

