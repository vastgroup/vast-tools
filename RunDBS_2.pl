#!/usr/bin/perl
# This pipeline takes PSI templates and adds PSIs from new samples.
use Cwd;
$cwd = getcwd;
($dir)=$cwd=~/(.+?\/AS_PIPE_S)/;

if (!$ARGV[0] || $ARGV[0]=~/\-h/){
    print "\nCommand:\nRunDBS_2.pl Species (Hsa or Mmu)\n\n";
    print ">> It will add all the samples that are in the running folder (=SAMPLES)\n";
    print ">>>> Root names for each sample must be identical (they will if coming from RunDBS_1)\n";
    die ">> Not yet for IR\n\n";
}

### Settings:
$sp=$ARGV[0];
die "Needs species 3-letter key\n" if !$sp;
@files=glob("$dir/$sp/SAMPLES/*exskX"); #gathers all exskX files (a priori, simple).
$N=$#files+1;

### Gets the PSIs for the events in the a posteriori pipeline
print "\nBuilding Table for COMBI (a posteriori pipeline)\n";
system "$dir/bin/Add_to_COMBI.pl $sp";

### Gets the PSIs for the a priori, SIMPLE
print "\nBuilding Table for EXSK (a priori pipeline, single)\n";
system "$dir/bin/Add_to_APR.pl $sp exskX";

### Gets the PSIs for the a priori, COMPLEX
print "\nBuilding Table for MULTI (a priori pipeline, multiexon)\n";
system "$dir/bin/Add_to_APR.pl $sp MULTI3X";

### Gets the PSIs for the MIC pipeline
print "\nBuilding Table for MIC (microexons)\n";
system "$dir/bin/Add_to_MIC.pl $sp";

### Adds those PSIs to the full database of PSIs (MERGE3m).
print "\nBuilding non-redundant PSI table (MERGE3m)\n";
system "$dir/bin/Add_to_MERGE3m.pl $dir/$sp/SAMPLES/INCLUSION_LEVELS_EXSK-$sp$N-n.tab $dir/$sp/SAMPLES/INCLUSION_LEVELS_MULTI-$sp$N-n.tab $dir/$sp/SAMPLES/INCLUSION_LEVELS_COMBI-$sp$N-n.tab $dir/$sp/SAMPLES/INCLUSION_LEVELS_MIC-$sp$N-n.tab";

### Gets PSIs for ALT5ss and adds them to the general database
print "\nBuilding Table for Alternative 5'ss choice events\n";
system "$dir/bin/Add_to_ALT5.pl $sp";

### Gets PSIs for ALT3ss and adds them to the general database
print "\nBuilding Table for Alternative 3'ss choice events\n";
system "$dir/bin/Add_to_ALT3.pl $sp";

