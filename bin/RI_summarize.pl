# This script gets, for each event, the normalized read counts for the three associated junctions in the sample

$file=$ARGV[0];
($species,$rle,$sample)=$file=~/(.{3})IJ\-(\d+?)\-(.+?)\_s\.out/;
$sp=$species;

$dir1="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/Hsa/FILES";
$dir2="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/Mmu/FILES";

my $maxcount = $rle - 15;

# Getting mappability inexit bash scriptformation (from output of uniquecount.IJ.pl)
my %ucount;

if ($sp eq 'Mmu') {
    open (UC, "$dir2/$sp.IntronJunctions.new.$rle.8.uniquecount.txt") || die "Can't find mappability for IJ ($dir2/$sp.IntronJunctions.new.$rle.8.uniquecount.txt)\n";
}
elsif ($sp eq 'Hsa') {
    open (UC, "$dir1/$sp.IntronJunctions.new.$rle.8.uniquecount.txt") || die "Can't find mappability for IJ (*$dir1/$sp.IntronJunctions.new.$rle.8.uniquecount.txt)\n";
}
else {
    die "Unkown species\n";
}

while(<UC>){
    chomp($_);
    my($junction0,$count) = split(/\t/,$_);
    my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
    my $junction = $gene.":".$trans.":".$en;
    $ucount{$junction} = $count;
}
close UC;

# Getting raw read counts
my %rcount;
open(RC,$file);
while(<RC>){
    chomp($_);
    my($junction0,$count,$pos) = split(/\t/,$_);
    my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
    my $junction = $gene.":".$trans.":".$en;
    $rcount{$junction} = $count;
}
close RC;

# Getting junction annotation (i.e. the IDs of the 3 junctions associated with each event) and generating output file
my %eventseen;
if ($sp eq "Mmu") {
    open (ANOT,"$dir2/$sp.IntronJunctions.new.annotation.txt") || die "Can't find annotations for IJ run\n";
}
elsif ($sp eq "Hsa") {
    open (ANOT, "$dir1/$sp.IntronJunctions.new.annotation.txt") || die "Can't find annotations for IJ run\n";
}

my $outfile = "/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sample.summary.txt";
open(OUT,">$outfile");
my $head = <ANOT>;
print OUT "Event\tEIJ1\tEIJ2\tEEJ\n";
while(<ANOT>){
    chomp($_);
    my ($event,$C1Aj,$AC2j,$C1C2j,@rest) = split(/\t/,$_);
    if($eventseen{$event} ne "Y" && ($rcount{$C1Aj} =~ /[0-9]/ || $rcount{$AC2j} =~ /[0-9]/ || $rcount{$C1C2j} =~ /[0-9]/)){
	my $C1A = 0;
	if($rcount{$C1Aj} =~ /[0-9]/ && $ucount{$C1Aj} =~ /[0-9]/ && $ucount{$C1Aj} > 0){
	    $C1A = $rcount{$C1Aj} / $ucount{$C1Aj} * $maxcount;
	}
	my $AC2 = 0;
	if($rcount{$AC2j} =~ /[0-9]/ && $ucount{$AC2j} =~ /[0-9]/ && $ucount{$AC2j} > 0){
	    $AC2 = $rcount{$AC2j} / $ucount{$AC2j} * $maxcount;
	}
	my $C1C2 = 0;
	if($rcount{$C1C2j} =~ /[0-9]/ && $ucount{$C1C2j} =~ /[0-9]/ && $ucount{$C1C2j} > 0){
	    $C1C2 = $rcount{$C1C2j} / $ucount{$C1C2j} * $maxcount;
	}
	print OUT "$event\t$C1A\t$AC2\t$C1C2\n";
	$eventseen{$event} = "Y";
    }
}
close ANOT;
close OUT;
