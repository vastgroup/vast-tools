# This script gets, for each event, the normalized read counts for the three associated junctions in the sample

$file=$ARGV[0];
($species,$rle,$sample)=$file=~/(.{3})INT\-(\d+?)\-(.+?)\_s\.out/;
$sp=$species;

$dir1="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/Hsa/FILES/IR";
$dir2="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/Mmu/FILES/IR";

my $maxcount = $rle - 15;

# Getting mappability information (from output of uniquecount.IJ.pl)
my %ucount;

if ($sp eq 'Mmu') {
    open (UC,"$dir2/MouseIntrons.sample.200.$rle.uniquecount.txt") || die "Can't find mappability for INT\n";
} 
elsif ($sp eq 'Hsa') {
    open (UC, "$dir1/$sp/$sp.Introns.sample.200.$rle.uniquecount.txt") || die "Can't find mappability for INT\n";
} 
else {
    die "Unkown species\n";
}

while(<UC>){
    chomp($_);
    my($junction0,$count,$pos) = split(/\t/,$_);
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
    my($junction0,$count) = split(/\t/,$_);
    my($gene,$trans,$en,$l1,$l2) = split(/\:/,$junction0);
    my $junction = $gene.":".$trans.":".$en;
    $rcount{$junction} = $count;
}
close RC;

# Getting junction annotation (i.e. the IDs of the 3 junctions associated with each event) and generating output file
my %eventseen;
open(ANOT,"/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sample.summary.txt");
my $outfile = "/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sample.cReadcount.txt";
open(OUT,">$outfile");
my $head = <ANOT>;
print OUT "Event\tEIJ1\tEIJ2\tEEJ\tI\n";
while(<ANOT>){
    chomp($_);
    my ($event,$C1A,$AC2,$C1C2) = split(/\t/,$_);
    my $I = 0;

    if($rcount{$event} =~ /[0-9]/ && $ucount{$event} =~ /[0-9]/ && $ucount{$event} > 0){
	    $I = $rcount{$event} / $ucount{$event} * $maxcount;
    }
    print OUT "$event\t$C1A\t$AC2\t$C1C2\t$I\n";
}
close ANOT;
close OUT;
