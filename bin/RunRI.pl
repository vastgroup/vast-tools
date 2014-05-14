#!/usr/bin/perl
#This is the command line:
#submitjob 200 -m 8 -c x perl RunRI_UBxxxxxx.pl Sample-50.fq Mmu/Hsa 1/0 [cores]
#Output fake.cReadcount.txt is final output with corrected read counts
# Manuel Irimia, 13/09/27
# Changes UB 16/10/13 etc.
# Changes: - instead of having each fq file in its own dir, have them all in the same one and create dir only temporarily
#          - using Sandy junction and intron libries as well as mappabilites fro mouse and Manuel's for human
#          - normalizing I reads to the same length as junction reads
# Assumes trimming to a length for which a library exists in the specified paths.

$sample1 = $ARGV[0];
$sp = $ARGV[1];
if (!$ARGV[3]) {
    $cores = 1;
} else { 
    $cores = $ARGV[3];
}
($sample,$readlength) = $sample1=~/(.+?)\-(.+?)\./;

$dir1="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/Hsa/FILES/IR";
$dir2="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/Mmu/FILES/IR";
$dir3="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/bin";
$outdir="/home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRcounts";

$map_done=$ARGV[2];
print "$sample1 $readlength $sp $cores\n\n";
if (!$map_done){ # check if there is a galigned file present
    system "mkdir TMP_$sample";
    if (-e "$sample-$readlength-galigned.fa.gz"){
	system "mv $sample-$readlength-galigned.fa.gz $sample-$readlength-galigned.fq.gz";
	system "gunzip $sample-$readlength-galigned.fq.gz";
    }
    elsif (-e "$sample-$readlength-galigned.fq.gz"){
	system "gunzip $sample-$readlength-galigned.fq.gz";
    }
    elsif (-e "$sample-$readlength-galigned.fq"){
	print "File found: $sample-$readlength-galigned.fq\n";
    }
    else {  # create galigned file (SAM file with only reads that map to junctions or intron windows)
	system "gunzip $sample1" if $sample1=~/\.gz/;
	if ($sp eq "Mmu" ) {
	    system "bowtie -m 1 -v 2 -p $cores --al $sample-$readlength-galigned.fq $dir2/mm9.EEJ.new.$readlength.8 $sample-$readlength.fq /dev/null";
	} 
	elsif ($sp eq "Hsa") {
	    system "bowtie -m 1 -v 2 -p $cores --al $sample-$readlength-galigned.fq $dir1/$sp/$sp.EEJ.new.$readlength.8 $sample-$readlength.fq /dev/null";
	}
	else {
	    die "Species not found\n";
	}
	system "gzip $sample-$readlength.fq";
    }

    # Map reads against junctions and intron internal windows
    if ($sp eq "Mmu" ) {
	system "bowtie -m 1 -v 2 -p $cores $dir2/MouseIntronJunctions.new.$readlength.8 $sample-$readlength-galigned.fq /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."IJ-$readlength-$sample.out";    
	system "bowtie -m 1 -v 2 -p $cores $dir2/MouseIntrons.sample.200 $sample-$readlength-galigned.fq /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."INT-$readlength-$sample.out";
    } 
    elsif ($sp eq "Hsa") {
	system "bowtie -m 1 -v 2 -p $cores $dir1/$sp/$sp.IntronJunctions.new.$readlength.8 $sample-$readlength-galigned.fq /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."IJ-$readlength-$sample.out";
	system "bowtie -m 1 -v 2 -p $cores $dir1/$sp/$sp.Introns.sample.200 $sample-$readlength-galigned.fq /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."INT-$readlength-$sample.out";
    }
    else {
	die "Species not found\n";
    }
    system "gzip $sample-$readlength-galigned.fq";
    
    # sort reads by name, remove duplicates, simplify
    system "perl $dir3/sort_outs.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."INT-$readlength-$sample.out";
    system "perl $dir3/sort_outs.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."IJ-$readlength-$sample.out";
    
    system "perl $dir3/MakeSummarySAM.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."INT-$readlength-$sample"."_s.out";
    system "perl $dir3/MakeSummarySAM.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."IJ-$readlength-$sample"."_s.out";
}


print "perl $dir3/RI_summarize.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."IJ-$readlength-$sample"."_s.outsum\n";
system "perl $dir3/RI_summarize.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."IJ-$readlength-$sample"."_s.outsum";
# --> produces a $sample/$sample.summary.txt as in the original 
print "perl $dir3/RI_summarize_intron.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."INT-$readlength-$sample"."_s.outsum\n";
system "perl $dir3/RI_summarize_intron.pl /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sp"."INT-$readlength-$sample"."_s.outsum";

# check if output dir exists and move result there
unless (-e $outdir) { 
    system "mkdir $outdir";
}
system "mv /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp/TMP_$sample/$sample.cReadcount.txt $outdir/";
#system "rm -R /home/blencowe/blencowe1/ulrich/AS_PIPE_S_dev/$sp/IRtemp";



