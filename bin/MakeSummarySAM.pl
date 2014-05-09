#!/usr/bin/perl

($r)=$ARGV[0]=~/(.+?)\.out/;# if $ARGV[0]!~/\//;
$sam=$&;

system "gunzip $ARGV[0]" if $ARGV[0]=~/\.gz/;

open (I, $sam);
open (O, ">$r.outsum");

while (<I>){
    $read="";
    @t=split(/\t/);
    
    ($read)=$t[0]=~/(.+) /;
    ($read)=$t[0]=~/(.+)\#/ if !$read;
    ($read)=$t[0]=~/(.+)\:/ if !$read;
    ($read)=$t[0]=~/(.+)\// if !$read;
    ($read)=$t[0]=~/(.+?)\-/ if !$read;
    $read=$t[0] if !$read;
    
    $hit=$t[2];
    
    if ($read ne $reB){
        $tally{$hit}++;
	$EJ{$hit}++;
        $POS{$hit}{$t[3]}++;
    }
    $reB=$read;
    $hitB=$hit;
}

foreach $hit (sort (keys %tally)){
    $p_p="";
    foreach $POS (sort {$a<=>$b}(keys %{$POS{$hit}})){
	$p_p.="$POS:$POS{$hit}{$POS},";
    }
    chop($p_p);

    print O "$hit\t$tally{$hit}\t$p_p\n";
}

system "rm $sam";
#system "gzip $r.outsum";
