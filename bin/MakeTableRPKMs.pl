#!/usr/bin/env perl
# This script puts all cRPKMs together in two tables (one alone, the other one with the raw read counts for each sample).

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use FuncBasics qw(:all);

use Cwd;
use Getopt::Long;

my $dbDir;
my $sp; # here, the internal sp key is provided only.
my $sp_assembly;
my $cRPKMCounts = 0;
my $normalize = 0;
my $install_limma = 0;

GetOptions("dbDir=s" => \$dbDir, "sp=s" => \$sp, "C" => \$cRPKMCounts, "norm" => \$normalize, "install_limma" => \$install_limma);

die "[vast combine cRPKM error] Needs Species\n" if !$sp;
my @files=glob("expr_out/*.cRPKM");
my $index=$#files+1;

### gets the species assembly
get_internal_sp_key($sp);

# creates output files where run
my $OUTPUT;
open ($OUTPUT, ">cRPKM_AND_COUNTS-$sp_assembly-$index.tab") if $cRPKMCounts;
open (RPKM, ">cRPKM-$sp_assembly-$index.tab");

my %names;
open (NAMES, "$dbDir/FILES/$sp.ID.names.txt");
while (<NAMES>){
    chomp;
    my @t=split(/\t/);
    $names{$t[0]}=$t[1] if (defined $t[0]);
}
close NAMES;

my %data;
my %RPKM;
my $head="ID\tNAME";
my $headRPKM=$head;
my $first_sample_count=0;
foreach my $f (@files){
    my ($root)=$f=~/([^\/]+).cRPKM/;
    $head.="\t$root-cRPKM\t$root-Counts";
    $headRPKM.="\t$root";
    my $sample_count=0;
    
    open (INPUT, $f);
    while (<INPUT>){
        chomp;
        my @t=split(/\t/);
        my $cRPKM = sprintf("%.2f", 0);
        my $raw_count = 0;
	$sample_count++;

        if ($t[1] eq 'NA') {
            $cRPKM = 'NA';
            $raw_count = 'NA';
        } elsif ($t[1] != 0) {
            $cRPKM=sprintf("%.2f",$t[1]);
            $raw_count=$t[2];
        }
        
        $data{$t[0]}.="\t$cRPKM\t$raw_count";
        $RPKM{$t[0]}.="\t$cRPKM";
    }
    close INPUT;

    if ($sample_count >0 && $first_sample_count==0){
	$first_sample_count=$sample_count;
    }
    elsif ($sample_count > 0 && $sample_count != $first_sample_count){
	die "[vast combine cRPKM error] Count files do not have the same number of genes (incorrect versions?)\n";
    }
}

print $OUTPUT "$head\n" if $cRPKMCounts;
print RPKM "$headRPKM\n";

foreach my $gene (sort (keys %data)){
    $names{$gene}="NA" if (!defined $names{$gene});
    print $OUTPUT "$gene\t$names{$gene}$data{$gene}\n" if $cRPKMCounts;
    print RPKM "$gene\t$names{$gene}$RPKM{$gene}\n";
}
close $OUTPUT if $cRPKMCounts;
close RPKM;


### Makes normalized table:
my %norm_cRPKMs;
if ($normalize){
    my $input_file = "cRPKM-$sp_assembly-$index.tab";
    open (GE_2, $input_file) || die "[vast combine cRPKM error] Needs a cRPKM table\n";
    my ($input_path,$root_input);
    if ($input_file=~/\//){
	($input_path,$root_input) = $input_file =~/(.+)\/(.+)\./;
    }
    else {
	$input_path=".";
	($root_input) = $input_file =~/(.+)\./;
    }
    open (TEMP, ">$input_path/temp_cRPKMs.tab");
    while (<GE_2>){
	chomp($_);
	my @t=split(/\t/,$_);
	print TEMP "$t[0]";
	foreach my $i (2..$#t){ # not count table
	    print TEMP "\t$t[$i]";
	}
	print TEMP "\n";
    }
    close TEMP;
    close GE_2;
    
    open (Temp_R, ">$input_path/temp.R");
    print Temp_R "source(\"https://bioconductor.org/biocLite.R\")
biocLite(\"limma\")\n" if ($install_limma);
    
    print Temp_R "
library(limma)
setwd(\"$input_path/\")
matrix=as.matrix(read.table(\"temp_cRPKMs.tab\", row.names=1, header=TRUE,sep=\"\\t\"))
Nmatrix=normalizeBetweenArrays(as.matrix(matrix))
NmatrixF=cbind(Names=row.names(matrix),Nmatrix)
write.table(NmatrixF,\"$root_input-NORM.tab\",
            sep=\"\\t\",col.names=T,row.names=F,quote=F)";
    close Temp_R;
    
    system "Rscript $input_path/temp.R";
    system "rm $input_path/temp*";
}

sub get_internal_sp_key {
    my @temp_assembly = @_;

    my %assembly_to_species;
    $assembly_to_species{hg19}="Hsa"; $assembly_to_species{hg38}="Hs2"; $assembly_to_species{panTro4}="Ptr";
    $assembly_to_species{rheMac2}="Mma"; $assembly_to_species{mm9}="Mmu"; $assembly_to_species{mm10}="Mm2";
    $assembly_to_species{bosTau6}="Bta"; $assembly_to_species{bosTau9}="Bt2"; $assembly_to_species{monDom5}="Mdo";
    $assembly_to_species{galGal3}="Gg3"; $assembly_to_species{galGal4}="Gg4"; $assembly_to_species{galGal6}="Gga";
    $assembly_to_species{xenTro3}="Xt1"; $assembly_to_species{xenTro9}="Xtr"; 
    $assembly_to_species{danRer10}="Dre"; $assembly_to_species{danRer11}="Dr2";
    $assembly_to_species{lepOcu1}="Loc"; $assembly_to_species{esoLuc2}="Elu"; $assembly_to_species{eshark1}="Cm1";
    $assembly_to_species{calMil1}="Cmi"; $assembly_to_species{braLan2}="Bl1"; $assembly_to_species{braLan3}="Bla";
    $assembly_to_species{strPur4}="Spu"; $assembly_to_species{dm6}="Dme"; $assembly_to_species{AaegL5}="Aae";
    $assembly_to_species{bomMor1}="Bmo"; $assembly_to_species{triCas5}="Tca"; $assembly_to_species{apiMel4}="Ame";
    $assembly_to_species{blaGer1}="Bge"; $assembly_to_species{cloDip2}="Cdi"; $assembly_to_species{strMar1}="Sma";
    $assembly_to_species{WBcel235}="Cel"; $assembly_to_species{octBim1}="Obi"; $assembly_to_species{octMin1}="Omi";
    $assembly_to_species{schMed31}="Sme"; $assembly_to_species{nemVec1}="Nve"; $assembly_to_species{TAIR10}="Ath";
    my %species_to_assembly;
    $species_to_assembly{Hsa}="hg19"; $species_to_assembly{Hs2}="hg38"; $species_to_assembly{Ptr}="panTro4";
    $species_to_assembly{Mma}="rheMac2"; $species_to_assembly{Mmu}="mm9"; $species_to_assembly{Mm2}="mm10";
    $species_to_assembly{Bta}="bosTau6"; $species_to_assembly{Bt2}="bosTau9"; $species_to_assembly{Mdo}="monDom5";
    $species_to_assembly{Gg3}="galGal3"; $species_to_assembly{Gg4}="galGal4"; $species_to_assembly{Gga}="galGal6";
    $species_to_assembly{Xt1}="xenTro3"; $species_to_assembly{Xtr}="xenTro9";
    $species_to_assembly{Dre}="danRer10"; $species_to_assembly{Dr2}="danRer11";
    $species_to_assembly{Loc}="lepOcu1"; $species_to_assembly{Elu}="esoLuc2"; $species_to_assembly{Cm1}="eshark1";
    $species_to_assembly{Cmi}="calMil1"; $species_to_assembly{Bl1}="braLan2"; $species_to_assembly{Bla}="braLan3";
    $species_to_assembly{Spu}="strPur4"; $species_to_assembly{Dme}="dm6"; $species_to_assembly{Aae}="AaegL5";
    $species_to_assembly{Bmo}="bomMor1"; $species_to_assembly{Tca}="triCas5"; $species_to_assembly{Ame}="apiMel4";
    $species_to_assembly{Bge}="blaGer1"; $species_to_assembly{Cdi}="cloDip2"; $species_to_assembly{Sma}="strMar1";
    $species_to_assembly{Cel}="WBcel235"; $species_to_assembly{Obi}="octBim1"; $species_to_assembly{Omi}="octMin1";
    $species_to_assembly{Sme}="schMed31"; $species_to_assembly{Nve}="nemVec1"; $species_to_assembly{Ath}="TAIR10";

    if (defined $species_to_assembly{$temp_assembly[0]}){
	$sp_assembly = $species_to_assembly{$temp_assembly[0]};
    }
    else {
	die "$temp_assembly[0] is not a valid species key\n";
    }
}
