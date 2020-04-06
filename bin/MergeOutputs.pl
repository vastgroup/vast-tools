#!/usr/bin/perl
# Script to merge vast-tools outputs from different sub-samples into samples
# Format of the table with groupings must be:
# Subsample1\tSampleA
# Subsample2\tSampleA
# Subsample3\tSampleB
# ...

use warnings;
use strict;
use Cwd qw(abs_path);
use Getopt::Long;

# INITIALIZE PATH AND FLAGS
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

my $helpFlag = 0;
my $verboseFlag = 1; 
my $dbDir; # directory of VASTDB
my $species; # needed to run expr merge automatically
my $sp_assembly; # the new input from v2.4.0
my $groups; # list with the groupings: sample1_rep1\tgroup_1\n sample1_rep2\tgroup_1\n...
my $folder; # actual folder where the vast-tools outputs are (the to_combine folder!)
my $mapcorr_fileswitch="";
my $effective; #effective file for expression. Obtained automatically from VASTDB
my $expr; #activates merging of cRPKMs
my $exprONLY; # if you want to do expr only, just write anything.
my $move_to_PARTS; # to move the merged subfiles into the PARTS folder
my $IR_version = 2; # version of IR pipeline [new default 01/04/16]
my $noIR; # to avoid warnings

Getopt::Long::Configure("no_auto_abbrev");
GetOptions("groups=s" => \$groups,
	"g=s" => \$groups,
	"dbDir=s" => \$dbDir,
	"outDir=s" => \$folder,
	"o=s" => \$folder,
	"IR_version=i" => \$IR_version,
	"sp=s" => \$sp_assembly,
	"expr" => \$expr,
	"exprONLY" => \$exprONLY,
	"help" => \$helpFlag,
	"noIR" => \$noIR,
	"move_to_PARTS" => \$move_to_PARTS);


our $EXIT_STATUS = 0;

sub errPrint {
    my $errMsg = shift;
    print STDERR "[vast merge error]: $errMsg\n";
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
	print STDERR "[vast merge]: $verbMsg\n";
    }
}

sub time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
    $year += 1900;
    $mon += 1;
    my $datetime = sprintf "%04d-%02d-%02d (%02d:%02d)", $year, $mday, $mon, $hour, $min;
    return $datetime;
}

### Gets the version
my $version;
open (VERSION, "$binPath/../VERSION");
$version=<VERSION>;
chomp($version);
$version="No version found" if !$version;

if (!defined($groups) || $helpFlag){
    die "
VAST-TOOLS v$version

Usage: vast-tools merge -g path/groups_file [-o align_output] [options]

Merges vast-tools outputs from multiple subsamples into grouped samples

OPTIONS: 
        -g, --groups file        File with groupings (subsample1\\tsampleA\\nsubsample2\\tsampleA...)
        -o, --outDir             Path to output folder of vast-tools align (default vast_out)
                                 Must contain sub-folders to_combine or expr_out from align steps.
        --sp hg19/mm10/etc       Assembly of the species (e.g. hg19, mm10).
                                   The legacy three letter code for the database can also be used.
        --dbDir db               Database directory (default VASTDB)
        --IR_version 1/2         Version of the Intron Retention pipeline (1 or 2) (default 2)
        --expr                   Merges cRPKM files (default OFF)
        --exprONLY               Merges only cRPKM files (default OFF)
        --move_to_PARTS          Moves the subsample files to PARTS\/ within output folders (default OFF)
        --noIR                   Skips IR files
        --help                   Prints this help message


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

### Redefines species and sp_assembly
get_internal_sp_key($sp_assembly);

# Check database directory
if (defined $expr or defined $exprONLY){
    errPrintDie "Needs to provide species (\-\-sp)\n" if (!defined $sp_assembly);

    unless(defined($dbDir)) {
	$dbDir = "$binPath/../VASTDB";
    }
    $dbDir = abs_path($dbDir);
    $dbDir .= "/$species";
    errPrint "The database directory $dbDir does not exist" unless (-e $dbDir or $helpFlag);
}

# If exprONLY, activates expr
$expr=1 if (defined $exprONLY);

# Sanity checks
errPrintDie "Needs a file with the groupings\n" if (!defined $groups);
errPrintDie "IR version must be either 1 or 2\n" if ($IR_version != 1 && $IR_version != 2 && !defined $noIR);

# get full path to group definition file
# should come before we "leave" the working directory
# from where user has called vast-tools merge
my $groups_fullpath=abs_path($groups);

# prints version (05/05/19)                                                                                                                                             
verbPrint "VAST-TOOLS v$version";


verbPrint "Using VASTDB -> $dbDir" if (defined $expr);
# change directories
errPrintDie "The output directory \"$folder/to_combine\" does not exist. Check path specified with argument -o." unless (-e "$folder/to_combine");
chdir($folder) or errPrint "Unable to change directories into output" and die;
verbPrint "Setting output directory to $folder";

### Creates the LOG
open (LOG, ">>$folder/VTS_LOG_commands.txt");
my $all_args="-o $folder -groups $groups -IR_version $IR_version";
$all_args.=" -move_to_PARTS" if $move_to_PARTS;
$all_args.=" -sp $sp_assembly" if defined $sp_assembly;
$all_args.=" -expr" if $expr;
$all_args.=" -exprONLY" if $exprONLY;
$all_args.=" -noIR" if $noIR;

print LOG "[VAST-TOOLS v$version, ".&time."] vast-tools merge $all_args\n";


my %groups;
my %files_2b_merged;
my %file_2_groups;
my %file_grpchk;
my %file_is_ss;
my %group_is_ss;


if (defined $move_to_PARTS){
    system "mkdir to_combine/PARTS" unless (-e "to_combine/PARTS");
    system "mkdir expr_out/PARTS" unless (-e "expr_out/PARTS") || (!defined $expr);    
}


### Loading group info
my %info_files_of_subsamples=();
open (GROUPS, $groups_fullpath) || errPrintDie "Cannot open groupings: $groups_fullpath\n";
while (<GROUPS>){
    # cleaning in case they were made in Mac's excel
    $_ =~ s/\r//g;
    $_ =~ s/\"//g;
    chomp($_);
    
    my $chk_line=$_;
    $chk_line =~ s/\s//g; # delete all white spaces
    if(length($chk_line)==0){next;} 
    
    my @temp = split(/\t/,$_);
    
    if(!defined($file_2_groups{$temp[0]})){$file_2_groups{$temp[0]}=[];}  # for each file stores the groups it should get merged into 
    push(@{$file_2_groups{$temp[0]}},$temp[1]);                           # if a file should get merged several times into the same group, we store this group several times
    							                  # the same file might get merged into different groups
    if(!defined($file_grpchk{$temp[0]})){$file_grpchk{$temp[0]}={};}
    $file_grpchk{$temp[0]}->{$temp[1]}++;

    $groups{$temp[1]}=1;
    $files_2b_merged{$temp[0]}=1;
    my $tmp_file="to_combine/${temp[0]}.info";
    
    unless(-e "to_combine/${temp[0]}.info"){ verbPrint "$temp[0]: do not find to_combine/${temp[0]}.info. Sample will be treated as being not strand-specific.";
    }else{
    	open(my $fh_info,"to_combine/${temp[0]}.info") or die "$!"; my $line=<$fh_info>; close($fh_info);
    	$info_files_of_subsamples{"${temp[0]}.info"}=1;
    	my @fs=split("\t",$line);
    	if($fs[@fs-2] eq "-SS"){
    		$file_is_ss{$temp[0]}=1;
    		verbPrint "$temp[0]: found to_combine/${temp[0]}.info. Sample will be treated as being strand-specific."
    	}else{
    		verbPrint "$temp[0]: found to_combine/${temp[0]}.info. Sample will be treated as being not strand-specific."
    	}
    }

    # check if all necessary files exists; if not the groups file might contain incorrect subsample names
    $tmp_file="expr_out/${temp[0]}.cRPKM";
    if($expr && !(-e $tmp_file)){errPrintDie "File $tmp_file does not exist. Probable reason: wrong subsample names in group-definition file OR vast-tools align was run without mapping for expression analysis.";}
    unless(defined $exprONLY){
    	 $tmp_file="to_combine/${temp[0]}.IR";if($IR_version==2){$tmp_file="to_combine/${temp[0]}.IR2"};
	 if(!defined $noIR && !(-e $tmp_file)){goto STOPSCRIPT;}
    	 $tmp_file="to_combine/${temp[0]}.IR.summary_v2.txt";
    	 if(!defined $noIR && $IR_version==2 && !(-e $tmp_file)){goto STOPSCRIPT;}
	 $tmp_file="to_combine/${temp[0]}.micX";unless(-e $tmp_file){goto STOPSCRIPT;}
	 $tmp_file="to_combine/${temp[0]}.eej2";unless(-e $tmp_file){goto STOPSCRIPT;}
	 $tmp_file="to_combine/${temp[0]}.exskX";unless(-e $tmp_file){goto STOPSCRIPT;}
	 $tmp_file="to_combine/${temp[0]}.MULTI3X";unless(-e $tmp_file){goto STOPSCRIPT;}
    }
    $tmp_file=0;
    STOPSCRIPT: if($tmp_file){errPrintDie "File $tmp_file does not exist. Probable reason: wrong subsample names in group-definition file.";}
}
close GROUPS;

# check if all files within each group are either all strand-specific or strand-unspecific
# and output info files for groups which are needed in RunDBS_2.pl for deciding which mappability correction should be applied
foreach my $grp (keys %groups){
	open(my $fh,">to_combine/${grp}.info");
	my ($N_ss,$N_nss)=(0,0);
	foreach my $f (keys %file_2_groups){
		foreach my $grp2 (@{$file_2_groups{$f}}){
			if($grp2 eq $grp){
				if($file_is_ss{$f}){$N_ss++;}else{$N_nss++};
			}
		}
	}
	if($N_ss>0 && $N_nss>0){
		verbPrint("Attention: group $grp consists of $N_ss strand-specific and $N_nss strand-unspecific samples. Mappability correction for this group will be done in strand-unspecific mode.\n");
		print $fh "\t\tdone"; # means strand-unspecific
	}
	if($N_ss>0 && $N_nss==0){
		$group_is_ss{$grp}=1;
		print $fh "--norc\t-SS\tdone"; # means strand-specific
	}
	if($N_ss==0 && $N_nss>0){
		print $fh "\t\tdone"; # means strand-unspecific
	}
	close($fh);
}

# output warning if one file will be merged into the same group several times; maybe this is a mistake!
foreach my $fn (keys %file_grpchk){foreach my $gr (keys %{$file_grpchk{$fn}}){
	my $num=$file_grpchk{$fn}->{$gr};
	if($num>1){verbPrint("Attention: sample $fn gets merged $num times into group $gr.\n");}
}}


### variables for merging:
my $N_expr = 0; # number of merged expression files
my $N_IR = 0;
my $N_IRsum = 0;
my $N_MIC = 0;
my $N_EEJ = 0;
my $N_MULTI = 0;
my $N_EXSK = 0;

my %eff_nss;
my %eff_ss;
my %READS_EXPR;
my %TOTAL_READS_EXPR;
my %IR;
my %IRsum;
my $IR_head;
my $IRsum_head;
my $MIC_head;
my %MIC;
my %dataMIC;
my %EEJ;
my %EEJpositions;
my %EXSK;
my $EXSK_head;
my %dataEXSK_post;
my %dataEXSK_pre;
my $MULTI_head;
my %dataMULTI_pre;
my %dataMULTI_mid;
my %dataMULTI_post;
my %MULTIa;
my %MULTIb;

### For expression
verbPrint "Doing merging for Expression files only\n" if (defined $exprONLY);

if (defined $expr){

	verbPrint "Loading Effective data\n";
	open (EFF, "$dbDir/EXPRESSION/$species"."_mRNA-50.eff") or errPrintDie "$!";
	while (<EFF>){ chomp; my @temp=split(/\t/,$_); $eff_nss{$temp[0]}=$temp[1]; } close EFF;
	open (EFF, "$dbDir/EXPRESSION/$species"."_mRNA-50-SS.eff") or errPrintDie "$!";
	while (<EFF>){ chomp; my @temp=split(/\t/,$_); $eff_ss{$temp[0]}=$temp[1]; } close EFF;
		
	verbPrint "Loading Expression files\n";
	my @files=glob("expr_out/*.cRPKM");
	foreach my $file (@files){
	    my ($root)=$file=~/.+\/(.+?)\.cRPKM/;
	    next if !$files_2b_merged{$root};
	    $N_expr++;
	    
	    verbPrint "  Processing $file\n";
	    open (I, $file) or errPrint "Can't open $file";
	    while (<I>){
		chomp;
		my @temp=split(/\t/,$_);
		my $gene=$temp[0];
		foreach my $group (@{$file_2_groups{$root}}){
			if ($temp[1] eq "NA" || $temp[1] eq "ne"){
		    		$READS_EXPR{$group}{$gene}="NA";   # XXX check with Manu's answer
			}
			else {
		    		$temp[2] = 0 if (!defined $temp[2]);
		    		$READS_EXPR{$group}{$gene}+=$temp[2];
		    		$TOTAL_READS_EXPR{$group}+=$temp[2];
			}
	    	}
	    }
	    close I;
	    system "mv $file expr_out/PARTS/" if (defined $move_to_PARTS);
	}

}

unless (defined $exprONLY){
### For IR (v1 and v2)
    if (!defined $noIR){
	verbPrint "Loading IR files (version $IR_version)\n";
	if ($IR_version == 1){
	    my @files=glob("to_combine/*.IR"); 
	    foreach my $file (@files){
		my ($root)=$file=~/.+\/(.+?)\.IR/;
		next if !$files_2b_merged{$root};
		$N_IR++; 
		
		verbPrint "  Processing $file\n";
		open (I, $file);
		$IR_head=<I>;
		while (<I>){
		    chomp;
		    my @temp=split(/\t/,$_);
		    my $event=$temp[0];
		    foreach my $group (@{$file_2_groups{$root}}){
		    	for my $i (1..$#temp){
				$IR{$group}{$event}[$i]+=$temp[$i];
		    	}
		    }
		}
		close I;
		system "mv $file to_combine/PARTS/" if (defined $move_to_PARTS);
	    }
	}
	elsif ($IR_version == 2){
	    my @files=glob("to_combine/*.IR2"); 
	    foreach my $file (@files){
		my ($root)=$file=~/.+\/(.+?)\.IR2/;
		next if !$files_2b_merged{$root};
		$N_IR++; 
		
		verbPrint "  Processing $file\n";
		open (I, $file);
		$IR_head=<I>;
		while (<I>){
		    chomp;
		    my @temp=split(/\t/,$_);
		    my $event=$temp[0];
		    foreach my $group (@{$file_2_groups{$root}}){
		    	for my $i (1..$#temp){
				$IR{$group}{$event}[$i]+=$temp[$i];
		    	}
		    }
		}
		close I;
		system "mv $file to_combine/PARTS/" if (defined $move_to_PARTS);
	    }
	    
	    verbPrint "Loading IR.summary_v2.txt files\n"; # only in v2
	    @files=glob("to_combine/*.IR.summary_v2.txt"); 
	    foreach my $file (@files){
		my ($root)=$file=~/.+\/(.+?)\.IR.summary_v2/;
		next if !$files_2b_merged{$root};
		$N_IRsum++;
		
		verbPrint "  Processing $file\n";
		open (I, $file) || errPrintDie "Can't open $file\n";
		$IRsum_head=<I>;
		while (<I>){
		    chomp($_);
		    my @temp=split(/\t/,$_);
		    my $event=$temp[0];
		    foreach my $group (@{$file_2_groups{$root}}){
		    	for my $i (1..$#temp){   # XXX check with Manu's answer
				$IRsum{$group}{$event}[$i]+=$temp[$i] unless ($temp[$i] eq "ne"); # in v2 it has 6 elements 1..6 (corr counts and raw counts)
				$IRsum{$group}{$event}[$i]=$temp[$i] if ($temp[$i] eq "ne"); # in v2 it has 6 elements 1..6 (corr counts and raw counts)
		    	}
		    }
		}
		close I;
		system "mv $file to_combine/PARTS/" if (defined $move_to_PARTS);
	    }
	}
    }
    
### For MIC
    verbPrint "Loading Microexon files\n";
    my @files=glob("to_combine/*.micX");
    foreach my $file (@files){
	my ($root)=$file=~/.+\/(.+?)\.micX/;
	next if !$files_2b_merged{$root};
	$N_MIC++;

	verbPrint "  Processing $file\n";
	open (I, $file) || errPrintDie "Can't open $file\n";
	$MIC_head=<I>;
	while (<I>){
	    chomp($_);
	    my @temp=split(/\t/,$_);
	    my $event=$temp[1];
	    $dataMIC{$event}=join("\t",@temp[0..5]);
	    foreach my $group (@{$file_2_groups{$root}}){
	    	for my $i (7..$#temp){ # Raw_reads_exc  Raw_reads_inc  Corr_reads_exc  Corr_reads_inc
			$MIC{$group}{$event}[$i]+=$temp[$i] if $temp[$i] ne "NA";
			$MIC{$group}{$event}[$i]="NA" if $temp[$i] eq "NA";    # XXX check with Manu's answer
	    	}
	    }
	}
	### Needs to recalculate PSI: from 9 and 10
	close I;
	system "mv $file to_combine/PARTS/" if (defined $move_to_PARTS);
    }
    
### For EEJ2
    verbPrint "Loading eej2 files\n";
    @files=glob("to_combine/*.eej2");
    foreach my $file (@files){
	my ($root)=$file=~/.+\/(.+?)\.eej2/;
	next if !$files_2b_merged{$root};
	$N_EEJ++;

	verbPrint "  Processing $file\n";
	open (I, $file) or errPrintDie "Can't open $file\n";
	while (<I>){
	    chomp($_);
	    my @temp=split(/\t/,$_);
	    my $event="$temp[0]\t$temp[1]";
	    foreach my $group (@{$file_2_groups{$root}}){ $EEJ{$group}{$event}+=$temp[2]; }
	    
	    my $tot_control=0;
	    my @positions=split(/\,/,$temp[4]);
	    foreach my $pos (@positions){
		my ($p,$n)=$pos=~/(\d+?)\:(\d+)/;
		foreach my $group (@{$file_2_groups{$root}}){ $EEJpositions{$group}{$event}[$p]+=$n; }
		$tot_control+=$n;
	    }
	    errPrint "Sum of positions ne total provided for $event in $file\n" if $tot_control ne $temp[2];
	}
	close I;
	system "mv $file to_combine/PARTS/" if (defined $move_to_PARTS);
    }
    
### For EXSK
    verbPrint "Loading EXSK files\n";
    @files=glob("to_combine/*.exskX");
    foreach my $file (@files){
	my ($root)=$file=~/.+\/(.+?)\.exskX/;
	next if !$files_2b_merged{$root};
	$N_EXSK++;
	
	verbPrint "  Processing $file\n";
	open (I, $file) or errPrintDie "Can't open $file\n";
	$EXSK_head=<I>;
	while (<I>){
	    chomp($_);
	    my @temp=split(/\t/,$_);
	    my $event=$temp[3];
	    $temp[25]="" if (!defined $temp[25]);
	    $dataEXSK_pre{$event}=join("\t",@temp[0..11]);
	    $dataEXSK_post{$event}=join("\t",@temp[22..25]);
	    
	    foreach my $group (@{$file_2_groups{$root}}){
	    	for my $i (13..16,19..21){ # PSI  Reads_exc  Reads_inc1  Reads_inc2  Sum_of_reads  .  Complexity  Corrected_Exc  Corrected_Inc1  Corrected_Inc2
			# only 13-16 and 19-21 really
			$EXSK{$group}{$event}[$i]+=$temp[$i] if $temp[$i] ne "NA";
			$EXSK{$group}{$event}[$i]="NA" if $temp[$i] eq "NA";
		}
	    }
	}
	### Needs to recalculate PSI: from 20+21 and 19
	close I;
	system "mv $file to_combine/PARTS/" if (defined $move_to_PARTS);
    }
    
### For MULTI
    verbPrint "Loading MULTI files\n";
    @files=glob("to_combine/*.MULTI3X");
    foreach my $file (@files){
	my ($root)=$file=~/.+\/(.+?)\.MULTI3X/;
	next if !$files_2b_merged{$root};
	$N_MULTI++;

	verbPrint "  Processing $file\n";
	open (I, $file) or errPrintDie "Can't open $file\n";
	$MULTI_head=<I>;
	while (<I>){
	    chomp($_);
	    my @temp=split(/\t/,$_);
	    my $event=$temp[3];
	    # fixed data for each event (it basically overwrites it each time)
	    $temp[25]="" if (!defined $temp[25]); # no name available
	    $dataMULTI_pre{$event}=join("\t",@temp[0..11]);
	    $dataMULTI_mid{$event}=join("\t",@temp[17..18]);
	    $dataMULTI_post{$event}=join("\t",@temp[23..25]);
	    
	    foreach my $group (@{$file_2_groups{$root}}){
	    	for my $i (13..16){ # only sum of raw reads 
		    unless (defined($MULTIa{$group}{$event}[$i]) && $MULTIa{$group}{$event}[$i] eq "NA"){
			$MULTIa{$group}{$event}[$i]+=$temp[$i] if $temp[$i] ne "NA" && (defined $temp[$i]);
			$MULTIa{$group}{$event}[$i]="NA" if $temp[$i] eq "NA";
		    }
	    	}
	    	for my $i (19..21){ 
		    my ($a,$b,$c)=$temp[$i]=~/(.*?)\=(.*?)\=(.*)/;
		    unless (defined($MULTIb{$group}{$event}[$i][0]) && $MULTIb{$group}{$event}[$i][0] eq "NA"){
			$MULTIb{$group}{$event}[$i][0]+=$a if $a ne "NA" && $a=~/\d/;
			$MULTIb{$group}{$event}[$i][0]="NA" if $a eq "NA";
		    }
		    unless (defined($MULTIb{$group}{$event}[$i][1]) && $MULTIb{$group}{$event}[$i][1] eq "NA"){
			$MULTIb{$group}{$event}[$i][1]+=$b if $b ne "NA" && $b=~/\d/;
			$MULTIb{$group}{$event}[$i][1]="NA" if $b eq "NA";
		    }
		    unless (defined($MULTIb{$group}{$event}[$i][2]) && $MULTIb{$group}{$event}[$i][2] eq "NA"){
			$MULTIb{$group}{$event}[$i][2]+=$c if $c ne "NA" && $c=~/\d/;
			$MULTIb{$group}{$event}[$i][2]="NA" if $c eq "NA";
		    }
	    	}
	    }
	}
	### Needs to recalculate PSI: from 20a+21a and 19a
	close I;
	system "mv $file to_combine/PARTS/" if (defined $move_to_PARTS);
    }

### Doing sample number check
    errPrintDie "Different number of samples in each Cassette module\n" if $N_EXSK != $N_MULTI || $N_EXSK != $N_MIC || $N_EXSK != $N_EEJ;
    errPrintDie "Different number of samples in each IR file\n" if ($N_IR != $N_IRsum && $IR_version == 2 && !defined $noIR);
    verbPrint "Warning: Number of IR samples doesn't match those of other events\n" if ($N_IR != $N_EXSK && !defined $noIR);
    verbPrint "Warning: Number of EXPR samples ($N_expr) doesn't match those of other events ($N_EXSK)\n" if $N_expr != $N_EXSK && (defined $expr);
} 


# move info files of subsamples which have been merged into at least one group into subfolder PARTS
if (defined $move_to_PARTS){
	foreach my $infof (keys %info_files_of_subsamples){
    		system "mv to_combine/$infof to_combine/PARTS/";
	}
}

### Print output files
verbPrint "Printing group files\n";
foreach my $group (sort keys %groups){
    verbPrint ">>> $group\n";
    ### EXPR
    if (defined $expr){
	unless (-e "expr_out/$group.cRPKM"){
	    my $eff_href;
	    open (EXPR, ">expr_out/$group.cRPKM") || errPrintDie "Cannot open output file";
	    if($group_is_ss{$group}){$eff_href=\%eff_ss}else{$eff_href=\%eff_nss}
	    foreach my $g (sort keys %{$READS_EXPR{$group}}){
		my $cRPKM = "";
		if ($READS_EXPR{$group}{$g} eq "NA" || $eff_href->{$g}==0){
		    $cRPKM="NA";
		}
		else {
		    $cRPKM=sprintf("%.2f",1000000*(1000*$READS_EXPR{$group}{$g}/$eff_href->{$g})/$TOTAL_READS_EXPR{$group}); 
		}
		print EXPR "$g\t$cRPKM\t$READS_EXPR{$group}{$g}\n";
	    }
	    close EXPR;
	}
    }
    unless (defined $exprONLY){
	### IR
	if (!defined $noIR){
	    if ($IR_version == 1){
	    	unless (-e "to_combine/$group.IR"){
		    open (IR, ">to_combine/$group.IR") || errPrintDie "Cannot open IR output file";
		    print IR "$IR_head";
		    foreach my $ev (sort keys %{$IR{$group}}){
			print IR "$ev\t$IR{$group}{$ev}[1]\t$IR{$group}{$ev}[2]\t$IR{$group}{$ev}[3]\t$IR{$group}{$ev}[4]\n";
		    }
		    close IR;
	    	}
	    }
	    elsif ($IR_version == 2){
		unless (-e "to_combine/$group.IR2"){
		    open (IR, ">to_combine/$group.IR2") || errPrintDie "Cannot open IR output file";
		    print IR "$IR_head";
		    foreach my $ev (sort keys %{$IR{$group}}){
			print IR "$ev\t$IR{$group}{$ev}[1]\t$IR{$group}{$ev}[2]\t$IR{$group}{$ev}[3]\t$IR{$group}{$ev}[4]\n";
		    }
		    close IR;
		}
		### IRsum (only v2)
		unless (-e "to_combine/$group.IR.summary_v2.txt"){
		    open (IRsum, ">to_combine/$group.IR.summary_v2.txt") || errPrintDie "Cannot open IRsum output file";
		    print IRsum "$IRsum_head";
		    foreach my $ev (sort keys %{$IRsum{$group}}){
			print IRsum "$ev\t$IRsum{$group}{$ev}[1]\t$IRsum{$group}{$ev}[2]\t$IRsum{$group}{$ev}[3]\t$IRsum{$group}{$ev}[4]\t$IRsum{$group}{$ev}[5]\t$IRsum{$group}{$ev}[6]\n";
		    }
		    close IRsum;
		}
	    }
	}
	### MIC
	unless (-e "to_combine/$group.micX"){
	    open (MIC, ">to_combine/$group.micX") || errPrintDie "Cannot open micX output file";
	    print MIC "$MIC_head";
	    foreach my $ev (sort keys %{$MIC{$group}}){
		my $PSI_MIC_new = "";
		if ($MIC{$group}{$ev}[10] ne "NA" && $MIC{$group}{$ev}[9] ne "NA"){
		    if (($MIC{$group}{$ev}[10]+$MIC{$group}{$ev}[9])>0){
			$PSI_MIC_new=sprintf("%.2f",100*$MIC{$group}{$ev}[10]/($MIC{$group}{$ev}[10]+$MIC{$group}{$ev}[9]));
		    }
		    else {
			$PSI_MIC_new="NA";
		    }
		}
		else {
		    $PSI_MIC_new="NA";
		}
		print MIC "$dataMIC{$ev}\t$PSI_MIC_new\t$MIC{$group}{$ev}[7]\t$MIC{$group}{$ev}[8]\t$MIC{$group}{$ev}[9]\t$MIC{$group}{$ev}[10]\n";
	    }
	    close MIC;
	}
	### EXSK
	unless (-e "to_combine/$group.exskX"){
	    open (EXSK, ">to_combine/$group.exskX") || errPrintDie "Cannot open exskX output file";
	    print EXSK "$EXSK_head";
	    foreach my $ev (sort keys %{$EXSK{$group}}){
		my $PSI_EXSK_new = "";
		if (($EXSK{$group}{$ev}[19]+($EXSK{$group}{$ev}[20]+$EXSK{$group}{$ev}[21])/2)>0){
		    $PSI_EXSK_new=sprintf("%.2f",100*($EXSK{$group}{$ev}[20]+$EXSK{$group}{$ev}[21])/(($EXSK{$group}{$ev}[20]+$EXSK{$group}{$ev}[21])+2*$EXSK{$group}{$ev}[19]));
		}
		else {
		    $PSI_EXSK_new="NA";
		}
		print EXSK "$dataEXSK_pre{$ev}\t$PSI_EXSK_new\t$EXSK{$group}{$ev}[13]\t$EXSK{$group}{$ev}[14]".
		    "\t$EXSK{$group}{$ev}[15]\t$EXSK{$group}{$ev}[16]\t.\tS\t$EXSK{$group}{$ev}[19]\t$EXSK{$group}{$ev}[20]".
		    "\t$EXSK{$group}{$ev}[21]\t$dataEXSK_post{$ev}\n";
	    }
	    close EXSK;
	}
	### MULTI
	unless (-e "to_combine/$group.MULTI3X"){
	    open (MULTI, ">to_combine/$group.MULTI3X") || errPrintDie "Cannot open MULTI3X output file";
	    print MULTI "$MULTI_head";
	    foreach my $ev (sort keys %{$MULTIa{$group}}){
		my $PSI_MULTI_new = "";

		# if any junction has an NA (i.e. mappability == 0), then the PSI is killed
		for my $loop1 (19..21){ # 19 => excl, 20 => inc1, 21 => inc2
		    for my $loop2 (0..2){ # 0 => corr_all, 1 => raw_ref, 2 => corr_ref
			$MULTIb{$group}{$ev}[$loop1][$loop2] = 0 if (!defined $MULTIb{$group}{$ev}[$loop1][$loop2]);
			$PSI_MULTI_new = "NA" if $MULTIb{$group}{$ev}[$loop1][$loop2] eq "NA";
		    }
		}
		
		if ($PSI_MULTI_new eq "NA"){}
		elsif ((($MULTIb{$group}{$ev}[20][0]+$MULTIb{$group}{$ev}[21][0])+2*$MULTIb{$group}{$ev}[19][0])>0){
		    $PSI_MULTI_new=sprintf("%.2f",100*($MULTIb{$group}{$ev}[20][0]+$MULTIb{$group}{$ev}[21][0])/(($MULTIb{$group}{$ev}[20][0]+$MULTIb{$group}{$ev}[21][0])+2*$MULTIb{$group}{$ev}[19][0]));
		}
		else {
		    $PSI_MULTI_new="NA";
		}
		
		### Recalculates complexity
		$MULTIb{$group}{$ev}[19][2]=0 if (!defined $MULTIb{$group}{$ev}[19][2]);
		$MULTIb{$group}{$ev}[20][2]=0 if (!defined $MULTIb{$group}{$ev}[20][2]);
		$MULTIb{$group}{$ev}[21][2]=0 if (!defined $MULTIb{$group}{$ev}[21][2]);	    
		
		my $Q;
		if ($PSI_MULTI_new ne "NA"){
		    my $from_S=$MULTIb{$group}{$ev}[19][2]+$MULTIb{$group}{$ev}[20][2]+$MULTIb{$group}{$ev}[21][2]; # reads coming only from the reference EEJs (refI1, refI2 and refE)
		    my $from_C=($MULTIb{$group}{$ev}[19][0]+$MULTIb{$group}{$ev}[20][0]+$MULTIb{$group}{$ev}[21][0])-$from_S; # all other reads
		    if ($from_C > ($from_C+$from_S)/2) {$Q="C3";}
		    elsif ($from_C > ($from_C+$from_S)/5 && $from_C <= ($from_C+$from_S)/2){$Q="C2";}
		    elsif ($from_C > ($from_C+$from_S)/20 && $from_C <= ($from_C+$from_S)/5){$Q="C1";}
		    else {$Q="S";}
		}
		else {
		    $Q="NA";
		}
		
		$MULTIa{$group}{$ev}[13]="" if (!defined $MULTIa{$group}{$ev}[13]);
		$MULTIa{$group}{$ev}[14]="" if (!defined $MULTIa{$group}{$ev}[14]);
		$MULTIa{$group}{$ev}[15]="" if (!defined $MULTIa{$group}{$ev}[15]);	    
		$MULTIa{$group}{$ev}[16]="" if (!defined $MULTIa{$group}{$ev}[16]);
		$dataMULTI_pre{$ev}="" if (!defined $dataMULTI_pre{$ev});
		$dataMULTI_mid{$ev}="" if (!defined $dataMULTI_mid{$ev});
		$MULTIb{$group}{$ev}[19][1]=0 if (!defined $MULTIb{$group}{$ev}[19][1]);
		$MULTIb{$group}{$ev}[20][1]=0 if (!defined $MULTIb{$group}{$ev}[20][1]);
		$MULTIb{$group}{$ev}[21][1]=0 if (!defined $MULTIb{$group}{$ev}[21][1]);
		
		print MULTI "$dataMULTI_pre{$ev}\t$PSI_MULTI_new\t$MULTIa{$group}{$ev}[13]\t$MULTIa{$group}{$ev}[14]\t$MULTIa{$group}{$ev}[15]\t$MULTIa{$group}{$ev}[16]\t$dataMULTI_mid{$ev}\t".
		    "$MULTIb{$group}{$ev}[19][0]=$MULTIb{$group}{$ev}[19][1]=$MULTIb{$group}{$ev}[19][2]\t".
		    "$MULTIb{$group}{$ev}[20][0]=$MULTIb{$group}{$ev}[20][1]=$MULTIb{$group}{$ev}[20][2]\t".
		    "$MULTIb{$group}{$ev}[21][0]=$MULTIb{$group}{$ev}[21][1]=$MULTIb{$group}{$ev}[21][2]\t".
		    "$Q\t$dataMULTI_post{$ev}\n";
	    }
	    close MULTI;
	}
	### EEJ2
	unless (-e "to_combine/$group.eej2"){
	    open (EEJ2, ">to_combine/$group.eej2") || errPrintDie "Cannot open eej2 output file";
	    foreach my $ev (sort keys %{$EEJ{$group}}){
		my $pos="";
		for my $i (0..$#{$EEJpositions{$group}{$ev}}){
		    $pos.="$i:$EEJpositions{$group}{$ev}[$i]," if $EEJpositions{$group}{$ev}[$i];
		}
		chop($pos);
		print EEJ2 "$ev\t$EEJ{$group}{$ev}\tNA\t$pos\n";
	    }
	}
    }
}

verbPrint "Merge finished successfully\n";


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
    $assembly_to_species{strPur4}="Spu"; $assembly_to_species{dm6}="Dme"; $assembly_to_species{aedAeg5}="Aae";
    $assembly_to_species{bomMor1}="Bmo"; $assembly_to_species{triCas5}="Tca"; $assembly_to_species{apiMel4}="Ame";
    $assembly_to_species{blaGer1}="Bge"; $assembly_to_species{cloDip2}="Cdi"; $assembly_to_species{strMar1}="Sma";
    $assembly_to_species{ce11}="Cel"; $assembly_to_species{octBim1}="Obi"; $assembly_to_species{octMin1}="Omi";
    $assembly_to_species{schMed31}="Sme"; $assembly_to_species{nemVec1}="Nve"; $assembly_to_species{araTha10}="Ath";
    my %species_to_assembly;
    $species_to_assembly{Hsa}="hg19"; $species_to_assembly{Hs2}="hg38"; $species_to_assembly{Ptr}="panTro4";
    $species_to_assembly{Mma}="rheMac2"; $species_to_assembly{Mmu}="mm9"; $species_to_assembly{Mm2}="mm10";
    $species_to_assembly{Bta}="bosTau6"; $species_to_assembly{Bt2}="bosTau9"; $species_to_assembly{Mdo}="monDom5";
    $species_to_assembly{Gg3}="galGal3"; $species_to_assembly{Gg4}="galGal4"; $species_to_assembly{Gga}="galGal6";
    $species_to_assembly{Xt1}="xenTro3"; $species_to_assembly{Xtr}="xenTro9";
    $species_to_assembly{Dre}="danRer10"; $species_to_assembly{Dr2}="danRer11";
    $species_to_assembly{Loc}="lepOcu1"; $species_to_assembly{Elu}="esoLuc2"; $species_to_assembly{Cm1}="eshark1";
    $species_to_assembly{Cmi}="calMil1"; $species_to_assembly{Bl1}="braLan2"; $species_to_assembly{Bla}="braLan3";
    $species_to_assembly{Spu}="strPur4"; $species_to_assembly{Dme}="dm6"; $species_to_assembly{Aae}="aedAeg5";
    $species_to_assembly{Bmo}="bomMor1"; $species_to_assembly{Tca}="triCas5"; $species_to_assembly{Ame}="apiMel4";
    $species_to_assembly{Bge}="blaGer1"; $species_to_assembly{Cdi}="cloDip2"; $species_to_assembly{Sma}="strMar1";
    $species_to_assembly{Cel}="ce11"; $species_to_assembly{Obi}="octBim1"; $species_to_assembly{Omi}="octMin1";
    $species_to_assembly{Sme}="schMed31"; $species_to_assembly{Nve}="nemVec1"; $species_to_assembly{Ath}="araTha10";

    if (defined $assembly_to_species{$temp_assembly[0]}){ # it's a proper assembly
	$species = $assembly_to_species{$temp_assembly[0]};
    }
    else {
	if (defined $species_to_assembly{$temp_assembly[0]}){
	    $sp_assembly = $species_to_assembly{$temp_assembly[0]};
	    $species = $temp_assembly[0];
	}
	else {
	    errPrintDie "$temp_assembly[0] is not a valid species\n";
	}
    }
}
