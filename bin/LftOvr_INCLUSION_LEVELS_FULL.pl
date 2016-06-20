#!/usr/bin/perl
use strict;
use warnings;
use File::Temp qw(tempfile);
use Storable;  #store \%table, 'file';  $hashref = retrieve('file');
use File::Path 'rmtree';



#######
## GLOBAL VARS
#######
my $output_new_exon_borders=1;  # if 1, all exon borders of the translated events during construction of dictionary will be written to STDOUT
my $elim_coord_replicates=1;

#######
## SUBS
#######
sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

sub sort_num_up_a{
	my $aref=$_[0];
	if(@$aref<1){return([]);}
	#foreach my $v (@$aref){if(!is_integer($v)){die;}}
	my @ret = sort { $a <=> $b } @{ $aref };
return(\@ret);	
}

sub sort_num_down_a{
	my $aref=$_[0];
	if(@$aref<1){return([]);}
	my @ret = sort { $b <=> $a } @{ $aref };
return(\@ret);	
}

sub get_min_a{
	my $aref=$_[0];
	my $min_v=9**9**9;
	foreach my $v (@{$aref}){if($v<$min_v){$min_v=$v;}}
return($min_v);
}

sub get_max_a{
	my $aref=$_[0];
	my $max_v=-9**9**9;
	foreach my $v (@{$aref}){if($v>$max_v){$max_v=$v;}}
return($max_v);
}

########
## CALL
########

# lft_vts_annots.pl build_dict vast_out.tab liftOver_chain_map dict_file translations_file lost_events_file all_new_exon_border_file
# lft_vts_annots.pl translate vast_out.tab dict_file vast_out_translated.tab




########
## MAIN
########
my ($cmd,$dict_f,$vts_tab,$lftOvr_chain_map,$lost_events_file,$vts_tab_translated,$events_file,$new_ex_borders);

my ($tmp_str,$chr,$event_id,$strand,$coord,$start,$end,$coord_id,$coord_old_1,$coord_old_2,$coord_new_1,$coord_new_2,$fh,$fh_lost,$fh_out,$new_fields,$rest,$as_type,$coord_new,$coord_col,$fullco_col,$fullco_old);

my (@fs,@fs2,@fs2_new,@fs3,@fs3_new,@fs4,@fs4_new,@fs5,@fs5_new,@fs6,@fs6_new);

$cmd=$ARGV[0];
if($cmd eq "build_dict"){
	$vts_tab=$ARGV[1];
	$lftOvr_chain_map=$ARGV[2];
	$dict_f=$ARGV[3];
	$lost_events_file=$ARGV[5];
	$events_file=$ARGV[4];
	$new_ex_borders=$ARGV[6];
	

	# hash with all coordinates which should be translated
	# key=COORD_ID: CHROMOSOME:COORDINATE value: COORDINATE
	my %coords=();

	# hash for column COORD
	# key: rowkey value: string with three sub-strings CHR-COORD_ID-COORD_ID
	my %COORD=();
	my %FULLCOORD=();

	# get all coordinates from vast-tools table and write them into a gtf file
	open($fh,"<".$vts_tab) or die $!;  <$fh>;
	while(<$fh>){
		@fs=split("\t");
		$event_id=$fs[1];
		# Coord column
		($chr,$start,$end) = $fs[2] =~ /(.*):(.*)-(.*)/;
		unless($start eq ""){$coords{"$chr:$start"}=$start;}unless($end eq ""){$coords{"$chr:$end"}=$end;}
		$COORD{$event_id}="$chr-$start-$end";
		#print "|$chr| |$start| |$end|\n";

		# FullCo column chr19:36243988-36244094+36244154,36244918  chr17:3991972-3992187=3990728-3990828:- (IR)
		# we cut off the last 2 characters in case of IR to have all coordinates more similar
		$tmp_str=$fs[4];
		if(substr($tmp_str,length($tmp_str)-2,1) eq ":"){$tmp_str=substr($tmp_str,0,length($tmp_str)-2);}
		($chr,$rest)=split(":",$tmp_str);
		my @all_coords = $rest =~ /(\d+)/g;
		foreach $coord (@all_coords){$coords{"$chr:$coord"}=$coord;}
		$FULLCOORD{$event_id}=$fs[5].":".$tmp_str;
	}
	close($fh);

	# write temporary file with coordinates to be translated, call liftOver, read-in translated coordinates
	# create temporary file; will be deleted automatically at the end of script or if error occurs
	unless(-d "./tmp_lftOvr"){mkdir("./tmp_lftOvr");}
	my ($fh_tmp,$tmpf_nm) = tempfile( DIR => "./tmp_lftOvr", UNLINK => 1);
	while ( my ($coord_id, $coord) = each %coords ){($chr,$coord)=split(":",$coord_id);print $fh_tmp "$chr\t$chr\t$coord_id\t$coord\t$coord\t.\t.\t.\n";}close($fh_tmp);
	my $lft_ovr_output=`liftOver -gff $tmpf_nm $lftOvr_chain_map ${tmpf_nm}_transl ${tmpf_nm}_unmapped`;

	%coords=();
	open($fh,"<"."${tmpf_nm}_transl") or die $!;
	while(<$fh>){
		my @fs=split("\t");
		# chromosome has changed -> neglect this case
		if($fs[0] ne $fs[1]){next;}
		$coords{$fs[2]}=$fs[3];
	};close($fh);

	# only used if new exon borders should be output
	my %new_ex_borders;

	open(my $fh_trans,">".$events_file) or die "$!";
	
	# collect event ids of event which could not be translated
	my %lost_events=();
	my %DICT=();
	my $skip_this_event;
	# good to now: $s="4-5"; split("-",$s) produces (4,5), $s="4-" produces (4), $s="-5" produces ("",5)
	# go over COORD column and translate coordinate

	my $coord_counter=0;
	EVENTLOOP: foreach my $event_id (keys %COORD){
		$skip_this_event=0;$coord_col="";$fullco_col="";

		# translate COORD column
		# empty coordinates are eq "" like in chr3:-472628  or in chr7:34987-
		($chr,$coord_old_1,$coord_old_2)=split("-",$COORD{$event_id});
		
		($coord_new_1,$coord_new_2)=("",""); unless($coord_old_1 eq ""){$coord_new_1=$coords{"$chr:$coord_old_1"};} unless($coord_old_2 eq ""){$coord_new_2=$coords{"$chr:$coord_old_2"};}

		# old coordinates could not be translated into new coordinates
		if(!defined($coord_new_1) || !defined($coord_new_2) ){$skip_this_event=1;goto END_OF_EVENTLOOP;}

		# check if region still has same length
		if($coord_old_1 ne "" && $coord_old_2 ne ""){
			# length of region should be the same with old and new coordinates
			if($coord_old_2-$coord_old_1 != $coord_new_2-$coord_new_1){$skip_this_event=1;goto END_OF_EVENTLOOP;}
		}

		# event id and COORD column translated into new coordinates
		$coord_col="$chr:$coord_new_1-$coord_new_2";
		#################


		# translate FULLCO column
		$fullco_old=$FULLCOORD{$event_id};
		($as_type,$chr,$rest)=split(":",$fullco_old);		

		# TAKE CARE: the chromosome annotation of the examples actually does not exists in the $rest variable, e.g., for chr17:56424894+56424890-56424945,56424602 $rest=56424894+56424890-56424945,56424602
		if($as_type =~ /Alt5/){ # chr17:56424894+56424890-56424945,56424602  chr17:3+2-4,1 -> case 1 (- strand)
					            # chr7:50571195-50571267+50571352,50571803   chr7:1-2+3,4 -> case 2 (+ strand)
					            # we also have  Alt5:chr1:247267499+247267241-,247265415 
			@fs=split(",",$rest);
			$coord_new=$coords{"$chr:$fs[1]"};if(!defined($coord_new)){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			$fs[1]=$coord_new;

			@fs2=split("-",$fs[0]);
		
			@fs3_new=();
			@fs4_new=();
			if(@fs2==2){
				# -50571267+50571352
				unless($fs2[0] eq ""){
					@fs3=split("\\+",$fs2[0]);
					for(my $i=0;$i<@fs3;$i++){$coord_new=$coords{"$chr:$fs3[$i]"};if(defined($coord_new)){push(@fs3_new,$coord_new);}}
					if(@fs3==1 && @fs3_new<1 || @fs3>1 && @fs3_new<2){;$skip_this_event=1;goto END_OF_EVENTLOOP;}
				}
				
				@fs4=split("\\+",$fs2[1]);
				if(@fs4==0){$fs4_new[0]="";}
				for(my $i=0;$i<@fs4;$i++){$coord_new=$coords{"$chr:$fs4[$i]"};if(defined($coord_new)){push(@fs4_new,$coord_new);}}
				if(@fs4==1 && @fs4_new<1 || @fs4>1 && @fs4_new<2){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			}
			
			# 50571267+50571352-
			if(@fs2==1){
				@fs3=split("\\+",$fs2[0]);
				for(my $i=0;$i<@fs3;$i++){$coord_new=$coords{"$chr:$fs3[$i]"};if(defined($coord_new)){push(@fs3_new,$coord_new);}}
				if(@fs3==1 && @fs3_new<1 || @fs3>1 && @fs3_new<2){print "skip2\n";$skip_this_event=1;goto END_OF_EVENTLOOP;}
			}
			
			$fullco_col="$chr:".join("+",@{sort_num_up_a(\@fs3_new)})."-".join("+",@{sort_num_up_a(\@fs4_new)}).",".$fs[1];
			
			if($output_new_exon_borders){
				foreach my $v (@fs3_new,@fs4_new,$fs[1]){$new_ex_borders{"alt5:".($coord_counter++).":$chr:$v"}=1;}  
			}

		}elsif($as_type =~ /Alt3/){ # chr8:1844615,1846602+1846599-1846691                ->  1,3+2-4    (case 1)
					                # chr8:17752734,17749190-17749338+17749316+17749309   ->  5,1-4+3+2  (case 2)
					                # we also have  Alt3:chr15:85200407,-85198640+85198635
					                
			@fs=split(",",$rest);

			$coord_new=$coords{"$chr:$fs[0]"};if(!defined($coord_new)){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			$fs[0]=$coord_new;

			@fs2=split("-",$fs[1]);
			
			@fs3_new=();
			@fs4_new=();
			if(@fs2==2){
				# -50571267+50571352
				unless($fs2[0] eq ""){
					@fs3=split("\\+",$fs2[0]);
					for(my $i=0;$i<@fs3;$i++){$coord_new=$coords{"$chr:$fs3[$i]"};if(defined($coord_new)){push(@fs3_new,$coord_new);}}
					if(@fs3==1 && @fs3_new<1 || @fs3>1 && @fs3_new<2){$skip_this_event=1;goto END_OF_EVENTLOOP;}
				}
				
				@fs4=split("\\+",$fs2[1]);
				if(@fs4==0){$fs4_new[0]="";}
				for(my $i=0;$i<@fs4;$i++){$coord_new=$coords{"$chr:$fs4[$i]"};if(defined($coord_new)){push(@fs4_new,$coord_new);}}
				if(@fs4==1 && @fs4_new<1 || @fs4>1 && @fs4_new<2){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			}
			
			# 50571267+50571352-
			if(@fs2==1){
				@fs3=split("\\+",$fs2[0]);
				for(my $i=0;$i<@fs3;$i++){$coord_new=$coords{"$chr:$fs3[$i]"};if(defined($coord_new)){push(@fs3_new,$coord_new);}}
				if(@fs3==1 && @fs3_new<1 || @fs3>1 && @fs3_new<2){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			}
			
			$fullco_col="$chr:".$fs[0].",".join("+",@{sort_num_up_a(\@fs3_new)})."-".join("+",@{sort_num_up_a(\@fs4_new)});
			
			if($output_new_exon_borders){
				foreach my $v (@fs3_new,@fs4_new,$fs[0]){$new_ex_borders{"alt3:".($coord_counter++).":$chr:$v"}=1;}  
			}

		}elsif($as_type =~ /IR/){  # 108243001-108243113=108234570-108234631    ->   3-4=1-2:-
					               # 169706028-169706147=169710382-169716161    ->   1-2=3-4:+
					               # rember: chr and strand have been cut off before
			@fs=split("=",$rest);
			@fs2=split("-",$fs[0]);

			@fs2_new=();
			for(my $i=0;$i<@fs2;$i++){$coord_new=$coords{"$chr:$fs2[$i]"};if(defined($coord_new)){push(@fs2_new,$coord_new);}}
			# at least two coordinates must survive
			if(@fs2_new<2){$skip_this_event=1;goto END_OF_EVENTLOOP;}

			@fs3=split("-",$fs[1]);
			@fs3_new=();
			for(my $i=0;$i<@fs3;$i++){$coord_new=$coords{"$chr:$fs3[$i]"};if(defined($coord_new)){push(@fs3_new,$coord_new);}}
			# at least two coordinates must survive
			if(@fs3_new<2){$skip_this_event=1;goto END_OF_EVENTLOOP;}

			$strand="+";if($fs2[0]>$fs3[0]){$strand="-";}

			# consistency check
			if($strand eq "-"){unless($fs3_new[0]<$fs3_new[1] && $fs3_new[1]<$fs2_new[0] && $fs2_new[0]<$fs2_new[1]){$skip_this_event=1;goto END_OF_EVENTLOOP;}}
			if($strand eq "+"){unless($fs2_new[0]<$fs2_new[1] && $fs2_new[1]<$fs3_new[0] && $fs3_new[0]<$fs3_new[1]){$skip_this_event=1;goto END_OF_EVENTLOOP;}}

			$fullco_col="$chr:".join("-",@fs2_new)."=".join("-",@fs3_new).":".$strand;

			if($output_new_exon_borders){
				foreach my $v (@fs2_new,@fs3_new){$new_ex_borders{"ir:".($coord_counter++).":$chr:$v"}=1;}  
			}

		}else{# all C1...C3,S,MIC
			# chr1:100843177,100844743+100849022-100849196,100856288   ->  1,2+3-4,5   (case 1, + strand)
			# chr3:171330069,171326091-171326156+171329424,171323206   ->  5,2-3+4,1   (case 2, - strand)
			@fs=split(",",$rest);
			@fs2=split("\\+",$fs[0]);
			@fs2_new=();
			for(my $i=0;$i<@fs2;$i++){$coord_new=$coords{"$chr:$fs2[$i]"};if(defined($coord_new)){push(@fs2_new,$coord_new);}}
			if(@fs2_new<@fs2 && @fs2_new<1){$skip_this_event=1;goto END_OF_EVENTLOOP;}

			@fs3=split("\\+",$fs[2]);
			@fs3_new=();
			for(my $i=0;$i<@fs3;$i++){$coord_new=$coords{"$chr:$fs3[$i]"};if(defined($coord_new)){push(@fs3_new,$coord_new);}}
			if(@fs3_new<@fs3 && @fs3_new<1){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			
			# get strand information
			if($fs2[0]<$fs3[0]){$strand="+";}else{$strand="-";}
			
			if($strand eq "+" && get_max_a(\@fs2_new)>get_min_a(\@fs3_new)){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			if($strand eq "-" && get_min_a(\@fs2_new)<get_max_a(\@fs3_new)){$skip_this_event=1;goto END_OF_EVENTLOOP;}

			@fs4=split("-",$fs[1]);
			@fs5=split("\\+",$fs4[0]);
			@fs5_new=();
			for(my $i=0;$i<@fs5;$i++){$coord_new=$coords{"$chr:$fs5[$i]"};if(defined($coord_new)){push(@fs5_new,$coord_new);}}
			if(@fs5_new<@fs5 && @fs5_new<1){$skip_this_event=1;goto END_OF_EVENTLOOP;}

			@fs6=split("\\+",$fs4[1]);
			@fs6_new=();
			for(my $i=0;$i<@fs6;$i++){$coord_new=$coords{"$chr:$fs6[$i]"};if(defined($coord_new)){push(@fs6_new,$coord_new);}}
			if(@fs6_new<@fs6 && @fs6_new<1){$skip_this_event=1;goto END_OF_EVENTLOOP;}

			if(get_max_a(\@fs5) < get_min_a(\@fs6) && get_max_a(\@fs5_new) > get_min_a(\@fs6_new)){$skip_this_event=1;goto END_OF_EVENTLOOP;}
			if(get_min_a(\@fs5) > get_max_a(\@fs6) && get_max_a(\@fs5_new) < get_min_a(\@fs6_new)){$skip_this_event=1;goto END_OF_EVENTLOOP;}

			$fullco_col="$chr:".join("+",@{sort_num_up_a(\@fs2_new)}).",".join("+",@{sort_num_up_a(\@fs5_new)})."-".join("+",@{sort_num_up_a(\@fs6_new)}).",".join("+",@{sort_num_up_a(\@fs3_new)});
			
			if($output_new_exon_borders){
				foreach my $v (@fs2_new,@fs3_new,@fs5_new,@fs6_new){$new_ex_borders{"se:".($coord_counter++).":$chr:$v"}=1;}  
			}
				
		}#if($as_type)
		##############

		END_OF_EVENTLOOP:
		#delete($FULLCOORD{$event_id});delete($COORD{$event_id}); #XXX
		if($skip_this_event){
			$lost_events{$event_id}=1;
			next EVENTLOOP;
		}

		print $fh_trans "$fullco_old\t$fullco_col\n";

		# store translated coord and fullco
		$DICT{$event_id}=$coord_col."|".$fullco_col;
		
		# test output
		#print "$COORD{$event_id} -> $coord_col\n";
		#print "$FULLCOORD{$event_id} -> $fullco_col\n";
	}
	
	store \%DICT, "$dict_f";
	
	open($fh_lost,">".$lost_events_file) or die $!;
	open($fh,"<".$vts_tab) or die $!;
	# header
	$tmp_str=<$fh>; print $fh_lost $tmp_str;
	while(<$fh>){
		@fs=split("\t");
		if($lost_events{$fs[1]}){print $fh_lost $_;}
	}
	close($fh_lost);close($fh);
	
	open($fh,">".$new_ex_borders) or die "$!";
	foreach my $new_ex_border (keys %new_ex_borders){print $fh $new_ex_border."\n";}
	close($fh);
	
	open($fh,">".$new_ex_borders) or die "$!";
	foreach my $new_ex_border (keys %new_ex_borders){print $fh $new_ex_border."\n";}
	close($fh);
	close($fh_trans);
	
	rmtree("./tmp_lftOvr");
}


if($cmd eq "translate"){
	$vts_tab=$ARGV[1];
	$dict_f=$ARGV[2];
	$vts_tab_translated=$ARGV[3];

	my $DICT_href = retrieve($dict_f);
	
	open($fh,"<".$vts_tab) or die $!;
	open($fh_out,">".$vts_tab_translated) or die $!;

	# header
	my $header=<$fh>;
	print $fh_out $header;
	while(<$fh>){
		@fs=split("\t");
		my $tmp_str=$DICT_href->{$fs[1]};
		if(!defined($tmp_str)){next;}
		
		($fs[2],$fs[4])=split("\\|",$tmp_str);

		print $fh_out join("\t",@fs);
	}

	close($fh);close($fh_out);
}
