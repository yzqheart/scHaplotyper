#!/usr/bin/perl

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use SVG;

##get options from command line into variables and set default values
my ($Verbose,$Help, $Window_size);
GetOptions(
	"window:i"=>\$Window_size,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Window_size ||= 500000;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $chr_len_file = "$Bin/../db/hg38.len.pure";
my $chr_gap_file = "$Bin/../db/hg38.genome_gap_500k";
my $chr_centromere_file = "$Bin/../db/hg38.centromere";

my $pairwise_crossovers_file = shift;
my $pairwise_hiddenseq_file = shift;
my $reference_crossover_file = shift; #yan, cluster file
my $Aneuploidy_chrs_file = shift;
my $draw_pos = shift;

my $Min_dist_two_crossover = 1000000;
my $ChrID = $1 if($pairwise_crossovers_file =~ /(chr\w+)/);

my %ChrLen;
my %ChrGap;
my %Centro;
my %CO;
my %HS;
my @AllSamples;
my $Ref_sample;
my $Sample_num;
my @RefCO;
my %Aneuploidy;

read_chrlen($chr_len_file, \%ChrLen);
read_chrgap($chr_gap_file, \%ChrGap);
read_centromere($chr_centromere_file, \%Centro);


if (-f $Aneuploidy_chrs_file) {
	open IN, $Aneuploidy_chrs_file || die "fail $Aneuploidy_chrs_file";
	while (<IN>) {
		chomp;
		my ($sample,$chr,$copy) = split /\s+/;
		$Aneuploidy{$sample}{$chr} = $copy;
	}
	close IN;
}

open IN, $pairwise_crossovers_file || die "fail $pairwise_crossovers_file";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	my $sample = $t[1];
	next if(exists $Aneuploidy{$sample}{$ChrID});
	push @{$CO{$t[1]}}, [$t[3], 1]; #yan,CO{other_sample} = [crossover_pos, 1]
}
close IN;


open IN, $pairwise_hiddenseq_file || die "fail $pairwise_hiddenseq_file";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	if (!$Ref_sample) {
		$Ref_sample = $t[0];
	}
	my $sample_id = $t[1];
	next if(exists $Aneuploidy{$sample_id}{$ChrID});
	my $pos = $t[4]; #yan, snp pos
	my $win_id = int($pos / $Window_size);
	my $hs = $t[6]; #yan, hiddenseq (F or M)
	$HS{$sample_id}{$win_id}{$hs} ++;
}
close IN;


open IN, $reference_crossover_file || die "fail $reference_crossover_file";
while (<IN>) {
	
	if(/Confident Reference Crossover:\s+(\S+)\s+(\S+)/){
		my $crossover_pos = $1;
		my $ratio = $2;
		push @RefCO, $crossover_pos;
		#print STDERR "$crossover_pos\t$ratio\n";
	}
	
}
close IN;


$Sample_num = keys %HS;
@AllSamples = sort keys %HS;
unshift @AllSamples, $Ref_sample;
my $left_edge_pxi = 50;
my $right_edge_pxi = 20;
my $top_edge_pxi = 150;
my $bottom_edge_pxi = 30;

my $Max_win_num = int($ChrLen{$ChrID} / $Window_size) + 1; ##chr1 is the largest chromosome

my $chr_x_shift = 30;
my $win_y_shift = 1;
my $figure_width = $left_edge_pxi + $right_edge_pxi + $Sample_num*$chr_x_shift;
my $figure_height = $top_edge_pxi + $bottom_edge_pxi + $win_y_shift*$Max_win_num;
my $chr_width = $chr_x_shift * 0.3;


my $svg = SVG->new('width',$figure_width,'height',$figure_height);
$svg->rect('x',0, 'y',0,'width',$figure_width,'height',$figure_height,'fill',"white");
$svg->text('x',$left_edge_pxi-40,'y',$figure_height/2,'fill','black', 'font-family','Arial','font-size',16, '-cdata',$ChrID);

my $c = 0;
my $flag=1;
my $deletion_end_pos=1;
foreach my $sample_id (@AllSamples) {
	my $chr = $ChrID;
	my $x = $left_edge_pxi + $c*$chr_x_shift;
	my $y;
	my $text_x = $x+10;
	my $text_y = $top_edge_pxi-10;
	my $transform_format = "rotate(-90 $text_x,$text_y)";
	$svg->text('x',$text_x ,'y',$text_y,'fill','black', 'font-family','Arial','font-size',12, '-cdata',$sample_id, 'transform',$transform_format);
	my $Window_num = int($ChrLen{$chr} / $Window_size) + 1;

	#yan, draw mutation site
	if ($flag == 1){		
		
		my $mut_pos = $draw_pos;
        $y = $top_edge_pxi + int($mut_pos/$Window_size)*$win_y_shift;
        my $x1 = $left_edge_pxi;
        my $x2 = $figure_width - $right_edge_pxi;
        $svg->line('x1',$x1,'y1',$y,'x2',$x2,'y2',$y,'stroke','blue','stroke-width',2);

	}
	
	
	##draw centromere (draw on the top of figure)
	my ($centro_start,$centro_end) = ($1,$2) if($Centro{$chr} =~ /(\d+),(\d+)/);
	my $centro_pos = int(($centro_start + $centro_end) / 2);
	$y = $top_edge_pxi + int($centro_pos/$Window_size)*$win_y_shift;
	$svg->circle('cx',$x+$chr_width/2,'cy',$y,'r',$chr_width*0.6,'stroke','black','fill','white','stroke-width',3);
	
		
	##draw haplotype color
	my $hidden_p = $HS{$sample_id};
	for (my $win_id = 0; $win_id < $Window_num; $win_id ++) {
		my $F_num = $hidden_p->{$win_id}{"F"} + 1;
		my $M_num = $hidden_p->{$win_id}{"M"} + 1;
		my $total_F_M_num = $F_num + $M_num;
		my $color;
		if ( ($F_num >= $M_num && $F_num >= 2) || ($M_num > $F_num && $M_num >= 2) ) {
			my $green_intensity = int(255 * ($M_num / $total_F_M_num));
			my $red_intensity = int(255 * ($F_num / $total_F_M_num));

			#my $red_intensity = int(255 * ($M_num / $total_F_M_num)); # for FigS2
			#my $green_intensity = int(255 * ($F_num / $total_F_M_num)); # for FigS2
			$color = "rgb($red_intensity,$green_intensity,0)";
		}elsif ($win_id*$Window_size >= $deletion_end_pos && $flag == 1 ){
			$color = "rgb(255,0,0)";

		}else{
			$color = "rgb(127,127,0)";
			
		}
		
		$y = $top_edge_pxi + $win_id*$win_y_shift;
		$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$win_y_shift,'stroke',"none", "fill", $color);
	}
	$flag++;


	##draw gap region in the chromosome
	my $gap_p = $ChrGap{$chr};
	for (my $win_id = 0; $win_id < $Window_num; $win_id++) {
		my $gap_ratio = 1 - $gap_p->[$win_id] / $Window_size;
		my $is_gap = 0;
		if($gap_ratio > 0.9){  ##超过90%为gap的认为整个window为gap,无信号
			$y = $top_edge_pxi + $win_id*$win_y_shift;
			$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$win_y_shift,'stroke',"none", "fill", "black");
			$is_gap = 1;
		}
		my $ratio_percent = int($gap_ratio * 100);
	}


	##draw chromosome frame
	$y = $top_edge_pxi;
	$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$Window_num*$win_y_shift,'stroke',"black", "fill", "none");

	$c ++;
}


##output the svg and png figures
my $svg_file = "$pairwise_crossovers_file.svg";
print STDERR $svg_file."\n";

open OUT, ">$svg_file" || die "fail create $svg_file\n";
print OUT $svg->xmlify();
close OUT;


####################################################
################### Sub Routines ###################
####################################################

sub read_chrlen{
	my $file = shift;
	my $ChrLen_p = shift;
	
	#print STDERR "reading $file\n";
	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my ($chr,$len) = ($1,$2) if(/^(\S+)\s+(\d+)/);
		$ChrLen_p->{$chr} = $len;
	}
	close IN;
}

sub read_chrgap{
	my $file = shift;
	my $ChrGap_p = shift;

	#print STDERR "reading $file\n";
	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\s+/;
		my $chr = $t[0];
		next if($chr ne $ChrID);
		my $pos = $t[1];
		my $bases_num = $t[4];
		my $win_id = int($pos / $Window_size);

		$ChrGap_p->{$chr}[$win_id] += $bases_num;
	}
	close IN;
}


##23      chr1    121535434       124535434       1270    N       3000000 centromere      no
sub read_centromere{
	my $file = shift;
	my $Centro_p = shift;

	#print STDERR "reading $file\n";
	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\s+/;
		my $chr = $t[0];
		my $start = $t[1];
		my $end = $t[2];
		$Centro_p->{$chr} = "$start,$end";
	}
	close IN;
}
