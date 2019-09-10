#!/usr/bin/perl

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Verbose,$Help,$Likelihood_cutoff);
GetOptions(
	"likelihood:f"=>\$Likelihood_cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Likelihood_cutoff ||= 0.95;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $hetsnp_filter_file = shift;
my $probandfile = shift;
my $draw_pos = shift;
my $Aneuploidy_chrs_file = "$Bin/Aneuploidy_chrs_file";

my $hetsnp_filter_dir = dirname($hetsnp_filter_file);


my $pairwisehmm = "$Bin/pairwisehmm -n 20 ";
my $hmm_para = "$Bin/pairwise_2states.hmm";

my $ChrID = $1 if($hetsnp_filter_file =~ /(chr\w+)/);


print STDERR "chrid is $ChrID\n";

my @Heads;
my @Pos;
my @Data;
my %Acell;

my %Aneuploidy;

if (-f $Aneuploidy_chrs_file) {
	open IN, $Aneuploidy_chrs_file || die "fail $Aneuploidy_chrs_file";
	while (<IN>) {
		chomp;
		my ($sample,$chr,$copy) = split /\s+/;
		$Aneuploidy{$sample}{$chr} = $copy;
		#print STDERR $sample,$chr,$copy;
	}
	close IN;
}


open IN, $hetsnp_filter_file || die "fail $hetsnp_filter_file";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	if (/^\#/) {
		for (my $i=9; $i<@t; $i+=2) {
			if (!exists $Aneuploidy{$t[$i]}{$ChrID}) {
				push @Heads, $t[$i]; #yan, @Heads存储了所有精子的名称
				#print STDERR "sample $t[$i]\n";
			}else{
				$Acell{$i} = "A"; #yan, @Acell 存储了非整倍体精子的id，不是名称。
				#print STDERR "fail $t[$i]\n";
			}
		}
		next;
	}

	push @Pos, [$t[0],$t[1]]; #yan, @Pos=[[chr, pos]...]
	my @temp;
	for (my $i=9; $i<@t; $i+=2) {
		if ($t[$i+1] < $Likelihood_cutoff) {
			$t[$i] = "N";
		}
		push @temp, $t[$i] if(!exists $Acell{$i});
	}
	push @Data, \@temp; #yan, @Data 每个元素表示一行，存储了所有精子在该位点的碱基

}
close IN;

open IN, $probandfile || die "fail $probandfile";
my $ancestor_id;
while (<IN>) {
	chomp;
	$ancestor_id = $_;	
}
close IN;

my ($ancestor)=grep{$Heads[$_] eq $ancestor_id} 0..$#Heads;
#print STDERR "ancestor\t".$ancestor_id."\t".$ancestor."\n";

##有aneuploidy的sample/chromosome, 不参与计算, 从源头保证数据的纯洁性
next if(exists $Aneuploidy{$ancestor_id}{$ChrID});
my $Sample_total = 0;
my $crossover_file = $hetsnp_filter_file.".$ancestor_id.crossover"; #yan, 对于每一个样本，生成一个hmm处理的文件
my $hiddenseq_file = $hetsnp_filter_file.".$ancestor_id.hiddenseq";
open CO, ">$crossover_file" || die "fail $crossover_file";
open HS, ">$hiddenseq_file" || die "fail $hiddenseq_file";
for (my $offspring=0; $offspring<@Heads; $offspring++) {
	next if($Heads[$offspring] eq $ancestor_id);
	
	my $offspring_id = $Heads[$offspring];
	
	next if(exists $Aneuploidy{$offspring_id}{$ChrID});
	
	$Sample_total ++;

	my $temp_file = $hetsnp_filter_file.".temp"; #yan,对于每一个ancestor和offspring的组合，生成一个temp file
	open OUT, ">$temp_file" || die "fail $temp_file";
	print OUT "#chr\tposition\t$ancestor_id\t$offspring_id\n";
	for (my $pos=0; $pos<@Pos; $pos++) {
		my $chr = $Pos[$pos][0];
		my $position = $Pos[$pos][1];
		my $ancestor_base = $Data[$pos][$ancestor];
		my $offspring_base = $Data[$pos][$offspring];
		if ($ancestor_base ne "N" && $offspring_base ne "N") {
			print OUT "$chr\t$position\t$ancestor_base\t$offspring_base\n";
		}
	}
	close OUT;
	my $SNP_pair_num = `cat $temp_file | wc -l`;
	chomp $SNP_pair_num;
	#print STDERR "\n##SNP pair number between $ancestor_id and $offspring_id is:  $SNP_pair_num\n\n";
	
	`$pairwisehmm  $hmm_para $temp_file`; #yan, 由于temp文件里只存放两个精子样本，因而每次只使用hmm处理两个样本
	
	my $str = `cat $temp_file.crossover`;
	my $str2 = `cat $temp_file.hiddenseq`;
	`rm $temp_file  $temp_file.crossover  $temp_file.hiddenseq`;
	
	my @lines = split /\n/, $str;
	my @lines2 = split /\n/, $str2;
	foreach my $line (@lines) { #yan, 对于每一个 CO 位置
		chomp;
		next if($line =~ /^\#/);
		print CO "$ancestor_id\t$offspring_id\t$line\n";
	}
	foreach my $line2 (@lines2) {
		chomp;
		next if($line2 =~ /^ChrID/);
		print HS "$ancestor_id\t$offspring_id\t$line2\n"; #yan, 注意，此处已经将两个样本的id输出来了
	}
}
close CO;
close HS;

`perl $Bin/filter_crossover.pl $crossover_file > $crossover_file.filter`;

`perl $Bin/infer_crossover_by_merging.pl $crossover_file.filter $Aneuploidy_chrs_file $Sample_total > $crossover_file.filter.cluster`;
##手动更改cluster文献，可以实现人工矫正的功能
`perl $Bin/draw_pairwise_crossovers_for_curation.green.red.pl $crossover_file $hiddenseq_file $crossover_file.filter.cluster  $Aneuploidy_chrs_file $draw_pos`;

