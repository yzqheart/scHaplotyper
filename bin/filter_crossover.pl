use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

my $chr_len_file = "$Bin/../db/hg38.len.pure";
my $human_gap_file = "$Bin/../db/hg38.gap.centrome.telemore";

my $pairwise_crossovers_file = shift;

my $ChrID = $1 if($pairwise_crossovers_file =~ /(chr\w+)/);

my %ChrLen;
my %CO;
my %MidCent;
my %TerCent;

my $Min_dist_two_crossover = 1000000;

read_chrlen($chr_len_file, \%ChrLen); #yan, ChrLen is a hash, ChrLen{chrID}=chrIDLen

read_chrgap($human_gap_file, \%MidCent, \%TerCent); #yan, 近端着丝粒 和 近中着丝粒 %MidCent={$chr}:[$start, $end, $len]

##print STDERR Dumper \%MidCent;
##print STDERR Dumper \%TerCent;

open IN, $pairwise_crossovers_file || die "fail $pairwise_crossovers_file";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	my $sample = $t[1]; #yan, each offspring
	my $pos = $t[3]; #yan, each crossover site in chr
	$CO{$sample}{$pos} = \@t;

}
close IN;


foreach my $sample_id (sort keys %CO) {
	my $sample_p = $CO{$sample_id};
	
	##含有两个或以上crossover时，检查相邻间距，小于某个cutoff,则拔出这两个crossover, 每循环一轮只拔一次
	while (1) {
		my @Pos = sort {$a<=>$b} keys %$sample_p;
		last if(@Pos < 2);
		my $has_delete = 0;
		for (my $i=1; $i<@Pos; $i++) {
			if ($Pos[$i] - $Pos[$i-1] < $Min_dist_two_crossover) {
				delete $sample_p->{$Pos[$i]};
				delete $sample_p->{$Pos[$i-1]};
				$has_delete = 1;
				last;
			}elsif(exists $MidCent{$ChrID}){
				my $centro_start = $MidCent{$ChrID}[0];
				my $centro_end = $MidCent{$ChrID}[1];
				my $centro_len = $centro_end - $centro_start;
				if ( ($Pos[$i-1] < $centro_start && $Pos[$i] > $centro_end && $Pos[$i] - $Pos[$i-1] - $centro_len < $Min_dist_two_crossover) || ($Pos[$i-1] > $centro_start && $Pos[$i] < $centro_end) ) {
					delete $sample_p->{$Pos[$i]};
					delete $sample_p->{$Pos[$i-1]};
					$has_delete = 1;
					last;
				}
			
			}
		}
		last if($has_delete == 0);
	}

	##当剩下大于等于1个crossover时，检查其是否位于染色体两端
	my @Pos = sort {$a<=>$b} keys %$sample_p;
	if (@Pos >= 1) {
		delete $sample_p->{$Pos[0]} if($Pos[0] < $Min_dist_two_crossover);
		delete $sample_p->{$Pos[-1]} if($Pos[-1] > $ChrLen{$ChrID} - $Min_dist_two_crossover);
	}

	##对于末端centromere染色体，删除所有距离gap近的crossovers
	if (exists $TerCent{$ChrID}) {
		foreach my $pos (sort {$a<=>$b} keys %$sample_p) {
			my $gap_end = $TerCent{$ChrID}[1];
			if ($pos < $gap_end + $Min_dist_two_crossover) {
				delete $sample_p->{$pos};
			}
		}
	}
	
	##Output合格的crossovers
	foreach my $pos (sort {$a<=>$b} keys %$sample_p) {
		my $pos_p = $sample_p->{$pos};
		$pos_p->[7] = 1;
		my $line = join("\t", @$pos_p);
		print $line."\n";
	}

} 


####################################################################################

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


##read centromere and telomere gaps
sub read_chrgap{
	my $file = shift;
	my $mid_p = shift;
	my $ter_p = shift;

	#print STDERR "reading $file\n";
	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my ($chr, $start, $end, $len, $type) = split /\s+/;
		if ($type eq "MiddleCentromere") {
			$mid_p->{$chr} = [$start, $end, $len];
		}
		if ($type eq "TerminalCentromere") {
			$ter_p->{$chr} = [$start, $end, $len];
		}
	}
	close IN;
}

