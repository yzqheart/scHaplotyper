use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory


my $pairwise_crossovers_file = shift;
my $Aneuploidy_chrs_file = shift;
my $Sample_total = shift; #yan, offspring num

my $Distance_cutoff = 1000000;
my $Reference_crossover_ratio = 0.6;

my $ChrID = $1 if($pairwise_crossovers_file =~ /(chr\w+)/);

my @CO;

my %Aneuploidy;

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
	my $sample = $t[1]; #yan, offspring_id
	next if(exists $Aneuploidy{$sample}{$ChrID});
	my $pos = $t[3]; #yan, crossover_pos
	push @CO, \@t; #yan, @CO 包含了原始文件所有信息
}
close IN;

#print STDERR "Sample_total: $Sample_total\n";

##sort all crossovers according to the genomic positions
@CO = sort {$a->[3] <=> $b->[3]} @CO;

##print Dumper \@CO;

##predict the reference（就是ancestor） crossover by clustering method
if (@CO >= 1) {
	my @cluster;
	my $Confident_crossover_id = 1;
	push @cluster, $CO[0]; #yan, 将crossover file的第一行push到@cluster中
	for (my $i=1; $i<@CO; $i++) { #yan, 剩余所有行
		my $this_pos = $CO[$i][3];
		my $cluster_end = $cluster[-1][3];
		if ($this_pos - $cluster_end < $Distance_cutoff) {
			push @cluster, $CO[$i];
		}else{ 
			my %Count;
			for (my $j=0; $j<@cluster; $j++) {
				my $sample = $cluster[$j][1];
				$Count{$sample} = $cluster[$j];
			}
			my $number = keys %Count;
			my $pos_centered;
			foreach my $sample (keys %Count) {
				$pos_centered += $Count{$sample}[3];
			}
			$pos_centered /= $number;
			$pos_centered = int($pos_centered);
			my $ratio = $number / $Sample_total;
			if ($ratio >= $Reference_crossover_ratio) {#yan, 大部分crossover发生在offspring中，表明实际的crossover在reference中
				print "\n**$Confident_crossover_id  Confident Reference Crossover:  $pos_centered  $ratio\n";
				$Confident_crossover_id ++;
			}else{
				print "\n## Reference crossover:  $pos_centered  $ratio\n";
			}

			foreach my $p (@cluster) {
				my $line = join("\t",@$p);
				print $line."\n";
			}
			
			@cluster = ();
			push @cluster, $CO[$i];
		}
	}

	if (@cluster) {  ##output the last cluster's content
		my %Count;
		for (my $j=0; $j<@cluster; $j++) {
			my $sample = $cluster[$j][1];
			$Count{$sample} = $cluster[$j];
		}
		my $number = keys %Count;
		my $pos_centered;
		foreach my $sample (keys %Count) {
			$pos_centered += $Count{$sample}[3];
		}
		$pos_centered /= $number;
		$pos_centered = int($pos_centered);
		my $ratio = $number / $Sample_total;
		if ($ratio >= $Reference_crossover_ratio) {
			print "\n**$Confident_crossover_id  Confident Reference Crossover:  $pos_centered  $ratio\n";
			$Confident_crossover_id ++;
		}else{
			print "\n## Reference Crossover:  $pos_centered  $ratio\n";
		}

		foreach my $p (@cluster) {
			my $line = join("\t",@$p);
			print $line."\n";
		}
	}

}

