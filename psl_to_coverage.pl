use strict;
use warnings;

my($psl,$outfile) = @ARGV;

my $usage = "USAGE:\nperl $0 <psl file> <output file>\n";

die $usage unless(@ARGV == 2);

my %hash_t;

open(IN,"<$psl") or die $!;
while(<IN>){
	chomp;
	my $inline = $_;
	next unless ($inline =~ /^\d+/);
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$inline;
	unless(exists $hash_t{$tname}){
		$hash_t{$tname}{length} = $tsize;
	}
	my @arr_blocksizes = split/,/,$blocksizes;
	my @arr_tbstarts = split/,/,$tbstarts;
	
	unless(@arr_blocksizes == @arr_tbstarts){
		die "# ERROR: number of block sizes and tbstarts not identical at line: \n$inline\n";
	}
	
	for(my $i = 0; $i < @arr_blocksizes; $i++){
		my $start = $arr_tbstarts[$i] + 1;
		my $end = $arr_tbstarts[$i] + $arr_blocksizes[$i];
		if(exists $hash_t{$tname}{covered}{$start}){
			next unless($hash_t{$tname}{covered}{$start} < $end);
		}
		$hash_t{$tname}{covered}{$start} = $end;
	}
}
close IN;

open(OUT,">$outfile");

foreach my $tname(sort keys %hash_t){
	my @starts = sort {$a <=> $b} keys %{$hash_t{$tname}{covered}};
	my $tlen = $hash_t{$tname}{length};
	my $tcovered = 0;
	
	for(my $i = 0; $i < @starts; $i++){
		my $starti = $starts[$i];
		my $endi = $hash_t{$tname}{covered}{$starti};
		for(my $j = $i+1; $j < @starts; $j++){
			my $startj = $starts[$j];
			my $endj = $hash_t{$tname}{covered}{$startj};
			if($startj > $endi + 1){
				splice(@starts, $i, 1);
				$i--;
				last;
			}
			# merge two blocks into one
			if($endj > $endi){
				$endi = $endj;
				$hash_t{$tname}{covered}{$starti} = $endj;
			}
			delete($hash_t{$tname}{covered}{$startj});
			splice(@starts, $j, 1);
			$j--;
		}
	}
	
	@starts = sort {$a <=> $b} keys %{$hash_t{$tname}{covered}};
	for(my $i = 0; $i < @starts; $i++){
		my $start = $starts[$i];
		my $end = $hash_t{$tname}{covered}{$start};
		my $len = $end - $start + 1;
		$tcovered += $len;
	}
	
	my $covered_ratio = $tcovered/$tlen;
	print OUT "$tname\t$tlen\t$tcovered\t$covered_ratio\n";
}

close OUT;