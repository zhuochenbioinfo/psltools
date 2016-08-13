use strict;
use warnings;

my($psl,$outfile) = @ARGV;

my $usage = "USAGE:\nperl $0 <psl file> <output file>\n";

die $usage unless(@ARGV == 2);

my %hash_q;

open(IN,"<$psl") or die $!;
while(<IN>){
	chomp;
	my $inline = $_;
	next unless ($inline =~ /^\d+/);
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$inline;
	unless(exists $hash_q{$qname}){
		$hash_q{$qname}{length} = $qsize;
	}
	my @arr_blocksizes = split/,/,$blocksizes;
	my @arr_qbstarts = split/,/,$qbstarts;
	
	unless(@arr_blocksizes == @arr_qbstarts){
		die "# ERROR: number of block sizes and qbstarts not identical at line: \n$inline\n";
	}
	
	for(my $i = 0; $i < @arr_blocksizes; $i++){
		my $start = "";
		my $end = "";
		if($strand eq "+"){
			$start = $arr_qbstarts[$i] + 1;
			$end = $arr_qbstarts[$i] + $arr_blocksizes[$i];
		}elsif($strand eq "-"){
			$end = $qsize - $arr_qbstarts[$i];
			$start = $end - $arr_blocksizes[$i] + 1;
			# 1 2 3 4 5 6 7 8 9
			# 9 8 7 6 5 4 3 2 1
		}else{
			die "# ERROR: strand error.\n";
		}
		
		if(exists $hash_q{$qname}{covered}{$start}){
			next unless($hash_q{$qname}{covered}{$start} < $end);
		}
		$hash_q{$qname}{covered}{$start} = $end;
	}
}
close IN;

open(OUT,">$outfile");

foreach my $name(sort keys %hash_q){
	my @starts = sort {$a <=> $b} keys %{$hash_q{$name}{covered}};
	my $qlen = $hash_q{$name}{length};
	my $qcovered = 0;
	
	for(my $i = 0; $i < @starts; $i++){
		my $starti = $starts[$i];
		my $endi = $hash_q{$name}{covered}{$starti};
		for(my $j = $i+1; $j < @starts; $j++){
			my $startj = $starts[$j];
			my $endj = $hash_q{$name}{covered}{$startj};
			if($startj > $endi + 1){
				splice(@starts, $i, 1);
				$i--;
				last;
			}
			# merge two blocks into one
			if($endj > $endi){
				$endi = $endj;
				$hash_q{$name}{covered}{$starti} = $endj;
			}
			delete($hash_q{$name}{covered}{$startj});
			splice(@starts, $j, 1);
			$j--;
		}
	}
	
	@starts = sort {$a <=> $b} keys %{$hash_q{$name}{covered}};
	for(my $i = 0; $i < @starts; $i++){
		my $start = $starts[$i];
		my $end = $hash_q{$name}{covered}{$start};
		my $len = $end - $start + 1;
		$qcovered += $len;
	}
	
	my $covered_ratio = $qcovered/$qlen;
	print OUT "$name\t$qlen\t$qcovered\t$covered_ratio\n";
}

close OUT;