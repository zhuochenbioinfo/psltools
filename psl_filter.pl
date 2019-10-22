use strict;
use warnings;

# 1. mark duplicate
# 2. for each query, sort by target
# 3. find deletions

my($pslfile,$outfile) = @ARGV;

my $usage = "USAGE:\nperl $0 <psl output> <outfile>\n";

die $usage unless(@ARGV == 2);

# 1. Check duplicate in each query

my %hash_query;

open(IN,"<$pslfile") or die $!;
while(<IN>){
	chomp;
	#my($query,$qlen,$qstart,$qend,$strand,$target,$tlen,$tstart,$tend,$match,$blocklen,@others) = split/\t/;
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$query,$qlen,$qstart,$qend,$target,$tlen,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$_;
	
	my $blocklen = $match + $mismatch;
	
	if($strand eq "-"){
		my $qstart_new = $qlen - $qend;
		my $qend_new = $qlen - $qstart;
		$qstart = $qstart_new;
		$qend = $qend_new;
	}
	
	if(exists $hash_query{$query}{$qstart}{$qend}){
		$hash_query{$query}{$qstart}{$qend}{dup} ++;
		if($hash_query{$query}{$qstart}{$qend}{ratio} < $match/$blocklen){
			$hash_query{$query}{$qstart}{$qend}{ratio} = $match/$blocklen;
		}
	}else{
		$hash_query{$query}{$qstart}{$qend}{dup} = 1;
		$hash_query{$query}{$qstart}{$qend}{ratio} = $match/$blocklen;
	}
}
close IN;

foreach my $query(keys %hash_query){
	my @starts = sort {$a <=> $b} keys %{$hash_query{$query}};
	my @pairs = ();
	foreach my $start(@starts){
		foreach my $end(sort {$b <=> $a} keys %{$hash_query{$query}{$start}}){
			push @pairs, "$start,$end";
		}
	}
	
	for(my $i = 0; $i < @pairs; $i++){
		my($starti,$endi) = split/,/,$pairs[$i];
		my $lengthi = $endi - $starti + 1;
		for(my $j = $i + 1; $j < @pairs; $j++){
			my($startj,$endj) = split/,/,$pairs[$j];
			my $lengthj = $endj - $startj + 1;
			if($startj > $endi){
				splice(@pairs, $i, 1);
				$i--;
				$j--;
				last;
			}
			my($overlap,$type) = overlapper($starti,$endi,$startj,$endj);
			if($overlap > $lengthi * 0.9 and $overlap > $lengthj * 0.9){
				if($hash_query{$query}{$starti}{$endi}{ratio} > $hash_query{$query}{$startj}{$endj}{ratio}){
					$hash_query{$query}{$startj}{$endj}{dup} ++;
				}else{
					$hash_query{$query}{$starti}{$endi}{dup} ++;
				}
			}elsif($type == 1){
				$hash_query{$query}{$startj}{$endj}{contained} = "";
			}elsif($type == 2){
				$hash_query{$query}{$starti}{$endi}{contained} = "";
			}
		}
	}
}

open(OUT,">$outfile");
open(IN,"<$pslfile") or die $!;
while(<IN>){
	chomp;
	#my($query,$qlen,$qstart,$qend,$strand,$target,$tlen,$tstart,$tend,$match,$blocklen,@others) = split/\t/;
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$query,$qlen,$qstart,$qend,$target,$tlen,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$_;
	
	my $blocklen = $match + $mismatch;
	
	if($strand eq "-"){
		my $qstart_new = $qlen - $qend;
		my $qend_new = $qlen - $qstart;
		$qstart = $qstart_new;
		$qend = $qend_new;
	}
	
	next if(exists $hash_query{$query}{$qstart}{$qend}{contained});
	next if($hash_query{$query}{$qstart}{$qend}{dup} > 1);
	next if($match/$blocklen < 0.8);
	next if($match < 200);
	next if($hash_query{$query}{$qstart}{$qend}{ratio} > $match/$blocklen);
	print OUT "$_\n";
}
close IN;
close OUT;

sub overlapper{
	my($start1,$end1,$start2,$end2) = @_;
	my $overlap = 0;
	my $type = 0;
	# -1 for no overlap; 0 for overlap ; 1 for the first region contains the second region; 2 for the second region contains the first region
	if($end1 < $start2 or $end2 < $start1){
		$type = -1;
	}elsif($end1 <= $end2 and $start1 >= $start2){
		$type = 2;
		$overlap = $end1 - $start1 + 1;
	}elsif($end2 <= $end1 and $start2 >= $start1){
		$type = 1;
		$overlap = $end2 - $start2 + 1;
	}elsif($end1 >= $end2){
		$type = 0;
		$overlap = $end2 - $start1 + 1;
	}elsif($end2 >= $end1){
		$type = 0;
		$overlap = $end1 - $start2 + 1;
	}
	return($overlap,$type);
}
