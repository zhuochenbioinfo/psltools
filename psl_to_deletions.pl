use strict;
use warnings;

my($psl,$outfile) = @ARGV;

my $usage = "USAGE:\nperl $0 <psl file> <output file>\n";

die $usage unless(@ARGV == 2);

my %hash_q;
my %hash_t;

open(IN,"cat $psl|sort -k1rn -k2n|") or die $!;
open(OUT,">$outfile");

print OUT "#REF\tSTART\tEND\tGAPLEN\tQUERY\tQSTART\tQEND\tQGAPLEN\tSTRAND\tBLOCKLEFT\tBLOCKRIGHT\tMATCHRATIO\n";

while(<IN>){
	chomp;
	my $inline = $_;
	next unless ($inline =~ /^\d+/);
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$inline;
	
	next if(exists $hash_q{$qname});
	$hash_q{$qname} = "";
	
	my @arr_blocksizes = split/,/,$blocksizes;
	my @arr_tbstarts = split/,/,$tbstarts;
	my @arr_qbstarts = split/,/,$qbstarts;
	
	my $allblocksize = $match + $mismatch;
	my $match_ratio = $match/$allblocksize;
	
	my $blockleft = 0;
	my $blockright = $allblocksize;
	
	for(my $i = 0; $i < @arr_tbstarts - 1; $i++){
		my $tgapstart = $arr_tbstarts[$i] + $arr_blocksizes[$i] + 1;
		my $tgapend = $arr_tbstarts[$i+1];
		my $tgaplen = $tgapend - $tgapstart + 1;
		my $qgapstart = $arr_qbstarts[$i] + $arr_blocksizes[$i] + 1;
		my $qgapend = $arr_qbstarts[$i+1];
		my $qgaplen = $qgapend - $qgapstart + 1;
		
		$blockleft = $arr_blocksizes[$i];
		$blockright = $arr_blocksizes[$i+1];
		
		
		print OUT "$tname\t$tgapstart\t$tgapend\t$tgaplen\t$qname\t$qgapstart\t$qgapend\t$qgaplen\t$strand\t$blockleft\t$blockright\t$match_ratio\n";
	}
	
}
close IN;
close OUT;