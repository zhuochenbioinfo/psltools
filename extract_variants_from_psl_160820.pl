#########################################################
#	extract_variants_from_psl.pl
#	
#	1. Extracting snps and indels from BLAT psl output file 
#	2. Counting coverage of each base in query
#	
#	By Zhuo CHEN, IGDB, CAS
#	zhuochen@genetics.ac.cn
#	zhuochenbioinfo@gamil.com
#
#	Version 2.0 at 08.20.2016
#	
#########################################################

=note
1. Change the strategy to calculate coverage and reduce memory use.
	No longer reads every bases of the target into a hash for coverage and variants analysis. Allowing big target sequences.
2. Change the output format.
	Adding support and conflict query sequence information to each variant.
=cut

use strict;
use warnings;
use Getopt::Long;

my($psl,$qref,$tref,$outname,$maxgapsize,$q2bit,$t2bit);
my $usage = "USAGE:\nperl $0 --psl <psl file> --qref <query refence seq> --tref <target reference seq> --out <outname> --maxgap <max gap size>\n";
$usage .= "<psl file> is the BLAT output file. Better filtered and sorted by match length or other requirements.\n";
$usage .= "<qref> the query fasta file.\n";
$usage .= "<tref> the target fasta file.\n";
$usage .= "<outname> the prefix of the output files. There will be three output files named as <outname>.var <outname>.block and <outname>.gap\n";
$usage .= "<max gap size> is the up limit of the gap size to be retained in the output.\n";
$usage .= "--q2b <query 2bit list>: the 2bit list for query fasta.\n";
$usage .= "--t2b <target 2bit list>: the 2bit list for target fasta.\n";

GetOptions(
	"psl=s" => \$psl,
	"qref=s" => \$qref,
	"tref=s" => \$tref,
	"out=s" => \$outname,
	"maxgap=s" => \$maxgapsize,
	"q2b=s" => \$q2bit,
	"t2b=s" => \$t2bit,
) or die $usage;

unless(defined $psl and defined $qref and defined $tref and defined $outname){
	die $usage;
}

unless(defined $maxgapsize){
	$maxgapsize = 100;
}

# ---------------------- PART I: Reading the reference sequences -----------------------

my %hash_qref;
my %hash_tref;

my $qrefcount = 0;
my $trefcount = 0;

if(defined $q2bit){
	open(IN,"<$q2bit") or die $!;
	while(<IN>){
		chomp;
		my($twobit,$chr,$others) = split/\:/;
		$hash_qref{$chr}{twobit} = "";
	}
	close IN;
	$qrefcount = keys %hash_qref;
}
if(defined $t2bit){
	open(IN,"<$t2bit") or die $!;
	while(<IN>){
		chomp;
		my($twobit,$chr,$others) = split/\:/;
		$hash_tref{$chr}{twobit} = "";
	}
	close IN;
	$trefcount = keys %hash_tref;
}


print "Reading query reference...\n";
open(QREF,"<$qref") or die "$!";
local $/ = "\n>";
while(<QREF>){
	$_ =~ s/^>//g;
	chomp;
	my($head,$seq) = split/\n/,$_,2;
	my($seqname) = $head =~ /^(\S+)/;
	if(defined $q2bit){
		last if($qrefcount == 0);
		unless(exists $hash_qref{$seqname}){
			next;
		}else{
			$qrefcount--;
		}
	}
	$seq =~ s/\s+//g;
	$seq = uc($seq);
	
	my $seq_revcom = $seq;
	$seq_revcom =~ tr/ATCG/TAGC/;
	$seq_revcom = reverse($seq_revcom);
	
	$hash_qref{$seqname}{seq} = $seq;
	$hash_qref{$seqname}{revcom} = $seq_revcom;
	$hash_qref{$seqname}{len} = length($seq);
	#print "$seqname\t".length($seq)."\n";
}
local $/ = "\n";
close QREF;

print "Reading target reference...\n";

my %hash_pslcount;

open(TREF,"<$tref") or die "$!";
local $/ = "\n>";
while(<TREF>){
	$_ =~ s/^>//g;
	chomp;
	my($head,$seq) = split/\n/,$_,2;
	my($seqname) = $head =~ /^(\S+)/;
	if(defined $t2bit){
		last if($trefcount == 0);
		unless(exists $hash_tref{$seqname}){
			next;
		}else{
			$trefcount--;
		}
	}
	$seq =~ s/\s+//g;
	$seq = uc($seq);
	$hash_tref{$seqname}{seq} = $seq;
	$hash_tref{$seqname}{len} = length($seq);
	$hash_pslcount{$seqname} = 0;
}
local $/ = "\n";
close TREF;

# --------------------- PART I END -----------------------

# --------------------- PART II: Reading psl file -----------------------

print "Reading PSL file...\n";

open(PSL,"<$psl") or die "$!";

my %hash_var;
my %hash_query;

while(<PSL>){
	chomp;
	my $inline = $_;
	next unless ($inline =~ /^\d+/);
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$inline;
	
	# check if the query exists in more than one alignment
	next if(exists $hash_query{$qname});
	$hash_query{$qname}{match} = $match;
	$hash_query{$qname}{target} = $tname;
	$hash_query{$qname}{start} = $qstart;
	$hash_query{$qname}{end} = $qend;
	
	my @qbs = split/,/,$qbstarts;
	my @tbs = split/,/,$tbstarts;
	my @blks = split/,/,$blocksizes;
	
	#my($flankup,$flankdown) = split/,/,$flanks;
	
	# count the number of query alignments to cover the target
	$hash_pslcount{$tname}++;
	
	# set block cover
	# block region are covered both absolutely and relatively
	for(my $count = 0; $count < $blknum; $count++){
		# The gap and block start of block is 1-based
		my $tbstart = $tbs[$count] + 1;
		my $tbend = $tbs[$count] + $blks[$count];
		$hash_query{$qname}{block}{$tbstart} = $tbend;
	}
	
	# set gap cover
	# gap region are covered only relatively
	# gap number is one less than block number
	# the start pos of gap and block are both 1-based, there may be gaps with size equals 0, when gapend - gaptart + 1 = 0
	for(my $count = 0; $count < $blknum - 1; $count++){
		# important change here! 20160823
		my $tgapstart = $tbs[$count] + $blks[$count] + 1;
		my $tgapend = $tbs[$count+1];
		# gapend may be smaller than gapstart when an insertion occurs
		# query  ATGTACTCGTTCCAGAT
		# target ATGTAC-----CCAGAT
		$hash_query{$qname}{gap}{$tgapstart} = $tgapend;
	}
	
	# Mind that the cover and block border may change if there are blocks with wrong end base
	# which means the last bases of qblock and tblock are not the same
	
	for(my $count = 0; $count < $blknum; $count++){
		my $tseq = substr($hash_tref{$tname}{seq},$tbs[$count],$blks[$count]);
		my $qseq;
		
		# in blat psl results query block start positions are based on the negative strand if the alignment strand is negative
		# block start is 0-based
		if($strand eq "+"){
			$qseq = substr($hash_qref{$qname}{seq},$qbs[$count],$blks[$count]);
		}else{
			$qseq = substr($hash_qref{$qname}{revcom},$qbs[$count],$blks[$count]);
		}
		
		#dealing with mismatches -> snps
		for(my $i=0;$i<$blks[$count];$i++){
			my $qbase = substr($qseq,$i,1);
			my $tbase = substr($tseq,$i,1);
			my $tpos = $tbs[$count] + 1 + $i;
			
			if($qbase ne $tbase){
				# if the last bases of the two blocks are not the same, change the block and gap border
				if($i == $blks[$count]-1){
					print "# WRONG BLOCK END at query:$qname target:$tname blocknum:$count pos:$tpos\n";
					if(exists $hash_var{$tname}{$tpos}){
						if(exists $hash_var{$tname}{$tpos}{SNP}){
							if(exists $hash_var{$tname}{$tpos}{SNP}{ALT}{$qbase}){
								my @tmp = throw_item_from_array($qname,@{$hash_var{$tname}{$tpos}{SNP}{ALT}{$qbase}{query}});
								if(@tmp == 0){
									# if the query is the only one sample to introduce the variant, delete the variant when removing WBE
									delete($hash_var{$tname}{$tpos}{SNP}{ALT}{$qbase});
									delete($hash_query{$qname}{var}{$tname}{$tpos}{SNP});
								}else{
									@{$hash_var{$tname}{$tpos}{SNP}{ALT}{$qbase}{query}} = @tmp;
								}
							}
						}
					}
					# adjust the end of the block 
					$hash_query{$qname}{block}{$tbs[$count]+1}--;
					
					# change the border of the next gap
					if($count < $blknum - 1){
						my $nextgapstart = $tbs[$count] + $blks[$count] + 1;
						my $nextgapend = $tbs[$count+1];
						
						delete($hash_query{$qname}{gap}{$nextgapstart});
						
						# change the block size
						$blks[$count]--;
						
						$nextgapstart = $tbs[$count] + $blks[$count] + 1;
						$hash_query{$qname}{gap}{$nextgapstart} = $nextgapend;
					}else{
						# change the block size
						$blks[$count]--;
					}
					$i = $i - 2;
					next;
				}
				push @{$hash_var{$tname}{$tpos}{SNP}{ALT}{$qbase}{query}}, $qname;
				$hash_var{$tname}{$tpos}{SNP}{REF} = $tbase;
				$hash_query{$qname}{var}{$tname}{$tpos}{SNP} = $qbase;
			}
		}
		
		
		#dealing with gap -> indel
		if($count < $blknum - 1){
			my $tgapstart = $tbs[$count] + $blks[$count] + 1;
			my $tgapend = $tbs[$count + 1];
			my $qgapstart = $qbs[$count] + $blks[$count] + 1;
			my $qgapend = $qbs[$count + 1];
			
			my $tgaplen = $tgapend - $tgapstart + 1;
			my $qgaplen = $qgapend - $qgapstart + 1;
			
			my $tpos = $tbs[$count] + $blks[$count];
			
			
			# Set the conditions the abandon a gap
			my $throwgap = 0;
			if($tgaplen > $maxgapsize or $qgaplen > $maxgapsize){
				$throwgap = 1;
			}

			if($throwgap == 1){
				delete($hash_query{$qname}{gap}{$tgapstart});
				goto THROW2;
			}
			
			# get the last bases of the q and t block 
			my $qgap = substr($qseq,$blks[$count]-1,1);
			my $tgap = substr($tseq,$blks[$count]-1,1);
			
			my $tgapbases = substr($hash_tref{$tname}{seq},$tgapstart - 1, $tgapend - $tgapstart + 1);
			$tgap .= $tgapbases;
			
			my $qgapbases;
			if($strand eq "+"){
				$qgapbases = substr($hash_qref{$qname}{seq}, $qgapstart - 1, $qgapend - $qgapstart + 1);
			}else{
				$qgapbases = substr($hash_qref{$qname}{revcom}, $qgapstart - 1, $qgapend - $qgapstart + 1);
			}
			$qgap .= $qgapbases;
			
			push @{$hash_var{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{query}}, $qname;
			$hash_query{$qname}{var}{$tname}{$tpos}{INDEL}{REF} = $tgap;
			$hash_query{$qname}{var}{$tname}{$tpos}{INDEL}{ALT} = $qgap;
		}
		THROW2:
		
	}
}
close PSL;
print "Done!\n";
# --------------------- PART II END -----------------------

# --------------------- PART III: dissect coverage -----------------------

my %hash_cov;
my %hash_blk;
my %hash_gap;
my %hash_covre;

# collect covered regions of targets by queries
foreach my $qname(sort keys %hash_query){
	my $tname = $hash_query{$qname}{target};
	foreach my $tstart(sort {$a <=> $b} keys %{$hash_query{$qname}{block}}){
		my $tend = $hash_query{$qname}{block}{$tstart};
		push @{$hash_blk{$tname}{$tstart}{$tend}{query}}, $qname;
		if(exists $hash_cov{$tname}{block}{$tstart}){
			if($hash_cov{$tname}{block}{$tstart} < $tend){
				$hash_cov{$tname}{block}{$tstart} = $tend;
			}
		}else{
			$hash_cov{$tname}{block}{$tstart} = $tend;
		}
		# put the blocks into relatively covered region
		if(exists $hash_covre{$tname}{$tstart}){
			if($hash_covre{$tname}{$tstart} < $tend){
				$hash_covre{$tname}{$tstart} = $tend;
			}
		}else{
			$hash_covre{$tname}{$tstart} = $tend;
		}
	}
	foreach my $tstart(sort {$a <=> $b} keys %{$hash_query{$qname}{gap}}){
		my $tend = $hash_query{$qname}{gap}{$tstart};
		push @{$hash_gap{$tname}{$tstart}{$tend}{query}}, $qname;
		if(exists $hash_cov{$tname}{gap}{$tstart}){
			if($hash_cov{$tname}{gap}{$tstart} < $tend){
				$hash_cov{$tname}{gap}{$tstart} = $tend;
			}
		}else{
			$hash_cov{$tname}{gap}{$tstart} = $tend;
		}
		# put the gaps into relatively covered region
		if(exists $hash_covre{$tname}{$tstart}){
			if($hash_covre{$tname}{$tstart} < $tend){
				$hash_covre{$tname}{$tstart} = $tend;
			}
		}else{
			$hash_covre{$tname}{$tstart} = $tend;
		}
	}
}

# merge covered regions
foreach my $tname(sort keys %hash_cov){

	# merge absolutely covered region
	my @starts = sort {$a <=> $b} keys %{$hash_cov{$tname}{block}};
	for(my $i = 0; $i < @starts; $i++){
		my $start = $starts[$i];
		my $end = $hash_cov{$tname}{block}{$start};
		for(my $j = $i+1; $j < @starts; $j++){
			my $start_ = $starts[$j];
			my $end_ = $hash_cov{$tname}{block}{$start_};
			if($start_ <= $end + 1 and $end_ > $end){
				$end = $end_;
				$hash_cov{$tname}{block}{$start} = $end;
				delete($hash_cov{$tname}{block}{$start_});
				splice(@starts, $j, 1);
				$j--;
				next;
			}elsif($start_ <= $end + 1 and $end_ <= $end){
				delete($hash_cov{$tname}{block}{$start_});
				splice(@starts, $j, 1);
				$j--;
				next;
			}else{
				splice(@starts, $i, 1);
				$i--;
				last;
			}
		}
	}
	
	# merge relatively covered region
	# relatively covered region contains blocks and gaps
	@starts = sort {$a <=> $b} keys %{$hash_covre{$tname}};
	for(my $i = 0; $i < @starts; $i++){
		my $start = $starts[$i];
		my $end = $hash_covre{$tname}{$start};
		for(my $j = $i+1; $j < @starts; $j++){
			my $start_ = $starts[$j];
			my $end_ = $hash_covre{$tname}{$start_};
			if($start_ <= $end + 1 and $end_ > $end){
				$end = $end_;
				$hash_covre{$tname}{$start} = $end;
				delete($hash_covre{$tname}{$start_});
				splice(@starts, $j, 1);
				$j--;
				next;
			}elsif($start_ <= $end + 1 and $end_ <= $end){
				delete($hash_covre{$tname}{$start_});
				splice(@starts, $j, 1);
				$j--;
				next;
			}else{
				splice(@starts, $i, 1);
				$i--;
				last;
			}
		}
	}
}

# output coverage information
open(COV,">$outname.covab");
print COV "#CHROM\tSTART\tEND\n";
open(RCOV,">$outname.covre");
print RCOV "#CHROM\tSTART\tEND\n";

foreach my $tname (keys %hash_cov){
	my @starts = sort {$a <=> $b} keys %{$hash_cov{$tname}{block}};
	foreach my $start(@starts){
		my $end = $hash_cov{$tname}{block}{$start};
		print COV "$tname\t$start\t$end\n";
	}
	@starts = sort {$a <=> $b} keys %{$hash_covre{$tname}};
	foreach my $start(@starts){
		my $end = $hash_covre{$tname}{$start};
		print RCOV "$tname\t$start\t$end\n";
	}
}
close COV;
close RCOV;

# --------------------- PART III END -----------------------

# --------------------- PART III: dissect variants -----------------------

#CHROM POS   TYPE  REF  ALT  SUPPORT                         CONFLICT
#Chr1  5000  SNP   A    T    ctg1:SNP:2123M;ctg12:SNP:160M   ctg3:REF:637M;ctg8:GAP:117M

open(OUT,">$outname.var");

print OUT "#CHROM\tPOS\tTYPE\tREF\tALT\tSUPPORT\tCONFLICT\n";

foreach my $tname(sort keys %hash_var){
	my @positions = sort {$a <=> $b} keys %{$hash_var{$tname}};
	my @blkstarts = sort {$a <=> $b} keys %{$hash_blk{$tname}};
	my @gapstarts = sort {$a <=> $b} keys %{$hash_gap{$tname}};
	
	foreach my $pos(sort {$a <=> $b} keys %{$hash_var{$tname}}){
		
		# pick the query seqs that cover the position by block
		my @blk_qnames;
		for(my $i = 0; $i < @blkstarts;$i++){
			my $blkstart = $blkstarts[$i];
			if($blkstart > $pos){
				last;
			}
			my @blkends = sort {$a <=> $b} keys %{$hash_blk{$tname}{$blkstart}};
			if(@blkends == 0){
				splice(@blkstarts, $i ,1);
				$i--;
				next;
			}
			foreach my $blkend(@blkends){
				if($blkend < $pos){
					delete($hash_blk{$tname}{$blkstart}{$blkend});
					next;
				}
				if($blkend >= $pos){
					push @blk_qnames, @{$hash_blk{$tname}{$blkstart}{$blkend}{query}};
				}
			}
		}
		
		# pick the query seqs that cover the position by gap
		my @gap_qnames;
		for(my $i = 0; $i < @gapstarts;$i++){
			my $gapstart = $gapstarts[$i];
			if($gapstart > $pos){
				last;
			}
			my @gapends = sort {$a <=> $b} keys %{$hash_gap{$tname}{$gapstart}};
			if(@gapends == 0){
				splice(@gapstarts, $i ,1);
				$i--;
				next;
			}
			foreach my $gapend(@gapends){
				if($gapend < $pos){
					delete($hash_gap{$tname}{$gapstart}{$gapend});
					next;
				}
				if($gapend >= $pos){
					push @gap_qnames, @{$hash_gap{$tname}{$gapstart}{$gapend}{query}};
				}
			}
		}
		
		if(exists $hash_var{$tname}{$pos}{SNP}){
			my $tbase = $hash_var{$tname}{$pos}{SNP}{REF};
			my @ref_qnames;
			my @alt_qnames;
			foreach my $qname(@blk_qnames){
				if(exists $hash_query{$qname}{var}{$tname}{$pos}{SNP}){
					push @alt_qnames, $qname;
				}else{
					push @ref_qnames, $qname;
				}
			}
			foreach my $qbase(sort keys %{$hash_var{$tname}{$pos}{SNP}{ALT}}){
				my @supp_qnames = ();
				my @conf_qnames = ();
				foreach my $qname(@alt_qnames){
					if($hash_query{$qname}{var}{$tname}{$pos}{SNP} eq $qbase){
						push @supp_qnames, $qname;
					}else{
						push @conf_qnames, $qname;
					}
				}
				my $support_queries = "";
				my $conflict_queries = "";
				foreach my $qname(@supp_qnames){
					my $qmatch = $hash_query{$qname}{match};
					$support_queries .= "$qname:SNP:M$qmatch;";
				}
				$support_queries =~ s/;$//;
				if($support_queries eq ""){
					$support_queries = "-";
				}
				
				foreach my $qname(@ref_qnames){
					my $qmatch = $hash_query{$qname}{match};
					$conflict_queries .= "$qname:REF:M$qmatch;";
				}
				foreach my $qname(@gap_qnames){
					my $qmatch = $hash_query{$qname}{match};
					$conflict_queries .= "$qname:GAP:M$qmatch;";
				}
				foreach my $qname(@conf_qnames){
					my $qmatch = $hash_query{$qname}{match};
					$conflict_queries .= "$qname:SNP:M$qmatch;";
				}
				$conflict_queries =~ s/;$//;
				if($conflict_queries eq ""){
					$conflict_queries = "-";
				}
				print OUT "$tname\t$pos\tSNP\t$tbase\t$qbase\t$support_queries\t$conflict_queries\n";
			}
		}
		
		if(exists $hash_var{$tname}{$pos}{INDEL}){
			my @ref_qnames;
			my @alt_qnames;
			my @snp_qnames;
			foreach my $qname(@blk_qnames){
				if(exists $hash_query{$qname}{var}{$tname}{$pos}{INDEL}){
					push @alt_qnames, $qname;
				}elsif(exists $hash_query{$qname}{var}{$tname}{$pos}{SNP}){
					push @snp_qnames, $qname;
				}else{
					push @ref_qnames, $qname;
				}
			}
			foreach my $tgap(sort keys %{$hash_var{$tname}{$pos}{INDEL}{REF}}){
				foreach my $qgap(sort keys %{$hash_var{$tname}{$pos}{INDEL}{REF}{$tgap}{ALT}}){
					my @supp_qnames = ();
					my @conf_qnames = ();
					foreach my $qname(@alt_qnames){
						my $tgap_ = $hash_query{$qname}{var}{$tname}{$pos}{INDEL}{REF};
						my $qgap_ = $hash_query{$qname}{var}{$tname}{$pos}{INDEL}{ALT};
						if($tgap eq $tgap_ and $qgap_ eq $qgap){
							push @supp_qnames, $qname;
						}else{
							push @conf_qnames, $qname;
						}
					}
					my $support_queries = "";
					my $conflict_queries = "";
					foreach my $qname(@supp_qnames){
						my $qmatch = $hash_query{$qname}{match};
						$support_queries .= "$qname:INDEL:M$qmatch;";
					}
					$support_queries =~ s/;$//;
					if($support_queries eq ""){
						$support_queries = "-";
					}
					
					foreach my $qname(@ref_qnames){
						my $qmatch = $hash_query{$qname}{match};
						$conflict_queries .= "$qname:REF:M$qmatch;";
					}
					foreach my $qname(@gap_qnames){
						my $qmatch = $hash_query{$qname}{match};
						$conflict_queries .= "$qname:GAP:M$qmatch;";
					}
					foreach my $qname(@snp_qnames){
						my $qmatch = $hash_query{$qname}{match};
						$conflict_queries .= "$qname:SNP:M$qmatch;";
					}
					foreach my $qname(@conf_qnames){
						my $qmatch = $hash_query{$qname}{match};
						$conflict_queries .= "$qname:INDEL:M$qmatch;";
					}
					$conflict_queries =~ s/;$//;
					if($conflict_queries eq ""){
						$conflict_queries = "-";
					}
					print OUT "$tname\t$pos\tINDEL\t$tgap\t$qgap\t$support_queries\t$conflict_queries\n";
				}
			}
		}
	}
}

close OUT;

sub check_array{
	my($query,@targets) = @_;
	my %hash = map {$targets[$_], ""} 0..$#targets;
	my $check = 0;
	if(exists $hash{$query}){
		$check = 1;
	}
	return($check);
}

sub throw_item_from_array{
	my($query,@targets) = @_;
	for(my $i = 0; $i < @targets; $i++){
		if($targets[$i] eq $query){
			splice(@targets,$i,1);
			$i--;
		}
	}
	return(@targets);
}