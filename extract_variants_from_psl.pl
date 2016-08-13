# For spliced genome assembly
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
#	Version 1.0 at 2016.03.07
#	
#########################################################

=note
The first time I started my work on using BLAT to find variants in genes from assembly results, I aligned the reference sequences of the genes to the whole genome assembled contigs with BLAT.
So the former versions of program (psl2vcf.pl) to extract variants from BLAT results used the reference sequence (where you anchor the variants) as query sequence and the contigs as the target (or database).
But things have changed now when I need to find variants from a local assembly result against the reference chromosome, BLAT takes too long when using the large chromosome file as the query.
Here is a new program using the reference sequence as target (or database) and contigs as query.
=cut

use strict;
use warnings;
use Getopt::Long;

my($psl,$qref,$tref,$outname);
my $usage = "USAGE:\nperl $0 --psl <psl file> --qref <query refence seq> --tref <target reference seq> --out <outname>\n";
GetOptions(
	"psl=s" => \$psl,
	"qref=s" => \$qref,
	"tref=s" => \$tref,
	"out=s" => \$outname,
) or die $usage;

unless(defined $psl and defined $qref and defined $tref and defined $outname){
	die $usage;
}

# ---------------------- PART I: Reading the reference sequences -----------------------

my %hash_qref;
my %hash_tref;

print "Reading query reference...\n";
open(QREF,"<$qref") or die "$!";
local $/ = "\n>";
while(<QREF>){
	$_ =~ s/^>//g;
	chomp;
	my($head,$seq) = split/\n/,$_,2;
	my($seqname) = $head =~ /^(\S+)/;
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

my %hash_cov;
my %hash_pslcount;

open(TREF,"<$tref") or die "$!";
local $/ = "\n>";
while(<TREF>){
	$_ =~ s/^>//g;
	chomp;
	my($head,$seq) = split/\n/,$_,2;
	my($seqname) = $head =~ /^(\S+)/;
	$seq =~ s/\s+//g;
	$seq = uc($seq);
	$hash_tref{$seqname}{seq} = $seq;
	$hash_tref{$seqname}{len} = length($seq);
	for(my $i = 1; $i <= length($seq); $i++){
		$hash_cov{$seqname}{$i}{absolute} = 0;
		$hash_cov{$seqname}{$i}{relative} = 0;
	}
	$hash_pslcount{$seqname} = 0;
}
local $/ = "\n";
close TREF;

# --------------------- PART I END -----------------------

# --------------------- PART II: Reading psl file -----------------------

print "Reading PSL file...\n";

open(PSL,"<$psl") or die "$!";

my $tname_tmp = "";
my @covered_regions;
my %hash_vcf;
open(UNC,">unused.list");

while(<PSL>){
	chomp;
	my $inline = $_;
	next unless ($inline =~ /^\d+/);
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$inline;
	my $uncovered = 0;
	#my($flankup,$flankdown) = split/,/,$flanks;
	
	if($tname ne $tname_tmp){
		@covered_regions = ();
		$tname_tmp = $tname;
	}
	$hash_pslcount{$tname}++;
	
	for(my $i=$tstart+1;$i<=$tend;$i++){
		$hash_cov{$tname}{$i}{relative}++;
	}
	
	my @qbs = split/,/,$qbstarts;
	my @tbs = split/,/,$tbstarts;
	my @blks = split/,/,$blocksizes;
	
	for(my $count = 0; $count < $blknum; $count++){
		
		my $tseq = substr($hash_tref{$tname}{seq},$tbs[$count],$blks[$count]);
		my $qseq;
		if($strand eq "+"){
			$qseq = substr($hash_qref{$qname}{seq},$qbs[$count],$blks[$count]);
		}else{
			$qseq = substr($hash_qref{$qname}{revcom},$qbs[$count],$blks[$count]);
		}
		
		#coverage count
		for(my $i=$tbs[$count]+1;$i<=$tbs[$count]+$blks[$count];$i++){
			my $tpos =  $i;
			# Optional: ignore overlap region.
			foreach my $region (@covered_regions){
				my($rstart,$rend) = split/,/,$region;
				if($tpos > $rstart and $tpos <= $rend){
					goto THROW1;
				}
			}
			# Optional: ignore overlap region.
			$hash_cov{$tname}{$i}{absolute}++; 
			THROW1:
		}
		
		#dealing with mismatches -> snps
		for(my $i=0;$i<$blks[$count];$i++){
			my $qbase = substr($qseq,$i,1);
			my $tbase = substr($tseq,$i,1);
			my $tpos = $tbs[$count] + 1 + $i;
			
			
			# Optional: ignore overlap region.
			my $covered = 0;
			foreach my $region (@covered_regions){
				my($rstart,$rend) = split/,/,$region;
				if($tpos > $rstart and $tpos <= $rend){
					$covered = 1;
				}
			}
			next if($covered == 1);
			# Optional: ignore overlap region.
			$uncovered++;
			
			if($qbase ne $tbase){
				$hash_vcf{$tname}{$tpos}{REF} = $tbase; 
				$hash_vcf{$tname}{$tpos}{ALT} = $qbase; 
			}
=cu
			# find putative indel by flanking region:
			if($flankup > 0){
				if($count == 0 and $i == 0){
					$hash_vcf{$tname}{$tbs[$count]}{REF} = 'N';
					$hash_vcf{$tname}{$tbs[$count]}{ALT} = 'NNNN';
					$hash_cov{$tname}{$tbs[$count]}{relative}++;
				}
			}
			if($flankdown > 0){
				if($count == $blknum-1 and $i == $blks[$count]-1){
					$hash_vcf{$tname}{$tpos+1}{REF} = 'N';
					$hash_vcf{$tname}{$tpos+1}{ALT} = 'NNNN';
					$hash_cov{$tname}{$tpos+1}{relative}++;
				}
			}
=cut
		}
		
		
		#dealing with gap -> indel
		if($count < $blknum - 1){
			
			# changed 20150717 by zhuochen, for a extremely strange situation in which q and t block end bases is not identical
			# This problem first found when blating Pi21-Kasalath region sequence to Y256 unitig.
			WRONG_BLOCK_END:
			my $tgapstart = $tbs[$count] + $blks[$count] + 1;
			my $tgapend = $tbs[$count + 1];
			my $qgapstart = $qbs[$count] + $blks[$count] + 1;
			my $qgapend = $qbs[$count + 1];
			
			my $tpos = $tbs[$count] + $blks[$count];
			
			# Optional: ignore overlap region.
			foreach my $region(@covered_regions){
				my($rstart,$rend) = split/,/,$region;
				if(($tgapstart > $rstart and $tgapstart <= $rend)
					or ($tgapend > $rstart and $tgapend <= $rend)
					or ($tgapstart <= $rstart and $tgapend > $rend)
					or (exists $hash_vcf{$tname}{$tpos})){
					for(my $i = $tgapstart; $i <= $tgapend; $i++){
						$hash_cov{$tname}{$i}{relative}--;
					}
					goto THROW2;
				}
			}
			# Optional: ignore overlap region.
			
			# get the last bases of the q and t block 
			my $qgap = substr($qseq,$blks[$count]-1,1);
			my $tgap = substr($tseq,$blks[$count]-1,1);
			
			# check if the last bases are identical; if not, 
			if($qgap ne $tgap){
				delete($hash_vcf{$tname}{$tpos});
				$hash_cov{$tname}{$tpos}{absolute} --; 
				$blks[$count]--;
				goto WRONG_BLOCK_END;
			}
			
			my $tgapbases = substr($hash_tref{$tname}{seq},$tgapstart - 1, $tgapend - $tgapstart + 1);
			$tgap .= $tgapbases;
			
			my $qgapbases;
			if($strand eq "+"){
				$qgapbases = substr($hash_qref{$qname}{seq}, $qgapstart - 1, $qgapend - $qgapstart + 1);
			}else{
				$qgapbases = substr($hash_qref{$qname}{revcom}, $qgapstart - 1, $qgapend - $qgapstart + 1);
			}
			$qgap .= $qgapbases;
			
			$hash_vcf{$tname}{$tpos}{REF} = $tgap;
			$hash_vcf{$tname}{$tpos}{ALT} = $qgap;
		}
		THROW2:
		
	}
	if($uncovered == 0){
		print UNC "$qname\n";
	}
	push @covered_regions,"$tstart,$tend";
}
close PSL;
print "Done!\n";
# --------------------- PART II END -----------------------

open(COV,">$outname.coverage");
print COV "chr\tposition\tdepth_absolute\tdepth_relative\n";

foreach my $tname (keys %hash_cov){
	foreach my $pos (sort{$a<=>$b} keys %{$hash_cov{$tname}}){
		unless(exists $hash_cov{$tname}{$pos}{absolute}){
			print "$tname\t$pos\tno absolute\n";
			unless(exists $hash_cov{$tname}{$pos}{relative}){
				die "$tname\t$pos\tno relative\n";
			}else{
				die;
			}
		}
		print COV "$tname\t$pos\t$hash_cov{$tname}{$pos}{absolute}\t$hash_cov{$tname}{$pos}{relative}\n";
	}
}
close COV;

open(VCF,">$outname.vcf");
print VCF "#CHROM\tPOS\tREF\tALT\n";
foreach my $tname (keys %hash_vcf){
	foreach my $tpos (sort{$a<=>$b} keys %{$hash_vcf{$tname}}){
		my $qbase = $hash_vcf{$tname}{$tpos}{ALT};
		print VCF "$tname\t$tpos\t$hash_vcf{$tname}{$tpos}{REF}\t$qbase\n"; #20150318
	}
}
close VCF;
