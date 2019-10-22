#########################################################
#	extract_variants_from_psl.pl
#	
#	1. Extracting snps and indels from BLAT psl output file 
#	2. Counting coverage of each base in query
#	
#	By Zhuo CHEN, IGDB, CAS
#	zhuochen@genetics.ac.cn
#	chenomics@163.com
#	
#########################################################

=anno
The first time I started my work on using BLAT to find variants in genes from assembly results, I aligned the reference sequences of the genes to the whole genome assembled contigs with BLAT.
So the former versions of program (psl2vcf.pl) to extract variants from BLAT results used the reference sequence (where you anchor the variants) as query sequence and the contigs as the target (or database).
But things have changed now when I need to find variants from a local assembly result against the reference chromosome, BLAT takes too long when using the large chromosome file as the query.
So here is a new program using the reference sequence as target (or database) and contigs as query.
=cut

use strict;
use warnings;
use Getopt::Long;

my($psl,$qref,$tref,$outname,$gapmax);
my $usage = "USAGE:\nperl $0 --psl <psl file> --qref <query refence seq> --tref <target reference seq> --out <outname> --maxgap <ma gap size>\n";
$usage .= "--maxgap set the max gap size to call from the psl. Default=100\n";

GetOptions(
	"psl=s" => \$psl,
	"qref=s" => \$qref,
	"tref=s" => \$tref,
	"out=s" => \$outname,
	"maxgap=s" => \$gapmax,
) or die $usage;

unless(defined $psl and defined $qref and defined $tref and defined $outname){
	die $usage;
}

unless(defined $gapmax){
	$gapmax = 100;
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
my %hash_ctgs;

while(<PSL>){
	chomp;
	my $inline = $_;
	next unless ($inline =~ /^\d+/);
	my($match,$mismatch,undef,undef,undef,undef,undef,undef,$strand,
	$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,
	$blknum,$blocksizes,$qbstarts,$tbstarts,$flanks) = split/\t/,$inline;
	
	next if(exists $hash_ctgs{$qname});
	$hash_ctgs{$qname}{match} = $match;
	
	my @qbs = split/,/,$qbstarts;
	my @tbs = split/,/,$tbstarts;
	my @blks = split/,/,$blocksizes;
	
	#my($flankup,$flankdown) = split/,/,$flanks;
	
	$hash_pslcount{$tname}++;
	
	# set absolute cover
	for(my $count = 0; $count < $blknum; $count++){
		for(my $i=$tbs[$count]+1;$i<=$tbs[$count]+$blks[$count];$i++){
			my $tpos =  $i;
			#next if(exists $hash_cov{$tname}{$tpos}{final});
			$hash_cov{$tname}{$i}{absolute}++; 
		}
	}
	
	# set relative cover
	for(my $i=$tstart+1;$i<=$tend;$i++){
		#next if(exists $hash_cov{$tname}{$i}{final});
		$hash_cov{$tname}{$i}{relative}++;
	}
	
	for(my $count = 0; $count < $blknum; $count++){
		
		my $tseq = substr($hash_tref{$tname}{seq},$tbs[$count],$blks[$count]);
		my $qseq;
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
			
			#next if(exists $hash_cov{$tname}{$tpos}{final});
			if(exists $hash_vcf{$tname}{$tpos}{REF}){
				my $base = $hash_vcf{$tname}{$tpos}{REF}{base};
				if($tbase ne $base){
					die "# ERROR: multiple tbase type at ref:$tname pos:$tpos\n";
				}
			}
			$hash_vcf{$tname}{$tpos}{REF}{base} = $tbase;
			
			if($qbase ne $tbase){
				#$hash_vcf{$tname}{$tpos}{ALT} = $qbase; 
				$hash_vcf{$tname}{$tpos}{SNP}{$qbase}{CTG}{$qname} = "";
			}else{
				# ref type contig
				$hash_vcf{$tname}{$tpos}{REF}{CTG}{$qname} = "";
				# ref type contig with out gap
				$hash_vcf{$tname}{$tpos}{REF}{GCTG}{$qname} = "";
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
			
			my $tgaplen = $tgapend - $tgapstart + 1;
			my $qgaplen = $qgapend - $qgapstart + 1;
			
			my $tpos = $tbs[$count] + $blks[$count];
			
			
			# Set the conditions the abandon a gap
			my $throwgap = 0;
			if($tgaplen > $match or $qgaplen > $match or $tgaplen > $gapmax or $qgaplen > $gapmax){
				$throwgap = 1;
			}
=cu
			if(exists $hash_cov{$tname}{$tgapstart - 1}{final}){
				$throwgap = 1;
			}else{
				for(my $i = $tgapstart; $i <= $tgapend; $i++){
					if(exists $hash_cov{$tname}{$i}{final}){
						$throwgap = 1;
						last;
					}
				}
			}
=cut
			if($throwgap == 1){
				for(my $i = $tgapstart; $i <= $tgapend; $i++){
					#unless(exists $hash_cov{$tname}{$i}{final}){
						$hash_cov{$tname}{$i}{relative} --;
					#}
				}
				goto THROW2;
			}
			
			# get the last bases of the q and t block 
			my $qgap = substr($qseq,$blks[$count]-1,1);
			my $tgap = substr($tseq,$blks[$count]-1,1);
			
			# check if the last bases are identical; if not, 
			if($qgap ne $tgap){
				delete($hash_vcf{$tname}{$tpos}{SNP}{$qgap}{CTG}{$qname});
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
			
			delete($hash_vcf{$tname}{$tpos}{REF}{GCTG}{$qname});
			$hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{CTG}{$qname} = "";
		}
		THROW2:
		
	}
	
	# set final cover
	for(my $i=$tstart+1;$i<=$tend;$i++){
		$hash_cov{$tname}{$i}{final} = 1;
	}
}
close PSL;
print "Done!\n";
# --------------------- PART II END -----------------------

open(COV,">$outname.covab");
print COV "#CHROM\tTYPE\tSTART\tEND\n";
open(RCOV,">$outname.covre");
print RCOV "#CHROM\tTYPE\tSTART\tEND\n";

foreach my $tname (keys %hash_cov){
	my $gapstart = "";
	my $covstart = "";
	my $gapstart_re = "";
	my $covstart_re = "";
	my $endpos = "";
	foreach my $pos (sort{$a<=>$b} keys %{$hash_cov{$tname}}){
		unless(exists $hash_cov{$tname}{$pos}{absolute}){
			print "$tname\t$pos\tno absolute\n";
			unless(exists $hash_cov{$tname}{$pos}{relative}){
				die "$tname\t$pos\tno relative\n";
			}else{
				die;
			}
		}
		my $depth = $hash_cov{$tname}{$pos}{absolute};
		my $depth_re = $hash_cov{$tname}{$pos}{relative};
		
		if($depth == 0){
			if($gapstart eq "" and $covstart eq ""){
				$gapstart = $pos;
			}elsif($gapstart eq ""){
				$gapstart = $pos;
				my $covend = $pos - 1;
				print COV "$tname\tcov\t$covstart\t$covend\n";
				$covstart = "";
			}
		}else{
			if($gapstart eq "" and $covstart eq ""){
				$covstart = $pos;
			}elsif($covstart eq ""){
				$covstart = $pos;
				my $gapend = $pos - 1;
				print COV "$tname\tgap\t$gapstart\t$gapend\n";
				$gapstart = "";
			}
		}
		
		if($depth_re == 0){
			if($gapstart_re eq "" and $covstart_re eq ""){
				$gapstart_re = $pos;
			}elsif($gapstart_re eq ""){
				$gapstart_re = $pos;
				my $covend_re = $pos - 1;
				print RCOV "$tname\tcov\t$covstart_re\t$covend_re\n";
				$covstart_re = "";
			}
		}else{
			if($gapstart_re eq "" and $covstart_re eq ""){
				$covstart_re = $pos;
			}elsif($covstart_re eq ""){
				$covstart_re = $pos;
				my $gapend_re = $pos - 1;
				print RCOV "$tname\tgap\t$gapstart_re\t$gapend_re\n";
				$gapstart_re = "";
			}
		}
		
		$endpos = $pos;
		#print COV "$tname\t$pos\t$hash_cov{$tname}{$pos}{absolute}\t$hash_cov{$tname}{$pos}{relative}\n";
	}
	if($gapstart eq "" and $covstart eq ""){
		print STDERR "# ERROR: both gapstart and covstart are empty at the end of chr:$tname\n";
	}elsif($gapstart eq ""){
		print COV "$tname\tcov\t$covstart\t$endpos\n";
	}else{
		print COV "$tname\tgap\t$gapstart\t$endpos\n";
	}
	
	if($gapstart_re eq "" and $covstart_re eq ""){
		print STDERR "# ERROR: both gapstart_re and covstart_re are empty at the end of chr:$tname\n";
	}elsif($gapstart_re eq ""){
		print RCOV "$tname\tcov\t$covstart_re\t$endpos\n";
	}else{
		print RCOV "$tname\tgap\t$gapstart_re\t$endpos\n";
	}
}
close COV;
close RCOV;

open(VCF,">$outname.vcf");
print VCF "#CHROM\tPOS\tTYPE\tREF\tALT\tREFLEN\tALTLEN\tALTCTG\n";
foreach my $tname (keys %hash_vcf){
	foreach my $tpos (sort{$a<=>$b} keys %{$hash_vcf{$tname}}){
		if(exists $hash_vcf{$tname}{$tpos}{SNP}){
			my $refbase = $hash_vcf{$tname}{$tpos}{REF}{base};
			my $reflen = "-";
			
			if(exists $hash_vcf{$tname}{$tpos}{REF}{CTG}){
				foreach my $contig(keys %{$hash_vcf{$tname}{$tpos}{REF}{CTG}}){
					my $match = $hash_ctgs{$contig}{match};
					if($reflen ne "-"){
						next unless($match > $reflen);
					}
					$reflen = $match;
				}
			}
			
			foreach my $altbase(keys %{$hash_vcf{$tname}{$tpos}{SNP}}){
				foreach my $contig(keys %{$hash_vcf{$tname}{$tpos}{SNP}{$altbase}{CTG}}){
					my $match = $hash_ctgs{$contig}{match};
					if(exists $hash_vcf{$tname}{$tpos}{SNP}{$altbase}{MLEN}){
						my $mlen = $hash_vcf{$tname}{$tpos}{SNP}{$altbase}{MLEN};
						next unless($match > $mlen);
					}
					$hash_vcf{$tname}{$tpos}{SNP}{$altbase}{MLEN} = $match;
					$hash_vcf{$tname}{$tpos}{SNP}{$altbase}{LCTG} = $contig;
				}
				unless(exists $hash_vcf{$tname}{$tpos}{SNP}{$altbase}{LCTG}){
					delete($hash_vcf{$tname}{$tpos}{SNP}{$altbase});
				}
			}

			my @altbases = sort {$hash_vcf{$tname}{$tpos}{SNP}{$b}{MLEN} <=> $hash_vcf{$tname}{$tpos}{SNP}{$a}{MLEN}} keys %{$hash_vcf{$tname}{$tpos}{SNP}};
			next if(@altbases == 0);
			my @altlens = ();
			my @altctgs = ();
			foreach my $altbase(@altbases){
				my $mlen = $hash_vcf{$tname}{$tpos}{SNP}{$altbase}{MLEN};
				my $contig = $hash_vcf{$tname}{$tpos}{SNP}{$altbase}{LCTG};
				push @altlens, $mlen;
				push @altctgs, $contig;
			}
			print VCF "$tname\t$tpos\tSNP\t$refbase\t".join(",",@altbases)."\t$reflen\t".join(",",@altlens)."\t".join(",",@altctgs)."\n";
		}
		next unless(exists $hash_vcf{$tname}{$tpos}{INDEL});
		
		my $reflen = "-";
		if(exists $hash_vcf{$tname}{$tpos}{REF}{GCTG}){
			foreach my $contig(keys %{$hash_vcf{$tname}{$tpos}{REF}{GCTG}}){
				my $match = $hash_ctgs{$contig}{match};
				if($reflen ne "-"){
					next unless($match > $reflen);
				}
				$reflen = $match;
			}
		}
		foreach my $tgap(keys %{$hash_vcf{$tname}{$tpos}{INDEL}{REF}}){
			foreach my $qgap(keys %{$hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}}){
				foreach my $contig(keys %{$hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{CTG}}){
					my $match = $hash_ctgs{$contig}{match};
					if(exists $hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{MLEN}){
						my $mlen = $hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{MLEN};
						next unless($match > $mlen);
					}
					$hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{MLEN} = $match;
					$hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{LCTG} = $contig;
				}
			}
			
			my @qgaps = sort {$hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$b}{MLEN} <=> $hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$a}{MLEN}} keys %{$hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}};
			my @altlens = ();
			my @altctgs = ();
			foreach my $qgap(@qgaps){
				my $mlen = $hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{MLEN};
				my $contig = $hash_vcf{$tname}{$tpos}{INDEL}{REF}{$tgap}{ALT}{$qgap}{LCTG};
				push @altlens, $mlen;
				push @altctgs, $contig;
			}
			print VCF "$tname\t$tpos\tInDel\t$tgap\t".join(",",@qgaps)."\t$reflen\t".join(",",@altlens)."\t".join(",",@altctgs)."\n";
		}
	}
}
close VCF;
