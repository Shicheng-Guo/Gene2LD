#!	/usr/bin/perl	-w

=head1	header

	###########################################################
	#  
	#  This script extracts UCSC genes (transcripts) included in the pre-defined LD regions of the SNPs
	#  
	#  perl AssignGeneinLDblock.pl LDrange.txt UCSC_RefGene_Transcript_Position.txt
	#  
	#  Input : LDrange.txt UCSC_RefGene_Transcript_Position.txt
	#  Output : LDrange-UCSC_RefGene_Transcript_Position.txt
	#  
	#  LDrange.txt could be generated from output file of LDRangeAssign.R
	#  Col1 : SNP ID
	#  Col2 : SNP position
	#  Col3 : LD region start position
	#  Col4 : LD region end position
	#  Col5 : LD region length
	#  
	#  Position information is based on ChrNum*100000000000+BasePaierPosition
	#  
	#  Any questions to Yukinori Okada (http://plaza.umin.ac.jp/~yokada/datasource/software.htm   yokada@broadinstitute.org)

=cut

use	strict;


&main(@ARGV);
exit;

=head1	main
=cut

sub	main(@){
	my ($SNPPosifile, $GenePosifile)=@_;
	my $prefix1 = $SNPPosifile;
	my $prefix2 = $GenePosifile;
	for (my $i=0;$i<4;$i++)	{
		chop($prefix1);
		chop($prefix2);
	}
	my $outfile = $prefix1."-".$prefix2.".txt";
	
	my $SNP;
	my $SNPChr;
	my $SNPPosi;
	my $SNPPosiStr;
	my $SNPPosiEnd;
	my $numtmpGene;
	my $tmpGene;
	my $tmpSNPPosiStr;
	my $tmpSNPPosiEnd;
	
	my @Gene;
	my @GeneChr;
	my @GenePosiStr;
	my @GenePosiEnd;
	
	my $counter = 0;
	my $out = "";
	my $base = 100000000000;

	open (INPUT1, "$GenePosifile");
	while(<INPUT1>){
		chomp;
		my $last = substr $_,-1;
		if ($last eq "\r") {
			substr ($_,-1) = "";
		}
		my @inline = split(/\t/);
		$Gene[$counter] = $inline[0]."-".$inline[1];
		$GeneChr[$counter] = $inline[2];
		$GenePosiStr[$counter] = $inline[3];
		$GenePosiEnd[$counter] = $inline[4];
		$counter++;
	}
	close INPUT1;

	open (OUT, "> $outfile");
	$out = "SNPID\tOrigianlLDregionStart\tOriginalLDregionEnd\tlength\tLD+GeneregionStart\tLD+GeneregionEnd\tlength\tNo.Genes(transcripts)inRegion\tGenes(transcripts)inRegion\tNearestGenes(transcripts)inRegion\n";
	print OUT $out;
	
	open (INPUT2, "$SNPPosifile");
	while(<INPUT2>){
		chomp;
		my $last = substr $_,-1;
		if ($last eq "\r") {
			substr ($_,-1) = "";
		}
		my @inline = split(/\t/);
		$SNP = $inline[0];
		$SNPPosi = $inline[1];
		$SNPPosiStr = $inline[2];
		$SNPPosiEnd = $inline[3];
		$SNPChr = ($inline[2]-$inline[2]%$base)/$base;
		
		$numtmpGene = 0;
		$tmpGene ="";
		$tmpSNPPosiStr = $SNPPosiStr;
		$tmpSNPPosiEnd = $SNPPosiEnd;
		
		my $Dis = 999999999;
		my @NearGene;
		my $numNear;
		
		for (my $i=0;$i<$counter;$i++) {
				if (!($GenePosiStr[$i]<$SNPPosiStr && $GenePosiEnd[$i]<$SNPPosiStr) && !($GenePosiStr[$i]>$SNPPosiEnd && $GenePosiEnd[$i]>$SNPPosiEnd)) {
					$numtmpGene++;
					$tmpGene .= $Gene[$i].",";
					if ($tmpSNPPosiStr>$GenePosiStr[$i]) {
						$tmpSNPPosiStr=$GenePosiStr[$i];
					}
					if ($tmpSNPPosiEnd<$GenePosiEnd[$i]) {
						$tmpSNPPosiEnd=$GenePosiEnd[$i];
					}
				}
				
				my $tmpDis;
				my $tmpDisStr = abs($GenePosiStr[$i]-$SNPPosi);
				my $tmpDisEnd = abs($GenePosiEnd[$i]-$SNPPosi);
				if ($tmpDisStr <= $tmpDisEnd) {
					$tmpDis = $tmpDisStr;
				} else {
					$tmpDis = $tmpDisEnd;
				}
				if ($GenePosiStr[$i]<$SNPPosi && $GenePosiEnd[$i]>$SNPPosi) {
					$tmpDis = 0;
				}
				if ($tmpDis<$Dis) {
					$Dis = $tmpDis;
					@NearGene = ();
					$numNear = 1;
					$NearGene[0] = $Gene[$i];
				} elsif ($tmpDis==$Dis) {
					$NearGene[$numNear] = $Gene[$i];
					$numNear++;
				}
		}
		
		if ($tmpGene eq "") {
			$tmpGene = "-";
		} else {
			chop($tmpGene);
		}
		
		my $outNear = "";
		if ($numNear>=1) {
			for (my $j=0;$j<@NearGene;$j++) {
				$outNear .= $NearGene[$j].",";
			}
		} else {
			$outNear ="-";
		}
		
		$out = $SNP."\t".$SNPPosiStr."\t".$SNPPosiEnd."\t".($SNPPosiEnd-$SNPPosiStr)."\t".$tmpSNPPosiStr."\t".$tmpSNPPosiEnd."\t".($tmpSNPPosiEnd-$tmpSNPPosiStr)."\t".$numtmpGene."\t".$tmpGene."\t".$outNear."\n";
		print OUT $out;
		
	}
	close INPUT2;
	close OUT;
	
	
}
