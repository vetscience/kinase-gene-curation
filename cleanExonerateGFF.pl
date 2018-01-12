#!/usr/local/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use POSIX;
my %gff;
my %exonCnt;
use Storable 'dclone';

open (GFFIN, '<', $ARGV[0]) or die "Cannot open file $ARGV[0]: $!";

while (<GFFIN>)
{

	if(!($_ =~ /^#/))
	{
		chomp;
		my @inline=split(/\t/,$_);
		my $id;
		my $scaffold=$inline[0];
		#my $source=".";
		my $source=$inline[1];
		my $type=$inline[2];
		my $start=$inline[3];
		my $end=$inline[4];
		my $score=$inline[5];
		my $strand;

		if($ARGV[2] eq "f")
		{	
			my %flip = ("+","-","-","+");
			$strand=$flip{$inline[6]};
			#print $strand."\n";	
		}
		else
		{
			$strand=$inline[6];
		}
		
		my $phase=$inline[7];

		my %attrH;
		my @attrA=split(/;/,$inline[8]);

		foreach my $attr (@attrA)
		{
			#print "Attr: ".$attr."\n";
			my ($keyAttr, $valueAttr) = split(/=/,$attr);
			$attrH{$keyAttr."="} = $valueAttr;
			#print $keyAttr."="."\t".$valueAttr."\n";
		}

		my %line = (scaffold => $scaffold, source => $source, type => $type, start => $start, end => $end, score => $score, strand => $strand, phase => $phase, attr => \%attrH);
		
		if(exists($line{attr}->{"ID="}))
		{
			$id = $line{attr}->{"ID="};
		}
		else
		{
			$id = $line{attr}->{"Parent="};
			$line{attr}->{"ID="} = $id;
		}
		

		#for each found mRNA entry create a CDS and gene entry
		if($type eq "mRNA")
		{
			
			$line{attr}{"Parent="} = "gene:".$id;
		 	#print "before ".$line{attr}{"ID="}."\n";
		 	$line{attr}->{"ID="} = "transcript:".$id;
		 	#print "after ".$line{attr}{"ID="}."\n\n";
		 	my $mRNA = dclone \%line;

		 	$line{type} = "gene";
		 	$line{attr}->{"ID="} = "gene:".$id;
			delete($line{attr}->{"Parent="});
		 	my $gene = dclone \%line;

		 	# $line{type} = "CDS";
		 	# $line{attr}{"Parent="} = "transcript:".$id;
		 	# $line{attr}->{"ID="} = "cds:".$id;
		 	# $line{phase} = $ARGV[1];
		 	# my $CDS = dclone \%line;

		 	my %features = (gene => $gene, mRNA => $mRNA);
			$gff{$id}{features} = \%features;
	 	}
	 	else
	 	{
		 	if($type eq "exon")
		 	{
		 		if(!(exists($exonCnt{$id})))
		 		{
		 			$exonCnt{$id} = 0;
		 			my @exonA;
		 			$gff{$id}{exons} = \@exonA;	
		 		}

		 		$exonCnt{$id}++;
		 		$line{attr}{"Parent="} = "transcript:".$id;
		 		$line{attr}->{"ID="} = $id;
		 		#$line{attr}->{"ID="} = "exon:".$id.".".$exonCnt{$id};

		 		push(@{$gff{$id}{exons}},\%line);
		 	}
	 	}
	}
}
close GFFIN;

print "##gff-version 3\n";

foreach my $ids (keys(%gff))
{

	print "# Gene gene:".$ids."\n";

	foreach my $keys ("gene","mRNA")
	{
		my $attr;

		if(exists($gff{$ids}{features}->{$keys}{attr}->{"Parent="}))
		{
		$attr = join(";", "ID=".$gff{$ids}{features}->{$keys}{attr}->{"ID="}, "Parent=".$gff{$ids}{features}->{$keys}{attr}->{"Parent="});
		}
		else
		{
		$attr = "ID=".$gff{$ids}{features}->{$keys}{attr}->{"ID="};
		}

		#print "test".$gff{$ids}{features}->{$keys}{scaffold}."\n";

		print	join("\t", $gff{$ids}{features}->{$keys}{scaffold},
					$gff{$ids}{features}->{$keys}{source},
					$gff{$ids}{features}->{$keys}{type},
					$gff{$ids}{features}->{$keys}{start},
					$gff{$ids}{features}->{$keys}{end},
					$gff{$ids}{features}->{$keys}{score},
					$gff{$ids}{features}->{$keys}{strand},
					$gff{$ids}{features}->{$keys}{phase},$attr);
		print "\n";
			
	}


	if($ARGV[2] eq "r")
	{
		@{$gff{$ids}{exons}} = reverse(@{$gff{$ids}{exons}});
	}

	my $count = 0;
	#foreach my $exon (reverse(@{$gff{$ids}{exons}}))
	foreach my $exon (@{$gff{$ids}{exons}})
	{
		$count++;
		my $attr;
		
		$attr = join(";", "ID=exon:".$exon->{attr}->{"ID="}."."."$count", "Parent=".$exon->{attr}->{"Parent="});

		print	join("\t", $exon->{scaffold},
					$exon->{source},
					$exon->{type},
					$exon->{start},
					$exon->{end},
					$exon->{score},
					$exon->{strand},
					$exon->{phase},$attr);
		print "\n";
	}

	$count = 0;
	my $nextPhase;

	foreach my $exon (@{$gff{$ids}{exons}})
	{
		$count++;
		my $attr;
		my $phase;

		if($count == 1)
		{
			$phase = $ARGV[1];
			$nextPhase = (abs($exon->{end}-$exon->{start})-$phase)%3
		}
		else
		{
			$phase = 3-$nextPhase;
			if ($phase == 3)
			{
				$phase = 0;
			}
			$nextPhase = (abs($exon->{end}-$exon->{start})-$phase)%3
		}
		
		$attr = join(";", "ID=cds:".$exon->{attr}->{"ID="}, "Parent=".$exon->{attr}->{"Parent="});

		print	join("\t", $exon->{scaffold},
					$exon->{source},
					"CDS",
					$exon->{start},
					$exon->{end},
					$exon->{score},
					$exon->{strand},
					$phase,$attr);
		print "\n";
	}
print "###\n";
}

#print Dumper %gff;