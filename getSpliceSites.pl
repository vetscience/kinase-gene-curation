#!/usr/local/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use POSIX;
my %opts;
my %scaffolds;

getopts('a:g:t:f:', \%opts);


if(!defined($opts{a}) || !defined($opts{g}) || !defined($opts{t}) || !defined($opts{f}))
{
	print STDERR "Usage: getSpliceSites.pl
	-a <properly formatted GFF file>
	-g <associated genome scaffold file>
	-t <prefix for 3' output file (Exonerate .pssm format)>
	-f <prefix for 5' output file (Exonerate .pssm format)>\n";
	exit 1;
}

my $current;

open (GENIN, '<', $opts{g}) or die "Cannot open file $opts{g}: $!";
while (<GENIN>)
{
	chomp;
	if($_ =~ m/^>/)
	{
		my @header=split(/[>\s\t]/,$_);
		my %features = (start => 0, end => 0, seq => "");
		$current = $header[1];
		$scaffolds{$current} = \%features;
	}
	else
	{
		$scaffolds{$current}->{seq} = $scaffolds{$current}->{seq}.$_;
		$scaffolds{$current}->{end} = (length($scaffolds{$current}->{seq})-1);
	}
}
close GENIN;

#Read in GFF3 file (assume properly sorted)
#number of exons within current mRNA stream
my $exoncnt = 0;
#strand of previous exon
my $prevstr = '';
#stack of 5nt strings that represent 3' splice (mainly AG, beginning of exon) sites
my @stack3;
#stack of 5nt strings that represent 5' splice (mainly GT, end of exon) sites
my @stack5;

open (GFFIN, '<', $opts{a}) or die "Cannot open file $opts{a}: $!";
while (<GFFIN>)
{
	chomp;
	my @inline = split(/\t/,$_);
	if (scalar(@inline) == 9)
	{
		my $scaff = $inline[0];
		my $strand = $inline[6];
		my $type = $inline[2];
		my $start = $inline[3]-1;
		#print "start gff: ".$start."\n";
		my $end = $inline[4]-1;
		#rint "end gff: ".$end."\n";

		if($type eq "mRNA")
		{
			if($exoncnt == 1)
			{
				pop(@stack5);
				pop(@stack3);
			}
			$exoncnt = 0;

			if($prevstr eq "-")
			{
				pop(@stack5);	
			}
			elsif($prevstr eq "+")
			{
				pop(@stack3);	
			}	
		}

		if($type eq "exon")
		{

			# if(($end-$start)<=19)
			# {
			# 	#print "Super short exon\n";
			# }

			$exoncnt++;
			my $startseq;
			my $endseq;

			if(($start-5) < $scaffolds{$scaff}->{start})
			{
				$startseq = substr($scaffolds{$scaff}->{seq},0,($start+5));

				if(length($startseq) == 11){print "Warning start ".$start." ".($start-5)." ".$scaffolds{$scaff}->{start}." ".$startseq.": ".length($startseq)."\n"}
			}
			else
			{
				$startseq = substr($scaffolds{$scaff}->{seq},($start-5),11);
			}

			if(($end+5) > $scaffolds{$scaff}->{end})
			{
				$endseq = substr($scaffolds{$scaff}->{seq},($end-5),($scaffolds{$scaff}->{end}-$end));

				if(length($endseq) == 11){print "Warning end ".$end." ".($end+5)." ".$scaffolds{$scaff}->{end}." ".$endseq.": ".length($endseq)."\n"}
			}
			else
			{
				$endseq = substr($scaffolds{$scaff}->{seq},($end-5),11);
			}

			
			my $splice5;
			my $splice3;

			if($strand eq "-")
			{
				#print "end before: $endseq\n\n";
				$endseq =~ tr/ACGTacgt/TGCAtgca/;
				$startseq =~ tr/ACGTacgt/TGCAtgca/;
				$splice3 = reverse($endseq);
				$splice5 = reverse($startseq);
				#print "end after: $splice3\n\n";

				# if(length($splice5) == 11)
				# {
				# 	print substr($splice5,6,2)."\n"; #TODO GT
				# }


				# if(length($splice3) == 11)
				# {
				# 	print substr($splice3,3,2)."\n"; # TODO AG
				# }
				
			}
			elsif($strand eq "+")
			{
				$splice3 = $startseq;
				$splice5 = $endseq;
				# if(length($splice5) == 11)
				# {
				# 	print substr($splice5,6,2)."\n"; #TODO GT
				# }

				# if(length($splice3) == 11)
				# {
				# 	print substr($splice3,3,2)."\n"; #TODO AG
				# }
			}
			else
			{
				die "Wrong strand identifier (should be either + or -)! Exiting...\n";
			}
			
			push(@stack5, $splice5);
			push(@stack3, $splice3);
			$prevstr = $strand;
		}
	}
}
close GFFIN;

if($prevstr eq "-")
{
	pop(@stack5);	
}
elsif($prevstr eq "+")
{
	pop(@stack3);	
}

my %out5;
my %out3;

calcFreqs(\@stack5,\%out5,"$opts{f}.pssm",6);
calcFreqs(\@stack3,\%out3,"$opts{t}.pssm",5);

sub calcFreqs
{
	my $array = shift;
	my $outHash = shift;
	my $printname = shift;
	my $splice = shift;
	my %freqs;
	open (OUT, '>', $printname) or die "Cannot open file $printname for writing: $!";
	print OUT "#A\tC\tG\tT\n";

	foreach my $bases (@{$array})
	{
		if(length($bases) == 11)
		{
			my @pos = split(//,$bases);
			foreach my $i (0..$#pos)
			{
				unless($pos[$i] eq 'N')
				{
					$outHash->{$i}{$pos[$i]}++;
					$outHash->{$i}{total}++;
				}
			}
		}
	}

	foreach my $pos (sort keys %{$outHash})
	{
		if($pos == $splice)
		{
			print OUT "splice\n";
		}

		my $total = 0;
		my $largest = "N";

		foreach my $nt ("A","C","G","T")
		{
			my $freq = ($outHash->{$pos}->{$nt}/$outHash->{$pos}->{total})*100;
			$freq = int($freq + 0.5);
			$freqs{$nt} = $freq;
			if($largest eq "N" or $freq > $freqs{$largest})
			{
				$largest = $nt;
			}
			$total = $total+$freq;
		}

		$freqs{$largest}=$freqs{$largest}+(100-$total);
		$total = 0;
		my @final;

		foreach my $nt ("A","C","G","T")
		{
			push (@final, $freqs{$nt});	
			$total = $total + $freqs{$nt};
		}

		print OUT join("\t",@final)."\n";
	}

	close OUT;
}
