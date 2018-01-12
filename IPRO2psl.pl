#!/usr/local/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use POSIX;
use Storable 'dclone';
my %opts;

getopts('g:i:r:o:', \%opts);
my $DOMAININFO = "/media/Space2/home/andreas/data/domain_lookup_05_16.txt";
print "Domain information will be retrieved from ".$DOMAININFO.".
If the description for domains (tab-separated text file containing domain ID in column 1
and free text description in column 2) is at another location,
please modify the constant \$DOMAININFO in line 11 of this script.\n";

#This program takes a GFF3 file (-g option) and an InterProScan output file (.tsv, -i option)
#and creates a .psl file that can be read into a genome browser such as IGV, which displays all domains and their positions within a gene
#AJ Stroehlein, November 2016, The University of Melbourne

#Read in fasta file containing scaffolds
my %fasta;
open (FASTA, '<', $opts{r}) or die "Cannot open file $opts{g}: $!";
my $current;
while (<FASTA>)
{
	chomp;
	if($_ =~ /^>/)
	{
		my @header = split(/[>\s]/,$_);
		$fasta{$header[1]} = 0;
		$current = $header[1];
	}
	else
	{
		$_ =~ s/\s//g;
		$fasta{$current} = $fasta{$current} + length($_);
	}
}
close FASTA;

#Read in GFF file
my %cds;
open (GFFIN, '<', $opts{g}) or die "Cannot open file $opts{g}: $!";

while (<GFFIN>)
{

	if(!($_ =~ /^#/))
	{
		chomp;
		my @inline=split(/\t/,$_);

		if($inline[2] eq "CDS")
		{
			my @attr = split(/[:;]/,$inline[8]);
			my $id = $attr[1];
			my $scaffold = $inline[0];
			my @startEndPhase = ($inline[3] , $inline[4] , $inline[7]);
			my $strand = $inline[6];

			if(exists($cds{$id}))
			{
				if($cds{$id}->{strand} ne $strand)
				{
					die "FATAL ERROR: Strand of current CDS\n.".$_."\nis not equal to strand of previous CDSs of this feature (should be ".$cds{$id}->{strand}.")\n";
				}

				push(@{$cds{$id}->{exons}}, \@startEndPhase);
			}
			else
			{
				my @posArray = (\@startEndPhase);
				my %infoHash = (scaffold => $scaffold, strand => $strand, exons => \@posArray);
				$cds{$id} = \%infoHash;
			}
		}
	}
}
close GFFIN;

#Read in lookup table
my %look;
open (DOMS, '<', $DOMAININFO) or die "Cannot open file $DOMAININFO: $!";

while (<DOMS>)
{
	chomp;
	my @inline = split(/\t/,$_);
	$look{$inline[0]} = $inline[1];	
}
close DOMS;

#Read in InterProScan file
my %ipro;
open (IPRO, '<', $opts{i}) or die "Cannot open file $opts{i}: $!";

while (<IPRO>)
{
	chomp;
	my @inline = split(/\t/,$_);
	my $id = $inline[0];
	my $domID;
	my $domDescr;
	my $start;
	my $length;

	if($inline[3]=~/^PANTHER$|^Pfam$|^SUPERFAMILY$|^ProSitePatterns$/)
	{

		my @splitSubFams = split(/:/,$inline[4]);
		$inline[4] = $splitSubFams[0];

		if($inline[5] ne "")
		{
			$inline[5] =~ s/\s/_/g;
			$domDescr = $inline[5];
		}
		else
		{
			if(exists($look{$inline[4]}))
			{
				$domDescr = $look{$inline[4]};
				$domDescr =~ s/\s/_/g;
			}
			else
			{
				$domDescr = "ID_NA";
			}
		}

		$domID = $inline[4];
		$start = ((($inline[6]-1)*3)+1);
		$length = (($inline[7]-$inline[6]+1)*3);

		if(exists($ipro{$id}))
		{
			my %infoHash = (description => $domID."|".$domDescr, start => $start, length => $length);
			push(@{$ipro{$id}->{domStack}}, \%infoHash);
		}
		else
		{
			$ipro{$id}{domStack} = [];	
			my %infoHash = (description => $domID."|".$domDescr, start => $start, length => $length);
			push(@{$ipro{$id}->{domStack}}, \%infoHash);
		}
	}	
}
close IPRO;


open (OUT, '>', $opts{o}) or die "Cannot open file $opts{o} for writing: $!";
#Go through all features that have an IPRO domain associated with them
foreach my $id (keys %ipro)
{
	my $startpos;
	my $length;
	my @scaffArray = ();
	my $check = 0;

	if(defined $cds{$id}->{scaffold})
	{
		if(defined $fasta{$cds{$id}->{scaffold}})
		{
			$length = $fasta{$cds{$id}->{scaffold}};
			@scaffArray = (0) x $length;
			$check = 1;
		}
		else
		{
			print "No fasta entry for $cds{$id}->{scaffold}\n";
		}
	}
	else
	{
		print "WARNING: No scaffold entry for IPRO $id. This sequence is not represented in the given GFF file\n";
	}
	
	if($check)
	{
		#sort all CDSs of this feature according to their position
		if($cds{$id}->{strand} eq "-")
		{
			@{$cds{$id}->{exons}} = reverse(sort { $a->[1] <=> $b->[1] || $a->[0] <=> $b->[0] } @{$cds{$id}->{exons}});
		}
		else
		{
			@{$cds{$id}->{exons}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$cds{$id}->{exons}};
		}

		#populate the array areas that contain a cds region
		foreach my $pos (@{$cds{$id}->{exons}})
		{

			for my $walkpos ($pos->[0]-1..$pos->[1]-1)
			{
				$scaffArray[$walkpos] = 1;
			}
		}

		#define offset
		my $phase = $cds{$id}->{exons}->[0]->[2];
		print "Phase: $phase\n";

		my $domCount = 0;

		foreach my $domain (@{$ipro{$id}->{domStack}})
		{

			$domCount++;
			print "Domain number: $domCount\n";
			print "$domain->{description}
			Length: $domain->{length}
			Start: $domain->{start}\n";
			#Initialise the line with all the stuff that is already know without mapping the domain onto the exons

			#Negative
			#1062    0       0       0       0       0       1       466     -       transcript:feature_scaffold116_2
			#1062    0       1062    scaffold128     404211  34099   35627   2       1021,41,        0,1021, 34099,35586,

			#Positive
			#1059    0       0       0       0       0       3       179     +       transcript:feature_scaffold142_1
			#1059    0       1059    scaffold175     119545  56702   57940   4       180,281,289,309,        0,180,461,750,  56702,56934,57273,57631,

			my $sumIntronPos = 0;
			my $tStart;
			my $tStop;
			my $numBlocks = 0;
			my @blockSize = ();
			my @qStarts = ();
			my @tStarts = ();

			# this is how many 1s in the scaffold array have to be skipped at the start before you starting "putting down" the domain
			my $offset = $domain->{start} + $phase;
			print "Offset: $offset\n";
			my $startIDX;
			my $endIDX;

						if($cds{$id}->{strand} eq "+")
						{
							$startIDX = 0;
							$endIDX = $#scaffArray;

							my $currPos = 0;
							my $inside = 0;
							my $currBlock = 0;

							# if($cds{$id}->{strand} eq "-")
							# {
							# 	print $cds{$id}->{strand}."\n";	
							# }
							print "For this domain $domCount we will run through the array from $startIDX to $endIDX\n";

							foreach my $scaffpos (0..$#scaffArray)
							{	
								if($currPos == $domain->{length})
								{
									print "Reached end of domain for $id domain $domCount\n";
									$tStop = $scaffpos-1;
									$currPos = 0;
									last;
								}
								else
								{
									if($scaffArray[$scaffpos] == 1)
									{	
										#if already inside an cds just increase $currPos
										if(!$inside)
										#cds start
										{
											if($offset == 0){
											print "Offset is used up: End up in block\n";
											$inside = 1;
											$numBlocks++;
											#$currBlock++;

											push(@qStarts, $currPos);
											push(@tStarts, $scaffpos);

											if(!defined($tStart)){print "\$tStart will be defined now for $id domain $domCount\n"; $tStart = $scaffpos};
											}
											else
											{
												$offset--;
												#print "New Offset: $offset\n";
											}	
										}

										#FIXME this seems to be off
										if($offset == 0){$currPos++; $currBlock++;};
									}
									else
									{
										#end of cds
										if($inside)
										{	
											$inside = 0;
											push (@blockSize, $currBlock);
											$currBlock = 0;	
										}
										#still in intron
										else
										{
											if($numBlocks > 0){$sumIntronPos++;}
										}
									}	
								}
							}

							#final push if last position of domain is inside a cds
							#this should always be the case
							if($inside){push (@blockSize, $currBlock)};	


							my @line = ($domain->{length},0,0,0,0,0,0,$sumIntronPos,$cds{$id}->{strand},$domain->{description},$domain->{length},0,$domain->{length},$cds{$id}->{scaffold},$length,$tStart,$tStop,$numBlocks,join(",",@blockSize).",",join(",",@qStarts).",",join(",",@tStarts).",");

							print OUT join("\t", @line)."\n";
						}
			else
			{
				$startIDX = $#scaffArray;
				$endIDX = 0;

				my $currPos = 0;
				my $inside = 0;
				my $currBlock = 0;

			# if($cds{$id}->{strand} eq "-")
			# {
			# 	print $cds{$id}->{strand}."\n";	
			# }
			print "For this domain $domCount we will run through the array from $startIDX to $endIDX\n";

			foreach (my $scaffpos = $#scaffArray; $scaffpos >= 0; $scaffpos--)
			{	
				# if($cds{$id}->{strand} eq "-")
				# {
				# 	print "I can get into this loop\n";	
				# }
				

				#if($scaffpos >= $length){die "ERROR: Loop runs into the wrong direction for $cds{$id}->{strand}\n"}
				#print "\$scaffpos: $scaffpos\n";
				if($currPos == $domain->{length})
				{
					print "Reached end of domain for $id domain $domCount\n";
					$tStop = $scaffpos-1;
					$currPos = 0;
					last;
				}
				else
				{
					if($scaffArray[$scaffpos] == 1)
					{	

						#if already inside an cds just increase $currPos
						if(!$inside)
						#cds start
						{
							if($offset == 0){
							print "Offset is used up: End up in block\n";
							$inside = 1;
							$numBlocks++;

							push(@qStarts, $currPos);
							push(@tStarts, $scaffpos);

							if(!defined($tStart)){print "\$tStart will be defined now for $id domain $domCount\n"; $tStart = $scaffpos};
							}
							else
							{
								$offset--;
								#print "New Offset: $offset\n";
							}	
						}

						if($offset == 0){$currPos++; $currBlock++;};
					}
					else
					{
						#end of cds
						if($inside)
						{	
							$inside = 0;
							push (@blockSize, $currBlock);
							$currBlock = 0;	
						}
						#still in intron
						else
						{
							if($numBlocks > 0){$sumIntronPos++;}
						}
					}	
				}
			}

			#final push if last position of domain is inside a cds
			#this should always be the case
			if($inside){push (@blockSize, $currBlock)};	


			my @line = ($domain->{length},0,0,0,0,0,0,$sumIntronPos,$cds{$id}->{strand},$domain->{description},$domain->{length},0,$domain->{length},$cds{$id}->{scaffold},$length,$tStart,$tStop,$numBlocks,join(",",@blockSize).",",join(",",@qStarts).",",join(",",@tStarts).",");

			print OUT join("\t", @line)."\n";
			}		
		}
	}
}

close OUT;

#print Dumper %ipro;
#print Dumper %cds;

exit 0;