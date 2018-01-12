#!/usr/local/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use POSIX;

my %opts;		#Input parameters and files
my %psl_in; 	#Kinase domains (or other domains of interest) mapped to the genome
my %all_psl; 	#All mapped de novo transcript information
my $incount;	#Counter variable for the number of read kinase transcripts
my $allcount;	#Counter variable for the number of all read transcripts
my %psl_len;	#Hash containing scaffold length
my %outHash; 	#Hash containing exon output info coordinates 
my %outLen;		#Hash containing scaffold length (for output)
my %outScaff;	#Hash containing scaffold name (for output)
my %kinCand;

if(defined($opts{h}) || @ARGV == 0)
{
	die "Usage: psl2cap3.pl
	\n-i <BLAT psl-file with sequences mapped to areas of interest (e.g., all transcripts with a kinase domain)>
	\n-a <BLAT psl-file with all mapped transcripts>
	\n-f <transcript sequences>
	\n-g <target genome scaffolds>
	\n-o <output list of de novo transcripts to assemble>
	\n-l <log file>
	\n";
}

#Example:
#psl2cap3.pl -i best_hits_kinaseAA_ORFsFE.psl -a best_hits_all_dn_FE.psl -o curateFE.ls -l log_curateFE.txt

getopts('i:o:l:a:', \%opts);
open(LOGFH, '>', $opts{l}) or die "Could not open file $opts{l}: $!\n";
open(OUTF, '>', $opts{o}) or die "Could not open file $opts{o}: $!\n";

print LOGFH (time - $^T)." Start program\n";

open (INFILE, '<', $opts{i}) or die "Cannot open file $opts{i}: $!\n";
	while (<INFILE>)
	{
		chomp;
		my @inline=split(/\t/,$_);
		my %pos_hash;

		#check if bestBLAT was run and only one best hit exist for a transcript/scaffold pair in the .psl file
		if(!exists($psl_in{$inline[13]}{$inline[9]}))
		{
			#remove trailing comma and split exon sizes into array
			$inline[18] = substr($inline[18],0,-1);
			my @exSize = split(/,/,$inline[18]);
			#since the kinase mapping was protein to genome using BLAT, the blocksizes have to be converted to nucleotide length
			#foreach my $size (@exSize) { $size = $size * 3; };

			#remove trailing comma and split exon position on scaffolds into array
			$inline[20] = substr($inline[20],0,-1);
			my @exTPos = split(/,/,$inline[20]);

			%pos_hash = (start => $inline[15], stop => $inline[16], strand => $inline[8], exSize => \@exSize, exTPos => \@exTPos);
			#hash reference containing information on start, stop, strand, exon size, and target position of exons on scaffold
			#for each transcript/scaffold pair
			$psl_in{$inline[13]}{$inline[9]} = \%pos_hash;
			#write scaffold length
			$psl_len{$inline[13]} = $inline[14];
			$kinCand{$inline[9]}=$inline[13];
			$incount++;
		}
		else
		{
			die "The scaffold/id combination ".$inline[13]." and ".$inline[9]." has been seen before!\nPlease run bestBLAT first to obtain only the lines with the best BLAT hit";
		}
	}

print LOGFH (time - $^T)." Done reading in $incount sequences from $opts{i}\n";
close INFILE;

open (ALLPSL, '<', $opts{a}) or die "Cannot open file $opts{a}: $!";
	while (<ALLPSL>)
	{
		chomp;
		my @inline=split(/\t/,$_);
		my %pos_hash;
		
		#check if bestBLAT was run and only one best hit exist for a transcript/scaffold pair in the .psl file
		if(!exists($all_psl{$inline[13]}{$inline[9]}))
		{
			#remove trailing comma and split exon position on scaffolds into array
			$inline[20] = substr($inline[20],0,-1);
			my @exons = split(/,/,$inline[20]);
			%pos_hash = (start => $inline[15], stop => $inline[16], exons => \@exons);
			$all_psl{$inline[13]}{$inline[9]} = \%pos_hash;
			$allcount++;
			$psl_len{$inline[13]} = $inline[14];
		}
		else
		{
			die "The scaffold/id combination ".$inline[13]." and ".$inline[9]." has been seen before!\nPlease run bestBLAT first to obtain only the lines with the best BLAT hit";
		}
	}
print LOGFH (time - $^T)." Done reading in $allcount sequences from $opts{a}\n";
close ALLPSL;

##############################Expand matched kinase domain hits based on all matched transcripts########################

#Loop through all scaffold that contain at least one kinase candidate transcript
foreach my $scaffs (keys(%psl_in))
{
	my $feature = 1;
	#initialise scaffold-sized array with zeros 
	my @scaffArr = (0) x $psl_len{$scaffs};
	
	#increment each position that is covered by a kinase candidate transcript by 1
	foreach my $ids (keys(%{$psl_in{$scaffs}}))
	{
		for(my $i=$psl_in{$scaffs}{$ids}->{start}; $i <= ($psl_in{$scaffs}{$ids}->{stop}-1); $i++)
		{
			$scaffArr[$i]++;
		}
	}

	my $inside = 0; 				#bool to check if inside a feature
	my $qName; 						#feature name
	my $qSize = 0;					#size of the feature
	my $tSize = scalar @scaffArr;	#size of the scaffold
	my $tStart; 					#Alignment start position in query.
	my $tEnd; 						#Alignment end position in query.

	#run through entire scaffold array
	foreach my $pos (0 .. $#scaffArr)
	{
		#if position is EMPTY (=0) and we ARE NOT inside a feature, continue
		if(!$scaffArr[$pos] && !$inside)
		{
			next;
		}

		#if position is NOT EMPTY (!=0) and we ARE inside a feature, increase the size of the feature by one and continue
		if($scaffArr[$pos] && $inside)
		{
			$qSize++;
			next;
		}
		
		#marks the START of a feature
		if($scaffArr[$pos] && !$inside)
		{
			$inside = 1;
			$qName = "feature_".$scaffs."_".$feature++;
			$qSize++;
			$tStart = $pos;
		}

		#marks the END of a feature
		if(!$scaffArr[$pos] && $inside)
		{
			$tEnd = $pos-1;
			$inside = 0;

			## foreach my $scaffs (keys(%all_psl))
			## {
			## 	if($scaffs eq $scaffs)
			## 	{

					#feature has ended, now add all overlapping non kinase transcripts on this scaffold
					#loop through all (non-kinase AND kinase) transcripts 
					foreach my $trans (keys(%{$all_psl{$scaffs}}))
					{	
						#check if transcript boundaries are outside of feature boundaries
						unless($all_psl{$scaffs}{$trans}{start} > $tEnd || $all_psl{$scaffs}{$trans}{stop} < $tStart)
						{
							#Check all exons of this transcript, if at least one of them lies within the region of the original gene prediction
							#add this transcript to the list and exit the loop (last).
							my $written = 0;

							foreach my $exon (@{$all_psl{$scaffs}{$trans}{exons}})
							{
								if($exon >= $tStart && $exon <= $tEnd)
								{

									my $exSize = "";
									my $exTPos = "";
									
									#If the transcript that is about to be written exists in %psl_in add exon information, otherwise leave blank
									if(exists($psl_in{$scaffs}{$trans}))
									{
										my $blockCnt = 0;

										foreach my $block (@{$psl_in{$scaffs}{$trans}->{exSize}})
										{

											my $start;
											my $end;

											# if($psl_in{$scaffs}{$trans}->{strand} eq "+-")
											# {
											# 	$start = $tSize - ($psl_in{$scaffs}{$trans}->{exTPos}[$blockCnt]+$block);
											# 	$end = $tSize - $psl_in{$scaffs}{$trans}->{exTPos}[$blockCnt];
											# }
											# else 
											# {
												$start = $psl_in{$scaffs}{$trans}->{exTPos}[$blockCnt];
												$end = $psl_in{$scaffs}{$trans}->{exTPos}[$blockCnt]+$block;	
											# }
											

											my @startEnd = ($start , $end);
											push(@{$outHash{$qName}}, \@startEnd);
											$outLen{$qName} = $tSize;
											$outScaff{$qName} = $scaffs;
											$blockCnt++;
										}
									
										print LOGFH "Transcript $scaffs $qName $trans is a mapped kinases (psl_in), adding to exonfile\n";
									}
									else
									{
										if(exists($kinCand{$trans}))
										{
											print LOGFH "### Transcript $trans in psl_in but not mapped to $scaffs, instead mapped to $kinCand{$trans}\n";
										}
									}
									#Feature name 				transcripts     scaffold 		start  end    length
									# feature_scaffold139_1   >dntrans012496  scaffold139     105821  107633  499514
									# feature_scaffold139_1   >dntrans111631  scaffold139     102839  120587  499514
									# feature_scaffold139_2   >dntrans526278  scaffold139     398568  399316  499514
									# feature_scaffold139_2   >dntrans311860  scaffold139     399478  403022  499514
									print OUTF $qName."\t>".$trans."\t".$scaffs."\t".$all_psl{$scaffs}{$trans}{start}."\t".$all_psl{$scaffs}{$trans}{stop}."\t".$tSize."\n";
									#This loop exits as soon as one exon was found that lies within the mapping region of the original kinase transcript
									$written = 1;
									last;
								}#End if($exon >= $tStart && $exon <= $tEnd)					
							}#End foreach my $exon (@{$all_psl{$scaffs}{$trans}{exons}})	

							if (!$written)
							{
								print LOGFH "$trans $scaffs : No exons within boundaries. $trans $scaffs will not be written\n";
							}
						}#End unless
					}#foreach my $trans (keys(%{$all_psl{$scaffs}}))
			##}#End if($scaffs eq $scaffs)
			##}#foreach my $scaffs (keys(%all_psl))
		}#End of if (if(!$scaffArr[$pos] && $inside)) -> End of feature
	}#End of foreach my $pos (0 .. $#scaffArr)
}#End foreach my $scaffs (keys(%psl_in))

#Write exon position of old kinase candidate transcripts for comparison with new repredicted sequences created in dnCap3.sh script

my @tcsv = split(/\./,$opts{o});
my $exonfile = $tcsv[0].".tcsv";
open(EXLOG, '>', $exonfile) or die "Could not open file $exonfile: $!\n";

my $numout = scalar keys(%outHash);
print STDERR "Number of elements in outhash: $numout\n";
foreach my $features (keys(%outHash))
{

	#Sort the array that is holding arrays of two values (start, stop) according to their start position then according to their stop position
	@{$outHash{$features}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$outHash{$features}};
	print EXLOG $features;
	print EXLOG "\t".$outScaff{$features};
	print EXLOG "\t".$outLen{$features};
	foreach my $pos (@{$outHash{$features}})
	{
		print EXLOG "\t".@{$pos}[0].",".@{$pos}[1];
	}
	print EXLOG "\n";
}

close EXLOG;
close OUTF;
print LOGFH (time - $^T)." End program\n";
close LOGFH;
#END