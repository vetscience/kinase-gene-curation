#!/usr/local/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use POSIX;
my %opts;
my @exonStats;
my $oldOnly = 0;
my $newOnly = 0;
my $exonMatch = 0;
#This programs compares two sets of files containing info of exon positions on scaffolds
#and returns the percentage of old (o) exons covered by new (n) exons as well as the percentage of sequence covered by old exons in the new exon
#file format:
# feature_scaffold147_1	668321	272971,273145	272971,273145	272977,273145	272977,273145	272977,273145	273007,273145	273007,273145
# feature_scaffold158_1	268310	93518,93560	93616,93796	93845,93938	93986,94184	94238,94427	94477,94588	94641,95226	95279,95402	95452,95602
my $test;
open($test, ">>checkexons.log") or die "Could not open file checkexons.log: $!\n";


if(defined($opts{h}) || @ARGV == 0){die "Usage: compExons.pl
\n-o <positions of old exons>
\n-n <positions of newly predicted exons>
\n-f <feature name to assess>
\n";}

getopts('o:n:f:', \%opts);

#print "compExon $opts{f}\n";

#add 1s for old exon positions

my $oPtr = readFile( $opts{o} , 1 , $opts{f} );
#print "\$oPtr".$oPtr."\n";
my @oArray = @{$oPtr};

#add 2s for new exon postions
my $nPtr = readFile( $opts{n} , 2 , $opts{f} );
#print "\$nPtr".$nPtr."\n";
my @nArray = @{$nPtr};

if (scalar(@oArray) != scalar(@nArray))
{
	print $test "old: ".scalar(@oArray)."\n";
	print $test "new: ".scalar(@nArray)."\n";
	print $test "old: ".$opts{o}."\n";
	print $test "new: ".$opts{n}."\n";
	print $test "Fatal error: Scaffolds of $opts{o} and $opts{n} are not the same!\n";
	die;
} else
	{
		#print $test "old: ".scalar(@oArray)."\n";
		#print $test "new: ".scalar(@nArray)."\n";
		#print $test "Scaffolds of $opts{o} and $opts{n} are equal!\n";

		for my $pos (0 .. $#oArray) 
		{

			if ($oArray[$pos] > 0 && $nArray[$pos] > 0)
			{
				#print $test "$oArray[$pos]\n";
				#print $test "$nArray[$pos]\n\n";	
			}
			my $newVal = ($oArray[$pos] + $nArray[$pos]);
			#if ($newVal == 3){print $test "newVal\n";};			 
			push(@exonStats, $newVal);
		}
		#print $test "\n\n";
	}

foreach my $stat (@exonStats)
{
	if ($stat == 1){$oldOnly++;};
	if ($stat == 2){$newOnly++;};
	if ($stat == 3)
	{
		$exonMatch++;
		#print $test "stat3 $stat";
	}
}


#print $test "1: $oldOnly\n";
#print $test "2: $newOnly\n";
#print $test "3: $exonMatch\n";


my $coverage = $exonMatch/($oldOnly+$exonMatch)*100;
$coverage = int($coverage + 0.5);
print $test "Coverage comp: $coverage\n";
print $test "feature: ".$opts{f}."\n\n";
print $coverage;

sub readFile
{

	my $infile = shift;
	my $increm = shift;
	my $feature = shift;
	my @posAr = ();

	open (EXIN, '<', $infile) or die "Cannot open file $infile: $!";

	while (<EXIN>)
	{
		chomp;
		my @inline=split(/\t/,$_);
		my $featname = shift(@inline);
		my $scaffname = shift(@inline);
		my $length = shift(@inline);
		$length++;
		@posAr = (0) x $length;

		if($featname eq $feature)
		{
		#print "$featname equals $feature\n";
			foreach my $exon (@inline)
			{

				#print $test "$exon\n";
				my ($start, $end) = split(/,/,$exon);
				#print $test "$increm\n";
				#print $test "$start\n";
				#print $test "$end\n\n";

				for(my $i=$start; $i <= $end; $i++)
				{
					$posAr[$i]=$increm;
					#print "increm: $increm\n";
				}
			}
			#print join("\t",@posAr)."\n";
			last;
		}
		#else {print "$featname does not equal $feature\n";}
	}

	
	close EXIN;
	return \@posAr;
}

close $test;