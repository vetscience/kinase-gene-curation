#!/usr/local/bin/perl -w

#Andreas Stroehlein, 2016, The University of Melbourne
use strict;
use warnings;
use Getopt::Std;
my %opts;
my %scaffolds;
my $current;
my $out;

getopts('g:p:o:h', \%opts);

if(defined($opts{h}) || !defined($opts{g}) || !defined($opts{p}))
{
	print STDERR "Usage: getCdsFromProtPsl
	-g <genome scaffold file>
	-p <prot/dnax .psl file>
	-o <output fasta file> (if not defined, print to STDOUT)
	-h print this help\n";
	exit 1;
}

if(defined($opts{o}))
{
	open ($out, '>>', $opts{o}) or die "Cannot open file $opts{o}: $!";
}
else
{
	$out = *STDOUT;
}

open (GENIN, '<', $opts{g}) or die "Cannot open file $opts{g}: $!";
while (<GENIN>)
{
	chomp;
	if($_ =~ m/^>/)
	{
		my @header=split(/[>\s]/,$_);
		$scaffolds{$header[1]}{fw}="";
		$current = $header[1];
	}
	else
	{
		$scaffolds{$current}{fw} = $scaffolds{$current}{fw}.$_;
	}
}
close GENIN;

foreach my $keys (%scaffolds)
{
	$scaffolds{$keys}{rv} = reverse($scaffolds{$keys}{fw});
	$scaffolds{$keys}{rv} =~ tr/ACGTacgt/TGCAtgca/;
}

open (PSLIN, '<', $opts{p}) or die "Cannot open file $opts{p}: $!";
while (<PSLIN>)
{
	chomp;
	my @psl=split(/\t/,$_);
	my $header = $psl[9];
	my $scaffold = $psl[13];
	my $str = $psl[8];

	if($str =~ m/-/)
	{
		$str = "rv";
	}
	else
	{
		$str = "fw";
	}

	if(defined($scaffolds{$scaffold}{$str}))
	{
		my $seq = $scaffolds{$scaffold}{$str};
		
		#get block sizes
		$psl[18] = substr($psl[18],0,length($psl[18])-1);
		#get block starts
		$psl[20] = substr($psl[20],0,length($psl[20])-1);
		
		my @blockSize = split(/,/,$psl[18]);
		my @blockStarts = split(/,/,$psl[20]);
		my $cds = "";

		foreach my $block (0..$#blockSize)
		{
			my $blockSeq = substr($seq,$blockStarts[$block],($blockSize[$block]*3));
			$cds = $cds.$blockSeq;
		}

		print $out ">".$header."\n".$cds."\n";
	}
	else
	{
		print STDERR "WARNING: Target sequence $scaffold not found in $opts{g}! Skipping record...\n";
	}
}
close PSLIN;

if(defined($opts{o}))
{
	close $out;
}
exit 0;