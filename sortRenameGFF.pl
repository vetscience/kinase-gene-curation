#!/usr/local/bin/perl -w

#Andreas Stroehlein, 2016, The University of Melbourne
#This program takes a properly formatted GFF3 file as input (-g option),
#sorts the GFF entries according to the scaffold name and start position within a scaffold,
#and renames all genes using a given prefix (-p option) and an increasing number
#or a tab-separated lookup table (-l option). Genes that are not found in the lookup_table will be renamed "NOTFOUND_<number>"
#Output (-o defines the output prefix) will be a <prefix>.gff and a tab-separated lookup file (<prefix>.lkp)
#TODO To test if this works, you can simply feed the <prefix>.lkp file back into the script and it should recreate the original file from the first run

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use POSIX;
use Storable 'dclone';
my %opts;
my $linecnt = 0;
my $genecnt = 0;
my %lookup;
my %genes;

print STDERR "WARNING: This script currently only handles lines of the types \"gene\", \"mRNA\", \"exon\" and \"CDS\".\n";

getopts('g:p:s:o:l:h', \%opts);

if(defined($opts{h}) || !defined($opts{g}) || (!defined($opts{p}) && !defined($opts{l}) && !defined($opts{s})) || !defined($opts{o}))
{
	print STDERR "Usage: sortRenameGFF.pl
	-g <GFF3 input file>
	-o <output prefix>

	Either -p or -l or -s need to be defined:
	-p <name prefix> (this will create gene names in the format <prefix>_00X)
	-s <name suffix> (this will create gene names in the format 00X_<suffix>)
	-l <tab-separated lookup file> (in the format \"oldname\tnewname\")

	-h print this help\n";
	exit 1;
}

if(defined($opts{l}))
{
	open (LOOK, '<', $opts{l}) or die "Cannot open file $opts{l}: $!";
	while (<LOOK>)
	{
		chomp;
		my @line = split(/\t/,$_);
		$lookup{$line[0]} = $line[1];
	}
	close LOOK;
}

#Do first quick run, just to get all gene type lines, count the number of genes, sort them according to their position and scaffold and store all of this in a hash
open (GFFIN, '<', $opts{g}) or die "Cannot open file $opts{g}: $!";

while (<GFFIN>)
{
	chomp;
	my @inline=split(/\t/,$_);
	if(!($inline[0] =~ /^#/) && $inline[2] eq "gene")
	{
		$genecnt++;
		my $scaffold = $inline[0];
		my $start = $inline[3];
		my $attr = $inline[8];
		my @attrA = split(/[=:]/,$attr);
		my $id = $attrA[2];
		my @lines;
		my $pt = "";

		if($id =~ m/_pt\d+$/)
		{
			$pt = $id;
		}
		
		my %geneinfo = (scaffold => $scaffold, start => $start, pt => $pt, lines => \@lines, newName => "");

		$genes{$id} = \%geneinfo;
	}
}
close GFFIN;

#Sort hash and add newIDs. The hash needs to be sorted again for the output later
my $ptcnt = 0;
my $current = 0;
foreach my $oldid (sort {$genes{$a}->{pt} cmp $genes{$b}->{pt} || ($a =~ /([a-zA-Z_]+)\d+/)[0] cmp ($b =~ /([a-zA-Z_]+)\d+/)[0] || ($a =~ /[a-zA-Z_]+(\d+)/)[0] <=> ($b =~ /[a-zA-Z_]+(\d+)/)[0] || $genes{$a}->{start} <=> $genes{$b}->{start} } keys %genes)
{	
	if(defined($opts{l}))
	{
		my $newid = $lookup{$oldid};
		$genes{$oldid}->{newName} = $newid;
	}
	else
	{
		my $part = "";
		$current++;
		if($genes{$oldid}->{pt} ne "")
		{
			$ptcnt++;
			$part = "_pt".$ptcnt;
			if($ptcnt > 1)
			{
				$current--;
			}
		}
		else
		{
			$ptcnt = 0;
		}

		my $newid;
		my $padlen = length($genecnt)-length($current);
		my $padding = "";

		for(1..$padlen)
		{
			$padding=$padding."0";
		}

		$newid = $opts{p}."_".$padding.$current."_".$opts{s}.$part;
		$genes{$oldid}->{newName} = $newid;
	}
}

#Now run through all lines and change gene names
open (GFFIN, '<', $opts{g}) or die "Cannot open file $opts{g}: $!";

while (<GFFIN>)
{

	if(!($_ =~ /^#/))
	{
		chomp;
		my @inline=split(/\t/,$_);
		my $id;
		my $scaffold=$inline[0];
		my $source=$inline[1];
		my $type=$inline[2];
		my $start=$inline[3];
		my $end=$inline[4];
		my $score=".";
		my $strand=$inline[6];
		my $phase=$inline[7];
		my $attr = $inline[8];
		
		if($type eq "gene")
		{
			my @attrA = split(/[=:]/,$attr);
			$id = $attrA[2];
			$attr = "ID=gene:".$genes{$id}->{newName};
		}
		elsif($type eq "mRNA")
		{
			my @attrA = split(/[=:;]/,$attr);
			$id = $attrA[2];

			$attr = "ID=transcript:".$genes{$id}->{newName}.";Parent=gene:".$genes{$id}->{newName};
	 	}
	 	elsif($type eq "exon")
	 	{
	 		my @attrA = split(/[=:;\.]/,$attr);
			$id = $attrA[2];

			$attr = "ID=exon:".$genes{$id}->{newName}.".".$attrA[3].";Parent=transcript:".$genes{$id}->{newName};
	 	}
	 	elsif($type eq "CDS")
	 	{
	 		my @attrA = split(/[=:;]/,$attr);
			$id = $attrA[2];

			$attr = "ID=cds:".$genes{$id}->{newName}.";Parent=transcript:".$genes{$id}->{newName};
	 	}
	 	else
	 	{
	 		print STDERR "Skipping line $linecnt:\n$_\n";
	 		next;
	 	}

	 	my $cleanline = join("\t", $scaffold, $source, $type, $start, $end, $score, $strand, $phase, $attr);
		push(@{$genes{$id}->{lines}}, $cleanline);
	}
}
close GFFIN;

my $outlkp = $opts{o}.".lkp";
my $outgff = $opts{o}.".gff";

open (OUTLKP, '>>', $outlkp) or die "Cannot open file $outlkp: $!";
open (OUTGFF, '>>', $outgff) or die "Cannot open file $outgff: $!";
print OUTGFF "##gff-version 3\n";

foreach my $oldid (sort {$genes{$a}->{pt} cmp $genes{$b}->{pt} || ($a =~ /([a-zA-Z_]+)\d+/)[0] cmp ($b =~ /([a-zA-Z_]+)\d+/)[0] || ($a =~ /[a-zA-Z_]+(\d+)/)[0] <=> ($b =~ /[a-zA-Z_]+(\d+)/)[0] || $genes{$a}->{start} <=> $genes{$b}->{start} } keys %genes)
{
	print OUTLKP $oldid."\t".$genes{$oldid}->{newName}."\n";

	print OUTGFF "# Gene gene:".$genes{$oldid}->{newName}."\n";
	foreach(@{$genes{$oldid}->{lines}})
	{
		print OUTGFF $_."\n";
	}
	print OUTGFF "###\n";
}
exit 0;