#!/usr/local/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use POSIX;

my %scaff_in;
my %opts;
my $incount; #counter for read in sequences 
my $gencount; #counter for read in scaffolds
my $facount; #counter for fasta (protein) sequences
my %outseqs;
my %typecounts;
my %gff_in;
my $date = `date`;
my $look=0;
my %scaff_seq;
my $curr_scaff;
my %coding;
my $header = "";
my $seq;
my $scaffs;
my $ids;
my $firstcall = 1; #when gff sub is first called insert GFF header
my $count = 0; #count for number of gff entries

if(defined($opts{h}) || @ARGV == 0){
	die "Usage: ./BLAT2Exonerate2GFF -i <BLAT psl-file query vs target> -f <cds sequences query> -g <target genome sequences> -r <GFF file target> -o <output folder> -l <name logfile>\n";
}

getopts('hi:o:g:f:l:r:', \%opts);


#Create a folder for this run
mkdir $opts{o},0755;
my $path="./$opts{o}/";
my $log_fh; #log file file handler
my $overlap_genes = $path."overlap_genes.tsv";

if($opts{l}){
    open($log_fh, '>', $path.$opts{l}) or die "Could not open file '$opts{l}: ' $!\n";
    print $log_fh (time - $^T)." Created folder for run: ./$opts{o}/\n";
	print $log_fh (time - $^T)." Created logfile for run: $path"."$opts{l}\n"; 
}
	else{
		$log_fh = *STDERR;
		print $log_fh (time - $^T)." Created folder for run: ./$opts{o}/\n";
		print $log_fh (time - $^T)." No location for logfile specified (-l option). Log will be output to STDERR\n";
	}

print $log_fh (time - $^T)." Run started: $date";

#STEP 1: Get information from the best BLAT hit .tsv file (no .psl header, each scaffold id combination only occurs once)
#Assuming input file is a BLAT file

open (INFILE, "<$opts{i}") or die "Cannot open file $!";
while (<INFILE>){
	chomp;
	my @inline=split(/\t/,$_);
	my %pos_hash;
	$inline[8] =~ s/^.//;

	if(!exists($scaff_in{$inline[13]}{$inline[9]})){
		%pos_hash = (start => $inline[15], stop => $inline[16], strand => $inline[8], genseq => "", length => $inline[14]);
		$scaff_in{$inline[13]}{$inline[9]} = \%pos_hash;
		$incount++;
	}
		else{
			die "The scaffold/id combination ".$inline[13]." and ".$inline[9]." has been seen before!\nPlease run get_best_BLAT_hits.sh first to obtain only the lines with the best BLAT hit";
		}
}

print $log_fh (time - $^T)." Reading in $incount sequences from $opts{i}\n";
close INFILE;

#Read in GFF file of the target species (i.e. the scaffolds against which the template AA sequences have been mapped)

open (GFFIN, "<$opts{r}") or die "Cannot open file $!";
print $log_fh (time - $^T)." Reading in GFF file $opts{r}\n";
	while (<GFFIN>){
		chomp;
		next if($_ =~ m/^#/);
		my @gffline=split(/\t/,$_);
		my %gffhash;


        if($gffline[2] eq "gene"){
			#print "start: $gffline[3]\nstop: $gffline[4]\n";
        	$gffline[8] = (split(/;/,$gffline[8]))[0];
        	$gffline[8] =~ s/^ID=//;
        	%gffhash = (start => $gffline[3], stop => $gffline[4], strand => $gffline[6]);
        	$gff_in{$gffline[0]}{$gffline[8]} = \%gffhash;
        }
	}
close GFFIN;

#STEP 2: Read in the scaffold sequences from the target file based on all scaffolds that had a hit. This file is provided by the user via the g (genome) option

open (GENOME, "<$opts{g}") or die "Cannot open file $!";
	while(<GENOME>){
		if($_ =~ m/^>/){
			my @scaff_cand = split(/[\s>]/,$_);
				if(exists($scaff_in{$scaff_cand[1]})){
					#print $log_fh (time - $^T)." \$scaff_in{\$scaff_cand[1]} exists. Setting \$look to 1. Start reading in sequence in the next round.\n";
					$look = 1; $curr_scaff=$scaff_cand[1]; $gencount++;
				}else{$look=0;}
		}else{
        		if($look){
					chomp $_;
       					if(exists($scaff_seq{$curr_scaff})){
							$scaff_seq{$curr_scaff}=$scaff_seq{$curr_scaff}.$_;
						}
							else{$scaff_seq{$curr_scaff}=$_;}
        		}
			}
	}
print $log_fh (time - $^T)." Reading in $gencount sequences from $opts{g}\n";
close GENOME;

#STEP 3: Extract all regions were a hit was found, one by one, and feed this sequence into Exonerate together with the AA fasta file of the original hit. This file is provided by the user using the -f (fasta) option
#Read in the cds sequences we are interested in first

open (FASTA, "<$opts{f}") or die "Cannot open file $!";
	while(my $line = <FASTA>){
       chomp($line);
       		if ($line =~ m/^>/) {
               if ($header ne "") {
                  	$coding{$header} = $seq;
					$facount++;
               	}

				my @head = split(/\s/,$line);
               	$header = substr($head[0],1);
               	$seq = "";

       		} else {
               		$seq .= $line;
       			}       
	}
$coding{$header} = $seq;
$facount++;
print $log_fh (time - $^T)." Reading in $facount CDS sequences from $opts{f}\n\n";
close FASTA;

#Go through a double loop, the first for the scaffold the second for the id and extract that piece of scaffold on which we would like to search.
#We save this and the respective cds file into two temporary files and then call Exonerate
#create and open file for gff check

open (OVER, ">$overlap_genes") or die "Cannot open file $overlap_genes for writing: $!";
print OVER "TemplID\tScaff\tScaffStart\tScaffEnd\tScaffStrand\tOverlapID\tOverlapStart\tOverlapEnd\tOverlapStrand\n";

#initialise loop variables and start loop

foreach $scaffs (keys(%scaff_in)){
	foreach $ids (keys(%{$scaff_in{$scaffs}})){
		print $log_fh (time - $^T)." Preparing run for $ids / $scaffs\n";
	
		if(defined($scaff_seq{$scaffs})){

			#create an array that holds the sequence string one base per array field
			my @sequence = split('',$scaff_seq{$scaffs});
			#reference this array to the respective scaffold/id pair
			$scaff_in{$scaffs}{$ids}->{genseq}= \@sequence;
			#this prints the last index of the array, not the array size, so to get the size we need to add 1

			#The -r option specificies a gff file (r for region or reference). Loop through all gene models on the current scaffold (%gff_in)  and report any genes that are within the region were a template mapped on that scaffold.
			#Then print out the id of scaffolds and sequence and the respective original gene ID into a file

			print $log_fh (time - $^T)." Checking if hit regions of $ids on $scaffs contains previously predicted gene model\n";

			my $orig_gene;
				
				foreach my $orig_gene (keys(%{$gff_in{$scaffs}})){

					if($gff_in{$scaffs}{$orig_gene}->{strand} eq $scaff_in{$scaffs}{$ids}->{strand}
						&& !(
						($gff_in{$scaffs}{$orig_gene}->{start} <= $scaff_in{$scaffs}{$ids}->{start} && $gff_in{$scaffs}{$orig_gene}->{stop} <= $scaff_in{$scaffs}{$ids}->{start})
						||
						($gff_in{$scaffs}{$orig_gene}->{start} >= $scaff_in{$scaffs}{$ids}->{stop} && $gff_in{$scaffs}{$orig_gene}->{stop} >= $scaff_in{$scaffs}{$ids}->{stop})
						)
					){

						print OVER "$ids\t$scaffs\t$scaff_in{$scaffs}{$ids}->{start}\t$scaff_in{$scaffs}{$ids}->{stop}\t$scaff_in{$scaffs}{$ids}->{strand}\t"
						."$orig_gene\t$gff_in{$scaffs}{$orig_gene}->{start}\t$gff_in{$scaffs}{$orig_gene}->{stop}\t$gff_in{$scaffs}{$orig_gene}->{strand}\n";
		
						print $log_fh (time - $^T)." Overlapping gene ($orig_gene) found for $ids on scaffold $scaffs\n";	

					}

				}


			my $start_pos;
			my $end_pos;
			my $final_seq = "";

				if(($scaff_in{$scaffs}{$ids}->{start}-500) <= 0){
					$start_pos=0;
				}
					else{
						$start_pos = ($scaff_in{$scaffs}{$ids}->{start}-501);
					}

			if(($scaff_in{$scaffs}{$ids}->{stop}+500) >= $scaff_in{$scaffs}{$ids}->{length}){
				$end_pos = ($scaff_in{$scaffs}{$ids}->{length}-1);
			}
				else{
					$end_pos = ($scaff_in{$scaffs}{$ids}->{stop}+499);
				}

			for(my $i = $start_pos; $i <= $end_pos; $i++){
				$final_seq = $final_seq.$scaff_in{$scaffs}{$ids}->{genseq}->[$i];
			}

			my $fasta_file = "$path$ids".".fasta";
			my $region_file = "$path$scaffs"."_"."$ids"."_region.fa";


			if(exists($coding{$ids})){
				open (CDS, ">$fasta_file") or die "Cannot open file $fasta_file for writing: $!";
				print CDS ">$ids\n";
				print CDS "$coding{$ids}\n";
				close CDS;
			}else
			{print $log_fh "No coding sequence found for $ids. Skipping sequence.\n"}

			open (REGION, ">$region_file") or die "Cannot open file $region_file for writing: $!";
			print REGION ">".$scaffs."::".$ids."::".$scaff_in{$scaffs}{$ids}->{start}."::".$scaff_in{$scaffs}{$ids}->{stop}."\n";
			print REGION "$final_seq\n";
			close REGION;

			#FIXME make small test set for 3-4 samples with edge cases
			#FIXME test if works by extracting gffread from both files using the truncated scaffold and the original scaffold
			#FIXME GFFs should be exactly the same using the complete sequence vs truncated sequence
			#FIXME Get rid of single comment lines. keep double comment lines for first entry get rid of them for others as well
			#FIXME Now add everything to the same file, do not create multiple single files.
			my $tmpfile = $path."tmp.cds";
			my $tmppep = $path."tmp.pep";
			my $cmd = "exonerate -E --showvulgar FALSE --showalignment no --singlepass FALSE --model coding2genome --showtargetgff no --bestn 1 --ryo \">%ti\\n%tcs\\n\" --target $region_file --query $fasta_file > $tmpfile 2> /dev/null";
			my $cleancmd = "sed -i -e \'/^Command line:/d\' -e \'/^Hostname:/d\' -e \'/^-- completed exonerate analysis/d\' -e \'/^\$/d\' $tmpfile";


			#my $cmd = "/media/Space2/home/andreas/bin/exonerate -E --showvulgar FALSE --showalignment no --singlepass FALSE --model coding2genome --showtargetgff yes --bestn 1 --target $region_file --query $fasta_file --gff3 > $region_file.gff 2>> /dev/null";
			print $log_fh (time - $^T)." Running Exonerate for $ids / $scaffs\n";
			system($cmd);
			print $log_fh (time - $^T)." Running postprocessing for $ids / $scaffs\n";
			system($cleancmd);

			my $transcheck = "translate_seqs.pl $tmpfile 1";

			system($transcheck);

			my $checkPep = 1;

				if(-e "$tmppep"){ #file exists
					open (PEP, "<$tmppep") or die "Cannot open file $!";
						while(my $line = <PEP>){
							if($line =~ m/\*/ && !eof){
								$checkPep = 0;
								last;
							}
						}
					close PEP;

					#TODO log if a pep sequence gets rejected due to stop codons
					if($checkPep){
						my $finalpep = $path."output.pep";
						my $finalcds = $path."output.cds";
						system("cat $tmppep >> $finalpep");
						system("cat $tmpfile >> $finalcds");
						system("rm $tmppep");
						system("rm $tmpfile");
					}else{
						system("rm $tmppep");
						system("rm $tmpfile");
					}
				}			
			# my $finalGFF = processGFF("$region_file.gff", $count++);
			# open (FGFF, ">>".$path.$opts{o}.".gff") or die "Cannot open file ".$path.$opts{o}.".gff for writing: $!";
			# foreach my $line (@{$finalGFF}){
			# 	my $printline = join("\t", @{$line});
			# 	print FGFF $printline;
			# }
			#FIXME once it works clean up the working directory by removing the current cds plus truncated scaffold file
			#TODO -v verbose switch -d debug switch?
		}
	}
}

#FIXME integrate this part into the pipeline
# blat -t=dna -q=dna -fine -noHead ../trichuris_suis.MA.WBPS5.genomic.fa output.cds mapBackMA.psl

# bestBLAT mapBackMA.psl

# blat2gff.pl best_hits_mapBackMA.psl > mapBackMA.gff3

# sed 's/Parent/ID/1' mapBackMA.gff3 | sed 's/_mid1//1' | sed 's/Target/Parent/1' | awk 'OFS="\t"{print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' | grep -v "^#" | grep -v "mid*" |gffread - -o - | sed 's/transcript/mRNA/1' > reformat.gff3

# gffread reformat.gff3 -g ../trichuris_suis.MA.WBPS5.genomic.fa -w mapBackMA.cds

# translate_seqs.pl mapBackMA.cds 1
#TODO log this
# —> then check if original (Exonerate gene prediction CDS and pep) are the same as the one derived from the new gff3 file
#########################END################



# my $transcmd = "translate_seqs.pl test_out.cds 1";
# print $log_fh (time - $^T)." Translating repredicted CDSs\n";
# system($transcmd);


# sub processGFF {
# my $file = shift;
# my $count = shift;
# open (procGFF, "<$file") or die "Cannot open file $file for reading: $!";
# my @finalGff;
# 	while (<procGFF>){
# 		if($_ =~ m/^##/ && $firstcall){
# 			my @header;
# 			$header[0] = $_;
# 			push (@finalGff,\@header);
# 		}else{
# 			unless($_ =~ m/^Command line:|^Hostname:|^-- completed exonerate analysis|^#|^$/){
# 				my @GffId = split(/::/,$_);
# 				my @GffCols = split(/\t/,$_);
# 					# if($GffCols[8] =~ m/^gene_id 0 ;/){
# 					# 	print $GffCols[8]."\n";
# 					# 	$GffCols[8] =~ s/^gene_id 0 ;/gene_id $count ;/;
# 					# 	print $GffCols[8]."\n";
# 					# 	print "\n";
# 					# }
# 				unless($GffCols[2] eq "similarity"){
# 					$GffCols[3] = $GffCols[3] + $GffId[2];
# 					$GffCols[4] = $GffCols[4] + $GffId[2];
# 					$GffCols[0] = $GffId[0];
# 					push (@finalGff,\@GffCols);
# 				}
# 			}
# 		}
# 	}
# $firstcall = 0;
# return \@finalGff;
# }
