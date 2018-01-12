#!/bin/bash
echo "Start dnCap3.sh: $(date)"

shopt -s expand_aliases
source ~/.bashrc

#$1 is a file that contains the feature names in the first column and the FASTA headers
#for respective de novo transcripts in the second column
#Example:
#Feature name 				transcripts     scaffold 		start  end    length
# feature_scaffold139_1   >dntrans012496  scaffold139     105821  107633  499514
# feature_scaffold139_1   >dntrans111631  scaffold139     102839  120587  499514
# feature_scaffold139_2   >dntrans526278  scaffold139     398568  399316  499514
# feature_scaffold139_2   >dntrans311860  scaffold139     399478  403022  499514

#$2 is the fasta file containing these de novo transcripts


GEN=$3 #Genome scaffold fasta file

GFF=$4 #old gff file (replace repredicted genes in this file?)

if [[ "$GEN" = "" || "$GFF" = "" ]]; then
echo "No path for genome and/or GFF file specified (Arguments 3 and 4). Exiting..."
exit 1
fi

#feature_list.ls example
# feature_C404737_1
# feature_C411119_1
# feature_scaffold104_1
# feature_scaffold104_2

cut -f1 $1 | sort -u > feature_list.ls

COUNT=0

cat feature_list.ls | while read FEAT; do

ORFCNT=0
GOODORF=0

awk -v feat=$FEAT -F"\t" '{if($1 == feat){print $2}}' $1 >> seq_list.ls
#seq_list.ls example
# >dntrans467434
# >dntrans315725
# >dntrans044075

SCAFF=$(awk -v feat=$FEAT -F"\t" '{if($1 == feat){scaff=$3}}END{print scaff}' $1)

SCAFFLEN=$(awk -v feat=$FEAT -F"\t" '{if($1 == feat){scafflen=$6}}END{print scafflen}' $1)

wc -l seq_list.ls

fasta_formatter < $2 | grep -x -A1 -f seq_list.ls | sed '/--/d' > $FEAT".fa"

cap3 $FEAT".fa" > /dev/null

SINGLET=$(cat $FEAT".fa.cap.contigs" | wc -l)

#TODO add entries from singlet file
#TODO sort unique ORFs (save time running them all through the rest of the loop)
#TODO get mRNA assembled length vs ORF ratio (this value is applied for the check later)
#TODO add crosscheck with 
#TODO check if feature overlaps with existing features (maybe improve selection of transcripts for reassembly based on exon coverage/overlap, add only if it matches with first or last scaffold?)
#TODO define "consensus" among all overlapping kinase transcripts, the kinase catalytic domain should be in this area
#TODO make sure that in the initial mapping, the piece that has the best hit is also the piece containing the kinase domain, otherwise go down the list
#TODO If a mapped kinase candidate spans another mapped kinase candidate but without sharing any exons they have to be considered separately
#TODO Base Coverage value on the exons that are common to all mapped best hits kinase candidates in that area, so the one that gets selected is most likely to also contain the kinase domain

if [[ $SINGLET -gt 0 ]]; then
	getorf -find 2 -minsize 570 $FEAT".fa.cap.contigs" $FEAT".orfs" > /dev/null 2>&1
else
	getorf -find 2 -minsize 570 $FEAT".fa" $FEAT".orfs" > /dev/null 2>&1
fi

#TODO keep this file make sure no clashes regarding overwriting
rm $FEAT".fa.cap"*


NUMORF=$(grep ">" $FEAT".orfs" | wc -l)
echo "NUMORF: $NUMORF"
seqstat $FEAT".orfs"
ORFCNT=0

while [  $GOODORF -eq 0 -a $ORFCNT -lt $NUMORF ]; do

	ORFCNT=$(($ORFCNT + 1))
	echo "ORFCNT: $ORFCNT"

	#changed awk -v header=$FEAT"_"$ORFCNT 27.10. --> changed back to header=$FEAT
	fa2len $FEAT".orfs" | sort -t$'\t' -k2,2nr | head -$ORFCNT | tail -1 | cut -f1 | sed 's/\[/./;s/\]/./;s/\s/./g' | while read line; do fasta_formatter < $FEAT".orfs" | grep -A1 $line | awk -v header=$FEAT '{if($0 ~ "^>"){print ">"header}else{print $0}}' > "best_orf_"$FEAT"__"$ORFCNT".fa"; done

	if [ ! -f $SCAFF.fa ]; then
		echo $SCAFF".fa does not exist. Creating now..."
    	#10.08.17 changed to no space after scaffold id
	fasta_formatter < $3 | grep -A1 "^>"$SCAFF | awk '{print $1}' > $SCAFF.fa
	#fasta_formatter < $3 | grep -A1 "^>"$SCAFF" " | awk '{print $1}' > $SCAFF.fa
	fi
	
	echo "Starting BLAT"
	blat -t=dna -q=dna -fine -noHead $SCAFF.fa best_orf_"$FEAT"__"$ORFCNT".fa best_orf_"$FEAT"__"$ORFCNT".psl > /dev/null
	bestBLAT best_orf_"$FEAT"__"$ORFCNT".psl
	mv best_hits_best_orf_"$FEAT"__"$ORFCNT".psl best_orf_"$FEAT"__"$ORFCNT".psl
	echo "BLAT done"

	ls best_orf_"$FEAT"__"$ORFCNT".fa
	wc -l best_orf_"$FEAT"__"$ORFCNT".fa
	ls best_orf_"$FEAT"__"$ORFCNT".psl
	wc -l best_orf_"$FEAT"__"$ORFCNT".psl

	#get the length coordinates to call from the mapped current bes ORF, not from the curateX.ls file
	#(selected truncated region is way too long when considering all sequences that have at least one exon within the mapping region of kinase candidates)

	

	CDSstart=6000000000000
	echo "CDSstart before "$CDSstart
	CDSstart=$(awk -v feat=$FEAT -v sta=$CDSstart -F"\t" '{if($10 == feat){if($16 < sta){sta=$16}}}END{print sta}' best_orf_"$FEAT"__"$ORFCNT".psl)
	echo "CDSstart after "$CDSstart


	CDSstop=0
	echo "CDSstop before "$CDSstop
	CDSstop=$(awk -v feat=$FEAT -v sto=$CDSstop -F"\t" '{if($10 == feat){if($17 > sto){sto=$17}}}END{print sto}' best_orf_"$FEAT"__"$ORFCNT".psl)
	echo "CDSstop after "$CDSstop

	awk -F"\t" '{printf $10"\t"$14"\t"$15; sub(/,$/,"",$19); split($19,blocks,/,/); sub(/,$/,"",$21); split($21,blockPos,/,/); for(i=1; i<=length(blockPos); i++){printf "\t"blockPos[i]","(blockPos[i]+blocks[i])}; print ""}' best_orf_"$FEAT"__"$ORFCNT".psl > checkBLAT.tcsv
	#rm checkBLAT.psl

	#FIXME Check for strand specific mapping [20]
	#						genome
	# ----------------------------------------------------------
	#						mapped exons (kinase candidate)
	#	>>>>>>>>>>>			>>>>>>>				>>>>>>>>
	#
	#
	#					mapped new longest ORF (reverse strand, covering ~90 of the original exons.)
	#	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#TODO Translate sequence and run InterProScan (KinProScan) to check if kinase

	EXCHECK=$(ls $1 | awk '{split($0,prefix,/\./); print prefix[1]".tcsv";}')
	COVER=$(compExons.pl -o $EXCHECK -n checkBLAT.tcsv -f $FEAT 2> /dev/null)
	#echo "dnCAP $FEAT"
	#compExons.pl -o exonfile.tcsv -n checkBLAT.tcsv -f $FEAT
	echo "Cover: $COVER"

#changed 10.08.2017 don't check

#	if [[ $COVER -gt 50 ]]; then
#		GOODORF=1
#		echo "Good ORF (number $ORFCNT) found ($COVER coverage), done"
#		echo ""
#		#changed cp file to awk and changing header 27.10.
	awk -F"__" '{print $1}' best_orf_"$FEAT"__"$ORFCNT".fa > best_orf_"$FEAT".fa
#	else
#		if [[ $ORFCNT -eq $NUMORF ]]; then
#			echo "WARNING: None of the $ORFCNT predicted ORFs for $FEAT have been selected due to lack of coverage of original kinase prediction. Suggest manual curation for $FEAT!"
#		else
#			echo "No good ORF found ($COVER for $FEAT ORF# $ORFCNT), try next longest..."
#		fi
#	fi
done

GOODORF=1

if [[ $GOODORF -eq 1 ]]; then
	#Set boundaries for remapping with BLAT
	CDSstart=$(($CDSstart - 50))
	if [[ $CDSstart -lt 0 ]]; then
	CDSstart=0
	fi

	CDSstop=$(($CDSstop + 50))
	if [[ $CDSstop -gt $SCAFFLEN ]]; then
	CDSstop=$SCAFFLEN
	fi

	echo "Processing: $FEAT"
	echo "Scaffold: $SCAFF"
	echo "Start: $CDSstart"
	echo "Stop: $CDSstop"
	echo "Length: $(($CDSstop - $CDSstart + 1))"

	#Call script to map new transcript, run Exonerate to predict gene boundaries and convert to gff file
	#ExoBlatGFF.sh <CDS> <target scaffold name> <start of target region> <end of target region> <genome file> <old gff file (replace gene pred with new one)>
	echo "Starting ExoBlatGFF.sh:"
	ExoBlatGFF.sh "best_orf_"$FEAT".fa" $SCAFF $CDSstart $CDSstop $GEN $GFF $COUNT
	echo "Done ExoBlatGFF.sh:"
	COUNT=$(($COUNT + 1))
	echo "dnCap \$COUNT: $COUNT"
fi

rm seq_list.ls

done

echo "End dnCap3.sh: $(date)"
#rm feature_list.ls
