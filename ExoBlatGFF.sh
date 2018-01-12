#!/bin/bash

#AJS 07.07.2016
#Given a CDS sequence, a target area on a scaffold, repredict the CDS in that area of the genome and create a gff file
#TODO if user gives path for CDS this creates a problem when prefix is added to the file name
#NTH add CAP3 option, check if there is a gene prediction in that region, and if the cds covers this gene prediction, if not repredict cds from cap reassembly (all transcripts in this area +- 5 bp)
#then getorf -find 2
#get longest
#extract and clean up

if [[ $1 == "-h" ]]; then
	echo "Usage: ExoBlatGFF.sh <CDS> <target scaffold name> <start of target region> <end of target region> <genome file> <old gff file (replace gene pred with new one)>"
	exit 1;
fi

shopt -s expand_aliases
source ~/.bashrc

#echo "CDS=$1 TAR=$2 RSTART=$3 REND=$4 GEN=$5 GFF=$6"

CDS=$1
TAR=$2
RSTART=$3
REND=$4
GEN=$5
GFF=$6
FIRST=$7

HEAD=$(head -1 $CDS)

#logfile, will only be created if fatal error
#LOG="$(date | sed 's/\s\+/_/g;s/$/.diff/')".log

#create a truncated version
fa2tab $GEN | awk -v tar=$TAR '{if($1 == ">"tar){print $0}}' | awk -v start=$RSTART -v end=$REND -F"\t" '{print $1; print substr($2,start,(end-start)+1)}' | awk '{print $1}' > $TAR"_trunc.fa"

echo "$(($REND - $RSTART))"
seqstat $TAR"_trunc.fa"

#create a variable for the truncated scaffold
TARtrunc=$TAR"_trunc.fa"

#Run exonerate on this scaffold and the given CDS
exonerate -E --showvulgar FALSE --showalignment no --singlepass FALSE --model coding2genome --showtargetgff no --bestn 1 --ryo ">%ti\n%tcs\n" --target $TARtrunc --query $CDS > "repr_"$CDS

#TODO check that the resulting CDS is the same as the one as the query CDS or at least has a clean ORF
#TODO if not reverse complement the sequence before blatting
#TODO pass on + - from transcript then revcomp before

#Remove stuff from Exonerate output and replace fasta header by original fasta header
grep -v -e "^Command line:" -e "^Hostname" -e "--" "repr_"$CDS | sed '/^$/d' | awk -v head=$HEAD '{if(NR == 1){print head}else{print $0}}' > tmp

mv tmp "repr_"$CDS


#BLAT CDS produced by exonerate back to the original genome
blat -t=dna -q=dna -fine -noHead $GEN "repr_"$CDS "repr_"$CDS".psl"

#get best BLAT hit only
bestBLAT "repr_"$CDS".psl"

#clean up best_hits file
mv "best_hits_repr_"$CDS".psl" "repr_"$CDS".psl"

#Convert .psl file to GFF3
blat2gff.pl "repr_"$CDS".psl" > "mapBack_repr_"$CDS".gff3"

sed 's/Parent/ID/1' "mapBack_repr_"$CDS".gff3" | sed 's/_mid1//1' | sed 's/Target/Parent/1' | awk 'OFS="\t"{print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' | grep -v "^#" | grep -v "mid*" |gffread - -o - | sed 's/transcript/mRNA/1' > reformat.gff3
FIRST=$(($FIRST + 1))

echo "$CDS on $TAR reformat.gff3:"
echo "Lines in reformat.gff3: $(wc -l reformat.gff3)"

echo "Exo \$FIRST: "$FIRST

if [[ $FIRST -le 1 ]]; then
	echo "\$FIRST is less/equal than 1. Creating new GFF file"
cat reformat.gff3 > $GFF
else
	echo "\$FIRST is greater than 1. Appending to existing GFF file"
sed -i '/^#/d' reformat.gff3
cat reformat.gff3 >> $GFF
fi

rm "mapBack_repr_"$CDS".gff3"

#before merging into old GFF file check if there are multiple genes on that scaffold
#if so, print warning and exit
# gencheck=$(grep "$TAR" "$GFF" | awk -F"\t" '{if($3=="mRNA"){print $0}}' | wc -l)

# if [[ $gencheck -gt 1 ]]; then
# 	echo date >> $LOG
# 	echo "Multiple genes on $TAR! Cannot merge GFF files" >> $LOG
# 	echo "Output gff in reformat.gff3" >> $LOG
# 	echo $gencheck >> $LOG
# 	GFF="reformat.gff3"
# else
# 	#merge the old gff with the new one
# 	sed -i '/^#/d' reformat.gff3
# 	grep -v $TAR $GFF > tmp
# 	cat tmp reformat.gff3 > $GFF
# 	rm reformat.gff3
# 	rm tmp
# fi
# unset gencheck

#run some checks, if cds from gff is not equivalent to the one produced by exonerate produce warning and exit
fasta_formatter < "repr_"$CDS > diff1
grep $TAR reformat.gff3 | gffread -g $GEN -w - | fasta_formatter > diff2
diffcheck=$(diff diff1 diff2)


if [[ -n $diffcheck ]]; then
	echo $(date) # >> $LOG
	#echo "FATAL: diff not empty!"; >> $LOG
	echo "\$diffcheck not empty! Differences: $diffcheck" #>> $LOG
	#echo $diffcheck >> $LOG
	unset diffcheck
	#exit 1;
fi
unset diffcheck

rm diff1
rm diff2

echo $(date)" Done with $CDS on $TAR. Output in $GFF"
echo ""