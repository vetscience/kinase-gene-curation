#!/bin/bash


shopt -s expand_aliases
source ~/.bashrc

COUNT=0

GFF=$2
THREE=$3
FIVE=$4
LOOK=$5

grep -v '^#' $1 |

while read file
do
	echo "Start $file"
	
	GOOD=0
	HEAD=$(head -1 $file)
	FEAT=$(head -1 $file | sed 's/^>//')
	SCAFF=$(awk -v feat=$FEAT -F"\t" '{if($1 == feat){print $2}}' $LOOK)
	SCAFFLEN=$(fa2tab "$SCAFF".fa | awk -F"\t" '{print length($2)}')

	echo "Starting BLAT"
	blat -t=dna -q=dna -fine -noHead "$SCAFF".fa $file "$file".psl > /dev/null 2>&1
	bestBLAT "$file".psl
	mv best_hits_"$file".psl "$file".psl
	echo "BLAT done"

	#ls "$file".psl
	wc -l "$file".psl

	#get the length coordinates to call from the mapped current bes ORF, not from the curateX.ls file
	#(selected truncated region is way too long when considering all sequences that have at least one exon within the mapping region of kinase candidates)

	CDSstart=6000000000000
	#echo "CDSstart before "$CDSstart
	CDSstart=$(awk -v feat=$FEAT -v sta=$CDSstart -F"\t" '{if($10 == feat){if($16 < sta){sta=$16}}}END{print sta}' "$file".psl)
	#echo "CDSstart after "$CDSstart


	CDSstop=0
	#echo "CDSstop before "$CDSstop
	CDSstop=$(awk -v feat=$FEAT -v sto=$CDSstop -F"\t" '{if($10 == feat){if($17 > sto){sto=$17}}}END{print sto}' "$file".psl)
	#echo "CDSstop after "$CDSstop

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

	echo "Starting Exonerate"

	fa2tab "$SCAFF".fa | awk -v start=$CDSstart -v end=$CDSstop -F"\t" '{print $1; print substr($2,start,(end-start)+1)}' | awk '{print $1}' > "$SCAFF"_trunc.fa 
	#echo "$(($CDSstop - $CDSstart))"
	#seqstat "$SCAFF"_trunc.fa

	LENB=$(fasta_formatter < $file | grep -v ">" | awk '{print length($0)}')
#removed three and 5 splice parameters - AJS 11.08.2017
#	exonerate -E --showvulgar FALSE --showalignment no --singlepass FALSE --model coding2genome --showtargetgff yes --splice3 $THREE --splice5 $FIVE --bestn 1 --target "$SCAFF"_trunc.fa --query $file > repr_"$FEAT".gff 2> /dev/null
exonerate -E --showvulgar FALSE --showalignment no --singlepass FALSE --model coding2genome --showtargetgff yes --bestn 1 --target "$SCAFF"_trunc.fa --query $file > repr_"$FEAT".gff 2> /dev/null

	echo "Done Exonerate $file on $SCAFF"

	#check if gff file exists and is not empty
	if [[ -s repr_"$FEAT".gff ]]; then
		head -1 repr_"$FEAT".gff
		tail -4 repr_"$FEAT".gff | head -1
		#add coordinate offset depending on the truncated version of the scaffold vs. complete version
		#22 Nov: removed -1's in the awk statement, since GFF3 files are 1 based
		grep -v -e "^Command line:" -e "^Hostname:" -e "--" -e "^#" repr_"$FEAT".gff | sed '/^$/d' | grep "^$SCAFF" | awk -F"\t" -v feat=$FEAT -v start=$CDSstart 'BEGIN{print "##gff-version 3"}{$2="exonerate"; $4=($4+start); $5=($5+start); if($3 == "gene"){$3="mRNA"; $6="."; $9="ID="feat; printf $1; for(i=2; i<=NF; i++){printf "\t"$i}; print "";}; if($3 == "exon"){$9="Parent="feat; printf $1; for(i=2; i<=NF; i++){printf "\t"$i}; print "";}}' > tmp

		mv tmp repr_"$FEAT".gff

		#grep -v -e "^Command line:" -e "^Hostname" -e "--" "repr_"$file | sed '/^$/d' | awk -v head=$HEAD '{if(NR == 1){print head}else{print $0}}' > tmp
		#mv tmp "repr_"$file
		let COUNT=COUNT+1
		echo "Exonerate run # \$COUNT: "$COUNT

		PHASE=0
		FLIP="n"
		ORIG=0
		while [[ $GOOD -eq 0 && $PHASE -le 2 ]]
		do
			#check this for all frames cut off 0,1,2 nt from start to see if ORF can be found
			STOPS=$(cleanExonerateGFF.pl repr_"$FEAT".gff $PHASE $FLIP n | gffread - -g "$SCAFF".fa -y - | fasta_formatter | grep -v ">" | sed 's/\.$//' | fold -w1 | sort | uniq -c | grep "\." | awk '{print $1}')

			# translate_seqs.pl "repr_"$file 1 > /dev/null 2>&1
			# PEP=$(echo "repr_"$file | awk -F"." '{print $1}')
			# STOPS=$(fasta_formatter < $PEP".pep" | sed 's/\*$//' | grep -v ">" | fold -w1 | sort | uniq -c | grep "\*" | awk '{print $1}')

			if [[ -n $STOPS ]]; then
				echo $STOPS" stop codons found in translation of repr_"$FEAT".gff for phase "$PHASE"!"
				#echo "Attempting to simply map original ORF .fa file using BLAT"
				#echo "and create GFF entry from this mapping"
				if [[ $PHASE -lt 2 ]]; then
					let PHASE=PHASE+1
					echo "Trying next phase: "$PHASE". FLIP: "$FLIP
					#cp $file "repr_"$file
				elif [[ $FLIP == "n" ]]; then
					PHASE=0
					FLIP="f"
				elif [[ $ORIG -eq 0 ]]; then
					echo ""
					echo "EXONERATE could not produce a good gene model for "$file"!"
					echo "Trying all phases for original ORF ("$file") mapped using BLAT for "$FEAT
					ORIG=1
					FLIP="n"
					PHASE=0
					blat2gff.pl "$file".psl > "blat_"$FEAT".gff3"
					sed 's/Parent/ID/1' "blat_"$FEAT".gff3" | sed 's/_mid1//1' | sed 's/Target/Parent/1' | awk 'OFS="\t"{print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' | grep -v "^#" | grep -v "mid*" |gffread - -o - | sed 's/transcript/mRNA/1' > "rf_blat_"$FEAT".gff3"
					mv "rf_blat_"$FEAT".gff3" "blat_"$FEAT".gff3"
					cp "repr_"$FEAT".gff" "Exo_repr_"$FEAT".gff"
					mv "blat_"$FEAT".gff3" "repr_"$FEAT".gff"
				else
					PHASE=5
				fi
			else
				echo "No ("$STOPS") stop codons found in translation of repr_"$FEAT".gff for phase "$PHASE" + FLIP: "$FLIP
				cleanExonerateGFF.pl repr_"$FEAT".gff $PHASE $FLIP n > tmp
				mv tmp repr_"$FEAT".gff	
				GOOD=1;

				LENA=$(gffread repr_"$FEAT".gff -g "$SCAFF".fa -x - | fasta_formatter | grep -v ">" | awk '{print length($0)}')

				DIFF=$(awk -v orig=$LENB -v new=$LENA 'BEGIN{printf "%.0f\n", (new/orig)*100}')
				echo "diff: "$DIFF

				echo "Length original ORF: "$LENB
				echo "Length repredicted: "$LENA

				#compare originial ORF length ($FEAT) against repredicted. Should be at least 90% of original. Otherwise print WARNING
				if [[ $DIFF -lt 90 ]]; then
					echo "WARNING: Repredicted sequence only "$DIFF"% of original sequence for "$FEAT"!"
				fi	
			fi
		done

		#BLAT CDS produced by exonerate back to the original genome
		#blat -t=dna -q=dna -fine -noHead "$SCAFF".fa "repr_"$file "repr_"$file".psl"

		#get best BLAT hit only
		#bestBLAT "repr_"$file".psl"

		#clean up best_hits file
		#mv "best_hits_repr_"$file".psl" "repr_"$file".psl"

		#Convert .psl file to GFF3
		#blat2gff.pl "repr_"$file".psl" > "mapBack_repr_"$file".gff3"

		#sed 's/Parent/ID/1' "mapBack_repr_"$file".gff3" | sed 's/_mid1//1' | sed 's/Target/Parent/1' | awk 'OFS="\t"{print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' | grep -v "^#" | grep -v "mid*" |gffread - -o - | sed 's/transcript/mRNA/1' > reformat_"$file".gff3

		#echo "$file on $SCAFF reformat_"$file".gff3:"
		#echo "Lines in reformat_"$file".gff3: $(wc -l reformat_"$file".gff3)"

		# ORIG=$(diff $file "repr_"$file)
			
		# if [[ -z $ORIG ]]; then

		# 	echo "Testing original ORF for stop codons"
		# 	gffread reformat_"$file".gff3 -g "$SCAFF".fa -w - | fasta_formatter > to_"$file".cds
		# 	translate_seqs.pl to_"$file".cds 1 > /dev/null 2>&1
		# 	seqstat to_"$file".pep
		# 	STOPS=$(fasta_formatter < to_"$file".pep | sed 's/\*$//' | grep -v ">" | fold -w1 | sort | uniq -c | grep "\*" | awk '{print $1}')

		# 	if [[ -z $STOPS ]]; then
		# 		echo "No ("$STOPS") stop codons found in "$PEP".pep !"
		# 		GOOD=1
		# 	else
		# 		echo "Still "$STOPS" stop codons found in "$PEP".pep !"
		# 	fi

		# 	#rm to_"$file".cds
		# 	#rm to_"$file".pep
		# fi


		if [[ $GOOD -eq 1 ]]; then

			if [[ $COUNT -le 1 ]]; then
				echo "\$COUNT is less/equal than 1. Creating new GFF file"
				#cat reformat_"$file".gff3 > $GFF
				cat repr_"$FEAT".gff > $GFF
			else
				echo "\$COUNT is greater than 1. Appending to existing GFF file"
				#sed -i '/^#/d' reformat_"$file".gff3
				sed -i '/^#/d' repr_"$FEAT".gff
				echo "# Gene gene:"$FEAT >> $GFF
				cat repr_"$FEAT".gff >> $GFF
				echo "###" >> $GFF
			fi

					
			# fasta_formatter < "repr_"$file > diff1
			# grep "$SCAFF" reformat_"$file".gff3 | gffread -g "$SCAFF".fa -w - | fasta_formatter > diff2

			# diffcheck=$(diff diff1 diff2)

			# if [[ -n $diffcheck ]]; then
			# 	echo $(date)
			# 	echo "\$diffcheck not empty! Differences: $diffcheck"
			# fi
			# unset diffcheck

			# rm diff1
			# rm diff2

			echo "SUCCESS! Done "$file". New prediction in $GFF."
		else
			echo "BOO! Done "$file". No good ORF found."
		fi

		echo ""
		echo ""
		#rm reformat_"$file".gff3
		#rm "mapBack_repr_"$file".gff3"
	else
		echo "WARNING: No (or empty) GFF file produced for "$file" ("$FEAT")!"
	fi
done
