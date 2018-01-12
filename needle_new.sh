#!/bin/bash

while read line
do
tmp1=">$(echo "$line" | awk -F"\t" '{print $1}')"
tmp2=">$(echo "$line" | awk -F"\t" '{print $2}')"

prefix=$(echo "$line" | sed 's/\t/_/')


if [ "$tmp2" == ">NA" ]
then 
echo $tmp1"	N/A		"| sed 's/>//g' >> output_$1
continue
fi

if [ "$tmp2" != ">NA" ]
then
echo "prefix= $prefix"
echo "tmp1= $tmp1"
echo "tmp2= $tmp2"

grep -x -A1 $tmp1 $2 > tmp1.fasta
grep -x -A1 $tmp2 $2 > tmp2.fasta

needle -asequence tmp1.fasta -bsequence tmp2.fasta -gapopen 10 -gapextend 0.5 -outfile "$prefix".needle 2>&1 > /dev/null

# Length: 420
# Identity:     122/420 (29.0%)
# Similarity:   178/420 (42.4%)
# Gaps:         194/420 (46.2%)

#length=$(grep "# Length:" "$prefix".matcher)
ident=$(grep "# Identity:" "$prefix".needle)
simil=$(grep "# Similarity:" "$prefix".needle)
#gaps=$(grep "# Gaps:" "$prefix".matcher)

echo $tmp1"	"$tmp2"	"$ident"	"$simil | perl -n -e '/(^.+?\t.+?)\t.+?\((.+?)%\)\t.+?\((.+?)%\)/ and print "$1\t$2\t$3\n"' | sed 's/>//g' >> output_$1


rm "$prefix".needle
fi

done < $1

rm tmp1.fasta
rm tmp2.fasta


