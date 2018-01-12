awk -F"\t" '{id[$1]=1; anno[$4]=1; code[$5]=1; anno["IPRO"]=1; anno["GO"]=1; look[$1","$4","$5]=$5; if($12!=""){code[$12]=1; look[$1",IPRO,"$12]=$12}; if($14!=""){code[$14]=1; look[$1",GO,"$14]=$14};}



END{	for(i in id){
	printf i"\t"; 
		for(k in anno){
			for(n in code){if(look[i","k","n]){printf look[i","k","n]"|"}};

		printf "\t";}
	printf "\n";}
}' $1 | sed 's/|/;/g' | sed 's/;\t/\t/g' | sed 's/\t$//' | awk -F"\t" '{split($0,line,/\t/); split(line[4], GO, /;/); for(i in GO){numGO[GO[i]]=GO[i]}; for(i in numGO){n++; uniqGO[n]=i}; n=0; line[4]=uniqGO[1]; for(i=2; i<=length(uniqGO); i++){line[4]=line[4]";"uniqGO[i]}; out=line[1]; for(i=2; i<=length(line); i++){out=out"\t"line[i]}; print out; delete line; delete GO; delete uniqGO; delete numGO; out="";}'
