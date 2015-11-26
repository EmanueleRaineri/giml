BEGIN{OFS="\t"}
{
	split($5,pos,":");
	split($6,nc,":");
	split($7,c,":");
	for(i=1;i<=length(pos);i++){
		met[i]=nc[i]/(nc[i]+c[i])
	};
	printf("%s\t%d\t%d\t%d\t%d\n", $1,$2,$3,$4,$3-$2,length(pos));
}

#Rscript -e "v<-as.numeric(readLines(file('stdin')));print(v)"
