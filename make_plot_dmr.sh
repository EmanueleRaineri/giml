#C000S5A1bs.DMR.chr5.28927763.28928071.txt
#monocytes
echo "unset key"
echo -n "plot "
for m in C000S5A1bs C0010KA2bs C004SQ51 C005PS51 S000RD54 S007G756;  do
	i=${m}.DMR.$1.$2.$3.txt
	echo "\"$i\" u 2:(\$6/(\$6+\$7)) pt 7 lc rgb \"blue\",\\" 
done
#macrophages
for m in C005VG51 S001S751 S0039051 S00BHQ51; do
	i=${m}.DMR.$1.$2.$3.txt
	echo "\"$i\" u 2:(\$6/(\$6+\$7)) pt 7 lc rgb \"red\",\\" 
done
m=S00DVR51
i=${m}.DMR.$1.$2.$3.txt
echo "\"$i\" u 2:(\$6/(\$6+\$7)) pt 7 lc rgb \"red\"" 
echo "set style arrow 1 nohead back nofilled lw 2 lc rgb \"green\""
echo "set arrow 1 from $2,0 to $2,1 arrowstyle 1"
echo "set arrow 1 from $3,0 to $3,1 arrowstyle 1"
