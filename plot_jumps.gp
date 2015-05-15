set xra[1655618-4000:1656083+4000]
set yrange[-0.1:1.1]
unset key
set grid
set border 2
set xtics nomirror
set ytics nomirror
plot "C000S5A1bs_cpg.chr1.1655618.1656083.txt" u 2:4 pt 7,\
	"C0010KA2bs_cpg.chr1.1655618.1656083.txt" u 2:4 pt 7,\
	"C001UYA3bs_cpg.chr1.1655618.1656083.txt" u 2:4 pt 7,\
	"C004SQ51_cpg.chr1.1655618.1656083.txt" u 2:4 pt 7
