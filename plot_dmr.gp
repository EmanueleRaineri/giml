unset key
plot "C000S5A1bs.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "blue",\
"C0010KA2bs.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "blue",\
"C004SQ51.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "blue",\
"C005PS51.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "blue",\
"S000RD54.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "blue",\
"S007G756.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "blue",\
"C005VG51.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "red",\
"S001S751.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "red",\
"S0039051.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "red",\
"S00BHQ51.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "red",\
"S00DVR51.DMR.chr1.27994303.27994356.txt" u 2:($6/($6+$7)) pt 7 lc rgb "red"
set style arrow 1 nohead back nofilled lw 2 lc rgb "green"
set arrow 1 from 27994303,0 to 27994303,1 arrowstyle 1
set arrow 2 from 27994356,0 to 27994356,1 arrowstyle 1
