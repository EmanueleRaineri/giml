set term postscript eps enhanced solid color
set output "fig2.eps"
unset key
set tics out
set xlab "read depth"
set ylab "delta lik"
set grid
plot "delta.vs.cov.mean.txt" u 1:2:3 w err lw 1.5 lc rgb "blue"
