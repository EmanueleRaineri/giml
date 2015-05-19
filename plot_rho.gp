set term postscript eps solid color enhanced
set xra[0:50000]
set xtics nomirror
set ytics nomirror
set grid
d0=10000
d1=100
set output "figrho.eps"
plot d1/(1+exp(-x/d0))-d1/2+1 lw 2
