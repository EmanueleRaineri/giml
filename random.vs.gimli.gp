set term postscript eps solid enhanced color
set tics out
set xtics nomirror
set ytics nomirror
set grid
set log y
set output "fig3.eps"
plot "./gimli_mean_var_le15.txt" u 3:4 , "./random_mean_var_le15.txt" u 3:4
