set term postscript eps solid enhanced color
set tics out
set xtics nomirror
set ytics nomirror
set grid
set log y
set ylabel "Variance"
set xlabel "Average methylation"
set output "fig_variance.eps"
plot "./gimli_mean_var_le15.txt" u 3:4 pt 7 t "gimli" ,\
"./random_mean_var_le15.txt" u 3:4 pt 7 t "random"
