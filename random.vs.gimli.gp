set term postscript eps solid enhanced color
set tics out
set output "fig3.eps"
set log y
plot "./gimli_mean_var_le15.txt" u 3:4 , "./random_mean_var_le15.txt" u 3:4
