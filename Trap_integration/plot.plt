set term png

set output "mean_vs_rootn.png"

set title " Error vs Inverse Square"
set xlabel "1/n**2"

set ylabel "error"

plot "error.dat" with lines


set output "error_vs_n.png"
set title "Error vs n"
set xlabel "n"

set ylabel "error"

plot "error_vs_n.dat"
