set terminal png
set output "scatterplot.png"

set title "Scatter Plot of uniform random number"


plot "scatterdata.dat" with points

set output "mean_vs_inverserootn.png"
set xlabel "1/sqrt(n)"

set title "Mean plot"
plot "meandata.dat" with lines

set output "correlation.png" 
set title "Correlation of random number"

set xlabel "k"
set ylabel "correlation"
plot "correlationdata.dat"


set output "correlation_small.png" 

set title "correlation with small value of k"
set xlabel "k"
set ylabel "correlation"
plot "correlationdata_small.dat" with lines
