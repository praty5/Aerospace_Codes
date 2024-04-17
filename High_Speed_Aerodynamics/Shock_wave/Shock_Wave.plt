set output 'plot.png'

set terminal pngcairo enhanced font "Arial,12" size 1200,800

set title 'Data Plot'
set xlabel 'M_1'
set xrange [1:5.1]

set ylabel 'rho, P, T'
set yrange [0:10.2]

set datafile separator ','

plot 'data.csv' using 1:2 with lines linewidth 2 title 'M_2', \
     'data.csv' using 1:3 with lines linewidth 2 title '{/Symbol r}2/{/Symbol r}1', \
     'data.csv' using 1:4 with lines linewidth 2 title 'p2/p1', \
     'data.csv' using 1:5 with lines linewidth 2 title 'T2/T1', \
     'data.csv' using 1:6 with lines linewidth 2 title 'Po2/Po1'
