set datafile separator ','
set terminal svg
set key autotitle columnhead
set xlabel "Frequency (Hz)"
set ylabel "Magnitude"
set key off
set yrange[0:1]

set title "W3HH Region of Interest (wsprd)"
set output "W3HH_region_of_interest.svg"
plot "W3HH_region_of_interest.csv" using 1:2 with lines

set title "W5BIT Region of Interest (wsprd)"
set output "W5BIT_region_of_interest.svg"
plot "W5BIT_region_of_interest.csv" using 1:2 with lines
