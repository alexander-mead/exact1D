set size square
set view map

set xlabel 'Position'
set xrange [-1:1]

set ylabel 'Velocity'
set yrange [-2:2]

set palette gray

splot 'data/phase_density.dat' w image
