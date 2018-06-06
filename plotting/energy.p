reset

set xlabel 'Time'

set ylabel 'Kinetic / Potential energy'
set yrange [0:2]

unset key

plot 0.5 w l lt -1 noti,\
     'data/energy.dat' u 1:5 w l
