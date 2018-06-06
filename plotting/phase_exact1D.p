reset

set size square

L=1.
set xrange [-L:L]
set xlabel 'Position'

v=2.*L
set yrange [-v:v]
set ylabel 'Velocity'

unset key

movie=1

n=500
if(movie==0){
#plot for [i=1:n] 'output.dat' u (column(2*i)):(column(2*i+1)) w p pt 1 lc -1
plot for [i=1:n] 'data/output.dat' u (column(2*i)):(column(2*i+1)) w l lw 3
pause 0.01
}

#pause -1

if(movie==1){
s=1001
do for [j=1:s] {
plot for [i=1:n] 'data/output.dat' every ::j::j u (column(2*i)):(column(2*i+1)) w p pt 1
}
}
