reset

set size square

print ''

L=1.
set xrange [-L:L]
set xlabel 'Position'

v=2.*L
set yrange [-v:v]
set ylabel 'Velocity'

unset key

if(!exists('movie')){movie=1}
if(!exists('n')){n=500}
if(!exists('s')){s=1001}

print 'Number of sheets: n: ', n
print 'Number of timesteps: s: ', s
print 'Movie option: movie: ', movie
print ''

if(movie==0){
plot for [i=1:n] 'data/output.dat' u (column(2*i)):(column(2*i+1)) w l lw 3
pause 0.01
}

#pause -1

if(movie==1){
do for [j=1:s] {
plot for [i=1:n] 'data/output.dat' every ::j::j u (column(2*i)):(column(2*i+1)) w p pt 1
}
}
