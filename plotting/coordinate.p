reset

set key top left

#T=50.
#set xrange [0:T]
set xlabel 'Time'

L=1.5
set yrange [-L:L]
set ylabel 'Position'

name(i)=sprintf('%d',i)

unset key

n=1000
#plot for [i=1:n] 'output.dat' u 1:(column(2*i)) w p pt 1 lc -1 ti name(i)
plot for [i=1:n] 'data/output.dat' u 1:(column(2*i)) w l lw 1 ti name(i)
