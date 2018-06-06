reset

file(i)=sprintf('particles/data_%04d.dat',i)

L=1.
set xrange [0:L]
set xlabel 'Mpc/h'

v=0.001
set yrange [-v:v]
set ylabel 'km s^{-1}'

n1=499
n2=499

do for [i=n1:n2]{
print i
plot file(i) u 1:2 w p noti
}
