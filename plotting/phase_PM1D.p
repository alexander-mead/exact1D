reset

if(print==0) {set term aqua}
if(print==1) {set term gif animate enh size 512,394 delay 1; set output 'phase.gif'}

phase(i)=sprintf('phase/data_%04d.dat',i)

set view map

L=1.
set xlabel 'Mpc/h'
set xrange [0:L]
set xtics

v=0.001
v=1e-5
set ylabel 'km s^{-1}'
set yrange [-v:v]
#set yrange [*:*]
set ytics

set cblabel '1+{/Symbol d}'
if(ilog==0){
dmin=0.
dmax=0.0025
unset log cb
}
if(ilog==1){
dmin=1e-1
dmax=1e2
set log cb
}
set cbrange [dmin:dmax]
set format cb '10^{%T}'
#set cbrange [*:*]
set colorbox

#set cbrange[0:100]

set palette defined (1 'black', 2 'white')

if(print==1) {
set tmargin at screen 1
set bmargin at screen 0
set lmargin at screen 0
set rmargin at screen 1
set xlabel ''
unset xtics
set ylabel ''
unset ytics
unset colorbox
}

n1=0
n2=128
do for [i=n1:n2] {
print i, ' ', phase(i)
splot phase(i) u 1:2:3 w image noti
}
