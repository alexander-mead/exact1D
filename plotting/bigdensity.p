reset

if(print==0) {set term aqua}
if(print==1) {set term png size 4096,3072; set output 'big.png'}

set view map

L=1.
set xlabel 'Mpc/h'
set xrange [0:L]

set ylabel 'log_{10}(a-a_{min})'
set yrange [0.:1.]

#set cblabel '1+{/Symbol d}'
if(ilog==0){
dmin=0.999
dmax=1.001
unset log cb
print 'Doing linear density plot'
print 'Minimum density: ', dmin
print 'Maximum density: ', dmax
}
if(ilog==1){
dmin=1e-1
dmax=1e1
set log cb
print 'Doing logarithmic density plot'
print 'Minimum density: ', dmin
print 'Maximum density: ', dmax
}
set palette defined (1 'black', 2 'white')
set cbrange [*:*]
set cbrange [dmin:dmax]

if(big==1){
unset xtics
unset ytics
unset colorbox
set tmargin at screen 1.
set bmargin at screen 0.
set lmargin at screen 0.
set rmargin at screen 1.
}

splot 'density/bigdensity.dat' u 1:2:(1+$3) w image noti
