reset

if(print==0) {set term aqua}

histo(i)=sprintf('histogram/data_%04d.dat',i)

ilog=0
if(ilog==0) {
xmin=0.
xmax=10.
unset log x
}
if(ilog==1) {
xmin=1e-3
xmax=1e3
set log x
}
set xrange [xmin:xmax]

ymin=0.
ymax=400.
#set yrange [ymin:ymax]
set yrange [0:1]

set cbrange [0:1]
set cblabel 'a'

n1=0
n2=4096
do for [i=n1:n2] {
plot histo(i) u 1:2:(real(i-1)/real(n-1)) w l lc palette lw 3 noti
}

#n=499
#plot for [i=1:n] histo(n-i+1) u 1:2:(real(i-1)/real(n-1)) w l lw 3 lc palette noti
