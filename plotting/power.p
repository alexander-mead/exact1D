reset

if(print==0) {set term aqua dashed}

power(i)=sprintf('power/data_%04d.dat',i)

set log x
set xlabel 'kL / 2{/Symbol p}'
set xrange [1:1e5]
set format x '10^{%T}'
set mxtics 10

set log y
set ylabel 'Half of {/Symbol D}^2_{1D}(k) = 2kL P(k) / 2{/Symbol p}'
set format y '10^{%T}'

set palette defined (1 'pink', 2 'red', 3 'black')
unset colorbox

#Initial spectrum amplitude and index
amp=1e-6
ind=2

#Number of particles (for shot noise)
np=2**20

#Timesteps
s=128

#Box size
L=1.

print 'Assumed initial amplitude: ', amp
print 'Assumed initial index: ', ind
print 'Assumed number of particles: ', np
print 'Assumed box size [Mpc/h]: ', L
print 'Timesteps: ', s

set key top left

plot for [i=0:s] power(i) u ($1*L/(2.*pi)):($3):((i-1.)/(s-1.)) w l lc palette noti,\
     amp*x**(ind+1) lw 3 lc -1 ti 'Initial spectrum',\
     (1./(np/L))*(2.*x*L/(2.*pi)) lw 3 dt 2 lc -1 ti 'Shot noise'

