unset multiplot
reset

# Initial white space
print ''

# Set the ranges for the position and velocity y-axes
X=1.
V=2.

# File names
name(i)=sprintf('%d',i)

# Remove key
unset key

# Number of sheets
# Setting this higher than the actual number is okay, but decreases plotting speed
n=50
print 'Number of sheets: n: ', n

# Set margins for multiplot
set lmargin 10
set rmargin 2

# Time axis label
set xlabel 'Time'

# Final white space
print ''

# Set muliplot
set multiplot layout 2,1

# Set the position axis range
set yrange [-X:X]
set ylabel 'Position'

# Plot position data
#plot for [i=1:n] 'output.dat' u 1:(column(2*i)) w p pt 1 lc -1 ti name(i)
plot for [i=1:n] 'data/output.dat' u 1:(column(2*i)) w l lw 1 ti name(i)

# Set the velocty axis range
set yrange [-V:V]
set ylabel 'Velocity'

# Plot velocity data
plot for [i=1:n] 'data/output.dat' u 1:(column(2*i+1)) w l lw 1 ti name(i)

unset multiplot
