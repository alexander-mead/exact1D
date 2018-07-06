reset

# Set the bin size
binwidth=0.01

# Parameters for size of x axis
L=1
min=-L
max=L
set xrange [min:max]

# Necessary
bin(x,width) = width*(floor((x-min)/width)+0.5) + Min
set boxwidth binwidth

plot 'data/end.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes

