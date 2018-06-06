#!/bin/bash

code=$1
if [ -z code ]; then
    echo 'Please specify a supported code to compile'
    exit 1
fi
if [ "$code" = "preIC1D" ]; then
    echo 'Compiling preIC1D'
elif [ "$code" = "IC1D" ]; then
    echo 'Compiling IC1D'
elif [ "$code" = "PM1D" ]; then
    echo 'Compiling PM1D'
elif [ "$code" = "exact1D" ]; then
    echo 'Compiling exact1D' 
else
    echo 'Please specify a supported code to compile'
    exit 1
fi

#Set the compiler
compiler=gfortran

#Normal compiler options
normal='-Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3 -L/usr/local/lib -lfftw3'

#Mead library location
mead='/Users/Mead/Physics/library'



#preIC1D
if [ "$code" = "preIC1D" ]; then
    precision='-fdefault-real-8'
    $compiler $mead/constants.f90 $mead/random_numbers.f90 src/simulations1D.f90 src/preIC1D.f90 $normal -o preIC1D.e
fi

#IC1D
if [ "$code" = "IC1D" ]; then
    precision='-fdefault-real-8'
    $compiler $mead/constants.f90 $mead/file_info.f90 $mead/string_operations.f90 $mead/special_functions.f90 $mead/table_integer.f90 $mead/fix_polynomial.f90 $mead/array_operations.f90 $mead/interpolate.f90 $mead/calculus_table.f90 $mead/random_numbers.f90 $mead/fft.f90 $mead/cosmology_functions.f90 src/simulations1D.f90 src/IC1D.f90 $normal $precision -o IC1D.e
fi

#PM1D
if [ "$code" = "PM1D" ]; then
    precision='-fdefault-real-8'
    $compiler $mead/constants.f90 $mead/file_info.f90 $mead/string_operations.f90 $mead/special_functions.f90 $mead/table_integer.f90 $mead/fix_polynomial.f90 $mead/array_operations.f90 $mead/interpolate.f90 $mead/calculus_table.f90 $mead/random_numbers.f90 $mead/fft.f90 $mead/cosmology_functions.f90 src/simulations1D.f90 src/PM1D.f90 $normal $precision -o PM1D.e
fi

#exact1D
if [ "$code" = "exact1D" ]; then
    precision=''
    #precision='-fdefault-real-16'
    $compiler $mead/constants.f90 $mead/fix_polynomial.f90 $mead/random_numbers.f90 $mead/array_operations.f90 src/exact1D.f90 $normal $precision -o exact1D.e
fi

#Clean up
rm *.mod
