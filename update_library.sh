#!/bin/bash

libs+=('constants')
libs+=('logical_operations')
libs+=('array_operations')
libs+=('fix_polynomial')
libs+=('random_numbers')
libs+=('numerology')
libs+=('file_info')
libs+=('sorting')
#libs+=('field_operations')
#libs+=('interpolate')
#libs+=('table_integer')
#libs+=('fft')

for code in "${libs[@]}"; do
    cp /Users/Mead/Physics/library/$code.f90 src/.
done    
