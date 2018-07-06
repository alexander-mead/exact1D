#!/bin/bash

libs+=('array_operations')
libs+=('constants')
libs+=('fix_polynomial')
libs+=('random_numbers')
libs+=('numerology')
libs+=('file_info')
libs+=('sorting')

for code in "${libs[@]}"; do
    cp /Users/Mead/Physics/library/$code.f90 .
done    
