#! /bin/csh

gfortran -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o calc_derived_parameters calc_derived_parameters.f90 \
-L/usr/local/lib/pgplot/ -lpgplot -L/usr/lib -lgcc  -L/usr/X11/lib -lX11

./calc_derived_parameters

rm *.mod

rm calc_derived_parameters





