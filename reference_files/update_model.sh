#! /bin/bash

python write_msa_metrology.py

python write_esa_msa_map.py

python read_siaf_file.py

# gfortran -ffree-line-length-none -fno-backslash -fno-range-check  -O2 -o fits_siaf_read_v2v3 fits_siaf_read_v2v3.f90 \
# -L/usr/local/lib/pgplot/ -lpgplot -L/usr/lib -lgcc -lcfitsio -L/usr/X11/lib
# 
# ./fits_siaf_read_v2v3

gfortran -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o calc_derived_parameters calc_derived_parameters.f90 \
-L/usr/local/lib/pgplot/ -lpgplot -L/usr/lib -lgcc  -L/usr/X11/lib -lX11

./calc_derived_parameters

rm *.mod

rm calc_derived_parameters
rm fits_siaf_read_v2v3





