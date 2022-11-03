#! /bin/bash
set -e

gfortran -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o ipa ipa_v2.4.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc   -L/usr/X11/lib -lX11 

gfortran  -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o k_make k_make_v1.7.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc   -L/usr/X11/lib -lX11

gfortran  -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o k_clean k_clean_v1.7.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc   -L/usr/X11/lib -lX11

gfortran  -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o m_make m_make_v2.0.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc   -L/usr/X11/lib -lX11

gfortran  -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o m_sort m_sort_v1.7.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc   -L/usr/X11/lib -lX11

gfortran  -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o m_pick m_pick_v2.4.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc   -L/usr/X11/lib -lX11

gfortran  -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o m_check m_check_v1.7.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc   -L/usr/X11/lib -lX11

gfortran  -std=legacy -ffree-line-length-none -fno-backslash -fno-range-check -O2 -o m_check_regions m_check_regions_v1.6.f90 \
-L/opt/local/lib/ -lpgplot -L/usr/local/lib -lgcc  -L/usr/X11/lib -lX11


rm *.mod


