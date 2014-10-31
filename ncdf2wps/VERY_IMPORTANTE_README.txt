#
# ncdf2wps by Markel García-Díez Barcelona 31-10-2014
#
# This program can read netCDF files and write the data into the intermediate 
# UNFORMATTED fortran files used by the WRF Preprocessor (WPS). It is meant to 
# replace ungrib.exe. UNFORMATTED means that these files do not have an actual format,
# as netCDF or GRIB files. Their format depends on the architecture, OS, 
# compiler and compiler options. Thus, these intermediate WPS files are NOT PORTABLE.
# Please DO NOT PORT THEM!!
#
# ncdf2wps_fortran.F90 needs to be compiled with EXACTLY the same compiler, compiler version,
# OS and architecture used to compile WPS (metgrid.exe). Ignoring this can result on
# errors or wrong simulations. Please be very careful with projection and grid specifications,
# as well as with the orientation of the arrays fed into ncdf2wps.
#
# ncdf2wps_fortran.F90 can be compiled with the following command, in a system with 
# ifort (IFORT) 14.0.1 20131008 installed. The flags --f90flags="-FR -convert big_endian"  
# are needed as they are present too in ungrib Makefile.

f2py --fcompiler=intelem -c ncdf2wps_fortran.F90 -m ncdf2wps_fortran --f90flags="-FR -convert big_endian"
