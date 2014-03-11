use intel
export NETCDFDIR="/software/meteo/netcdf/netcdf_4.1.3_intel11_netcdf4"
export HDF5DIR="/software/meteo/hdf5/hdf5_1.8.8_intel11"
ifort  -O3 -heap-arrays p_interp.F90 -o p_interp -I$NETCDFDIR/include -I$HDF5DIR/include -Bstatic -L$NETCDFDIR/lib -lnetcdf -lnetcdff -L$HDF5DIR/lib -lnetcdf -lhdf5_hl -lhdf5  -lz >& compile.log
