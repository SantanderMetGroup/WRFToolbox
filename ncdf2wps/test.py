#
# Minimal example of writing WPS format
#
import netCDF4 as ncdf
import numpy as np
from ncdf2wps_fortran import routines as f90
tr = np.transpose
cfname2wrf = {
"orog" : "SOILHGT"
}
cfvar = "orog"
ifile = "./test/orog_fx_HadGEM2-ES_historical_r1i1p1.nc"
print "Reading file {0}".format(ifile)
#
# Read the data
#
inc   = ncdf.Dataset(ifile, "r")
varobj = inc.variables[cfvar]
vardata = inc.variables[cfvar][:].squeeze()
lon = inc.variables["lon"][:]
lat = inc.variables["lat"][:]
#
# Somes definitions
#
field   = cfname2wrf[cfvar]
prefix = "OROG"
units = "m"
date  = "2001-01-01T00:00:00"
ofile = "./test/{0}_{1}".format(prefix, date)
#
# Arguments for writewps
#
writewps_args = {
"ofile"    : ofile,
"slab"     : tr(vardata),
"hdate"    : date,
"xfcst"    : 0.0,
"startloc" : "SWCORNER",
"field"    : field,
"units"    : units,
"desc"     : "Sea Surface Temperature",
"map_source" : "Met Office",
"xlvl"       : 200100.,
"iproj"      : 0,
"startlat"   : lat[0],
"startlon"   : lon[0],
"deltalat"   : lat[1] - lat[0],
"deltalon"   : lon[1] - lon[0],
"dx"         : 0,
"dy"         : 0,
"xlonc"      : 0,
"truelat1"   : 0,
"truelat2"   : 0,
"nlats"      : 0, 
}

f90.writewps(**writewps_args) 
writewps_args["field"] =  "OROG2"
f90.writewps(**writewps_args)
