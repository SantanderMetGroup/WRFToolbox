import netCDF4 as ncdf
import datetime
import sys
from ncdf2wps_fortran import routines as f90

class gridObj:
    def __init__(self, inc, proj="lonlat", lonlatnames = ("lon", "lat")):
        self.proj = proj
        self.lons = inc.variables[lonlatnames[0]][:]
        self.lats = inc.variables[lonlatnames[1]][:]
        
class varObj:
    def __init__(self, cfname, wrfname, units, desc, levels, source):
        self.cfname = cfname
        self.wrfname = wrfname
        self.units = units
        self.desc = desc
        self.levels = levels
        self.source = source

def get_timedelta_hours(tdelta):
    return tdelta.days*24 + tdelta.seconds/3600

def get_datelist(idate, fdate, freq):
    delta_hours = get_timedelta_hours(freq)
    nsteps = get_timedelta_hours(fdate - idate)/delta_hours + 1
    numdates = range(0, nsteps*delta_hours, delta_hours)
    datetimelist = ncdf.num2date(numdates, units="hours since %s" % (idate.isoformat(sep=" ")))
    return datetimelist

def writewps(data, ofile, varobj, units, date, gridobj):

    if gridobj.proj == "latlon":
        datestr  = date.isoformat("T")
        #
        # Arguments for writewps
        #
        writewps_args = {
        "ofile"    : ofile,
        "slab"     : tr(data),
        "hdate"    : datestr,
        "xfcst"    : 0.0,
        "startloc" : "SWCORNER",
        "field"    : varobj.wrfname,
        "units"    : varobj.units,
        "desc"     : varobj.desc,
        "map_source" : varobj.source,
        "xlvl"       : 200100.,
        "iproj"      : 0,
        "startlat"   : gridobj.lats[0],
        "startlon"   : gridobj.lons[0],
        "deltalat"   : gridobj.lats[1] - gridobj.lats[0],
        "deltalon"   : gridobj.lons[1] - gridobj.lons[0],
        "dx"         : 0,
        "dy"         : 0,
        "xlonc"      : 0,
        "truelat1"   : 0,
        "truelat2"   : 0,
        "nlats"      : 0, 
        }
        if len(varobj.levels) == 1:
            writewps_args["xlvl"] = varobj.levels[0]
            f90.writewps(**writewps_args) 
        else:
            for lev in varobj.levels:
                writewps_args["xlvl"] = lev
                f90.writewps(**writewps_args)
    else:
        print "Projection {0} not supported yet".format(gridobj.proj)
        sys.exit(1)