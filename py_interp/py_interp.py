########################################################################
#
# py_interp alpha by Markel Garcia Diez February 2014
#
# See the README for info and/or type python py_interp.py -h for help.
#
########################################################################
import numpy as np
import netCDF4 as ncdf
import os, sys
import py_interp_diags
from py_interp_fun import copy_n_filter_wrfout, BasicFields, add_pressure_axis, interp2plevs
from optparse import OptionParser
#
# Parse options
#
parser = OptionParser()
parser.set_defaults(quiet=False,singlerec=False)
parser.add_option("-i", dest="ifile", help="Input wrfout file")
parser.add_option("-v", dest="varlist", default="", help="Comma sepparated variable list")
parser.add_option("-p", "--plevs", dest="plevs", help="Comma sepparated pressure level list in hPa")
parser.add_option("--verbose", dest="verbose", action="store_true", default=False, help="Get a more detailed log")
(opt, args) = parser.parse_args()
#
# Some global variables
#
ifile = opt.ifile
ofile = ifile + "_PLEV"
varlist = opt.varlist.split(",")[:]
plevs = np.array(opt.plevs.split(",")[:], dtype=np.float32)
print "Starting py_interp"
if opt.verbose:
	print "Input file: %s" % ifile
	print "Output file: %s" % ofile
	print "Variables to be processed: %s" % varlist
	print "Levels to interpolate: %s" % plevs
if not os.path.exists(ifile):
	print "Input file %s does not exist" % ifile
	sys.exit(1)
#
# Reading input wrfout file
#
print "Reading file %s" % ifile
inc = ncdf.Dataset(ifile, "r")
#
# Separate the input variables in three lists: 2D variables (or soil variables)
# to copy, 3D variables to interpolate and diagnostics that need specific
# functions to be computed.
#
dims2D=("Time", "south_north", "west_east")
diagnostics = ["VIQC", "VIQI", "VIM", "MSLP", "CLT", "CLT_OLD", "CLH", "CLM", "CLL", "GHT", "PRES", "RH", "TT"]
copyvars = []
interpvars = []
diags = []
for var in varlist:
    if var in diagnostics:
        diags.append(var)
    else:
        if var not in inc.variables.keys():
            print "Error: Variable %s not found in file, and is not a diagnostic" % var
            sys.exit(1)
        vdims = inc.variables[var].dimensions
        if ("bottom_top" in vdims) or ("bottom_top_stag" in vdims):
            interpvars.append(var)
        else:
            copyvars.append(var)

if opt.verbose:
	print "Variables to copy %s" % copyvars
	print "3D variables to interpolate %s" % interpvars
	print "Diagnostics to compute %s" % diags
#
# Pressure levels from hPa to Pa.
#
plevs = plevs*100
#
# Get some basic fields: psfc, pres_field, ght, hgt, temp
#
bf = BasicFields(inc)
#
# Copy wrfnc structure with the 2D vars,  which are left untouched
#
if opt.verbose:
	print "Copying wrfout file structure and 2D fields"
onc = copy_n_filter_wrfout(inc, ofile, copyvars)
#
# Add the pressure axis and variable
#
if opt.verbose: print "Adding pressure axis to output file"
onc = add_pressure_axis(onc, plevs)
#
# Interpolate to pressure levels 3D variables
#
for var in interpvars:
	if opt.verbose: print "Interpolating %s to pressure levels" % var
	onc = interp2plevs(var, inc, onc, bf, plevs)
#
# Compute diagnostics
#
for var in diags:
    if opt.verbose: print "Computing diagnostic %s" % var
    #
    # Use getattr to get the function from py_interp_diags.py and then call it
    #
    compute_diag = getattr(py_interp_diags, "compute_%s" % var)
    onc = compute_diag(var, inc, onc, bf, plevs)
onc.sync()
print "SUCCESS: p_interp finished without errors"