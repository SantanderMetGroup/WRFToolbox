import numpy   as np
import netCDF4 as ncdf
import logging
from scipy.ndimage.filters import gaussian_filter
import py_interp.diags
from py_interp.fun import ( copy_n_filter_wrfout, BasicFields, 
                            add_pressure_axis, interp2plevs )
log = logging.getLogger(__name__)
#
# Main function
#
def py_interp_file(ifile, varlist, plevs, vars_to_filter=None, filter_std=1, overwrite=False, verbose=False):
    """
    Interpolate WRF raw files to pressure levels and get derived variables
    """
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    #
    # Some global variables
    #
    ofile = ifile + "_PLEV"
    plevs = np.array(plevs, dtype=np.float32)
    log.info("Starting py_interp")
    log.debug("Input file: {}".format(ifile))
    log.debug("Output file: {}".format(ofile))
    log.debug("Variables to be processed: {}".format(varlist))
    log.debug("Levels to interpolate: {}".format(plevs))
    #
    # Reading input wrfout file
    #
    log.info("Reading file {}".format(ifile))
    inc = ncdf.Dataset(ifile, "r")
    #
    # Classify the variables in those to copy, to interpolate and to diagnose
    #
    copyvars, interpvars, diags = classify_variables(varlist, inc)
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
    log.debug("Copying wrfout file structure and 2D fields")
    onc = copy_n_filter_wrfout(inc, ofile, copyvars)
    #
    # Add the pressure axis and variable
    #
    log.debug("Adding pressure axis to output file")
    onc = add_pressure_axis(onc, plevs)
    #
    # Interpolate to pressure levels 3D variables
    #
    for var in interpvars:
        log.debug("Interpolating {} to pressure levels".format(var))
        onc = interp2plevs(var, inc, onc, bf, plevs)
    #
    # Compute diagnostics
    #
    for var in diags:
        log.debug("Computing diagnostic {}".format(var))
        #
        # Use getattr to get the function from py_interp_diags.py and then call it
        #
        compute_diag = getattr(py_interp.diags, "compute_{}".format(var))
        onc = compute_diag(var, inc, onc, bf, plevs)
    onc.sync()
    #
    # Apply a gaussian filter to certain variables (only two dimensions!!)
    #
    if vars_to_filter is not None:
        onc = filter_variables(onc, vars_to_filter, filter_std)
    
    if overwrite:
        import shutil
        shutil.move(ofile,ifile)
    
    log.info("SUCCESS: p_interp finished without errors")
    
def classify_variables(varlist, inc):
    '''
    Separate the input variables in three lists: 2D variables 
    (or soil variables) to copy, 3D variables to interpolate and diagnostics 
    that need specific functions to be computed.
    '''
    
    diagnostics = ["VIQC", "VIQI", "VIM", "MSLP", "CLT", "CLT_OLD", "CLH", 
                   "CLM", "CLL", "GHT", "PRES", "RH", "TT", "PBLHh", "WSPD",
                   "WSPD10", "WGUST", "QCLOUDLL", "QICELL", "QSNOWLL", "CAPE", 
                   "TLAPSE", "TLAPSE200", "OMEGA", "MSLP_LEGACY"]
    copyvars = []
    interpvars = []
    diags = []
    for var in varlist:
        if var in diagnostics: 
            diags.append(var)
        else:
            if var not in inc.variables.keys():
                raise RuntimeError("Error: Variable {} not found in file, and is"
                                    " not a diagnostic".format(var))
            vdims = inc.variables[var].dimensions
            if ("bottom_top" in vdims) or ("bottom_top_stag" in vdims):
                interpvars.append(var)
            else: 
                copyvars.append(var)
    
    log.debug("Variables to copy {}".format(copyvars))
    log.debug("3D variables to interpolate {}".format(interpvars))
    log.debug("Diagnostics to compute {}".format(diags))
    return copyvars, interpvars, diags

def filter_variables(onc, vars_to_filter, filter_std):
    log.debug("Variables {} are going to be filtered".format(vars_to_filter))
    for varname in vars_to_filter:
        onc.variables[varname][:] = gaussian_filter(onc.variables[varname][:], [0, filter_std, filter_std])
    return onc
