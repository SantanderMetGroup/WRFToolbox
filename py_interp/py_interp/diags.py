#
# Diagnostics for py_interp
#
import numpy as np
from py_interp_fortran import routines as f90
from py_interp.fun import massvertint, de_stagger
tr = np.transpose

def compute_MSLP(ivar, inc, onc, bf, plevs):
    ovardata = f90.compute_mslp(tr(bf.pres_field), tr(bf.psfc), tr(bf.hgt), 
                                tr(bf.temp), tr(bf.qvapor))
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType  = 104
    ovarobj.MemoryOrder = "Z"
    ovarobj.description = "Pressure levels"
    ovarobj.units = "pa"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_GHT(ivar, inc, onc, bf, plevs):
    ovardata = f90.interp(tr(bf.ght), tr(bf.pres_field), plevs, tr(bf.psfc), 
                          tr(bf.hgt), tr(bf.temp), tr(bf.qvapor), linlog=1, 
                          extrapolate=1, geopt=True, missing=1.e36)
    dim_list = ["Time", "num_metgrid_levels", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType  = 104
    ovarobj.MemoryOrder = "XZY"
    ovarobj.description = "Geopotential Height"
    ovarobj.units = "m2s-2"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_PRES(ivar, inc, onc, bf, plevs):
    nj, ni = bf.pres_field.shape[2], bf.pres_field.shape[3]
    nplev = len(plevs)
    nt = bf.pres_field.shape[0]
    ovardata = np.repeat(np.nan, nt*nplev*nj*ni).reshape(nt, nplev, nj, ni)
    for n in range(nplev):
        ovardata[:, n, :, :] = plevs[n]
    
    dim_list = ["Time", "num_metgrid_levels", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)  
    ovarobj[:] = ovardata        
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XZY"
    ovarobj.description = "Pressure"
    ovarobj.units = "Pa"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_TT(ivar, inc, onc, bf, plevs):
    ovardata = f90.interp(tr(bf.temp), tr(bf.pres_field), plevs, tr(bf.psfc), 
                          tr(bf.hgt), tr(bf.temp), tr(bf.qvapor), linlog=1, 
                          extrapolate=1, geopt=False, missing=1.e36)
    
    dim_list = ["Time", "num_metgrid_levels", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XZY"
    ovarobj.description = "Temperature"
    ovarobj.units = "K"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_RH(ivar, inc, onc, bf, plevs):
    ovardata = f90.interp(tr(bf.rh), tr(bf.pres_field), plevs, tr(bf.psfc), 
                          tr(bf.hgt), tr(bf.temp), tr(bf.qvapor), linlog=1, 
                          extrapolate=1, geopt=False, missing=1.e36)
    
    dim_list = ["Time", "num_metgrid_levels", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)  
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType = 104 ;
    ovarobj.MemoryOrder = "XZY" ;
    ovarobj.description = "Relative Humidity" ;
    ovarobj.units = "%" ;
    ovarobj.stagger = "-" ;
    ovarobj.coordinates = "XLONG XLAT" ;
    return onc

def compute_PBLHh(ivar, inc, onc, bf, plevs):
    pt = inc.variables["T"][:] # Potential temperature
    ovardata = f90.pbl_height(tr(bf.ght), tr(pt), tr(bf.hgt))
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)  
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType = 104 ;
    ovarobj.MemoryOrder = "XZY" ;
    ovarobj.description = "PBL top height min theta + 1.5 method" ;
    ovarobj.units = "m" ;
    ovarobj.stagger = "-" ;
    ovarobj.coordinates = "XLONG XLAT" ;
    return onc

def compute_WSPD10(ivar, inc, onc, bf, plevs):
    ovardata = np.sqrt(inc.variables["U10"][:]**2 + inc.variables["V10"][:]**2)
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)  
    ovarobj[:] = ovardata
    ovarobj.FieldType = 104 ;
    ovarobj.MemoryOrder = "XZY" ;
    ovarobj.description = "Wind speed at 10 m" ;
    ovarobj.units = "m s-1" ;
    ovarobj.stagger = "-" ;
    ovarobj.coordinates = "XLONG XLAT" ;
    return onc

def compute_WSPD(ivar, inc, onc, bf, plevs):
    # 3D wind speed
    # Need to destagger 3D winds
    uobj = inc.variables["U"]
    u = de_stagger(uobj, uobj[:])
    vobj = inc.variables["V"]
    v = de_stagger(vobj, vobj[:])
    wspd = np.sqrt(u**2 + v**2)
    
    ovardata = f90.interp(tr(wspd), tr(bf.pres_field), plevs, tr(bf.psfc), 
                      tr(bf.hgt), tr(bf.temp), tr(bf.qvapor), linlog=1, 
                      extrapolate=1, geopt=False, missing=1.e36)
    
    dim_list = ["Time", "num_metgrid_levels", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)  
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType = 104 ;
    ovarobj.MemoryOrder = "XZY" ;
    ovarobj.description = "Wind speed" ;
    ovarobj.units = "m s-1" ;
    ovarobj.stagger = "-" ;
    ovarobj.coordinates = "XLONG XLAT" ;
    return onc

def compute_WGUST(ivar, inc, onc, bf, plevs):
    #
    # Wind gust following WRF-EMS 
    # method http://www.wrfforum.com/viewtopic.php?f=8&t=948
    # I prefer Brasseur (2008) but It requires TKE and some parameterizations
    # like YSU do not compute it.
    #
    wspd10 = np.sqrt(inc.variables["U10"][:]**2 + inc.variables["V10"][:]**2)
    # Need to destagger 3D winds
    uobj = inc.variables["U"]
    u = de_stagger(uobj, uobj[:])
    vobj = inc.variables["V"]
    v = de_stagger(vobj, vobj[:])
    
    wspd = np.sqrt(u**2 + v**2)
    
    if "PBLHh" not in onc.variables:
        msg = "PBLHh is required before WGUST in the variable list"
        raise RuntimeError(msg)
    
    pblh = onc.variables["PBLHh"][:]
    # Need to get the geopotential in meters above sea level in order to get
    # the level closest to the PBL top.
    ght_m = bf.ght/9.81
    ght_masl = ght_m - bf.hgt
    diff = np.abs(ght_masl - pblh[:, np.newaxis, :, :])
    pblh_levels = np.argmin(diff, axis=1)
    # Need to create a grid in order to index pblh correctly
    nt, nj, ni = pblh_levels.shape
    t, j, i = np.meshgrid(np.arange(nt), np.arange(nj), np.arange(ni), 
                          indexing='ij')
    
    wspd_pbl = wspd[t, pblh_levels, j, i]
    diff_wspd = wspd_pbl - wspd10
    # For PBLH > 1000 the difference between the wind gust and the 10m wind 
    # starts to be downweighted up to 0.5 when PBL > 2000 to avoid too large
    # speeds.
    diff_wspd = diff_wspd*(1.0 - np.where(pblh/2000. < 0.5, pblh/2000., 0.5))
    wgust = wspd10 + diff_wspd
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)  
    ovarobj[:] = wgust
    ovarobj.FieldType = 104 ;
    ovarobj.MemoryOrder = "XZY"
    ovarobj.description = "Wind gust WRF-EMS method"
    ovarobj.units = "m s-1"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_CLT_OLD(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_sundqvist(tr(cldfra))
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "total cloud fraction Sundqvist"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLT(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand(tr(cldfra))
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "total cloud fraction"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLL(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand_levels(tr(cldfra), tr(bf.pres_field), 
                                      maxpres=1000, minpres=680)
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Low cloud fraction (1000-680 hPa)"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLM(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand_levels(tr(cldfra), tr(bf.pres_field), 
                                      maxpres=680, minpres=440)
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Medium cloud fraction (680-440 hPa)"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLH(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand_levels(tr(cldfra), tr(bf.pres_field), 
                                      maxpres=440, minpres=10)
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "High cloud fraction (440-10 hPa)"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_VIM(ivar, inc, onc, bf, plevs):
    iarr = inc.variables["QVAPOR"][:]
    ovardata =  massvertint(iarr, inc)
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  ovardata
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Vertically integrated moisture"
    ovarobj.units = "Kg m-2"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_VIQC(ivar, inc, onc, bf, plevs):
    iarr = inc.variables["QCLOUD"][:]
    ovardata =  massvertint(iarr, inc)
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  ovardata
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Vertically integrated cloud water"
    ovarobj.units = "Kg m-2"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_VIQI(ivar, inc, onc, bf, plevs):
    iarr = inc.variables["QICE"][:]
    ovardata =  massvertint(iarr, inc)
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  ovardata
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Vertically integrated cloud ice"
    ovarobj.units = "Kg m-2"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_QCLOUDLL(ivar, inc, onc, bf, plevs):
    """
    Gets the lower hybrid level of the cloud water mixing ratio
    """
    qcloudll = inc.variables["QCLOUD"][:, 0, :, :]
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  qcloudll
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Cloud water mixing ratio in the lowest hybrid level"
    ovarobj.units = "Kg Kg-1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_QICELL(ivar, inc, onc, bf, plevs):
    """
    Gets the lower hybrid level of the ice mixing ratio
    """
    qicell = inc.variables["QICE"][:, 0, :, :]
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  qicell
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Ice mixing ratio in the lowest hybrid level"
    ovarobj.units = "Kg Kg-1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_QSNOWLL(ivar, inc, onc, bf, plevs):
    """
    Gets the lower hybrid level of the ice mixing ratio
    """
    qsnowll = inc.variables["QSNOW"][:, 0, :, :]
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  qsnowll
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Snow mixing ratio in the lowest hybrid level"
    ovarobj.units = "Kg Kg-1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_CAPE(ivar, inc, onc, bf, plevs):
    """
    Compute Convective Available Potential Energy
    """
    cape = f90.cape_3d(tr(bf.pres_field), tr(bf.temp), tr(bf.qvapor))
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  tr(cape)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Convective Available Potential Energy"
    ovarobj.units = "J Kg-1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_TLAPSE(ivar, inc, onc, bf, plevs):
    t2 = inc.variables["T2"][:]
    topo = inc.variables["HGT"][0, :]
    lmask = inc.variables["LANDMASK"][:].astype("bool")
    
    tlapse = f90.get_tlapser(tr(t2), tr(topo), tr(lmask[0]), 3)
    tlapse = np.where(lmask==0, -0.0098, tr(tlapse))
    tlapse = np.where(tlapse < -0.0098, -0.0098, tlapse)
    tlapse = np.where(tlapse > 0.1, 0.1, tlapse)
    tlapse = np.where(np.isnan(tlapse), -0.0098, tlapse)
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    
    ovarobj[:] =  tlapse
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Temperature lapse rate computed with linear regression"
    ovarobj.units = "K m-1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc 

def compute_TLAPSE200(ivar, inc, onc, bf, plevs):
    """
    Compute vertical lapse rate of temperature in the lowest 200 m.
    """
    ght_meters = bf.ght/9.8
    height_over_surface = ght_meters - bf.hgt
    maxheight = 200.
    indices_level = np.argmin(height_over_surface - maxheight, axis=1)
    nt, nj, ni = indices_level.shape
    t, j, i = np.meshgrid(np.arange(nt), np.arange(nj), np.arange(ni), 
                          indexing='ij')
    indices_level = np.where(indices_level == 0, 1, indices_level)
    ght_max = ght_meters[t, indices_level, j, i]
    ght_min = ght_meters[t, 0, j, i]
    ght_delta = ght_max - ght_min
    ta_max = bf.temp[t, indices_level, j, i]
    ta_min = bf.temp[t, 0, j, i]
    ta_delta = ta_max - ta_min
    lapse_rate = ta_delta/ght_delta # K m-1
    
    dim_list = ["Time", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, 
                                 complevel=4, shuffle=True)
    ovarobj[:] =  lapse_rate
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Temperature lapse rate"
    ovarobj.units = "K m-1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_OMEGA(ivar, inc, onc, bf, plevs):
    """
    Calculate approximate omega, based on vertical velocity w (dz/dt).
    It is approximate because it cannot take into account the vertical
    motion of pressure surfaces.
    """
    temp = bf.temp
    pres = bf.pres_field
    qvapor = bf.qvapor
    wobj = inc.variables["W"]
    vel_w = de_stagger(wobj, wobj[:])
    grav = 9.81
    rgas = 287.04
    eps = 0.622
    omega = -grav*pres/(rgas*((temp*(eps + qvapor))/(eps*(1.0 + qvapor))))*vel_w

    ovardata = f90.interp(tr(omega), tr(bf.pres_field), plevs, tr(bf.psfc),
                          tr(bf.hgt), tr(bf.temp), tr(bf.qvapor), linlog=1,
                          extrapolate=1, geopt=False, missing=1.e36)
    dim_list = ["Time","num_metgrid_levels", "south_north", "west_east"]
    ovarobj = onc.createVariable(ivar, 'float32', dim_list, zlib=True, complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XYZ"
    ovarobj.description = "Approximate omega"
    ovarobj.units = "Pa s-1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
 
