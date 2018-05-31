import xarray as xr
import netCDF4 as ncdf
import numpy as np
import scipy.stats
from itertools import product
from joblib import Parallel, delayed
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

def get_tlapse_regression(topo_tile, t2_tile, lmask_tile):
    if np.sum(lmask_tile ) < 25:
        print("Less than 25 earth points, returning nan")
        return np.nan
    topo_tile = topo_tile[lmask_tile]
    t2_tile = t2_tile[lmask_tile]
    ln = LinearRegression()
    X = topo_tile.flatten().reshape(-1, 1)
    y = t2_tile.flatten()
    ln.fit(X, y)
    #res = scipy.stats.linregress(topo_tile.flatten(), t2_tile.flatten())
    #tlapse = res.slope # degrees m-1
    return ln.coef_

def get_tlapse_tile(t2, topo, lmask, t, j, i, halo=3):
    print t, j, i
    imin = np.max([i - halo, 0])
    jmin = np.max([j - halo, 0])
    imax = np.min([i + halo, Ni - 1])
    jmax = np.min([j + halo, Nj - 1])
    t2_tile = t2[t, jmin:jmax, imin:imax]
    topo_tile = topo[0, jmin:jmax, imin:imax]
    lmask_tile = lmask[0, jmin:jmax, imin:imax]
    #tlapse[t, j, i] = get_tlapse_regression(topo_tile, t2_tile, lmask_tile)
    tlapse_point = get_tlapse_regression(topo_tile, t2_tile, lmask_tile)
    return tlapse_point

ifile = "/home/users/garciam/pruebas/wrfout_d04_20170531T120000Z_20170602T000000Z.nc"
inc = ncdf.Dataset(ifile, "r")
t2 = inc.variables["T2"][:]
topo = inc.variables["HGT"][:]
lmask = inc.variables["LANDMASK"][:].astype("bool")

tp = 15

Nt, Nj, Ni = t2.shape
#tlapse = np.empty([Nt, Nj, Ni])
halo = 3

tlapse = Parallel(n_jobs=8)(delayed(get_tlapse_tile)(t2, topo, lmask, t, j, i, halo=3) for t, j, i in product(xrange(Nt), xrange(Nj), xrange(Ni)))
tlapse = np.array(tlapse).reshape([Nt, Nj, Ni])
tlapse = np.where(np.isnan(tlapse), 0.0098, tlapse)
tlapse = np.where(tlapse > 0.0098, 0.0098, tlapse)
ds = xr.open_dataset(ifile)
da_tlapse = ds.T2.copy()
da_tlapse[:] = tlapse
da_tlapse.to_dataset(name="TLAPSE").to_netcdf("test_tlapse.nc")

inc.close()