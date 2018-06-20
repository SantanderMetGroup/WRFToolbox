import py_interp.mainfun as pm

ifile = "/home/users/garciam/pruebas/wrfout_d04_20170531T120000Z_20170602T000000Z.nc"
varlist = ["T2", "TSK", "Q2","PSFC", "U10", "V10", "CLT", "TLAPSEV", "TLAPSER"]
plevs = [1000, 925, 850, 700, 500, 300]
pm.py_interp_file(ifile, varlist, plevs, overwrite=False, verbose=True)