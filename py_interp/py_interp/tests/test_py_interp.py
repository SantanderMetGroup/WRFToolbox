import py_interp.mainfun as pm

ifile = "/home/users/garciam/pruebas/py_interp/wrfout_d03_20161001120000.nc"
varlist = ["T2", "TSK", "Q2","PSFC", "U10", "V10", "CLT", "TLAPSE200"]
vars_to_filter = ["CLT",]
plevs = [1000, 925, 850, 700, 500, 300]
pm.py_interp_file(ifile, varlist, plevs, vars_to_filter, overwrite=False,
                  verbose=True)