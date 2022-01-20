from psi_SOR_multigrid import *
from precision import *

import numpy as np
import xarray as xr
from wrf import getvar, destagger
from netCDF4 import Dataset
import metpy.constants as mpconst
import time

start_time = time.time()  

file = 'march_2020_wrfout.nc'
ds = xr.open_dataset(file)
wrfout = Dataset(file)
itime = 0

lats = hexize(ds.XLAT[itime,:])
rdx = 1/hexize_val(ds.RDX[itime])
rdy = 1/hexize_val(ds.RDY[itime])
msfm = hexize(ds.MAPFAC_M[itime,:])
msfu = hexize(ds.MAPFAC_U[itime,:])
msfv = hexize(ds.MAPFAC_V[itime,:])
a = int(mpconst.Re.magnitude)
f = 2*7.292e-5*np.sin(np.deg2rad(lats))
f0 = f.mean()
g = mpconst.earth_gravity.magnitude

(nt,) = ds.RDX.shape

streamfunc = np.zeros(ds.T.shape, dtype=np.float32)

for itime in range(nt):
    ### set initial precision to float64 so that streamfunction can be closed
    prc = np.float64

    u = ds.U.values[itime,:]
    v = ds.V.values[itime,:]

    # destagger for now, todo: work with staggered data
    u = destagger(u, stagger_dim=-1).astype(prc)
    v = destagger(v, stagger_dim=-2).astype(prc)

    Z = getvar(wrfout, "z", timeidx=itime)   # geopotential height
    nz,ny,nx = Z.shape

    psi = Z*g/f0  
    psi = psi_BC(psi, u, v, rdx, rdy, msfm, msfu, msfv)

    prc = np.float32

    psi = psi.astype(prc)
    avo = getvar(wrfout, 'avo', timeidx=itime, meta=False).astype(prc)
    rvo = avo*1e-5 - f[None,:,:].astype(prc)
    streamfunc = np.zeros((nz,ny,nx), dtype=prc)
    rdx = rdx.astype(np.int16)
    msfm = msfm.astype(prc)

    for k in range(nz):
        # use multigrid to quickly converge

        start_time_lv = time.time()
        itr = 1
        trypsi, eres_r = multigrid(psi[k,:], rvo[k,:], rdx, mf=msfm)
        peres = 1e10
        eres = abs(eres_r).mean()
        print(itr, eres)
        l_err = []
        while abs(peres - eres)/eres > 0.01 and itr<350:
            peres = eres
            trypsi, eres_r = multigrid(trypsi, rvo[k,:], rdx, mf=msfm)
            eres = abs(eres_r).mean()
            l_err.append(eres)
            itr += 1
            if itr%50==0:
                print(itr, eres)

        # use SOR to fine tune the boudary condition
        #   
        psi_sor = psi_BC(psi, u, v, rdx, rdy, msfm, msfu, msfv)[k,:]
        psi_sor[1:-1,1:-1] = trypsi[1:-1,1:-1]
        err_rec_list = []
        omg = 1.8
        while omg>=1:
            psi_sor, err_record, niter = vorticity_inversion(rvo[k,:], psi_sor, rdx, msfm, thrs=1, omega=omg)
            print(omg)
            if omg>=1.5 or (omg - 1)<1e-3:
                omg -= 0.1
            else:
                omg = 1
            err_rec_list.append(err_record[1:])
        
        streamfunc[itime, k, :] = psi_sor

np.save("mar2020_psi.npy", streamfunc)
print(f"-------------total time used----------- for a ({nt,nz,ny,nx}) file")
print("--- %s seconds ---" % (time.time() - start_time))