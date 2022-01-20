import numpy as np
import time

def psi_BC(psi, u, v, rdx, rdy, msfm, msfu, msfv, prc=np.float64):
    nz, ny, nx = psi.shape
    dsum = np.zeros(nz, dtype=prc)
    circum = 0
    for i in np.arange(nx-1):
        dsum -= rdx * (v[:,0,i]+v[:,0,i+1]) / 2 /msfu[0,i+1]    # southern boundary, outward, minus v
        circum += msfu[0,i+1] * rdx
        dsum += rdx * (v[:,-1,i]+v[:,-1,i+1]) / 2 /msfu[-1,i+1] # northern boundary, outward, plus v
        circum += msfu[-1,i+1] * rdx

    for j in np.arange(ny-1):
        dsum += rdy * (u[:,j,-1]+u[:,j+1,-1]) / 2 / msfv[j+1,-1]            # eastern boundary, outward, plus u 
        circum += msfv[j+1,-1] * rdy
        dsum -= rdy * (u[:,j,0]+u[:,j+1,0]) / 2 / msfv[j+1,0]            # western boundary, outward, minus u
        circum += msfv[j+1,0] * rdy

    avg_dsum = dsum / circum
    
    print('st:', psi[28,0,0]) # sanity check
    for i in range(1,nx):
        psi[:,0,i] =  psi[:,0,i-1] + rdx/msfu[0,i]*(v[:,0,i]+v[:,0,i-1])/2 + rdx*msfu[0,i]*avg_dsum[:]
    print(psi[28,0,-1])
    
    for j in range(1,ny):
        # start from lower right, going upward, eastern boundary
        psi[:,j,-1]= psi[:,j-1,-1] - rdy/msfv[j,-1]*(u[:,j,-1]+u[:,j-1,-1])/2 + rdy*msfv[j,-1]*avg_dsum[:]
    print(psi[28,-1,-1])
    
    for i in range(nx-2,-1,-1):
        # start from upper right, going westward, northern boundary
        psi[:,-1,i] = psi[:,-1,i+1] - rdx/msfu[-1,i+1]*(v[:,-1,i]+v[:,-1,i+1])/2 + rdx*msfu[-1,i+1]*avg_dsum[:]
    print(psi[28,-1,0])
    
    for j in range(ny-2,-1,-1):
        # start from upper left, going southward, western boundary
        psi[:,j,0] = psi[:,j+1,0] + rdy/msfv[j+1,0]*(u[:,j,0]+u[:,j+1,0])/2 + rdy*msfv[j+1,0]*avg_dsum[:]
    print('ed:', psi[28,0,0])
    
    return psi

def del2(psik, dx, prc=np.float32):
    # nyy, nxx = vor.shape
    lap = np.zeros(psik.shape, dtype=prc)
    dy = dx
    
    lap[1:-1,1:-1] = ((psik[:-2,1:-1] + psik[2:,1:-1] - 2*psik[1:-1,1:-1]) / dy**2 
                     + (psik[1:-1,:-2] + psik[1:-1,2:] - 2*psik[1:-1,1:-1]) / dx**2)
    
    return lap

def jacobi(psik, vor, dx, prc=np.float32):
    psikp1  = psik
    dy = dx
    nyy, nxx = psik.shape
    
    for j in range(1,nyy-1):
        for i in range(1,nxx-1):
            la = (psik[j,i+1] + psik[j,i-1])/dx**2 + (psik[j+1,i]+psik[j-1,i])/dy**2
            psikp1[j,i] = (vor[j,i] - la) / (-2/dx**2 - 2/dy**2)
        
    return psikp1

def multigrid(psi0, vor, dx, mf, prc=np.float32):
    for i in range(5):
        psi0 = psi0*0.333 + jacobi(psi0, vor, dx)*0.667
        
    res = vor/mf**2 - del2(psi0, dx)
    res = res*1.8
    # interpolate from fine to coarse
    res_coarse = res[::2,::2]
    err_coarse = np.zeros(res_coarse.shape, dtype=prc)
    mf_coarse = mf[::2,::2]

    nyy, nxx = vor.shape
    if min(nyy, nxx) > 3:
        err_coarse, _ = multigrid(err_coarse, res_coarse, dx, mf_coarse)

    # interpolation coarse to fine
    err = np.zeros(res.shape, dtype=prc)
    err[::2,::2] = err_coarse
    try:
        err[1::2, ::2] = (err[:-2:2, ::2] + err[2::2,::2])/2
    except ValueError:
        err[1:-1:2, ::2] = (err[:-2:2, ::2] + err[2::2,::2])/2
    try:
        err[:,1::2] = (err[:,:-2:2] + err[:,2::2])/2
    except ValueError:
        err[:,1:-1:2] = (err[:,:-2:2] + err[:,2::2])/2
    psik = psi0 + err
    
    return psik, err

def vorticity_inversion(vor, og, rdx, msfm, thrs, max_iter=5e3, omega=1.8, prc=np.float32):
    start_time = time.time()
    
    error_list = [1e10]
    unit_err = 999
    nitr = 0
    fny, fnx = vor.shape
    res = np.zeros((fny,fnx), dtype=np.float32)

    A1,A2,A3,A4,A5 = (1,1,-4,1,1)
    msfm2 = msfm**2
    msfm2 = msfm2.astype(np.float32)
    
    while unit_err>thrs:
        error = 0
        for j in range(1, fny-1):
            for i in range(1,fnx-1):
                lap = (A1 * og[j-1,i]
                     + A2 * og[j,i-1]
                     + A3 * og[j,i]
                     + A4 * og[j,i+1]
                     + A5 * og[j+1,i]) 
                
                res[j,i] = lap - vor[j,i] * (rdx*rdx)/(msfm2[j,i])
                og[j,i] = og[j,i] - omega*res[j,i]/A3 #1e9

        error = np.sum(abs(res))
        unit_err = error/(fnx-2)/(fny-2)
        err_dff = (unit_err - error_list[-1])/ unit_err
        error_list.append(unit_err)
        nitr += 1
        if nitr>max_iter or err_dff>0.01:
            break
        if omega==1 and err_dff>-5e-3:
            break
        
    print("--- %s seconds ---" % (time.time() - start_time))
    print('got here: ',unit_err, nitr)
    return og, error_list, nitr