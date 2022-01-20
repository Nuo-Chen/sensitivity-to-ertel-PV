import numpy as np

# https://met.nps.edu/~ldm/gempak/apxB2.pdf for mapfactor reference
def ddx(s, dx, mf):
    """
    calculate 1st order central finite differencing along x-axis with map factors 
    ------------
    s: array-like, 3d (nz,ny,nx)
        assume on the center grid of C-grid
    dx: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    ds/dx: array-like, 3d (nz,ny,nx)
    """
    c = (s[:,:,2:] - s[:,:,:-2])/2/dx
    l = (s[:,:,1] - s[:,:,0])/dx
    r = (s[:,:,-1] - s[:,:,-2])/dx
    
    ret = np.concatenate((l[:,:,None],c,r[:,:,None]), axis=-1)
    return ret*mf[None,:,:]

def ddy(s, dy, mf):
    """
    calculate 1st order central finite differencing along y-axis with map factors 
    ------------
    s: array-like, 3d (nz,ny,nx)
        assume on the center grid of C-grid
    dy: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    ds/dy: array-like, 3d (nz,ny,nx)
    """
    c = (s[:,2:,:] - s[:,:-2,:])/2/dy
    l = (s[:,1,:] - s[:,0,:])/dy
    r = (s[:,-1,:] - s[:,-2,:])/dy
    
    ret = np.concatenate((l[:,None,:],c,r[:,None,:]), axis=-2)
    return ret*mf[None,:,:]

def d2dx2(s, dx, mf):
    """
    calculate 2nd order central finite differencing along x-axis with map factors 
    ------------
    s: array-like, 3d (nz,ny,nx)
        assume on the center grid of C-grid
    dx: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    d2s/dx2: array-like, 3d (nz,ny,nx)
    """
    dx2 = dx**2
    c = (s[:,:,2:] + s[:,:,:-2] - 2*s[:,:,1:-1])/dx2
    l = (s[:,:,1] - s[:,:,0])/dx2
    r = (s[:,:,-1] - s[:,:,-2])/dx2
    
    ret = np.concatenate((l[:,:,None],c,r[:,:,None]), axis=-1)
    return ret*mf[None,:,:]**2

def d2dy2(s, dy, mf):
    """
    calculate 2nd order central finite differencing along y-axis with map factors 
    ------------
    s: array-like, 3d (nz,ny,nx)
        assume on the center grid of C-grid
    dy: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    d2s/dy2: array-like, 3d (nz,ny,nx)
    """
    dy2 = dy**2
    c = (s[:,2:,:] + s[:,:-2,:] - 2*s[:,1:-1,:])/dy2
    l = (s[:,1,:] - s[:,0,:])/dy2
    r = (s[:,-1,:] - s[:,-2,:])/dy2
    
    ret = np.concatenate((l[:,None,:],c,r[:,None,:]), axis=-2)
    return ret*mf[None,:,:]**2

def laplacian(s,dx,dy, mf):
    """
    calculate the laplacian of a scalar with map factors 
    ------------
    s: array-like, 3d (nz,ny,nx)
        assume on the center grid of C-grid
    dx, dy: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    lap(s): array-like, 3d (nz,ny,nx)
    """
    ret = np.zeros(s.shape)
    ret[:,1:-1,1:-1] = (s[:,2:,1:-1:] + s[:,:-2,1:-1] + s[:,1:-1,2:] + s[:,1:-1,:-2] - 4*s[:,1:-1,1:-1])/dx**2
    return ret*mf[None,:,:]**2

def grad(s, dx, dy, mf):
    """
    calculate gradient of a scalar with map factors 
    ------------
    s: array-like, 3d (nz,ny,nx)
        assume on the center grid of C-grid
    dx, dy: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    dsdx, dvdy: a tuple of 2 arrays, 3d (nz,ny,nx)
    """
    ret = (ddx(s, dx, mf), ddy(s, dy, mf))
    return ret 

def div(va, vb, dx, dy, mf):
    """
    calculate the divergence of a vector (va, vb) with map factors 
    ------------
    va, vb : array-like, 3d (nz,ny,nx)
        x and y components of a vector, assume on the center grid of C-grid
    dx, dy: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    div(v) : array-like, 3d (nz,ny,nx)
    """
    dudx = ddx(va, dx, mf)
    dvdy = ddy(vb, dy, mf)
    return dudx+dvdy

def dxdy(s, dx, dy, mf):
    """
    calculate 1st order central finite differencing along y-axis with map factors 
    ------------
    s: array-like, 3d (nz,ny,nx)
        assume on the center grid of C-grid
    dx,, dy: int
        grid spacing
    mf: array-like, 2d (ny,nx)
        map factor from wrfout
    ------------
    returns
    d2s/dxdy: array-like, 3d (nz,ny,nx)
    """
    ret = ddy(ddx(s, dx, mf), dy, mf)
    return ret 