import numpy as np

def laplacian(s, dx, dy, mf):
    dx2 = dx**2
    mf2 = mf**2 
    ret = np.zeros_like(s)
    ret[:,1:-1,1:-1] = (s[:,1:-1,2:] + s[:,1:-1,:-2] + s[:,2:,1:-1] + s[:,:-2,1:-1] - 4*s[:,1:-1,1:-1])/dx2 
    ret = ret * mf2
    return ret

def ddx(s, dx, mf):
    c = (s[:,:,2:] - s[:,:,:-2])/2/dx
    l = (s[:,:,1] - s[:,:,0])/dx
    r = (s[:,:,-1] - s[:,:,-2])/dx
    ret = np.concatenate((l[:,:,None],c,r[:,:,None]), axis=-1)
    return ret*mf

def ddy(s, dy, mf):
    c = (s[:,2:,:] - s[:,:-2,:])/2/dy
    l = (s[:,1,:] - s[:,0,:])/dy
    r = (s[:,-1,:] - s[:,-2,:])/dy
    ret = np.concatenate((l[:,None,:],c,r[:,None,:]), axis=-2)
    return ret*mf

def grad(s, dx, dy, mf):
    return ddx(s, dx, mf), ddy(s, dy, mf)

def div(va, vb, dx, dy, mf):
    dudx = ddx(va, dx, mf)
    dvdy = ddy(vb, dy, mf)
    return  dudx + dvdy

def d2dx2(s, dx, mf):
    dx2 = dx**2
    c = (s[:,:,2:]-s[:,:,:-2] - 2*s[:,:,1:-1])/dx2
    l = (s[:,:,1] - s[:,:,0])/dx2
    r = (s[:,:,-1] - s[:,:,-2])/dx2
    ret = np.concatenate((l[:,:,None],c,r[:,:,None]), axis=-1)
    return ret*mf**2

def d2dy2(s, dy, mf):
    dy2 = dy**2
    c = (s[:,2:,:] + s[:,:-2,:] - 2*s[:,1:-1,:])/dy2
    l = (s[:,1,:] + s[:,0,:])/dy2
    r = (s[:,-1,:] + s[:,-2,:])/dy2
    ret =  np.concatenate((l[:,None,:],c,r[:,None,:]), axis=-2)
    return  ret*mf**2 #???

def ddpi(s, h):
    him1 = h[:-1]
    hi = h[1:]
    si = s[1:-1,:,:]
    sip1 = s[2:,:,:]
    sim1 = s[:-2,:,:]
    c = (him1[:,None,None]**2 * sip1 - hi[:,None,None]**2 * sim1) / (hi[:,None,None]+him1[:,None,None]) - si * (him1[:,None,None] - hi[:,None,None])/hi[:,None,None]/him1[:,None,None]
    top = (s[-1,:,:] - s[-2,:,:])/h[-1]
    btn = (s[1,:,:] - s[0,:,:])/h[0]
    ret =  np.concatenate((btn[None,:,:],c,top[None,:,:]), axis=0)
    return ret

def d2dpi(s, h):
    him1 = h[:-1]
    hi = h[1:]
    si = s[1:-1,:,:]
    sip1 = s[2:,:,:]
    sim1 = s[:-2,:,:]
    c = ( (him1[:,None,None]*sip1 - hi[:,None,None]*sim1) / (hi[:,None,None]+him1[:,None,None]) - si) / hi[:,None,None]/ him1[:,None,None] * 2
    top = (s[-1,:,:] - s[-2,:,:])/h[-1]**2
    btn = (s[1,:,:] - s[0,:,:])/h[0]**2
    ret =  np.concatenate((btn[None,:,:],c,top[None,:,:]), axis=0)
    return ret

def dxdpi(s, dx, h, mf):
    return ddpi(ddx(s, dx, mf), h)

def dydpi(s, dy, h, mf):
    return ddpi(ddy(s, dy, mf), h)


def adj_R_phi(target, psi_bs, dx, dy, h, mf, A, f):
    rphi1 = (f+ laplacian(psi_bs, dx, dy, mf)) * d2dpi(target, h)
    rphi2 = dxdpi(psi_bs, dx, h, mf) * dxdpi(target, dx, h, mf)
    rphi3 = dydpi(psi_bs, dy, h, mf) * dydpi(target, dy, h, mf)
    rphi = A[:,None,None]*(rphi1 - rphi2 - rphi3)
    return rphi

def adj_R_psi(target, phi_bs, dx, dy, h, mf, A):
    rpsi1 = d2dpi(phi_bs, h) * laplacian(target, dx, dy, mf)
    rpsi2 = dxdpi(phi_bs, dx, h, mf) * dxdpi(target, dx, h, mf)
    rpsi3 = dydpi(phi_bs, dy, h, mf) * dydpi(target, dy, h, mf)
    return A[:,None,None] * (rpsi1 - rpsi2 - rpsi3)

def adj_T_phi(target, dx, dy, mf):
    tphi = - laplacian(target, dx, dy, mf)
    return tphi

def adj_T_psi(target, psi_bs, phi_bs, dx, dy, mf ,f):
    dpsidx, dpsidy = grad(psi_bs, dx, dy, mf)
    tpsi1 = div(f*dpsidx, f*dpsidy, dx, dy, mf)
    tpsi2 = d2dx2(psi_bs, dx, mf) * d2dy2(target, dy, mf)
    tpsi3 = d2dy2(phi_bs, dy, mf) * d2dx2(target, dx, mf)
    tpsi4 = 2 * psi_bs
    tpsi = tpsi1 + tpsi2 + tpsi3 - tpsi4

    return tpsi



def varphiP(target, psi_bs, phi_bs, h, dx, dy, mf, A, f):
    varphip1 = adj_T_phi(target, dx, dy, mf) * adj_R_psi(phi_bs, target, dx, dy, h, mf, A, f)
    varphip2 = adj_T_psi(target, psi_bs, phi_bs, dx, dy, mf, f) * adj_R_phi(target, psi_bs, dx, dy, h, mf, A, f)
    varphip = varphip1 - varphip2
    return varphip

def varphiB(target, psi_bs, phi_bs, h, dx, dy, mf, A, f):
    varphib1 = adj_R_phi(target, psi_bs, dx, dy, h, mf, A, f) * adj_T_psi(target, psi_bs, phi_bs, dx, dy, mf ,f)
    varphib2 = adj_R_psi(target, phi_bs, dx, dy, h, mf, A) * adj_T_phi(target, dx, dy, mf)
    varphib = varphib1 - varphib2
    return varphib

def UB(target, psi_bs, phi_bs, h, dx, dy, mf, A, f):
    ub = adj_R_phi(target, psi_bs, dx, dy, h, mf, A, f) - adj_R_psi(target, phi_bs, dx, dy, h, mf, A)
    return ub


def UP(target, psi_bs, phi_bs, dx, dy, mf, f):
    up = adj_T_phi(target, dx, dy, mf) - adj_T_psi(target, psi_bs, phi_bs, dx, dy, mf ,f)
    return up 



def varphiP_complete(target, psi_bs, phi_bs, h, dx, dy, mf, A, f):
    varphip1a = laplacian(A[:,None,None] * d2dpi(phi_bs, h), dx, dy, mf)
    varphip1b = laplacian(A[:,None,None] * dxdpi(phi_bs, dx, h, mf) * dxdpi(target, dx, h, mf), dx, dy, mf)
    varphip1c = laplacian(A[:,None,None] * dydpi(phi_bs, dy, h, mf) * dydpi(target, dy, h, mf), dx, dy, mf)
    varphip1 = varphip1a+varphip1b+varphip1c
    
    dRphi1dx, dRphi1dy = grad(A[:,None,None] *(f+laplacian(psi_bs, dx, dy, mf)) * d2dpi(target, h),
                              dx, dy, mf)
    varphip2_1a = div(f*dRphi1dx, f*dRphi1dy, dx, dy, mf)
    
    dRphi2dx, dRphi2dy = grad(A[:,None,None] * dxdpi(psi_bs, dx, h, mf) * dxdpi(target, dx, h, mf),
                              dx, dy, mf)
    varphip2_1b = div(f*dRphi2dx, f*dRphi2dy, dx, dy, mf)
    
    dRphi3dx, dRphi3dy = grad(A[:,None,None] *dydpi(psi_bs, dy, h, mf) * dydpi(target, dy, h, mf),
                              dx, dy, mf)
    varphip2_1c = div(f*dRphi3dx, f*dRphi3dy, dx, dy, mf)
    
    varphip2_1 = varphip2_1a - varphip2_1b - varphip2_1c
    
    
    varphip2_2a = d2dx2(psi_bs, dx, mf) * d2dy2(A[:,None,None]*(f+laplacian(psi_bs, dx, dy, mf)) * d2dpi(target, h),
                                                dy, mf)
    varphip2_2b = d2dx2(psi_bs, dx, mf) * d2dy2(A[:,None,None]*dxdpi(psi_bs, dx, h, mf) * dxdpi(target, dx, h, mf),
                                                dy, mf)
    varphip2_2c = d2dx2(psi_bs, dx, mf) * d2dy2(A[:,None,None]*dydpi(psi_bs, dy, h, mf) * dydpi(target, dy, h, mf),
                                                dy, mf)
    
    varphip2_2 = varphip2_2a - varphip2_2b - varphip2_2c
    
    
    varphip2_3a = d2dy2(phi_bs, dy, mf) * d2dx2(A[:,None,None]*(f+laplacian(psi_bs, dx, dy, mf)) * d2dpi(target, h),
                                                dx, mf)
    varphip2_3b = d2dy2(phi_bs, dy, mf) * d2dx2(A[:,None,None]*dxdpi(psi_bs, dx, h, mf) * dxdpi(target, dx, h, mf),
                                                dx, mf)
    varphip2_3c = d2dy2(phi_bs, dy, mf) * d2dx2(A[:,None,None]*dydpi(psi_bs, dy, h, mf) * dydpi(target, dy, h, mf),
                                                dx, mf)
    varphip2_3 = varphip2_3a - varphip2_3b - varphip2_3c
    
    varphip2_4 = 2*psi_bs * adj_R_phi(target, psi_bs, dx, dy, h, mf, A, f)
    varphip2 = varphip2_1 + varphip2_2 + varphip2_3 - varphip2_4
    varphip = varphip1 - varphip2
    return varphip