function ∂x(s, dx, mf, i, j, k) 
ret = (s[i+1,j,k] - s[i-1,j,k])/2/dx * mf[i,j]
end

function ∂y(s, dy, mf, i, j, k) 
    ret = (s[i,j+1,k] - s[i,j-1,k])/2/dy * mf[i,j]
end

function ∇(s, dx, dy, mf, i, j, k) 
    ret = (∂x(s, dx, mf, i, j, k), ∂y(s,dy,mf, i, j, k))
end

function ∇div(va, vb, dx, dy, mf, i, j, k) 
    dudx = ∂x(va, dx, mf, i, j, k)
    dvdy = ∂y(vb, dy, mf, i, j, k)
    return dudx + dvdy
end

function ∂π(s, h, i, j, k) 
    ret = (h[k-1]^2 * s[i,j,k+1] - h[k]^2 * s[i,j,k-1]) / (h[k]*h[k-1]* (h[k] + h[k-1])) - s[i,j,k]*(h[k-1]-h[k])/(h[k]*h[k-1])
end

# 2nd order derivatives

function ∂²x(s, dx2, mf2, i, j, k) 
    ret = (s[i+1,j,k] + s[i-1,j,k] - 2*s[i,j,k])/dx2 * mf2[i,j]
end

function ∂²y(s, dy2, mf2, i, j, k) 
    ret = (s[i,j+1,k] + s[i,j-1,k] - 2*s[i,j,k])/dy2 * mf2[i,j]
end

function ∇²(s, dx2, dy2, mf2, i, j, k) 
    ret = (s[i,j+1,k]+s[i,j-1,k]+s[i+1,j,k]+s[i-1,j,k]- 4*s[i,j,k])/dx2 * mf2[i,j]
end

function ∂²π(s, h, i, j, k) 
    ret = 2/h[k]/h[k-1]*( (h[k-1]*s[i,j,k+1] + h[k]*s[i,j,k-1]) / (h[k]+h[k-1]) - s[i,j,k])
end

function diff2nd(ap1, am1, da, mfij)
    ret = (ap1 - am1) /2/da * mfij 
end

function ∂xy(s, dx, dy, mf2, i, j, k) 
    ret = mf2[i,j]/(4*dx*dy) * ( s[i+1,j+1,k] + s[i-1,j-1,k] - s[i+1,j-1,k] - s[i-1,j+1,k])
end

function ∂xπ(s, dx, h, mf, i, j, k) 
    dπjm1 = ∂π(s, h, i-1, j, k) 
    dπjp1 = ∂π(s, h, i+1, j, k) 
    ret = diff2nd(dπjp1, dπjm1, dx, mf[i,j])
end

function ∂yπ(s, dy, h, mf, i, j, k) 
    dπjm1 = ∂π(s, h, i, j-1, k) 
    dπjp1 = ∂π(s, h, i, j+1, k) 
    ret = diff2nd(dπjp1, dπjm1, dy, mf[i,j])
end

function ∇f∇dis(s, dx, dy, mf, mf2, f, fx, fy, i, j, k) 
sx , sy = ∇(s, dx, dy, mf, i, j, k) 
∇f∇t = fx[i,j]*sx + fy[i,j]*sy + f[i,j]*∇²(s, dx^2, dy^2, mf2, i, j, k)
end

# 3rd order derivatives

"""i = 2:nx-1, j = 2:ny-1
∂³x(s, dx3, mf3, i, j, k)"""
function ∂³x(s, dx3, mf3, i, j, k)
    ret = (s[i+2,j,k] - 2*s[i+1,j,k] + 2*s[i-1,j,k] - s[i-2,j,k])/2/dx3 * mf3[i,j]
end

"""i = 2:nx-1, j = 2:ny-1
∂³y(s, dy3, mf3, i, j, k)"""
function ∂³y(s, dy3, mf3, i, j, k)
    ret = (s[i,j+2,k] - 2*s[i,j+1,k] + 2*s[i,j-1,k] - s[i,j-2,k])/2/dy3 * mf3[i,j]
end

function ∂³xy2(s, dx, dy2, mf3, i, j, k)
    ret = (s[i+1,j+1,k] + s[i+1,j-1,k] - 2*s[i+1,j,k] - s[i-1,j+1,k] - s[i-1,j-1,k] + 2*s[i-1,j,k])/2/dy2/dx * mf3[i,j]
end

function ∂³x2y(s, dx2, dy, mf3, i, j, k)
    ret = (s[i+1,j+1,k] + s[i-1,j+1,k] - 2*s[i,j+1,k] - s[i+1,j-1,k] - s[i-1,j-1,k] + 2*s[i,j-1,k])/2/dy/dx2 * mf3[i,j]
end

# 4th order derivatives

"""i = 2:nx-1, j = 2:ny-1
∇⁴(s, dx2, dy2, mf2, i, j, k) """
function ∇⁴(s, dx2, dy2, mf2, i, j, k) 
    ret = (s[i,j+2,k]+s[i,j-2,k]+s[i+2,j,k]+s[i-2,j,k] - 8*(s[i,j+1,k]+s[i,j-1,k]+s[i+1,j,k]+s[i-1,j,k]) + 2*(s[i+1,j+1,k]+s[i-1,j+1,k]+s[i+1,j-1,k]+s[i-1,j-1,k]) + 20*s[i,j,k])/(dx2^2) * mf2[i,j]^2
end

function ∂⁴x3π(s, dx3, h, mf3, i, j, k)
    dx3kp1 =  ∂³x(s, dx3, mf3, i, j, k+1)
    dx3k =  ∂³x(s, dx3, mf3, i, j, k)
    dx3km1 =  ∂³x(s, dx3, mf3, i, j, k-1)
    ret = (h[k-1]^2 * dx3kp1 - h[k]^2 * dx3km1) / (h[k]*h[k-1]* (h[k] + h[k-1])) - dx3k*(h[k-1]-h[k])/(h[k]*h[k-1])
end

function ∂⁴y3π(s, dy3, h, mf3, i, j, k)
    dy3kp1 =  ∂³y(s, dy3, mf3, i, j, k+1)
    dy3k =  ∂³y(s, dy3, mf3, i, j, k)
    dy3km1 =  ∂³y(s, dy3, mf3, i, j, k-1)
    ret = (h[k-1]^2 * dy3kp1 - h[k]^2 * dy3km1) / (h[k]*h[k-1]* (h[k] + h[k-1])) - dy3k*(h[k-1]-h[k])/(h[k]*h[k-1])
end

function ∂⁴x2π2(s, dx2, h, mf2, i, j, k)
    dx2kp1 = ∂²x(s, dx2, mf2, i, j, k+1) 
    dx2k = ∂²x(s, dx2, mf2, i, j, k) 
    dx2km1 = ∂²x(s, dx2, mf2, i, j, k-1) 
    ret = 2/h[k]/h[k-1]*( (h[k-1]*dx2kp1 + h[k]*dx2km1) / (h[k]+h[k-1]) - dx2k)
end

function ∂⁴y2π2(s, dy2, h, mf2, i, j, k)
    dy2kp1 = ∂²y(s, dy2, mf2, i, j, k+1) 
    dy2k = ∂²y(s, dy2, mf2, i, j, k) 
    dy2km1 = ∂²y(s, dy2, mf2, i, j, k-1) 
    ret = 2/h[k]/h[k-1]*( (h[k-1]*dy2kp1 + h[k]*dy2km1) / (h[k]+h[k-1]) - dy2k)
end

function ∂⁴xy2π(s, dx, dy2, h, mf3, i, j, k)
    dxy2kp1 = ∂³xy2(s, dx, dy2, mf3, i, j, k+1)
    dxy2k = ∂³xy2(s, dx, dy2, mf3, i, j, k)
    dxy2km1 = ∂³xy2(s, dx, dy2, mf3, i, j, k-1)
    ret = 2/h[k]/h[k-1]*( (h[k-1]*dxy2kp1 + h[k]*dxy2km1) / (h[k]+h[k-1]) - dxy2k)
end

function ∂⁴x2yπ(s, dx2, dy, h, mf3, i, j, k)
    dx2ykp1 = ∂³x2y(s, dx2, dy, mf3, i, j, k+1)
    dx2yk = ∂³x2y(s, dx2, dy, mf3, i, j, k)
    dx2ykm1 = ∂³x2y(s, dx2, dy, mf3, i, j, k-1)
    ret = 2/h[k]/h[k-1]*( (h[k-1]*dx2ykp1 + h[k]*dx2ykm1) / (h[k]+h[k-1]) - dx2yk)
end

"""i = 2:nx-1, j = 2:ny-1
∇f∇∂ππ(s, dx, dy, h, mf, mf2, f, fx, fy, i, j, k)"""
function ∇f∇∂ππ(s, dx, dx2, dy, h, mf, mf2, f, fx, fy, i, j, k) 
    pij = ∂²π(s, h, i, j, k) 
    pijp1 = ∂²π(s, h, i, j+1, k) 
    pijm1 = ∂²π(s, h, i, j-1, k) 
    pip1j = ∂²π(s, h, i+1, j, k) 
    pim1j = ∂²π(s, h, i-1, j, k) 
    
    px = (pip1j - pim1j)/2/dx * mf[i,j]
    py = (pijp1 - pijm1)/2/dy * mf[i,j]
    
    delp = (pijp1+pijm1+pip1j+pim1j-4*pij)/dx2 * mf2[i,j]
    
    ret = fx[i,j]*px + fy[i,j]*py + f[i,j]*delp
end

"""i = 3:nx-2, j = 3:ny-2
∇f∇∂xπ(s, dx, dy, h, mf, mf2, f, fx, fy, i, j, k) """
function ∇f∇∂xπ(s, dx, dx2, dy, h, mf, mf2, f, fx, fy, i, j, k) 
    pij = ∂xπ(s, dx, h, mf, i, j, k) 
    pijp1 = ∂xπ(s, dx, h, mf, i, j+1, k) 
    pijm1 = ∂xπ(s, dx, h, mf, i, j-1, k) 
    pip1j = ∂xπ(s, dx, h, mf, i+1, j, k) 
    pim1j = ∂xπ(s, dx, h, mf, i-1, j, k) 
    
    px = (pip1j - pim1j)/2/dx * mf[i,j]
    py = (pijp1 - pijm1)/2/dy * mf[i,j]
    
    delp = (pijp1+pijm1+pip1j+pim1j-4*pij)/dx2 * mf2[i,j]
    
    delfdelxπ = fx[i,j]*px + fy[i,j]*py + f[i,j]*delp
    return delfdelxπ, delp
end

"""i = 3:nx-2, j = 3:ny-2
∇f∇∂yπ(s, dx, dy, h, mf, mf2, f, fx, fy, i, j, k)"""
function ∇f∇∂yπ(s, dx, dy, dy2, h, mf, mf2, f, fx, fy, i, j, k) 
    pij = ∂yπ(s, dy, h, mf, i, j, k) 
    pijm1 = ∂yπ(s, dy, h, mf, i, j-1, k) 
    pijp1 = ∂yπ(s, dy, h, mf, i, j+1, k) 
    pip1j = ∂yπ(s, dy, h, mf, i+1, j, k) 
    pim1j = ∂yπ(s, dy, h, mf, i-1, j, k) 
    
    px = (pip1j - pim1j)/2/dx * mf[i,j]
    py = (pijp1 - pijm1)/2/dy * mf[i,j]
    
    delp = (pijp1+pijm1+pip1j+pim1j-4*pij)/dy2 * mf2[i,j]
    
    delfdelyπ = fx[i,j]*px + fy[i,j]*py + f[i,j]*delp
    return delfdelyπ, delp
end