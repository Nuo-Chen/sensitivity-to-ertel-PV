∇² = zeros(Float32, nxyz, nxyz)
for k in 1:nz, j in 3:ny-2, i in 3:nx-2
    tijk = ((k-1)*ny +j-1)*nx + i
    ∇² = logop(tijk, i, j-1, k, 1., ∇²)
    ∇² = logop(tijk, i-1, j, k, 1., ∇²)
    ∇² = logop(tijk, i, j, k, -4., ∇²)
    ∇² = logop(tijk, i, j+1, k, 1., ∇²)
    ∇² = logop(tijk, i+1, j, k, 1., ∇²)
end

∇⁴ = zeros(Float32, nxyz, nxyz)
for k in 1:nz, j in 3:ny-2, i in 3:nx-2
    tijk = ((k-1)*ny +j-1)*nx + i
    ∇⁴ = logop(tijk, i+2, j, k, 1., ∇⁴)
    ∇⁴ = logop(tijk, i-2, j, k, 1., ∇⁴)
    ∇⁴ = logop(tijk, i, j+2, k, 1., ∇⁴)
    ∇⁴ = logop(tijk, i, j-2, k, 1., ∇⁴)
    ∇⁴ = logop(tijk, i+1, j, k, -8., ∇⁴)
    ∇⁴ = logop(tijk, i-1, j, k, -8., ∇⁴)
    ∇⁴ = logop(tijk, i, j+1, k, -8., ∇⁴)
    ∇⁴ = logop(tijk, i, j-1, k, -8., ∇⁴)
    ∇⁴ = logop(tijk, i+1, j+1, k, 2., ∇⁴)
    ∇⁴ = logop(tijk, i-1, j-1, k, 2., ∇⁴)
    ∇⁴ = logop(tijk, i+1, j-1, k, 2., ∇⁴)
    ∇⁴ = logop(tijk, i-1, j+1, k, 2., ∇⁴)
    ∇⁴ = logop(tijk, i, j, k, 20., ∇⁴)
end

∂π = zeros(Float32, nxyz, nxyz)
h = dπ
for k in 2:nz-1, j in 1:ny, i in 1:nx
    tijk = ((k-1)*ny +j-1)*nx + i
    ∂π = logop(tijk, i, j, k+1, 1/(h[k]*(1+h[k]/h[k-1])), ∂π)
    ∂π = logop(tijk, i, j, k-1, -h[k]/(h[k-1]^2+h[k]*h[k-1]), ∂π)
    ∂π = logop(tijk, i, j, k, 1/h[k-1]-1/h[k], ∂π)
end

∂²π = zeros(Float32, nxyz, nxyz)
for k in 2:nz-1, j in 1:ny, i in 1:nx
    tijk = ((k-1)*ny +j-1)*nx + i
    ∂²π = logop(tijk, i, j, k+1, 2/(h[k]*(h[k]+h[k-1])), ∂²π)
    ∂²π = logop(tijk, i, j, k-1, 2/(h[k-1]*(h[k]+h[k-1])), ∂²π)
    ∂²π = logop(tijk, i, j, k, -2/h[k]/h[k-1], ∂²π)
end

xx = repeat(collect(Int32, 1:nx)/nx, outer=(1,ny))
yy = repeat(transpose(collect(Int32, 1:ny)/ny), outer=(nx,1))
dst = sqrt.((xx.-0.55).^2 + (yy.-0.55).^2)
gauss = repeat(exp.(-(dst.^2)/2), outer=(1,1,nz))
for i in 1:nx, j in 1:ny
    gauss[i,j,:] .*= sin.(collect(1:20))
end
supressout1 = 0

∇²gauss = reshape(∇²*vec(gauss), (nx, ny,nz))
∇⁴gauss = reshape(∇⁴*vec(gauss), (nx, ny,nz))
gaussₖ = reshape(∂π*vec(gauss), (nx, ny,nz))
gaussₖₖ = reshape(∂²π*vec(gauss), (nx, ny,nz))
supressout2 = 0