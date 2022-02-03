na = [CartesianIndex()]

"""calculate 1st order central finite differencing along x-axis with map factors 
------------
s: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
dx: Int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
ds/dx: array-like, 3d (nx, ny, nz)"""
function ∂x(s::Array{T,3}, dx::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    c = (s[3:end,:,:]-s[1:end-2,:,:])/2/dx
    l = (s[2,:,:]-s[1,:,:])/dx
    r = (s[end,:,:]-s[end-1,:,:])/dx 
    ret = cat(l[na,:,:],c,r[na,:,:], dims=1)
    ret = ret .* mf
end

"""calculate 1st order central finite differencing along y-axis with map factors 
------------
s: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
dy: Int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
ds/dy: array-like, 3d (nx, ny, nz)"""
function ∂y(s::Array{T,3}, dy::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    c = (s[:,3:end,:]-s[:,1:end-2,:])/2/dy
    l = (s[:,2,:]-s[:,1,:])/dy
    r = (s[:,end,:]-s[:,end-1,:])/dy
    ret = cat(l[:,na,:],c,r[:,na,:], dims=2)
    ret = ret .* mf
end

"""calculate 2nd order central finite differencing along x-axis with map factors 
------------
s: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
dx: int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
d2s/dx2: array-like, 3d (nx, ny, nz)"""
function ∂²x(s::Array{T,3}, dx::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    dx2 = dx^2
    c = (s[:,3:end,:]+s[:,1:end-2,:] - 2*s[:,2:end-1,:])/dx2
    l = (s[:,2,:]-s[:,1,:])/dx2
    r = (s[:,end,:]-s[:,end-1,:])/dx2
    ret = cat(l[:,na,:],c,r[:,na,:], dims=2)
    ret = ret .* mf
end

"""calculate 2nd order central finite differencing along y-axis with map factors 
------------
s: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
dy: int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
d2s/dy2: array-like, 3d (nx, ny, nz)"""
function ∂²y(s::Array{T,3}, dy::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    dy2 = dy^2
    c = (s[:,3:end,:]+s[:,1:end-2,:] - 2*s[:,2:end-1,:])/dy2
    l = (s[:,2,:]-s[:,1,:])/dy2
    r = (s[:,end,:]-s[:,end-1,:])/dy2
    ret = cat(l[:,na,:],c,r[:,na,:], dims=2)
    ret = ret .* mf
end

"""calculate the laplacian of a scalar with map factors 
------------
s: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
dx, dy: int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
lap(s): array-like, 3d (nx, ny, nz)"""
function ∇²(s::Array{T,3}, dx::Int, dy::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    dx2 = dx^2
    mf2 = mf.^2
    ret = zeros(Float64, 209, 143, 41)
    ret[2:end-1,2:end-1,:] = (s[2:end-1,3:end,:]+s[2:end-1,1:end-2,:]+s[3:end,2:end-1,:]+s[1:end-2,2:end-1,:]- 4*s[2:end-1,2:end-1,:])/dx2
    ret = ret .* mf2
end

"""calculate gradient of a scalar with map factors 
------------
s: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
dx, dy: int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
dsdx, dvdy: a tuple of 2 arrays, 3d (nx, ny, nz)"""
function ∇(s::Array{T,3}, dx::Int, dy::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    ret = (∂x(s, dx, mf), ∂y(s,dy,mf))
end

"""calculate the divergence of a vector (va, vb) with map factors 
------------
va, vb : array-like, 3d (nx, ny, nz)
    x and y components of a vector, assume on the center grid of C-grid
dx, dy: int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
div(v) : array-like, 3d (nx, ny, nz)"""
function ∇div(va::Array{T,3}, vb::Array{T,3}, dx::Int, dy::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    dudx = ∂x(va, dx, mf)
    dvdy = ∂y(vb, dy, mf)
    return dudx + dvdy
end

"""calculate 2nd order central finite differencing along x then y-axis with map factors 
------------
s: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
dx,, dy: int
    grid spacing
mf: array-like, 3d (nx, ny, None)
    map factor from wrfout
------------
returns
d2s/dxdy: array-like, 3d (nx, ny, nz)"""
function ∂xy(s::Array{T,3}, dx::Int, dy::Int, mf::Array{T,3}) where {T<:AbstractFloat}
    ret = ∂y(∂x(s, dx, mf), dy, mf)
end

function ∂xπ(s::Array{T,3}, dx::Int, h::Array{T,3}, mf::Array{T,3}) where {T<:AbstractFloat}
    ret = ∂π(∂x(s, dx, mf), h)
end

function ∂yπ(s::Array{T,3}, dy::Int, h::Array{T,3}, mf::Array{T,3}) where {T<:AbstractFloat}
    ret = ∂π(∂y(s, dy, mf), h)
end

"""calculate 1st order central finite differencing along vertical coordinate
Sundqvist and Veronis 1970
------------
f: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
h: array-like, 3d (None, None, nz)
    difference in vertical coordinate (Exner), 
------------
returns
ds/dpi: array-like, 3d (nx, ny, nz)"""
function ∂π(s::Array{T,3}, h::Array{T,3}) where {T<:AbstractFloat}
    him1 = h[:,:,1:end-1]
    hi = h[:,:,2:end]  #?
    si = s[:,:,2:end-1]
    sip1 = s[:,:,3:end]
    sim1 = s[:,:,1:end-2]
    c = (him1.^2 .* sip1 - hi.^2 .*sim1) ./ (hi .+ him1) - si .* (him1 .- hi) ./hi ./ him1
    top = (s[:,:,end] - s[:,:,end-1]) ./ h[:,:,end-1]
    btn = (s[:,:,2] - s[:,:,1]) ./ h[:,:,1]
    ret = cat(btn[:,:,na], c, top[:,:,na], dims=3)
end

"""calculate 2nd order central finite differencing along vertical coordinate
Sundqvist and Veronis 1970
------------
f: array-like, 3d (nx, ny, nz)
    assume on the center grid of C-grid
h: array-like, 3d (None, None, nz)
    difference in vertical coordinate (Exner), 
------------
returns
d2s/dpi2: array-like, 3d (nx, ny, nz)"""
function ∂²π(s::Array{T,3}, h::Array{T,3}) where {T<:AbstractFloat}
    him1 = h[:,:,1:end-1]
    hi = h[:,:,2:end]  #?
    si = s[:,:,2:end-1]
    sip1 = s[:,:,3:end]
    sim1 = s[:,:,1:end-2]
    c = ( (him1 .* sip1 - hi .* sim1) ./ (hi .+ him1) .- si) ./ hi ./him1 * 2
    top = (s[:,:,end] - s[:,:,end-1]) ./ h[:,:,end-1].^2
    btn = (s[:,:,2] - s[:,:,1]) ./ h[:,:,1].^2 
    ret = cat(btn[:,:,na], c, top[:,:,na], dims=3)
end