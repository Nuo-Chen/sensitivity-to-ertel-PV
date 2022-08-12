using Pkg
using NPZ#, NCDatasets
using JuMP, Ipopt

function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

grad = ingredients("/Volumes/Scratch/nchen67/ertel_inv/gradient1.jl")
adjt = ingredients("/Volumes/Scratch/nchen67/ertel_inv/adjoint1.jl")

path = "/Volumes/Scratch/nchen67/ertel_inv/"

ψ_bs = npzread(path*"psi_bs20lev.npy")
ϕ_bs = npzread(path*"phi_bs20lev.npy")
frompython = npzread(path*"forjulia.npz")
rdx = Float64.(frompython["rdx"])
rdy = Float64.(frompython["rdy"])
A = frompython["A"]
msfm = frompython["msfm"]
dπ = frompython["dexner"]
f = frompython["f"]
ψ̂ = npzread(path*"apsi20lev.npy");
guess = npzread(path*"2021041515z_inverted_20lev.npz")
initguess = reshape(guess["aq"], (20,143,209))
coeffext = npzread(path*"coeffext.npy")


ψ_bs_xyz = Float64.(permutedims(ψ_bs,[3,2,1]))
ϕ_bs_xyz = permutedims(ϕ_bs,[3,2,1])
coeffext_xyz  = permutedims(coeffext,[3,2,1])
msfm_xyz = repeat(permutedims(msfm,[2,1]), outer = [1,1, 20])
f_xyz = repeat(permutedims(f,[2,1]), outer = [1,1, 20])
ψ̂_xyz = permutedims(ψ̂,[3,2,1])	
initguess_xyz = permutedims(initguess,[3,2,1])
A_xyz = permutedims(repeat(A, outer=[1,209,143]), [2, 3, 1])
dπ_xyz = permutedims(repeat(dπ, outer=[1,209,143]), [2, 3, 1])
    

msfm_xy = permutedims(msfm,[2,1]); 
msfm2=msfm_xy.^2; 
f_xy = permutedims(f,[2,1]); 
rdx2 = rdx^2 ; rdy2 = rdy^2;

rhs = adjt.UP(ψ̂_xyz, ψ_bs_xyz, ϕ_bs_xyz, rdx, rdy, msfm_xyz, f_xyz);

nx, ny, nzz = size(ϕ_bs_xyz)
nz = nzz+2
nxy = nx*ny

rhs_ext = Array{Float64, 3}(undef, nx,ny,nz)
rhs_ext[:,:,2:nz-1] = rhs[:,:,1:nz-2]
rhs_ext[:,:,1] = rhs_ext[:,:,2]
rhs_ext[:,:,end] = rhs_ext[:,:,end-1]

ψ_bs_ext = Array{Float64, 3}(undef, nx,ny,nz)
ψ_bs_ext[:,:,2:nz-1] = ψ_bs_xyz[:,:,1:nz-2]
ψ_bs_ext[:,:,1] = ψ_bs_ext[:,:,2]
ψ_bs_ext[:,:,end] = ψ_bs_ext[:,:,end-1]

ϕ_bs_ext = Array{Float64, 3}(undef, nx,ny,nz)
ϕ_bs_ext[:,:,2:nz-1] = ϕ_bs_xyz[:,:,1:nz-2]
ϕ_bs_ext[:,:,1] = ϕ_bs_ext[:,:,2]
ϕ_bs_ext[:,:,end] = ϕ_bs_ext[:,:,end-1]

dπ_ext = Array{Float64, 1}(undef, nz-1)
dπ_ext[2:end-1] = dπ
dπ_ext[1] = -14
dπ_ext[end] = -100

initguess_ext = Array{Float64, 3}(undef, nx,ny,nz)
initguess_ext[:,:,2:nz-1] = initguess_xyz[:,:,1:nz-2]
initguess_ext[:,:,1] .= 0
initguess_ext[:,:,end] .= 0

A_ext = Array{Float64, 1}(undef, nz)
A_ext[2:nz-1] = A[1:nz-2]
A_ext[1] = 0.027
A_ext[end] = 0.25

coeffext_vec = vec(coeffext_xyz)
initguess_vec = vec(initguess_ext)
replace!(coeffext_vec, NaN=>0.01)
replace!(coeffext_vec, Inf=>9e10)
replace!(coeffext_vec, -Inf=>-9e10)

md = Model(Ipopt.Optimizer)
set_optimizer_attribute(md, "tol", 1e1)
@variable(md, a[i=1:nx*ny*nz], start=initguess_vec[i])
@variable(md, z)

@constraint(md, [k=1:nz], a[nxy*k-nx+1:nxy*k] .== a[nxy*k-2*nx+1:nxy*k-nx])
@constraint(md, [k=0:nz-1], a[nxy*k+1:nxy*k+nx] .== a[1+nx+nxy*k:2*nx+nxy*k])
@constraint(md, [k=0:nz-1], a[nxy*k+1:nx:nxy*(k+1)-nx+1] .== a[nxy*k+2:nx:nxy*(k+1)-nx+2])
@constraint(md, [k=0:nz-1], a[nxy*k+nx:nx:nxy*(k+1)] .== a[nxy*k+nx-1:nx:nxy*(k+1)-1])
@constraint(md, a[1:nxy] .== a[nxy+1:2*nxy])
@constraint(md, a[nxy*(nz-1)+1:nxy*nz] .== a[nxy*(nz-2)+1:nxy*(nz-1)])

myexp = @expression(md, 1/nx/ny/nz*sum(vec(rhs_ext) .- coeffext_vec.*a))

@constraint(md, z >= myexp)
@constraint(md, z >= -myexp)

@objective(md, Min, z)	
optimize!(md)

a_fin = value.(a)
npzwrite("matrixsen2ertelpv_fulllev.npy", a_fin)