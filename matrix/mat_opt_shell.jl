using NPZ
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

adjm = ingredients("/Volumes/Scratch/nchen67/ertel_inv/adjoint1.jl")

path = "/Volumes/Scratch/nchen67/ertel_inv/"
frompython = npzread(path*"forjulia.npz")
f_yx = frompython["f"]
rdx = Float64.(frompython["rdx"])
rdy = Float64.(frompython["rdy"])
msfm_yx = frompython["msfm"]

ψ_bs_zyx = npzread(path*"psi_bs20lev.npy")
ϕ_bs_zyx = npzread(path*"phi_bs20lev.npy")
ψ̂_zyx = npzread(path*"apsi20lev.npy")

guess = npzread(path*"2021041515z_inverted_20lev.npz")
initguess_zyx = reshape(guess["aq"], (20,143,209))
φP = npzread("./varphiP_mat.npy")

String("load original data")

ψ_bs = Float64.(permutedims(ψ_bs_zyx,[3,2,1]))
ϕ_bs = permutedims(ϕ_bs_zyx,[3,2,1])
ψ̂ = permutedims(ψ̂_zyx,[3,2,1])	
msfm = permutedims(msfm_yx,[2,1]); 
f = permutedims(f_yx,[2,1]); 

initguess = permutedims(initguess_zyx,[3,2,1])

String("permute dimensions from kji to ijk")

itv = 4
idx = [1,3, 11, 17, 20]

slicex = 45:itv:130
slicey = 40:itv:110
nzi = 5
nzs = 5+2

ψ5 = ψ_bs[slicex, slicey, idx]
ϕ5 = ϕ_bs[slicex, slicey, idx]
ψ̂5 = ψ̂[slicex, slicey, idx]
msfm5 = repeat(msfm[slicex, slicey], outer = [1,1,nzi])
f5 = repeat(f[slicex, slicey], outer = [1,1,nzi])

ig5 = initguess[slicex, slicey, [2,3, 11, 17, 19]]

nxs, nys, nzthrowaway = size(ψ5)
String("Extract necessary layers")

rhs = adjm.UP(ψ̂5, ψ5, ϕ5, rdx, rdy, msfm5, f5);


rhs_ext = Array{Float64, 3}(undef, nxs, nys, nzs)
rhs_ext[:,:,2:end-1] = rhs
rhs_ext[:,:,1] = rhs_ext[:,:,2] .+ (rand(Float64, (nxs, nys)) .- 0.5)*2e6
rhs_ext[:,:,end] = rhs_ext[:,:,end-1] .+ (rand(Float64, (nxs, nys)) .-0.5)*2e7

initguess_ext = Array{Float64, 3}(undef, nxs, nys, nzs)
initguess_ext[:,:,2:end-1] = ig5
initguess_ext[:,:,1] .= ig5[:,:,1] .+ (rand(Float64, (nxs, nys)) .- 0.5)*50000
initguess_ext[:,:,end] .= ig5[:,:,end] .+ (rand(Float64, (nxs, nys)) .- 0.5)*500
String("extend rhs and initial with buffer layers")


md = Model(Ipopt.Optimizer)
@variable(md, a[i=1:nxs, j=1:nys, k=1:nzs], start=initguess_ext[i,j,k])
# @constraint(md, [i=1:nxs, k=1:nzs], a[i,1,k] == a[i,2,k])
# @constraint(md, [i=1:nxs, k=1:nzs], a[i,nys,k] == a[i,nys-1,k])
# @constraint(md, [j=1:nys, k=1:nzs], a[1,j,k] == a[2,j,k])
# @constraint(md, [j=1:nys, k=1:nzs], a[nxs,j,k] == a[nxs-1,j,k])
# @constraint(md, [i=1:nxs, j=1:nys], a[i,j,1] == a[i,j,2])
# @constraint(md, [i=1:nxs, j=1:nys], a[i,j,nzs] == a[i,j,nzs-1])

rhs_vec = vec(rhs_ext)
nxyz = nxs*nys*nzs
myexp = @expression(md, 1/nxyz * sum((φP * vec(a) - rhs_vec).^2))
@objective(md, Min, myexp)
optimize!(md)
a_fin = value.(a)
npzwrite("sens2ertelpv_5lev_ipopt.npy", a_fin)