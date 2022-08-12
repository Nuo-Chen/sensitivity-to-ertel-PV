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
adjm = ingredients("/Volumes/Scratch/nchen67/ertel_inv/adjoint1.jl")
adjd = ingredients("/Volumes/Scratch/nchen67/ertel_inv/adjoint_d.jl")
st = ingredients("/Volumes/Scratch/nchen67/ertel_inv/constants.jl")

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


ψ_bs_xyz = Float64.(permutedims(ψ_bs,[3,2,1]))
ϕ_bs_xyz = permutedims(ϕ_bs,[3,2,1])
msfm_xyz = repeat(permutedims(msfm,[2,1]), outer = [1,1, 5])
f_xyz = repeat(permutedims(f,[2,1]), outer = [1,1, 5])
ψ̂_xyz = permutedims(ψ̂,[3,2,1])	
initguess_xyz = permutedims(initguess,[3,2,1])
A_xyz = permutedims(repeat(A, outer=[1,209,143]), [2, 3, 1])
dπ_xyz = permutedims(repeat(dπ, outer=[1,209,143]), [2, 3, 1])
    

msfm_xy = permutedims(msfm,[2,1]); 
msfm2=msfm_xy.^2; 
f_xy = permutedims(f,[2,1]); 
#rdx2 = rdx^2 ; rdy2 = rdy^2;

# 5 levels at 1000, 900, 500, 200, 50
idx = [1,3, 11, 17, 20]
vc = npzread(path*"varphi5extlev_inv_consts.npz")
gc = npzread(path*"general5lev_consts.npz")

ψ5 = ψ_bs_xyz[:,:,idx]
ϕ5 = ϕ_bs_xyz[:,:,idx]
initguess5 = initguess_xyz[:,:,idx]
A5 = A[idx]
dπ5 = [sum(dπ[1:2]), sum(dπ[3:10]), sum(dπ[11:16]), sum(dπ[17:end])]

V = st.φPconst(vc["ϕ̄ππ"] ,vc["ϕ̄xπ"],vc["ϕ̄yπ"] ,vc["∇²ϕ̄ππ"],vc["∇²ϕ̄xπ"] ,vc["∇²ϕ̄yπ"] ,vc["ζ"] ,vc["ψ̄xπ"],vc["ψ̄yπ"] ,vc["ψ̄xx"] ,vc["ψ̄yy"]  ,vc["∇f∇ζ"],vc["∇f∇ψ̄xπ"],vc["∇f∇ψ̄yπ"] ,vc["ζxx"] ,vc["ζyy"]  ,vc["ψ̄xπxx"] ,vc["ψ̄xπyy"]  ,vc["ψ̄yπxx"]  ,vc["ψ̄yπyy"], vc["ψ"],)
G = st.gconst(gc["fᵢ"], gc["fⱼ"], gc["A"], gc["rdx"], gc["rdx2"], gc["rdx3"], gc["rdy"], gc["rdy2"], gc["rdy3"], gc["m"], gc["m2"], gc["m3"], gc["dπ"], gc["f"])

nx, ny, nzz = size(ϕ_bs_xyz)
ψ̂5 = ψ̂_xyz[:,:,idx]
itv = 4
slicex = 1:3:nx
slicey = 1:3:ny

rhs = adjm.UP(ψ̂5[slicex, slicey, :], ψ5[slicex, slicey, :], ϕ5[slicex, slicey, :], rdx, rdy, msfm_xyz[slicex, slicey, :], f_xyz[slicex, slicey, :]);

nz = 5+2
nxs, nys, nzs = size(rhs)

rhs_ext = Array{Float64, 3}(undef, nxs, nys, nz)
rhs_ext[:,:,2:end-1] = rhs
rhs_ext[:,:,1] = rhs_ext[:,:,2]
rhs_ext[:,:,end] = rhs_ext[:,:,end-1]
    
ψ_ext = Array{Float64, 3}(undef, nxs, nys, nz)
ψ_ext[:,:,2:end-1] = ψ5[slicex,slicey,:]
ψ_ext[:,:,1] = ψ_ext[:,:,2]
ψ_ext[:,:,end] = ψ_ext[:,:,end-1]


ϕ_ext = Array{Float64, 3}(undef, nxs, nys, nz)
ϕ_ext[:,:,2:end-1] = ϕ5[slicex,slicey,:]
ϕ_ext[:,:,1] = ϕ_ext[:,:,2]
ϕ_ext[:,:,end] = ϕ_ext[:,:,end-1]
#size(ϕ_ext)

dπ_ext = Array{Float64, 1}(undef, nz-1)
dπ_ext[2:end-1] = dπ5
dπ_ext[1] = -14
dπ_ext[end] = -100

initguess_ext = Array{Float64, 3}(undef, nxs, nys, nz)
initguess_ext[:,:,2:end-1] = initguess5[slicex,slicey,:]
initguess_ext[:,:,1] .= 0
initguess_ext[:,:,end] .= 0

A_ext = Array{Float64, 1}(undef, nz)
A_ext[2:end-1] = A5
A_ext[1] = 0.027
A_ext[end] = 0.25

@time begin
  md = Model(Ipopt.Optimizer)
  # md = Model(Gurobi.Optimizer)
  # set_optimizer_attribute(md, "NonConvex", 2)
  set_optimizer_attribute(md, "tol", 1e1)
  @variable(md, a[i=1:nxs, j=1:nys, k=1:nz], start=initguess_ext[i,j,k])
  # @variable(md, z)

  # 	# set_start_value.(a[i=1:nx, j=1:ny, k=1], start=0). # for clp???
  # 	# set_start_value.(a[i=1:nx, j=1:ny, k=nz], start=0)
      
  # @constraint(md, [i=1:nxs, k=1:nz], a[i,1,k] == a[i,2,k])
  # @constraint(md, [i=1:nxs, k=1:nz], a[i,nys,k] == a[i,nys-1,k])
  # @constraint(md, [j=1:nys, k=1:nz], a[1,j,k] == a[2,j,k])
  # @constraint(md, [j=1:nys, k=1:nz], a[nxs,j,k] == a[nxs-1,j,k])
  # @constraint(md, [i=1:nxs, j=1:nys], a[i,j,1] == a[i,j,2])
  # @constraint(md, [i=1:nxs, j=1:nys], a[i,j,nz] == a[i,j,nz-1])

  phip = Array{AffExpr,3}(undef, nxs, nys, nz)
  phip .= 0
  # phip[i=3:nxs-2,j=3:nys-2,k=2:nz-1]
  for i=3:nxs-2,j=3:nys-2,k=2:nz-1
    add_to_expression!(phip[i,j,k], adjd.φP(a, i,j,k, A_ext, G.rdx*itv, G.rdx2*itv^2, G.rdx3*itv^3, G.rdy*itv, G.rdy2*itv^2, G.rdy3*itv^3, G.m[slicex,slicey], G.m2[slicex,slicey], G.m3[slicex,slicey], dπ_ext, f_xy[slicex,slicey], G.fᵢ[slicex,slicey], G.fⱼ[slicex,slicey], V.ϕ̄ππ[slicex,slicey,:], V.ϕ̄xπ[slicex,slicey,:], V.ϕ̄yπ[slicex,slicey,:], V.∇²ϕ̄ππ[slicex,slicey,:], V.∇²ϕ̄xπ[slicex,slicey,:], V.∇²ϕ̄yπ[slicex,slicey,:], V.ζ[slicex,slicey,:], V.ψ̄xπ[slicex,slicey,:], V.ψ̄yπ[slicex,slicey,:], V.ψ̄xx[slicex,slicey,:], V.ψ̄yy[slicex,slicey,:], V.∇f∇ζ[slicex,slicey,:], V.∇f∇ψ̄xπ[slicex,slicey,:], V.∇f∇ψ̄yπ[slicex,slicey,:], V.ζxx[slicex,slicey,:], V.ζyy[slicex,slicey,:], V.ψ̄xπxx[slicex,slicey,:], V.ψ̄xπyy[slicex,slicey,:], V.ψ̄yπxx[slicex,slicey,:], V.ψ̄yπyy[slicex,slicey,:], V.ψ[slicex,slicey,:]) ) 
  end
    
  rhs_inner = rhs_ext[3:nxs-2,3:nys-2,2:nz-1]
  #@constraint(md, residual[i,j,k] >= 0 for i in 3:nxs-2, j in 3:nys-2, k in 2:nz-1)
  @objective(md, Min, sum((phip[3:nxs-2,3:nys-2,2:nz-1] - rhs_inner).^2) )	
  optimize!(md)

  a_fin = value.(a)
  # # heatmap(a_fin)
  npzwrite("sens2ertelpv_5lev_ipopt.npy", a_fin)
end
