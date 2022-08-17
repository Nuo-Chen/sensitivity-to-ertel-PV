### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 66d74fba-1d88-11ed-3080-2f0eb57b24fa
begin
    import Pkg
   	Pkg.activate(".")
	using NPZ
	using RegularizedLeastSquares
	using Plots
end

# ╔═╡ b601bfd4-1aec-4f14-9462-a2858d673eaa
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

# ╔═╡ c28f2601-ea52-4bcb-9c40-06521988f0cd
adjm = ingredients("adjoint1.jl")

# ╔═╡ 912bb54b-69dd-4f8c-bce2-b656a2557131
begin
	φP = npzread("varphiP_mat.npy")
	
	path = "../"

	ψ_bs = npzread(path*"psi_bs20lev.npy")
	ϕ_bs = npzread(path*"phi_bs20lev.npy")
	frompython = npzread(path*"forjulia.npz")
	rdx = Float64.(frompython["rdx"])
	rdy = Float64.(frompython["rdy"])
	A = frompython["A"]
	msfm = frompython["msfm"]
	dπ = frompython["dexner"]
	f = frompython["f"]
	psiphihat = npzread("2021041515z_sens2psi_sens2phi_20lev.npz")
	ψ̂ = psiphihat["apsi"]
	ϕ̂ = psiphihat["aphi"]
	guess = npzread(path*"2021041515z_inverted_20lev.npz")
	initguess = reshape(guess["aq"], (20,143,209))

	ψ_bs_xyz = Float64.(permutedims(ψ_bs,[3,2,1]))
	ϕ_bs_xyz = permutedims(ϕ_bs,[3,2,1])
	msfm_xyz = repeat(permutedims(msfm,[2,1]), outer = [1,1, 7])
	f_xyz = repeat(permutedims(f,[2,1]), outer = [1,1, 7])
	ψ̂_xyz = permutedims(ψ̂,[3,2,1])	
	ϕ̂_xyz = permutedims(ϕ̂,[3,2,1])	
	initguess_xyz = permutedims(initguess,[3,2,1])
	A_xyz = permutedims(repeat(A, outer=[1,209,143]), [2, 3, 1])
	dπ_xyz = permutedims(repeat(dπ, outer=[1,209,143]), [2, 3, 1])
	    

	msfm_xy = permutedims(msfm,[2,1]); 
	msfm2=msfm_xy.^2; 
	f_xy = permutedims(f,[2,1]); 
	rdx2 = rdx^2 ; rdy2 = rdy^2;

	String("read data success")

end

# ╔═╡ eb48efb6-9fa2-4c09-8dbd-5027bcffe7cd
begin
	# # 5 levels at 1000, 900, 500, 200, 50
	idx = [1,3, 11, 17, 20]
	# vc = npzread(path*"varphi5extlev_inv_consts.npz")
	# gc = npzread(path*"general5lev_consts.npz")
	
	ψ5 = ψ_bs_xyz[:,:,idx]
	ϕ5 = ϕ_bs_xyz[:,:,idx]
	initguess5 = initguess_xyz[:,:,idx]
	A5 = A[idx]
	dπ5 = [sum(dπ[1:2]), sum(dπ[3:10]), sum(dπ[11:16]), sum(dπ[17:end])]

	itv = 2
	slicex = 45:itv:130
	slicey = 40:itv:110
	
	nxs, nys, nzz = size(ϕ_bs_xyz[slicex,slicey,idx])
	nz = nzz + 2
	
	# rhs_ext = Array{Float64, 3}(undef, nxs, nys, nz)
	# rhs_ext[:,:,2:end-1] = rhs
	# rhs_ext[:,:,1] = rhs_ext[:,:,2]
	# rhs_ext[:,:,end] = rhs_ext[:,:,end-1]

	ψ̂_ext = Array{Float64, 3}(undef, nxs, nys, nz)
	ψ̂_ext[:,:,2:end-1] = ψ̂_xyz[slicex,slicey,idx]
	ψ̂_ext[:,:,1] = ψ̂_ext[:,:,2]
	ψ̂_ext[:,:,end] = ψ̂_ext[:,:,end-1]
	
	
	ϕ̂_ext = Array{Float64, 3}(undef, nxs, nys, nz)
	ϕ̂_ext[:,:,2:end-1] = ϕ̂_xyz[slicex,slicey,idx]
	ϕ̂_ext[:,:,1] = ϕ̂_ext[:,:,2]
	ϕ̂_ext[:,:,end] = ϕ̂_ext[:,:,end-1]
	
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

	na = [CartesianIndex()]

	rhs = adjm.UP(ψ̂_ext, ϕ̂_ext, ψ_ext, ϕ_ext, rdx, rdy, msfm_xyz[slicex, slicey, :], f_xyz[slicex, slicey, :])./A_ext[na,na,:]

	# npzwrite("rhs.npy", rhs)
	# nz = 5+2
	# nxs, nys, nzs = size(rhs)

	String("UP success")
end

# ╔═╡ 30dc9b58-9094-45f6-abb1-b422415d675b
size(rhs), size(φP), 43*36*7, nxs, nys, nz

# ╔═╡ ff5a646b-e245-4d8d-856d-fbb2d64d8fe6
contour(φP[1000:1500,1000:1500]*1)

# ╔═╡ e2981b58-9d3a-4df0-94dc-282593d0541c
begin
	vrhs32= Float32.(vec(rhs*1e11))
	φP32 = Float32.(φP)

	nxyz = nxs*nys*nz
	reg = Regularization("L2", 0.000001; shape=(nxyz))
	solver = createLinearSolver("cgnr", φP, iterations=5000, regMatrix=reg)

	x_approx = solve(solver,vrhs32)
	ima = reshape(x_approx, nxs, nys, nz)
	String("inverse problem seccuss")
end

# ╔═╡ ae886008-7ef6-44be-8a19-4cb99a72ffde
npzwrite("sens2EPV.npy", ima)

# ╔═╡ 95764ddc-1d65-4cac-bcd0-4e6993081fca
begin
	ly = @layout [a b]
	p1 = contour(rhs[:,:,3], title = "rhs")
	p2 = contour(ima[:,:,3], title = "Eq̂")
	
	plot(p1, p2, layout=ly)
	plot!(size=(600,300))
end

# ╔═╡ 8455a83f-b412-43de-a2f2-05e1745903c5
# grad = ingredients("gradient1.jl")
# adjd = ingredients("adjoint_d.jl")
# st = ingredients("constants.jl")



# V = st.φPconst(vc["ϕ̄ππ"] ,vc["ϕ̄xπ"],vc["ϕ̄yπ"] ,vc["∇²ϕ̄ππ"],vc["∇²ϕ̄xπ"] ,vc["∇²ϕ̄yπ"] ,vc["ζ"] ,vc["ψ̄xπ"],vc["ψ̄yπ"] ,vc["ψ̄xx"] ,vc["ψ̄yy"]  ,vc["∇f∇ζ"],vc["∇f∇ψ̄xπ"],vc["∇f∇ψ̄yπ"] ,vc["ζxx"] ,vc["ζyy"]  ,vc["ψ̄xπxx"] ,vc["ψ̄xπyy"]  ,vc["ψ̄yπxx"]  ,vc["ψ̄yπyy"], vc["ψ"],)
# G = st.gconst(gc["fᵢ"], gc["fⱼ"], gc["A"], gc["rdx"], gc["rdx2"], gc["rdx3"], gc["rdy"], gc["rdy2"], gc["rdy3"], gc["m"], gc["m2"], gc["m3"], gc["dπ"], gc["f"])







# ╔═╡ Cell order:
# ╠═66d74fba-1d88-11ed-3080-2f0eb57b24fa
# ╠═b601bfd4-1aec-4f14-9462-a2858d673eaa
# ╠═c28f2601-ea52-4bcb-9c40-06521988f0cd
# ╠═912bb54b-69dd-4f8c-bce2-b656a2557131
# ╠═eb48efb6-9fa2-4c09-8dbd-5027bcffe7cd
# ╠═30dc9b58-9094-45f6-abb1-b422415d675b
# ╠═ff5a646b-e245-4d8d-856d-fbb2d64d8fe6
# ╠═e2981b58-9d3a-4df0-94dc-282593d0541c
# ╠═ae886008-7ef6-44be-8a19-4cb99a72ffde
# ╠═95764ddc-1d65-4cac-bcd0-4e6993081fca
# ╠═8455a83f-b412-43de-a2f2-05e1745903c5
