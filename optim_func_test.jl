### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ e5469a28-8456-11ec-2293-4b6e9a69d596
begin
	using Pkg
	using NPZ#, NCDatasets
	using JuMP, Clp, Ipopt, Gurobi
end

# ╔═╡ a7b4d828-f66b-4c35-a628-42565938c489
#using Plots

# ╔═╡ eb5064d5-b25a-45e7-9abb-569224a7f293
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

# ╔═╡ f34840c3-654d-417f-acb2-2369e851f450
begin
	grad = ingredients("/Volumes/Scratch/nchen67/ertel_inv/gradient1.jl")
	adjm = ingredients("/Volumes/Scratch/nchen67/ertel_inv/adjoint1.jl")
	adjd = ingredients("/Volumes/Scratch/nchen67/ertel_inv/adjoint_d.jl")
	st = ingredients("/Volumes/Scratch/nchen67/ertel_inv/constants.jl")
end

# ╔═╡ 9e01f440-0419-436f-910d-6766e6b4130f


# ╔═╡ 45051427-8f5a-4798-9603-ef4cf00682ab
begin
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
	slicex = 45:itv:130
	slicey = 40:itv:110
end

# ╔═╡ ffc108c0-daa9-42b2-9089-70e081732356
rhs = adjm.UP(ψ̂5[slicex, slicey, :], ψ5[slicex, slicey, :], ϕ5[slicex, slicey, :], rdx, rdy, msfm_xyz[slicex, slicey, :], f_xyz[slicex, slicey, :]);

# ╔═╡ 983121a8-34d2-42d6-a60e-3646a3f05ede
# md = Model(optimizer_with_attributes(Ipopt.Optimizer, "Threads"=>8)

# ╔═╡ 9917ad67-647e-4601-80d8-57b8fc19f402
size(rhs)

# ╔═╡ 671cadff-a122-4f9c-8686-17c878cf4e4f
npzwrite("rhs.npy", rhs)

# ╔═╡ 321bfb04-f450-4795-8c29-ef26254a9fed
begin
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
	
	# md = Model(Ipopt.Optimizer)
	# md = Model(Gurobi.Optimizer)
	md = Model(optimizer_with_attributes(Gurobi.Optimizer, "Threads"=>8))
	# set_optimizer_attribute(md, "NonConvex", 2)
	# set_optimizer_attribute(md, "tol", 1e1)
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
	# optimize!(md)
	
	a_fin = value.(a)
	# # heatmap(a_fin)
	npzwrite("sens2ertelpv_5lev_gurobi.npy", a_fin)
end

# ╔═╡ aaeaff04-56ac-4890-a108-c8a43eda06b5
# begin
# 	md = Model(Ipopt.Optimizer)
# 	nx, ny, nzz = size(ϕ_bs_xyz)
# 	nz = 5
	
# 	@variable(md, a[1:nx, 1:ny, 1:nz])
# 	@constraint(md, [i=1:nx, k=1:nz], a[i,1,k] == a[i,2,k])
# 	@constraint(md, [i=1:nx, k=1:nz], a[i,ny,k] == a[i,ny-1,k])
# 	@constraint(md, [j=1:ny, k=1:nz], a[1,j,k] == a[2,j,k])
# 	@constraint(md, [j=1:ny, k=1:nz], a[nx,j,k] == a[nx-1,j,k])
# 	@constraint(md, [i=1:nx, j=1:ny], a[i,j,1] == a[i,j,2])
# 	@constraint(md, [i=1:nx, j=1:ny], a[i,j,nz] == a[i,j,nz-1])

# 	@variable(md, t<=100)
# 	@constraint(md, [i=2:nx-1, j=2:ny-1, k=2:nz-1], t >= rhs[i,j,k] - ved.φP(a, ψ_bs_xyz[1:nx, 1:ny, 1:nz], ϕ_bs_xyz[1:nx, 1:ny, 1:nz], dπ[1:nz], rdx, rdy, rdx2, rdy2, msfm_xy[1:nx, 1:ny], msfm2[1:nx, 1:ny], A[1:nz], f_xy[1:nx, 1:ny], i,j,k))
# 	@constraint(md, [i=2:nx-1, j=2:ny-1, k=2:nz-1], t >= ved.φP(a, ψ_bs_xyz[1:nx, 1:ny, 1:nz], ϕ_bs_xyz[1:nx, 1:ny, 1:nz], dπ[1:nz], rdx, rdy, rdx2, rdy2, msfm_xy[1:nx, 1:ny], msfm2[1:nx, 1:ny], A[1:nz], f_xy[1:nx, 1:ny], i,j,k) - rhs[i,j,k])
	
# 	#myexp = @expression(md, sum(t[i,j,k]) for i in 2:nx-1, j in 2:ny-1, k in 2:nz-1)
# 	myexp = @expression(md, sum(rhs[i,j,k] - ved.φP(a, ψ_bs_xyz[1:nx, 1:ny, 1:nz], ϕ_bs_xyz[1:nx, 1:ny, 1:nz], dπ[1:nz], rdx, rdy, rdx2, rdy2, msfm_xy[1:nx, 1:ny], msfm2[1:nx, 1:ny], A[1:nz], f_xy[1:nx, 1:ny], i,j,k) for i in 2:nx-1, j in 2:ny-1, k in 2:nz-1))
	
# 	@objective(md, Min, myexp)	
# 	optimize!(md)
	
# 	a_fin = value.(a)
# 	heatmap(a_fin)
# end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Clp = "e2554f3b-3117-50c0-817c-e040a3ddf72d"
Gurobi = "2e9cd046-0924-5485-92f1-d5272153d98b"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
NPZ = "15e1cf62-19b3-5cfa-8e77-841668bca605"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[compat]
Clp = "~1.0.0"
Gurobi = "~0.11.1"
Ipopt = "~1.0.2"
JuMP = "~0.23.2"
NPZ = "~0.4.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6252039f98492252f9e47c312c8ffda0e3b9e78d"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Clp]]
deps = ["Clp_jll", "MathOptInterface"]
git-tree-sha1 = "8588186cc9b298ef08dda8732d96fb57ed993e1b"
uuid = "e2554f3b-3117-50c0-817c-e040a3ddf72d"
version = "1.0.0"

[[deps.Clp_jll]]
deps = ["Artifacts", "CoinUtils_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "MUMPS_seq_jll", "OpenBLAS32_jll", "Osi_jll", "Pkg"]
git-tree-sha1 = "b1031dcfbb44553194c9e650feb5ab65e372504f"
uuid = "06985876-5285-5a41-9fcb-8948a742cc53"
version = "100.1700.601+0"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.CoinUtils_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "44173e61256f32918c6c132fc41f772bab1fb6d1"
uuid = "be027038-0da8-5614-b30d-e42594cb92df"
version = "200.1100.400+0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.Gurobi]]
deps = ["Libdl", "MathOptInterface", "Pkg"]
git-tree-sha1 = "67d273949798c82781756b54af20d73b4dbe9e3f"
uuid = "2e9cd046-0924-5485-92f1-d5272153d98b"
version = "0.11.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.Ipopt]]
deps = ["Ipopt_jll", "MathOptInterface"]
git-tree-sha1 = "8b7b5fdbc71d8f88171865faa11d1c6669e96e32"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "1.0.2"

[[deps.Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "e3e202237d93f18856b6ff1016166b0f172a49a8"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "300.1400.400+0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuMP]]
deps = ["Calculus", "DataStructures", "ForwardDiff", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "NaNMath", "OrderedCollections", "Printf", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "c48de82c5440b34555cb60f3628ebfb9ab3dc5ef"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "0.23.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1d31872bb9c5e7ec1f618e8c4a56c8b0d9bddc7e"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.1+0"

[[deps.MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "29de2841fa5aefe615dea179fcde48bb87b58f57"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "5.4.1+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "Printf", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "a62df301482a41cb7b1db095a4e6949ba7eb3349"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.1.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NPZ]]
deps = ["Compat", "ZipFile"]
git-tree-sha1 = "fbfb3c151b0308236d854c555b43cdd84c1e5ebf"
uuid = "15e1cf62-19b3-5cfa-8e77-841668bca605"
version = "0.4.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c6c2ed4b7acd2137b878eb96c68e63b76199d0f"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.17+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Osi_jll]]
deps = ["Artifacts", "CoinUtils_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "28e0ddebd069f605ab1988ab396f239a3ac9b561"
uuid = "7da25872-d9ce-5375-a4d3-7a845f58efdd"
version = "0.10800.600+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "6976fab022fea2ffea3d945159317556e5dad87c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "3593e69e469d2111389a9bd06bac1f3d730ac6de"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.4"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═e5469a28-8456-11ec-2293-4b6e9a69d596
# ╠═a7b4d828-f66b-4c35-a628-42565938c489
# ╠═eb5064d5-b25a-45e7-9abb-569224a7f293
# ╠═f34840c3-654d-417f-acb2-2369e851f450
# ╠═9e01f440-0419-436f-910d-6766e6b4130f
# ╠═45051427-8f5a-4798-9603-ef4cf00682ab
# ╠═ffc108c0-daa9-42b2-9089-70e081732356
# ╠═983121a8-34d2-42d6-a60e-3646a3f05ede
# ╠═9917ad67-647e-4601-80d8-57b8fc19f402
# ╠═671cadff-a122-4f9c-8686-17c878cf4e4f
# ╠═321bfb04-f450-4795-8c29-ef26254a9fed
# ╠═aaeaff04-56ac-4890-a108-c8a43eda06b5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
