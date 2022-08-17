### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 087626d2-dd26-11ec-2371-cd85f62600fe
begin
    import Pkg
   	Pkg.activate()
	#using Plots
	using NPZ
	#using PlutoUI
end

# ╔═╡ 239a0ed4-d7ef-4429-a556-daca03b0b012
begin
	path = "../"
	
	ψ_bs = npzread(path*"psi_bs20lev.npy")
	ϕ_bs = npzread(path*"phi_bs20lev.npy")
	frompython = npzread(path*"forjulia.npz")
	rdx = Float64.(frompython["rdx"])
	rdy = Float64.(frompython["rdy"])
	A = frompython["A"]
	msfm_yx = frompython["msfm"]
	dπ = frompython["dexner"]
	f_yx = frompython["f"]
	ψ̂ = npzread(path*"apsi20lev.npy");
	guess = npzread(path*"2021041515z_inverted_20lev.npz")
	initguess = reshape(guess["aq"], (20,143,209))
	
	String("load origin field")
end

# ╔═╡ e8b841c7-5a9e-45c7-89a5-e9039d4e09f6
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

# ╔═╡ 407f03b8-8bcf-4425-84d3-91b5e542c41a
st = ingredients("../constants.jl")

# ╔═╡ c39217a8-f91f-4ff4-940c-a93068238856
begin
	# 5 levels at 1000, 900, 500, 200, 50
	idx = [1,3, 11, 17, 20]
    nzs = size(idx)[1] + 2
    
    itv = 2
	slicex = 45:itv:130
	slicey = 40:itv:110
	
	nxs = size(slicex)[1]
	nys = size(slicey)[1]
    nxyz = nxs*nys*nzs
    
    
	ψ_bs_3d = Float64.(permutedims(ψ_bs,[3,2,1]))
	ϕ_bs_3d = permutedims(ϕ_bs,[3,2,1])
	
	ψ̂_3d = permutedims(ψ̂,[3,2,1])	
	initguess_3d = permutedims(initguess,[3,2,1])
	A_3d = permutedims(repeat(A, outer=[1,nxs,nys]), [2, 3, 1])
	dπ_xyz = permutedims(repeat(dπ, outer=[1,nxs,nys]), [2, 3, 1])
		
    msfm = permutedims(msfm_yx,[2,1])[slicex, slicey]; 
	msfm_3d = repeat(msfm, outer = [1,1, nzs])
    msfm2 = msfm.^2; 
    
    f = permutedims(f_yx,[2,1])[slicex, slicey]; 
	f3d = repeat(f, outer = [1,1, nzs])

    
	vc = npzread(path*"varphi5extlev_inv_consts.npz")
	gc = npzread(path*"general5lev_consts.npz")
	
	ψ5 = ψ_bs_3d[slicex, slicey, idx]
	ϕ5 = ϕ_bs_3d[slicex, slicey, idx]
    ψ̂5 = ψ̂_3d[slicex, slicey, idx]
	ig5 = initguess_3d[slicex, slicey, idx]
	A5_3d = A_3d[:,:, idx]
	dπ5 = [sum(dπ[1:2]), sum(dπ[3:10]), sum(dπ[11:16]), sum(dπ[17:end])]
    dπ5_3d = permutedims(repeat(dπ5, outer=[1,nxs,nys]), [2, 3, 1])

	V = st.φPconst(vc["ϕ̄ππ"], vc["ϕ̄xπ"], vc["ϕ̄yπ"], vc["∇²ϕ̄ππ"], vc["∇²ϕ̄xπ"], vc["∇²ϕ̄yπ"], vc["ζ"], vc["ζx"], vc["ζy"], vc["ζxx"], vc["ζyy"], vc["∇²ζ"], vc["ζxy"], vc["ψ̄xπ"], vc["ψ̄yπ"], vc["ψ̄xx"], vc["ψ̄yy"], vc["ψ̄xy"], vc["∂xψ̄xπ"], vc["∂yψ̄xπ"], vc["∂yyψ̄xπ"], vc["∂xxψ̄xπ"], vc["∂xyψ̄xπ"], vc["∇²ψ̄xπ"], vc["∂xψ̄yπ"], vc["∂yψ̄yπ"], vc["∂yyψ̄yπ"], vc["∂xxψ̄yπ"], vc["∂xyψ̄yπ"], vc["∇²ψ̄yπ"], vc["∇f∇ζ"], vc["∇f∇ψ̄xπ"], vc["∇f∇ψ̄yπ"], vc["ψ"], vc["ϕ"], vc["fx"], vc["fy"], vc["f"])
		
	G = st.gconst(gc["fᵢ"], gc["fⱼ"], gc["A"], gc["rdx"], gc["rdx2"], gc["rdx3"], gc["rdy"], gc["rdy2"], gc["rdy3"], gc["m"], gc["m2"], gc["m3"], gc["dπ"], gc["f"])

	h = dπ
    dx = itv * rdx
    dy = itv * rdy
    rdx2 = rdx^2 ; rdy2 = rdy^2;
    
	String("permute axes and extend to 3d")
end

# ╔═╡ 509e5ee1-c3c7-47d5-ae74-8de5a40c43d9
begin
	ψ_ext = Array{Float64, 3}(undef, nxs, nys, nzs)
	ψ_ext[:,:,2:end-1] = ψ5
	ψ_ext[:,:,1] = ψ_ext[:,:,2]
	ψ_ext[:,:,end] = ψ_ext[:,:,end-1]
	
	
	ϕ_ext = Array{Float64, 3}(undef, nxs, nys, nzs)
	ϕ_ext[:,:,2:end-1] = ϕ5
	ϕ_ext[:,:,1] = ϕ_ext[:,:,2]
	ϕ_ext[:,:,end] = ϕ_ext[:,:,end-1]
	
	
	dπ_ext = Array{Float64, 3}(undef, nxs, nys, nzs-1)
	dπ_ext[:,:,2:end-1] = dπ5_3d
	dπ_ext[:,:,1] .= -14
	dπ_ext[:,:,end] .= -100
	
	initguess_ext = Array{Float64, 3}(undef, nxs, nys, nzs)
	initguess_ext[:,:,2:end-1] = ig5
	initguess_ext[:,:,1] .= 0
	initguess_ext[:,:,end] .= 0
	
	A_ext = Array{Float64, 3}(undef, nxs, nys, nzs)
	A_ext[:,:,2:end-1] = A5_3d
	A_ext[:,:,1] .= 0.027
	A_ext[:,:,end] .= 0.25
	
	String("Extend some variable matrix with top and bottom buffer layers")
end

# ╔═╡ 0dc4dbd5-a367-4c35-837d-3fc65205885e
function logop(row, m, n, l, op, a)
	col = ((l-1)*nys + n -1)*nxs + m
	a[col, row] = op
	return a
end

# ╔═╡ 3428bdd9-7405-4167-b760-ad5ef4a23a15
begin 
	∂x = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
	    tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	
	    coef = 1/dx * msfm[i,j]
	    ∂x = logop(tijk, ip1, j, k, 1/2*coef, ∂x)
	    ∂x = logop(tijk, im1, j, k, -1/2*coef, ∂x)
	end
	
	∂y = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
	    tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	
	    coef = 1/dy * msfm[i,j]
	    ∂y = logop(tijk, i, jp1, k, 1/2*coef, ∂y)
	    ∂y = logop(tijk, i, jm1, k, -1/2*coef, ∂y)
	end
	
	∂π = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
	    tijk = ((k-1)*nys +j-1)*nxs + i
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	
	    ∂π = logop(tijk, i, j, kp1, 1/(h[k]*(1+h[k]/h[km1])), ∂π)
	    ∂π = logop(tijk, i, j, km1, -h[k]/(h[km1]^2+h[k]*h[km1]), ∂π)
	    ∂π = logop(tijk, i, j, k, 1/h[km1]-1/h[k], ∂π)
	end

	∂²x = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	
	    coef = 1/dx^2 * msfm[i,j]^2
	    ∂²x = logop(tijk, ip1, j, k, coef, ∂²x)
	    ∂²x = logop(tijk, im1, j, k, coef, ∂²x)
	    ∂²x = logop(tijk, i, j, k, -2*coef, ∂²x)
	end

	∂²y = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	
	    coef = 1/dy^2 * msfm[i,j]^2
	    ∂²y = logop(tijk, i, jp1, k, coef, ∂²y)
	    ∂²y = logop(tijk, i, jm1, k, coef, ∂²y)
	    ∂²y = logop(tijk, i, j, k, -2*coef, ∂²y)
	end
	
	∂²xy = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
	
	    coef = 1/dx/dy * msfm[i,j]^2
	    ∂²xy = logop(tijk, ip1, jp1, k, 1/4*coef, ∂²xy)
	    ∂²xy = logop(tijk, im1, jm1, k, 1/4*coef, ∂²xy)
	    ∂²xy = logop(tijk, im1, jp1, k, -1/4*coef, ∂²xy)
	    ∂²xy = logop(tijk, ip1, jm1, k, -1/4*coef, ∂²xy)
	end
	
	∂²xπ  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dx * msfm[i,j]
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    ∂²xπ = logop(tijk, ip1, j, kp1, 1/2*coefh*coefkp1, ∂²xπ)
	    ∂²xπ = logop(tijk, im1, j, kp1, -1/2*coefh*coefkp1, ∂²xπ)
	    ∂²xπ = logop(tijk, ip1, j, k, 1/2*coefh*coefk, ∂²xπ)
	    ∂²xπ = logop(tijk, im1, j, k, -1/2*coefh*coefk, ∂²xπ)
	    ∂²xπ = logop(tijk, ip1, j, km1, 1/2*coefh*coefkm1, ∂²xπ)
	    ∂²xπ = logop(tijk, im1, j, km1, -1/2*coefh*coefkm1, ∂²xπ)
	end
	
	∂²yπ  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dy * msfm[i,j]
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    ∂²yπ = logop(tijk, i, jp1, kp1, 1/2*coefh*coefkp1, ∂²yπ)
	    ∂²yπ = logop(tijk, i, jm1, kp1, -1/2*coefh*coefkp1, ∂²yπ)
	    ∂²yπ = logop(tijk, i, jp1, k, 1/2*coefh*coefk, ∂²yπ)
	    ∂²yπ = logop(tijk, i, jm1, k, -1/2*coefh*coefk, ∂²yπ)
	    ∂²yπ = logop(tijk, i, jp1, km1, 1/2*coefh*coefkm1, ∂²yπ)
	    ∂²yπ = logop(tijk, i, jm1, km1, -1/2*coefh*coefkm1, ∂²yπ)
	end
	
	∂²π = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
		
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	
	    ∂²π = logop(tijk, i, j, kp1, 2/(h[k]*(h[k]+h[km1])), ∂²π)
	    ∂²π = logop(tijk, i, j, km1, 2/(h[km1]*(h[k]+h[km1])), ∂²π)
	    ∂²π = logop(tijk, i, j, k, -2/h[k]/h[km1], ∂²π)
	end
	
	∇² = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
		
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
		
		coefh = 1/dx^2 * msfm[i,j]^2
		
	    ∇² = logop(tijk, i, jm1, k, coefh, ∇²)
	    ∇² = logop(tijk, im1, j, k, coefh, ∇²)
	    ∇² = logop(tijk, i, j, k, -4*coefh, ∇²)
	    ∇² = logop(tijk, i, jp1, k, coefh, ∇²)
	    ∇² = logop(tijk, ip1, j, k, coefh, ∇²)
	end
	String("1st and 2nd order: ∂x, ∂y, ∂π, ∂²x, ∂²y, ∂²xy, ∂²xπ, ∂²yπ, ∂²π, ∇²")
end

# ╔═╡ 37f870cf-a9a7-42a0-8239-72b21270461c
begin
	∂³x = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    im2 = max(i-2,1)
	    ip2 = min(i+2,nxs)
	
	    coef = 1/dx^3 * msfm[i,j]^3
	    ∂³x = logop(tijk, ip2, j, k, 1/2*coef, ∂³x)
	    ∂³x = logop(tijk, ip1, j, k, -coef, ∂³x)
	    ∂³x = logop(tijk, im1, j, k, coef, ∂³x)
	    ∂³x = logop(tijk, im2, j, k, -1/2*coef, ∂³x)
	
	end
	
	∂³y = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	    jm2 = max(j-2,1)
	    jp2 = min(j+2,nys)
	
	    coef = 1/dy^3 * msfm[i,j]^3
	    ∂³y = logop(tijk, i, jp2, k, 1/2*coef, ∂³y)
	    ∂³y = logop(tijk, i, jp1, k, -coef, ∂³y)
	    ∂³y = logop(tijk, i, jm1, k, coef, ∂³y)
	    ∂³y = logop(tijk, i, jm2, k, -1/2*coef, ∂³y)
	
	end
	
	∂³xy2 = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	
	    coef = 1/dx/dy^2 * msfm[i,j]^3
	    ∂³xy2 = logop(tijk, ip1, jp1, k, 1/2*coef, ∂³xy2)
	    ∂³xy2 = logop(tijk, ip1, jm1, k, 1/2*coef, ∂³xy2)
	    ∂³xy2 = logop(tijk, ip1, j, k, -coef, ∂³xy2)
	    ∂³xy2 = logop(tijk, im1, jp1, k, -1/2*coef, ∂³xy2)
	    ∂³xy2 = logop(tijk, im1, jm1, k, -1/2*coef, ∂³xy2)
	    ∂³xy2 = logop(tijk, im1, j, k, coef, ∂³xy2)
	
	end
	
	∂³x2y = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	
	    coef = 1/dx^2/dy * msfm[i,j]^3
	    ∂³x2y = logop(tijk, ip1, jp1, k, 1/2*coef, ∂³x2y)
	    ∂³x2y = logop(tijk, im1, jp1, k, 1/2*coef, ∂³x2y)
	    ∂³x2y = logop(tijk, i, jp1, k, -coef, ∂³x2y)
	    ∂³x2y = logop(tijk, ip1, jm1, k, -1/2*coef, ∂³x2y)
	    ∂³x2y = logop(tijk, im1, jm1, k, -1/2*coef, ∂³x2y)
	    ∂³x2y = logop(tijk, i, jm1, k, coef, ∂³x2y)
	end

	∂³xπ2  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dx * msfm[i,j]
	    coefkp1 = 2/(h[k]*(h[k]+h[km1]))
	    coefk = -2/h[k]/h[km1]
	    coefkm1 = 2/(h[km1]*(h[k]+h[km1]))
	            
	    ∂³xπ2 = logop(tijk, ip1, j, kp1, 1/2*coefh*coefkp1, ∂³xπ2)
	    ∂³xπ2 = logop(tijk, im1, j, kp1, -1/2*coefh*coefkp1, ∂³xπ2)
	    ∂³xπ2 = logop(tijk, ip1, j, k, 1/2*coefh*coefk, ∂³xπ2)
	    ∂³xπ2 = logop(tijk, im1, j, k, -1/2*coefh*coefk, ∂³xπ2)
	    ∂³xπ2 = logop(tijk, ip1, j, km1, 1/2*coefh*coefkm1, ∂³xπ2)
	    ∂³xπ2 = logop(tijk, im1, j, km1, -1/2*coefh*coefkm1, ∂³xπ2)
	end
	
	∂³yπ2  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dy * msfm[i,j]
	    coefkp1 = 2/(h[k]*(h[k]+h[km1]))
	    coefk = -2/h[k]/h[km1]
	    coefkm1 = 2/(h[km1]*(h[k]+h[km1]))
	            
	    ∂³yπ2 = logop(tijk, i, jp1, kp1, 1/2*coefh*coefkp1, ∂³yπ2)
	    ∂³yπ2 = logop(tijk, i, jm1, kp1, -1/2*coefh*coefkp1, ∂³yπ2)
	    ∂³yπ2 = logop(tijk, i, jp1, k, 1/2*coefh*coefk, ∂³yπ2)
	    ∂³yπ2 = logop(tijk, i, jm1, k, -1/2*coefh*coefk, ∂³yπ2)
	    ∂³yπ2 = logop(tijk, i, jp1, km1, 1/2*coefh*coefkm1, ∂³yπ2)
	    ∂³yπ2 = logop(tijk, i, jm1, km1, -1/2*coefh*coefkm1, ∂³yπ2)
	end

	∂³x2π  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dx^2 * msfm[i,j]^2
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    ∂³x2π = logop(tijk, ip1, j, kp1, coefh*coefkp1, ∂³x2π)
		∂³x2π = logop(tijk, i, j, kp1, -2*coefh*coefkp1, ∂³x2π)
	    ∂³x2π = logop(tijk, im1, j, kp1, coefh*coefkp1, ∂³x2π)
	    ∂³x2π = logop(tijk, ip1, j, k, coefh*coefk, ∂³x2π)
		∂³x2π = logop(tijk, i, j, k, -2*coefh*coefk, ∂³x2π)
	    ∂³x2π = logop(tijk, im1, j, k, coefh*coefk, ∂³x2π)
	    ∂³x2π = logop(tijk, ip1, j, km1, coefh*coefkm1, ∂³x2π)
		∂³x2π = logop(tijk, i, j, km1, -2*coefh*coefkm1, ∂³x2π)
	    ∂³x2π = logop(tijk, im1, j, km1, coefh*coefkm1, ∂³x2π)
	end
	
	∂³y2π  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dy^2 * msfm[i,j]^2
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    ∂³y2π = logop(tijk, i, jp1, kp1, coefh*coefkp1, ∂³y2π)
		∂³y2π = logop(tijk, i, j, kp1, -2*coefh*coefkp1, ∂³y2π)
	    ∂³y2π = logop(tijk, i, jm1, kp1, coefh*coefkp1, ∂³y2π)
	    ∂³y2π = logop(tijk, i, jp1, k, coefh*coefk, ∂³y2π)
		∂³y2π = logop(tijk, i, j, k, -2*coefh*coefk, ∂³y2π)
	    ∂³y2π = logop(tijk, i, jm1, k, coefh*coefk, ∂³y2π)
	    ∂³y2π = logop(tijk, i, jp1, km1, coefh*coefkm1, ∂³y2π)
		∂³y2π = logop(tijk, i, j, km1, -2*coefh*coefkm1, ∂³y2π)
	    ∂³y2π = logop(tijk, i, jm1, km1, coefh*coefkm1, ∂³y2π)
	end


	∂³xyπ = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	
	    coefh = 1/dx/dy * msfm[i,j]^2
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    ∂³xyπ = logop(tijk, ip1, jp1, kp1, 1/4*coefh*coefkp1, ∂³xyπ)
	    ∂³xyπ = logop(tijk, im1, jm1, kp1, 1/4*coefh*coefkp1, ∂³xyπ)
	    ∂³xyπ = logop(tijk, im1, jp1, kp1, -1/4*coefh*coefkp1, ∂³xyπ)
	    ∂³xyπ = logop(tijk, ip1, jm1, kp1, -1/4*coefh*coefkp1, ∂³xyπ)
	    
	    ∂³xyπ = logop(tijk, ip1, jp1, k, 1/4*coefh*coefk, ∂³xyπ)
	    ∂³xyπ = logop(tijk, im1, jm1, k, 1/4*coefh*coefk, ∂³xyπ)
	    ∂³xyπ = logop(tijk, im1, jp1, k, -1/4*coefh*coefk, ∂³xyπ)
	    ∂³xyπ = logop(tijk, ip1, jm1, k, -1/4*coefh*coefk, ∂³xyπ)
	    
	    ∂³xyπ = logop(tijk, ip1, jp1, km1, 1/4*coefh*coefkm1, ∂³xyπ)
	    ∂³xyπ = logop(tijk, im1, jm1, km1, 1/4*coefh*coefkm1, ∂³xyπ)
	    ∂³xyπ = logop(tijk, im1, jp1, km1, -1/4*coefh*coefkm1, ∂³xyπ)
	    ∂³xyπ = logop(tijk, ip1, jm1, km1, -1/4*coefh*coefkm1, ∂³xyπ)
	end
	String("3rd order: ∂³x, ∂³y, ∂³x2y, ∂³xy2, ∂³xπ2, ∂³yπ2, ∂³x2π, ∂³y2π, ∂³xyπ")
end

# ╔═╡ 59586306-5cd6-4fc5-bb19-62d9153c077a
begin
	global ∇⁴ = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs, j in 3:nys-2, i in 3:nxs-2
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
	    jm2 = max(j-2,1)
	    im2 = max(i-2,1)
	    jp2 = min(j+2,nys)
	    ip2 = min(i+2,nxs)
	
	
	    global ∇⁴ = logop(tijk, ip2, j, k, 1., ∇⁴)
	    global ∇⁴ = logop(tijk, im2, j, k, 1., ∇⁴)
	    global ∇⁴ = logop(tijk, i, jp2, k, 1., ∇⁴)
	    global ∇⁴ = logop(tijk, i, jm2, k, 1., ∇⁴)
	    global ∇⁴ = logop(tijk, ip1, j, k, -8., ∇⁴)
	    global ∇⁴ = logop(tijk, im1, j, k, -8., ∇⁴)
	    global ∇⁴ = logop(tijk, i, jp1, k, -8., ∇⁴)
	    global ∇⁴ = logop(tijk, i, jm1, k, -8., ∇⁴)
	    global ∇⁴ = logop(tijk, ip1, jp1, k, 2., ∇⁴)
	    global ∇⁴ = logop(tijk, im1, jm1, k, 2., ∇⁴)
	    global ∇⁴ = logop(tijk, ip1, jm1, k, 2., ∇⁴)
	    global ∇⁴ = logop(tijk, im1, jp1, k, 2., ∇⁴)
	    global ∇⁴ = logop(tijk, i, j, k, 20., ∇⁴)
	end
	
	global ∂⁴x3π  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    im2 = max(i-2,1)
	    ip2 = min(i+2,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dx^3 * msfm[i,j]^3
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    global ∂⁴x3π = logop(tijk, ip2, j, kp1, 1/2*coefh*coefkp1, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, ip1, j, kp1, -coefh*coefkp1, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, im1, j, kp1, coefh*coefkp1, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, im2, j, kp1, -1/2*coefh*coefkp1, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, ip2, j, k, 1/2*coefh*coefk, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, ip1, j, k, -coefh*coefk, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, im1, j, k, coefh*coefk, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, im2, j, k, -1/2*coefh*coefk, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, ip2, j, km1, 1/2*coefh*coefkm1, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, ip1, j, km1, -coefh*coefkm1, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, im1, j, km1, coefh*coefkm1, ∂⁴x3π)
	    global ∂⁴x3π = logop(tijk, im2, j, km1, -1/2*coefh*coefkm1, ∂⁴x3π)
	end
	
	global ∂⁴y3π  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	    jm2 = max(j-2,1)
	    jp2 = min(j+2,nys)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dy^3 * msfm[i,j]^3
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    global ∂⁴y3π = logop(tijk, i, jp2, kp1, 1/2*coefh*coefkp1, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jp1, kp1, -coefh*coefkp1, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jm1, kp1, coefh*coefkp1, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jm2, kp1, -1/2*coefh*coefkp1, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jp2, k, 1/2*coefh*coefk, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jp1, k, -coefh*coefk, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jm1, k, coefh*coefk, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jm2, k, -1/2*coefh*coefk, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jp2, km1, 1/2*coefh*coefkm1, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jp1, km1, -coefh*coefkm1, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jm1, km1, coefh*coefkm1, ∂⁴y3π)
	    global ∂⁴y3π = logop(tijk, i, jm2, km1, -1/2*coefh*coefkm1, ∂⁴y3π)
	end
	
	global ∂⁴x2π2  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    im1 = max(i-1,1)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dx^2 * msfm[i,j]^2
	    coefkp1 = 2/(h[k]*(h[k]+h[km1]))
	    coefk = -2/h[k]/h[km1]
	    coefkm1 = 2/(h[km1]*(h[k]+h[km1]))
	            
	    global ∂⁴x2π2 = logop(tijk, ip1, j, kp1, coefh*coefkp1, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, i, j, kp1, -2*coefh*coefkp1, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, im1, j, kp1, coefh*coefkp1, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, ip1, j, k, coefh*coefk, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, i, j, k, -2*coefh*coefk, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, im1, j, k, coefh*coefk, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, ip1, j, km1, coefh*coefkm1, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, i, j, km1, -2*coefh*coefkm1, ∂⁴x2π2)
	    global ∂⁴x2π2 = logop(tijk, im1, j, km1, coefh*coefkm1, ∂⁴x2π2)
	end
	
	global ∂⁴y2π2  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    jp1 = min(j+1,nys)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dy^2 * msfm[i,j]^2
	    coefkp1 = 2/(h[k]*(h[k]+h[km1]))
	    coefk = -2/h[k]/h[km1]
	    coefkm1 = 2/(h[km1]*(h[k]+h[km1]))
	            
	    global ∂⁴y2π2 = logop(tijk, i, jp1, kp1, coefh*coefkp1, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, j, kp1, -2*coefh*coefkp1, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, jm1, kp1, coefh*coefkp1, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, jp1, k, coefh*coefk, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, j, k, -2*coefh*coefk, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, jm1, k, coefh*coefk, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, jp1, km1, coefh*coefkm1, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, j, km1, -2*coefh*coefkm1, ∂⁴y2π2)
	    global ∂⁴y2π2 = logop(tijk, i, jm1, km1, coefh*coefkm1, ∂⁴y2π2)
	end
	
	
	global ∂⁴xy2π = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	
	    coefh = 1/dx/dy^2 * msfm[i,j]^3
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    global ∂⁴xy2π = logop(tijk, ip1, jp1, kp1, 1/2*coefh*coefkp1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, ip1, jm1, kp1, 1/2*coefh*coefkp1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, ip1, j, kp1, -coefh*coefkp1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, jp1, kp1, -1/2*coefh*coefkp1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, jm1, kp1, -1/2*coefh*coefkp1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, j, kp1, coefh*coefkp1, ∂⁴xy2π)
	    
	    global ∂⁴xy2π = logop(tijk, ip1, jp1, k, 1/2*coefh*coefk, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, ip1, jm1, k, 1/2*coefh*coefk, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, ip1, j, k, -coefh*coefk, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, jp1, k, -1/2*coefh*coefk, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, jm1, k, -1/2*coefh*coefk, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, j, k, coefh*coefk, ∂⁴xy2π)
	    
	    global ∂⁴xy2π = logop(tijk, ip1, jp1, km1, 1/2*coefh*coefkm1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, ip1, jm1, km1, 1/2*coefh*coefkm1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, ip1, j, km1, -coefh*coefkm1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, jp1, km1, -1/2*coefh*coefkm1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, jm1, km1, -1/2*coefh*coefkm1, ∂⁴xy2π)
	    global ∂⁴xy2π = logop(tijk, im1, j, km1, coefh*coefkm1, ∂⁴xy2π)
	end
	
	global ∂⁴x2yπ = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	
	    coefh = 1/dx/dy^2 * msfm[i,j]^3
	    coefkp1 = 1/(h[k]*(1+h[k]/h[km1]))
	    coefk = 1/h[km1]-1/h[k]
	    coefkm1 = -h[k]/(h[km1]^2+h[k]*h[km1])
	    
	    global ∂⁴x2yπ = logop(tijk, ip1, jp1,  kp1, 1/2*coefh*coefkp1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, im1, jp1, kp1, 1/2*coefh*coefkp1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, i, jp1, kp1, -coefh*coefkp1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, ip1, jm1, kp1, -1/2*coefh*coefkp1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, im1, jm1, kp1, -1/2*coefh*coefkp1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, i, jm1, kp1, coefh*coefkp1, ∂⁴x2yπ)
	    
	    global ∂⁴x2yπ = logop(tijk, ip1, jp1,  k, 1/2*coefh*coefk, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, im1, jp1, k, 1/2*coefh*coefk, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, i, jp1, k, -coefh*coefk, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, ip1, jm1, k, -1/2*coefh*coefk, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, im1, jm1, k, -1/2*coefh*coefk, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, i, jm1, k, coefh*coefk, ∂⁴x2yπ)
	    
	    
	    global ∂⁴x2yπ = logop(tijk, ip1, jp1, km1, 1/2*coefh*coefkm1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, im1, jp1, km1, 1/2*coefh*coefkm1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, i, jp1, km1, -coefh*coefkm1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, ip1, jm1, km1, -1/2*coefh*coefkm1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, im1, jm1, km1, -1/2*coefh*coefkm1, ∂⁴x2yπ)
	    global ∂⁴x2yπ = logop(tijk, i, jm1, km1, coefh*coefkm1, ∂⁴x2yπ)
	end

	global ∇²∂²π2  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dx/dy * msfm[i,j]^2
	    coefkp1 = 2/(h[k]*(h[k]+h[km1]))
	    coefk = -2/h[k]/h[km1]
	    coefkm1 = 2/(h[km1]*(h[k]+h[km1]))
	            
	    global ∇²∂²π2 = logop(tijk, i, jm1, kp1, coefh*coefkp1, ∇²∂²π2)
	    global ∇²∂²π2 = logop(tijk, im1, j, kp1, coefh*coefkp1, ∇²∂²π2)
		global ∇²∂²π2 = logop(tijk, i, j, kp1, -4*coefh*coefkp1, ∇²∂²π2)
	    global ∇²∂²π2 = logop(tijk, i, jp1, kp1, coefh*coefkp1, ∇²∂²π2)
		global ∇²∂²π2 = logop(tijk, ip1, j, kp1, coefh*coefkp1, ∇²∂²π2)
		
	    global ∇²∂²π2 = logop(tijk, i, jm1, k, coefh*coefk, ∇²∂²π2)
	    global ∇²∂²π2 = logop(tijk, im1, j, k, coefh*coefk, ∇²∂²π2)
		global ∇²∂²π2 = logop(tijk, i, j, k, -4*coefh*coefk, ∇²∂²π2)
	    global ∇²∂²π2 = logop(tijk, i, jp1, k, coefh*coefk, ∇²∂²π2)
		global ∇²∂²π2 = logop(tijk, ip1, j, k, coefh*coefk, ∇²∂²π2)
		
	    global ∇²∂²π2 = logop(tijk, i, jm1, km1, coefh*coefkm1, ∇²∂²π2)
	    global ∇²∂²π2 = logop(tijk, im1, j, km1, coefh*coefkm1, ∇²∂²π2)
		global ∇²∂²π2 = logop(tijk, i, j, km1, -4*coefh*coefkm1, ∇²∂²π2)
	    global ∇²∂²π2 = logop(tijk, i, jp1, km1, coefh*coefkm1, ∇²∂²π2)
		global ∇²∂²π2 = logop(tijk, ip1, j, km1, coefh*coefkm1, ∇²∂²π2)
	end

	∇²∂²xπ  = ∇²*∂²xπ
	∇²∂²yπ  = ∇²*∂²yπ

	global ∂⁴xyπ2  = zeros(Float32, nxyz, nxyz)
	for k in 1:nzs-1, j in 1:nys, i in 1:nxs
		tijk = ((k-1)*nys +j-1)*nxs + i
	    jm1 = max(j-1,1)
	    im1 = max(i-1,1)
	    jp1 = min(j+1,nys)
	    ip1 = min(i+1,nxs)
	    km1 = max(k-1,1)
	    kp1 = min(k+1,nzs)
	    
	    coefh = 1/dx/dy * msfm[i,j]^2
	    coefkp1 = 2/(h[k]*(h[k]+h[km1]))
	    coefk = -2/h[k]/h[km1]
	    coefkm1 = 2/(h[km1]*(h[k]+h[km1]))
	            
	    global ∂⁴xyπ2 = logop(tijk, ip1, jp1, kp1, 1/4*coefh*coefkp1, ∂⁴xyπ2)
	    global ∂⁴xyπ2 = logop(tijk, im1, jm1, kp1, 1/4*coefh*coefkp1, ∂⁴xyπ2)
		global ∂⁴xyπ2 = logop(tijk, im1, jp1, kp1, -1/4*coefh*coefkp1, ∂⁴xyπ2)
		global ∂⁴xyπ2 = logop(tijk, ip1, jm1, kp1, -1/4*coefh*coefkp1, ∂⁴xyπ2)
		
	    global ∂⁴xyπ2 = logop(tijk, ip1, jp1, k, 1/4*coefh*coefk, ∂⁴xyπ2)
	    global ∂⁴xyπ2 = logop(tijk, im1, jm1, k, 1/4*coefh*coefk, ∂⁴xyπ2)
		global ∂⁴xyπ2 = logop(tijk, im1, jp1, k, -1/4*coefh*coefk, ∂⁴xyπ2)
		global ∂⁴xyπ2 = logop(tijk, ip1, jm1, k, -1/4*coefh*coefk, ∂⁴xyπ2)
		
	    global ∂⁴xyπ2 = logop(tijk, ip1, jp1, km1, 1/4*coefh*coefkm1, ∂⁴xyπ2)
	    global ∂⁴xyπ2 = logop(tijk, im1, jm1, km1, 1/4*coefh*coefkm1, ∂⁴xyπ2)
		global ∂⁴xyπ2 = logop(tijk, im1, jp1, km1, -1/4*coefh*coefkm1, ∂⁴xyπ2)
		global ∂⁴xyπ2 = logop(tijk, ip1, jm1, km1, -1/4*coefh*coefkm1, ∂⁴xyπ2)
	end
	String("4th order: ∇⁴, ∂⁴x3π, ∂⁴y3π, ∂⁴x2π2, ∂⁴y2π2, ∂⁴xy2π, ∂⁴x2yπ, ∇²∂²π2, ∇²∂²xπ, ∇²∂²yπ, ∂⁴xyπ2")
end

# ╔═╡ 050fcd6b-4d2b-45e5-b1d5-9b02ebdd97ab
# begin
# 	fx = ∂x*vec(f3d)
# 	fy = ∂y*vec(f3d)
# 	fx∂²π = zeros(Float32, nxyz, nxyz)
# 	fy∂²π = zeros(Float32, nxyz, nxyz)
# 	fx∂²xπ = zeros(Float32, nxyz, nxyz)
# 	fy∂²xπ = zeros(Float32, nxyz, nxyz)
# 	fx∂²yπ = zeros(Float32, nxyz, nxyz)
# 	fy∂²yπ = zeros(Float32, nxyz, nxyz)
	
# 	f∇² = zeros(Float32, nxyz, nxyz)
# 	for i in eachindex(fx)
# 	    fx∂²π[:,i] = fx.*∂²π[:,i]
# 	    fy∂²π[:,i] = fy.*∂²π[:,i]
# 	    fx∂²xπ[:,i] = fx.*∂²xπ[:,i]
# 	    fy∂²xπ[:,i] = fy.*∂²xπ[:,i]
# 	    fx∂²yπ[:,i] = fx.*∂²yπ[:,i]
# 	    fy∂²yπ[:,i] = fy.*∂²yπ[:,i]
# 	    f∇²[:,i] = vec(f3d).*∇²[:,i]
# 	end
# 	f∇²∂²π = f∇²*∂²π
# 	f∇²∂²xπ = f∇²*∂²xπ
# 	f∇²∂²yπ = f∇²*∂²yπ
	
# 	∇f∇∂²π = fx∂²π .+ fy∂²π .+ f∇²∂²π
# 	∇f∇∂²xπ = fx∂²xπ .+ fy∂²xπ .+ f∇²∂²xπ
# 	∇f∇∂²yπ = fx∂²yπ .+ fy∂²yπ .+ f∇²∂²yπ
# 	surp = 1
# 	String("∇f∇ terms: ∇f∇∂²π, ∇f∇∂²xπ, ∇f∇∂²yπ")
# end


# ╔═╡ b7a71daf-0f28-4bcc-8393-05d8dbd6e179
begin
	# ∇²
	C∇² = - vec(V.∇²ϕ̄ππ[slicex, slicey,:])
	# ∇⁴
	C∇⁴ = - vec(V.ϕ̄ππ[slicex, slicey,:])
	# ∂²xπ
	C∂²xπ = vec(V.∇²ϕ̄xπ[slicex, slicey,:] 
			+ V.fx[slicex, slicey,:].*V.∂xψ̄xπ[slicex, slicey,:] 
			+ V.fy[slicex, slicey,:].*V.∂yψ̄xπ[slicex, slicey,:] 
			+ f3d.*V.∇²ψ̄xπ[slicex, slicey,:] 
			+ 2*V.ψ̄xx[slicex, slicey,:].*V.∂yyψ̄xπ[slicex, slicey,:] 
			+ 2*V.ψ̄yy[slicex, slicey,:].*V.∂xxψ̄xπ[slicex, slicey,:] 
			- 4*V.ψ̄xy[slicex, slicey,:].*V.∂xyψ̄xπ[slicex, slicey,:])

	# ∇²∂²xπ
	C∇²∂²xπ = vec(V.ϕ̄xπ[slicex, slicey,:] + f3d.*V.ψ̄xπ[slicex, slicey,:])

	# ∂²yπ
	C∂²yπ = vec(V.∇²ϕ̄yπ[slicex, slicey,:] 
			+ V.fx[slicex, slicey,:].*V.∂xψ̄yπ[slicex, slicey,:] 
			+ V.fy[slicex, slicey,:].*V.∂yψ̄yπ[slicex, slicey,:] 
			+ f3d.*V.∇²ψ̄yπ[slicex, slicey,:] 
			+ 2*V.ψ̄xx[slicex, slicey,:].*V.∂yyψ̄yπ[slicex, slicey,:] 
			+ 2*V.ψ̄yy[slicex, slicey,:].*V.∂xxψ̄yπ[slicex, slicey,:] 
			- 4*V.ψ̄xy[slicex, slicey,:].*V.∂xyψ̄yπ[slicex, slicey,:])

	# ∇²∂²yπ
	C∇²∂²yπ = vec(V.ϕ̄yπ[slicex, slicey,:] + f3d.*V.ψ̄yπ[slicex, slicey,:])

	# ∂²π
	C∂²π = vec(-V.fx[slicex, slicey,:].*V.ζx[slicex, slicey,:] 
			- V.fy[slicex, slicey,:].*V.ζy[slicex, slicey,:] 
			- f3d.*V.∇²ζ[slicex, slicey,:] 
			- 2*V.ψ̄xx[slicex, slicey,:].*V.ζyy[slicex, slicey,:] 
			- 2*V.ψ̄yy[slicex, slicey,:].*V.ζxx[slicex, slicey,:] 
			+ 4*V.ψ̄xy[slicex, slicey,:].*V.ζxy[slicex, slicey,:])	

	# ∂³xπ2
	C∂³xπ2 = - vec(V.fx[slicex, slicey,:] .* V.ζ[slicex, slicey,:])
	# ∂³yπ2
	C∂³yπ2 = - vec(V.fy[slicex, slicey,:] .* V.ζ[slicex, slicey,:])
	# ∇²∂²π2
	C∇²∂²π2 = vec(- f3d .* V.ζ[slicex, slicey,:])
	#∂³x2π 
	C∂³x2π = vec(V.fx[slicex, slicey,:] .* V.ψ̄xπ[slicex, slicey,:])
	#∂³y2π 
	C∂³y2π = vec(V.fy[slicex, slicey,:] .* V.ψ̄yπ[slicex, slicey,:])
	#∂³xyπ
	C∂³xyπ = vec(V.fy[slicex, slicey,:] .* V.ψ̄xπ[slicex, slicey,:] 
			   + V.fx[slicex, slicey,:] .* V.ψ̄yπ[slicex, slicey,:])

	# ∂⁴x2π2
	C∂⁴x2π2 = vec(-2*V.ψ̄yy[slicex, slicey,:] .* V.ζ[slicex, slicey,:])
	# ∂⁴y2π2
	C∂⁴y2π2 = vec(-2*V.ψ̄xx[slicex, slicey,:] .* V.ζ[slicex, slicey,:])

	# ∂⁴x3π
	C∂⁴x3π = vec(2*V.ψ̄yy[slicex, slicey,:] .* V.ψ̄xπ[slicex, slicey,:])
	# ∂⁴y3π
	C∂⁴y3π = vec(2*V.ψ̄xx[slicex, slicey,:] .* V.ψ̄yπ[slicex, slicey,:])

	# ∂⁴xy2π
	C∂⁴xy2π = vec(2*V.ψ̄xx[slicex, slicey,:].*V.ψ̄xπ[slicex, slicey,:] 
				- 4*V.ψ̄xy[slicex, slicey,:].*V.ψ̄yπ[slicex, slicey,:])
	
	# ∂⁴x2yπ
	C∂⁴x2yπ = vec(2*V.ψ̄yy[slicex, slicey,:].*V.ψ̄yπ[slicex, slicey,:] 
				- 4*V.ψ̄xy[slicex, slicey,:].*V.ψ̄xπ[slicex, slicey,:])
	
	# ∂⁴xyπ2
	C∂⁴xyπ2 = vec(4*V.ψ̄xy[slicex, slicey,:].*V.ζ[slicex, slicey,:])

	####################################
	P∇² = zeros(Float32, nxyz, nxyz)
	P∇⁴ = zeros(Float32, nxyz, nxyz)
	P∂²xπ = zeros(Float32, nxyz, nxyz)
	P∇²∂²xπ = zeros(Float32, nxyz, nxyz)
	P∂²yπ = zeros(Float32, nxyz, nxyz)
	P∇²∂²yπ = zeros(Float32, nxyz, nxyz)
	P∂²π = zeros(Float32, nxyz, nxyz)
	P∂³xπ2 = zeros(Float32, nxyz, nxyz)
	P∂³yπ2 = zeros(Float32, nxyz, nxyz)
	P∇²∂²π2 = zeros(Float32, nxyz, nxyz)
	P∂³x2π  = zeros(Float32, nxyz, nxyz)
	P∂³y2π  = zeros(Float32, nxyz, nxyz)
	P∂³xyπ = zeros(Float32, nxyz, nxyz)
	P∂⁴x2π2 = zeros(Float32, nxyz, nxyz)
	P∂⁴y2π2 = zeros(Float32, nxyz, nxyz)
	P∂⁴x3π = zeros(Float32, nxyz, nxyz)
	P∂⁴y3π = zeros(Float32, nxyz, nxyz)
	P∂⁴xy2π = zeros(Float32, nxyz, nxyz)
	P∂⁴x2yπ = zeros(Float32, nxyz, nxyz)
	P∂⁴xyπ2 = zeros(Float32, nxyz, nxyz)

	for i in 1:nxyz
		P∇²[:,i] = C∇².* ∇²[:,i]
		P∇⁴[:,i] = C∇⁴.* ∇⁴[:,i]
		P∂²xπ[:,i] = C∂²xπ.* ∂²xπ[:,i]
		P∇²∂²xπ[:,i] = C∇²∂²xπ.* ∇²∂²xπ[:,i]
		P∂²yπ[:,i] = C∂²yπ.* ∂²yπ[:,i]
		P∇²∂²yπ[:,i] = C∇²∂²yπ.* ∇²∂²yπ[:,i]
		P∂²π[:,i] = C∂²π.* ∂²π[:,i]
		P∂³xπ2[:,i] = C∂³xπ2.* ∂³xπ2[:,i]
		P∂³yπ2[:,i] = C∂³yπ2.* ∂³yπ2[:,i]
		P∇²∂²π2[:,i] = C∇²∂²π2.* ∇²∂²π2[:,i]
		P∂³x2π[:,i] = C∂³x2π.* ∂³x2π[:,i]
		P∂³y2π[:,i] = C∂³y2π.* ∂³y2π[:,i] 
		P∂³xyπ[:,i] = C∂³xyπ.* ∂³xyπ[:,i]
		P∂⁴x2π2[:,i] = C∂⁴x2π2.* ∂⁴x2π2[:,i]
		P∂⁴y2π2[:,i] = C∂⁴y2π2.* ∂⁴y2π2[:,i]
		P∂⁴x3π[:,i] = C∂⁴x3π.* ∂⁴x3π[:,i]
		P∂⁴y3π[:,i] = C∂⁴y3π.* ∂⁴y3π[:,i]
		P∂⁴xy2π[:,i] = C∂⁴xy2π.* ∂⁴xy2π[:,i]
		P∂⁴x2yπ[:,i] = C∂⁴x2yπ.* ∂⁴x2yπ[:,i]
		P∂⁴xyπ2[:,i] = C∂⁴xyπ2.* ∂⁴xyπ2[:,i]
	end

	φP = (P∇² + P∇⁴ + P∂²xπ + P∇²∂²xπ + P∂²yπ + P∇²∂²yπ + P∂²π + P∂³xπ2 + P∂³yπ2 + P∇²∂²π2 + P∂³x2π + P∂³y2π  + P∂³xyπ + P∂⁴x2π2 + P∂⁴y2π2 + P∂⁴x3π + P∂⁴y3π + P∂⁴xy2π + P∂⁴x2yπ + P∂⁴xyπ2)
	
	String("Operator (TϕRψ + TψRΦ) on q̂")
end

# ╔═╡ 5ac112e9-707e-44f0-aca2-4f6d3816ee90
# begin
# 	# A1 = ∇²Φ̄ππ ∇² q̂ 
# 	# A2 = Φ̄ππ ∇⁴ q̂
# 	# A3 = ∇²Φ̄xπ ∂xπ q̂ 
# 	# A4 = Φ̄xπ ∇²∂xπ q̂ 
# 	# A5 = ∇²Φ̄yπ ∂yπ q̂ 
# 	# A6 = Φ̄yπ ∇²∂yπ q̂ 
# 	# TϕRψ = -A * (A1+A2-A3-A4-A5-A6)
	
# 	A1 = zeros(Float32, nxyz, nxyz)
# 	CA1 = vec(V.∇²ϕ̄ππ[slicex, slicey,:])
	
# 	A2 = zeros(Float32, nxyz, nxyz)
# 	CA2 = vec(V.ϕ̄ππ[slicex, slicey,:])
	
# 	A3 = zeros(Float32, nxyz, nxyz)
# 	CA3 = vec(V.∇²ϕ̄xπ[slicex, slicey,:])
	
# 	A4 = zeros(Float32, nxyz, nxyz)
# 	CA4 = vec(V.ϕ̄xπ[slicex, slicey,:])
	
# 	A5 = zeros(Float32, nxyz, nxyz)
# 	CA5 = vec(V.∇²ϕ̄yπ[slicex, slicey,:])
	
# 	A6 = zeros(Float32, nxyz, nxyz)
# 	CA6 = vec(V.ϕ̄yπ[slicex, slicey,:])
	
# 	for i in 1:nxyz
# 	    A1[:,i] = CA1.*∇²[:,i]
# 	    A2[:,i] = CA2.*∇⁴[:,i]
# 	    A3[:,i] = CA3.*∂²xπ[:,i]
# 	    A4[:,i] = CA4.*∇²∂²xπ[:,i]
# 	    A5[:,i] = CA5.*∂²yπ[:,i]
# 	    A6[:,i] = CA6.*∇²∂²yπ[:,i]
# 	end
	
# 	TϕRψ_noA = - (A1+A2-A3-A4-A5-A6)
# 	TϕRψ = zeros(Float32, nxyz, nxyz)
# 	vecA = vec(A_ext)
# 	for i in 1:nxyz
# 	    TϕRψ[:,i] = vecA.*TϕRψ_noA[:,i]
# 	end

# 	String("TϕRψ")
	
# 	                                    # B1 = + ∇(f∇)(f+∇²Ψ̄) [∂ππ⋅q̂]
# 	# B2 = + (f+∇²Ψ̄) [∇(f∇)∂ππ⋅q̂]
# 	                                    # B3 = - ∇(f∇)Ψ̄xπ [∂xπ⋅q̂]
# 	# B4 = - Ψ̄xπ [∇(f∇)∂xπ⋅q̂]
# 	                                    # B5 = - ∇(f∇)Ψ̄yπ [∂yπ⋅q̂]
# 	# B6 = - Ψ̄yπ [∇(f∇)∂yπ⋅q̂]
	
# 	                                    # B7 = + Ψ̄xx∂yy(f+∇²Ψ̄) [∂ππ⋅q̂]
# 	# B8 = + Ψ̄xx(f+∇²Ψ̄) [∂yy∂ππ⋅q̂]
# 	                                    # B9 = - Ψ̄xx∂yyΨ̄xπ [∂xπ⋅q̂]
# 	# B10 = - Ψ̄xxΨ̄xπ [∂yy∂xπ⋅q̂]
# 	                                    # B11 = - Ψ̄xx∂yyΨ̄yπ [∂yπ⋅q̂]
# 	# B12 = - Ψ̄xxΨ̄yπ [∂yy∂yπ⋅q̂]
	
# 	                                    # B13 = + Ψ̄yy∂xx(f+∇²Ψ̄) [∂ππ⋅q̂]
# 	# B14 = + Ψ̄yy⋅(f+∇²Ψ̄) [∂xx∂ππ⋅q̂]
# 	                                    # B15 = - Ψ̄yy∂xxΨ̄xπ [∂xπ⋅q̂]
# 	# B16 = - Ψ̄yyΨ̄xπ [∂xx∂xπ⋅q̂]
# 	                                    # B17 = - Ψ̄yy∂xxΨ̄yπ [∂yπ⋅q̂]
# 	# B18 = - Ψ̄yyΨ̄yπ [∂xx∂yπ⋅q̂]
	
# 	                                    # B19 = - 2Ψ̄(f+∇²Ψ̄) [∂ππ⋅q̂]
# 	                                    # B20 = + 2Ψ̄⋅Ψ̄xπ [∂xπ⋅q̂]
# 	                                    # B21 = + 2Ψ̄⋅Ψ̄yπ [∂yπ⋅q̂]
# 	# B = A*sum(B)
# 	# Bxπ = - ∇(f∇)Ψ̄xπ - Ψ̄xx∂yyΨ̄xπ - Ψ̄yy∂xxΨ̄xπ + 2Ψ̄⋅Ψ̄xπ # B3,9,15,20
# 	# Byπ = - ∇(f∇)Ψ̄yπ - Ψ̄xx∂yyΨ̄yπ - Ψ̄yy∂xxΨ̄yπ + 2Ψ̄⋅Ψ̄yπ # B5, 11, 17, 21
# 	# Bππ = ∇(f∇)(f+∇²Ψ̄) + Ψ̄xx∂yy(f+∇²Ψ̄) + Ψ̄yy∂xx(f+∇²Ψ̄) - 2Ψ̄(f+∇²Ψ̄) #B1,7,13,19
	
# 	B2 = zeros(Float32, nxyz, nxyz)
# 	CB2 = vec(V.ζ[slicex, slicey,:])
	
# 	B4 = zeros(Float32, nxyz, nxyz)
# 	CB4 = vec(-V.ψ̄xπ[slicex, slicey,:])
	
# 	B6 = zeros(Float32, nxyz, nxyz)
# 	CB6 = vec(-V.ψ̄yπ[slicex, slicey,:])
	
# 	B8 = zeros(Float32, nxyz, nxyz)
# 	CB8 = vec( (V.ψ̄xx .* V.ζ)[slicex, slicey,:])
	
# 	B10 = zeros(Float32, nxyz, nxyz)
# 	CB10 = vec(-(V.ψ̄xx .* V.ψ̄xπ)[slicex, slicey,:])
	
# 	B12 = zeros(Float32, nxyz, nxyz)
# 	CB12 = vec(-(V.ψ̄xx .* V.ψ̄yπ)[slicex, slicey,:])
	
# 	B14 = zeros(Float32, nxyz, nxyz)
# 	CB14 = vec((V.ψ̄yy .* V.ζ)[slicex, slicey,:])
	
# 	B16 = zeros(Float32, nxyz, nxyz)
# 	CB16 = vec(-(V.ψ̄yy .* V.ψ̄xπ)[slicex, slicey,:])
	
# 	B18 = zeros(Float32, nxyz, nxyz)
# 	CB18 = vec(-(V.ψ̄yy .* V.ψ̄yπ)[slicex, slicey,:])
	    
# 	# Bxπ = - ∇(f∇)Ψ̄xπ - Ψ̄xx∂yyΨ̄xπ - Ψ̄yy∂xxΨ̄xπ + 2Ψ̄⋅Ψ̄xπ # B3,9,15,20
# 	Bxπ = zeros(Float32, nxyz, nxyz)
# 	CBxπ = vec( (- V.∇f∇ψ̄xπ .* V.ψ̄xπ - V.ψ̄xx .* V.ψ̄xπyy - V.ψ̄yy .* V.ψ̄xπxx + 2*V.ψ .* V.ψ̄xπ)[slicex, slicey,:])
	
# 	# Byπ = - ∇(f∇)Ψ̄yπ - Ψ̄xx∂yyΨ̄yπ - Ψ̄yy∂xxΨ̄yπ + 2Ψ̄⋅Ψ̄yπ # B5, 11, 17, 21
# 	Byπ = zeros(Float32, nxyz, nxyz)
# 	CByπ = vec((- V.∇f∇ψ̄yπ .* V.ψ̄yπ - V.ψ̄xx .* V.ψ̄yπyy - V.ψ̄yy .* V.ψ̄yπxx + 2*V.ψ .* V.ψ̄yπ)[slicex, slicey,:])
	
# 	# Bππ = ∇(f∇)(f+∇²Ψ̄) + Ψ̄xx∂yy(f+∇²Ψ̄) + Ψ̄yy∂xx(f+∇²Ψ̄) - 2Ψ̄(f+∇²Ψ̄) #B1,7,13,19
# 	Bππ = zeros(Float32, nxyz, nxyz)
# 	CBππ = vec((V.∇f∇ζ + V.ψ̄xx .* V.ζyy + V.ψ̄yy .* V.ζxx - 2*V.ψ .* V.ζ)[slicex, slicey,:])
	
# 	for i in 1:nxyz
# 	    B2[:,i] = CB2.* ∇f∇∂²π[:,i]
# 	    B4[:,i] = CB4.* ∇f∇∂²xπ[:,i]
# 	    B6[:,i] = CB6.* ∇f∇∂²yπ[:,i]
# 	    B8[:,i] = CB8.* ∂⁴y2π2[:,i]
# 	    B10[:,i] = CB10.* ∂⁴xy2π[:,i]
# 	    B12[:,i] = CB12.* ∂⁴y3π[:,i]
# 	    B14[:,i] = CB14.* ∂⁴x2π2[:,i]
# 	    B16[:,i] = CB16.* ∂⁴x3π[:,i]
# 	    B18[:,i] = CB18.* ∂⁴x2yπ[:,i]
# 	    Bxπ[:,i] = CBxπ.* ∂²xπ[:,i]
# 	    Byπ[:,i] = CByπ.* ∂²yπ[:,i]
# 	    Bππ[:,i] = CBππ.* ∂²π[:,i]
# 	end
	
# 	TψRΦ_noA = B2+B4+B6+B8+B10+B12+B14+B16+B18+Bxπ+Byπ+Bππ
# 	TψRΦ = zeros(Float32, nxyz, nxyz)

# 	for i in 1:nxyz
# 	    TψRΦ[:,i] = vecA.*TψRΦ_noA[:,i]
# 	end
# 	String("TψRΦ")
	
# 	φP = TϕRψ - TψRΦ;

# 	String("TϕRψ - TψRΦ")
# end

# ╔═╡ 154e14a6-b8cd-4506-a2e3-ee4cd2a0b539
npzwrite("varphiP_mat.npy", φP)

# ╔═╡ Cell order:
# ╠═087626d2-dd26-11ec-2371-cd85f62600fe
# ╠═239a0ed4-d7ef-4429-a556-daca03b0b012
# ╠═e8b841c7-5a9e-45c7-89a5-e9039d4e09f6
# ╠═407f03b8-8bcf-4425-84d3-91b5e542c41a
# ╠═c39217a8-f91f-4ff4-940c-a93068238856
# ╠═509e5ee1-c3c7-47d5-ae74-8de5a40c43d9
# ╠═0dc4dbd5-a367-4c35-837d-3fc65205885e
# ╟─3428bdd9-7405-4167-b760-ad5ef4a23a15
# ╟─37f870cf-a9a7-42a0-8239-72b21270461c
# ╟─59586306-5cd6-4fc5-bb19-62d9153c077a
# ╠═050fcd6b-4d2b-45e5-b1d5-9b02ebdd97ab
# ╠═b7a71daf-0f28-4bcc-8393-05d8dbd6e179
# ╟─5ac112e9-707e-44f0-aca2-4f6d3816ee90
# ╠═154e14a6-b8cd-4506-a2e3-ee4cd2a0b539
