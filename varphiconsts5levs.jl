### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ cb296012-bf4a-11ec-3b05-432705488eaa
begin
	using NPZ#,Plots
	using Parameters	
end

# ╔═╡ 61b500b3-9a14-423f-9e4c-8592c855f368
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

# ╔═╡ 09954537-cb30-4094-96b8-01c6d02efd44
grad = ingredients("/Volumes/Scratch/nchen67/ertel_inv/gradient1.jl")

# ╔═╡ 945e1189-c20c-4f8a-b78d-4e9272c346fe
begin
	nz = 5
	idx = [1,3, 11, 17, 20]
	# all are in basic states
	phi = npzread("phi_bs20lev.npy") 
	psi = npzread("psi_bs20lev.npy")

	# python 3d files in [k,j,i]
	frompython = npzread("forjulia.npz")
	rdx = Float64.(frompython["rdx"])
	rdy = Float64.(frompython["rdy"])
	A_py = frompython["A"]
	msfm_py = frompython["msfm"]
	dπ_py = frompython["dexner"]
	f_py = frompython["f"]

	# permute to 3d [i,j,k]
	ψ = permutedims(psi,[3,2,1])
	ϕ = permutedims(phi,[3,2,1])
	msfm_3d = repeat(permutedims(msfm_py,[2,1]), outer = [1,1, nz])
	f_3d = repeat(permutedims(f_py,[2,1]), outer = [1,1, nz])
	A_3d = permutedims(repeat(A_py, outer=[1,209,143]), [2, 3, 1])
	

	# extracts to 5 levels
	ψ5 = ψ[:,:,idx]
	ϕ5 = ϕ[:,:,idx]
	A5 = A_py[idx]
	dπ5 = [sum(dπ_py[1:2]), sum(dπ_py[3:10]), sum(dπ_py[11:16]), sum(dπ_py[17:end])]
	dπ_3d = permutedims(repeat(dπ5, outer=[1,209,143]), [2, 3, 1])
	
	msfm_2d = permutedims(msfm_py,[2,1])
	f_2d = permutedims(f_py,[2,1])

	nx, ny, nzz = size(ϕ);
end

# ╔═╡ c8ef101d-520f-409e-87d9-68cfd55d01f3
dπ_py

# ╔═╡ 856e1504-827d-47a3-b031-cf89ec346829
begin
	# constants for TϕRψ
	ϕ̄ππ = grad.∂²π(ϕ5, dπ_3d);
	ϕ̄xπ = grad.∂xπ(ϕ5, rdx, dπ_3d, msfm_3d);
	ϕ̄yπ = grad.∂yπ(ϕ5, rdy, dπ_3d, msfm_3d);
	∇²ϕ̄ππ = grad.∇²(ϕ̄ππ, rdx, rdy, msfm_3d);
	∇²ϕ̄xπ = grad.∇²(ϕ̄xπ, rdx, rdy, msfm_3d);
	∇²ϕ̄yπ = grad.∇²(ϕ̄yπ, rdx, rdy, msfm_3d) 
	print("constants for TϕRψ")
end

# ╔═╡ befba63a-994f-4f28-8514-b88ac6d91ddb
begin
	m2_3d = msfm_3d.^2
	ζ = grad.∇²(ψ5, rdx, rdy, msfm_3d) .+ f_3d
	ψ̄xπ = grad.∂xπ(ψ5, rdx, dπ_3d, msfm_3d)
	ψ̄yπ = grad.∂yπ(ψ5, rdy, dπ_3d, msfm_3d)
	ψ̄xx = grad.∂²x(ψ5, rdx, m2_3d)
	ψ̄yy = grad.∂²y(ψ5, rdy, m2_3d)

	∇f∇ζ = grad.∇f∇(ζ, rdx, rdy, msfm_3d, f_3d)
	∇f∇ψ̄xπ = grad.∇f∇(ψ̄xπ, rdx, rdy, msfm_3d, f_3d)
	∇f∇ψ̄yπ = grad.∇f∇(ψ̄yπ, rdx, rdy, msfm_3d, f_3d)
	# delzeta = grad.∇²(ζ, rdx, rdy, msfm_3d)
	# maximum(filter(!isinf,(∇f∇ζ - f_3d.*delzeta)./∇f∇ζ))
	# having doubts on ∇f∇::function
	ζxx = grad.∂²x(ζ, rdx, m2_3d)
 	ζyy = grad.∂²y(ζ, rdy, m2_3d)

	ψ̄xπxx = grad.∂²x(ψ̄xπ, rdx, m2_3d)
	ψ̄xπyy = grad.∂²y(ψ̄xπ, rdy, m2_3d)
	ψ̄yπxx = grad.∂²x(ψ̄yπ, rdx, m2_3d)
	ψ̄yπyy = grad.∂²y(ψ̄yπ, rdy, m2_3d)

	fᵢ, fⱼ = grad.∇(f_3d, rdx, rdy, msfm_3d)
	print("constants for TψRϕ")
	
end

# ╔═╡ 68ccb535-595e-4309-bda1-cea3df0cc9d2
size(ϕ̄ππ)

# ╔═╡ 1d39f3e0-3b9b-402e-8cec-7c6cb4bd7b3b
begin
	nϕ̄ππ = Array{Float64, 3}(undef, nx, ny, nz+2)
	nϕ̄ππ[:,:,2:end-1] = ϕ̄ππ
	nϕ̄ππ[:,:,1] = nϕ̄ππ[:,:,2]
	nϕ̄ππ[:,:,end] = nϕ̄ππ[:,:,end-1]
	
	nϕ̄xπ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nϕ̄xπ[:,:,2:end-1] = ϕ̄xπ
	nϕ̄xπ[:,:,1] = nϕ̄xπ[:,:,2]
	nϕ̄xπ[:,:,end] = nϕ̄xπ[:,:,end-1]
	
	nϕ̄yπ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nϕ̄yπ[:,:,2:end-1] = ϕ̄yπ
	nϕ̄yπ[:,:,1] = nϕ̄yπ[:,:,2]
	nϕ̄yπ[:,:,end] = nϕ̄yπ[:,:,end-1]
	
	n∇²ϕ̄ππ = Array{Float64, 3}(undef, nx, ny, nz+2)
	n∇²ϕ̄ππ[:,:,2:end-1] = ∇²ϕ̄ππ
	n∇²ϕ̄ππ[:,:,1] = n∇²ϕ̄ππ[:,:,2]
	n∇²ϕ̄ππ[:,:,end] = n∇²ϕ̄ππ[:,:,end-1]
	
	n∇²ϕ̄xπ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	n∇²ϕ̄xπ[:,:,2:end-1] = ∇²ϕ̄xπ
	n∇²ϕ̄xπ[:,:,1] = n∇²ϕ̄xπ[:,:,2]
	n∇²ϕ̄xπ[:,:,end] = n∇²ϕ̄xπ[:,:,end-1]
	
	n∇²ϕ̄yπ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	n∇²ϕ̄yπ[:,:,2:end-1] = ∇²ϕ̄yπ
	n∇²ϕ̄yπ[:,:,1] = n∇²ϕ̄yπ[:,:,2]
	n∇²ϕ̄yπ[:,:,end] = n∇²ϕ̄yπ[:,:,end-1]
	
	nζ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nζ[:,:,2:end-1] = ζ 
	nζ[:,:,1] = nζ[:,:,2]
	nζ[:,:,end] = nζ[:,:,end-1]
	
	nψ̄xπ = Array{Float64, 3}(undef, nx, ny, nz+2)
	nψ̄xπ[:,:,2:end-1] = ψ̄xπ
	nψ̄xπ[:,:,1] = nψ̄xπ[:,:,2]
	nψ̄xπ[:,:,end] = nψ̄xπ[:,:,end-1]
	
	nψ̄yπ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nψ̄yπ[:,:,2:end-1] = ψ̄yπ
	nψ̄yπ[:,:,1] = nψ̄yπ[:,:,2]
	nψ̄yπ[:,:,end] = nψ̄yπ[:,:,end-1]
	
	nψ̄xx  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nψ̄xx[:,:,2:end-1] = ψ̄xx
	nψ̄xx[:,:,1] = nψ̄xx[:,:,2]
	nψ̄xx[:,:,end] = nψ̄xx[:,:,end-1]
	
	nψ̄yy   = Array{Float64, 3}(undef, nx, ny, nz+2)
	nψ̄yy[:,:,2:end-1] = ψ̄yy
	nψ̄yy[:,:,1] = nψ̄yy[:,:,2]
	nψ̄yy[:,:,end] = nψ̄yy[:,:,end-1]
	
	n∇f∇ζ = Array{Float64, 3}(undef, nx, ny, nz+2)
	n∇f∇ζ[:,:,2:end-1] = ∇f∇ζ
	n∇f∇ζ[:,:,1] = n∇f∇ζ[:,:,2]
	n∇f∇ζ[:,:,end] = n∇f∇ζ[:,:,end-1]
	
	n∇f∇ψ̄xπ = Array{Float64, 3}(undef, nx, ny, nz+2)
	n∇f∇ψ̄xπ[:,:,2:end-1] = ∇f∇ψ̄xπ
	n∇f∇ψ̄xπ[:,:,1] = n∇f∇ψ̄xπ[:,:,2]
	n∇f∇ψ̄xπ[:,:,end] = n∇f∇ψ̄xπ[:,:,end-1]
	
	n∇f∇ψ̄yπ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	n∇f∇ψ̄yπ[:,:,2:end-1] = ∇f∇ψ̄yπ
	n∇f∇ψ̄yπ[:,:,1] = n∇f∇ψ̄yπ[:,:,2]
	n∇f∇ψ̄yπ[:,:,end] = n∇f∇ψ̄yπ[:,:,end-1]
	
	nζxx   = Array{Float64, 3}(undef, nx, ny, nz+2)
	nζxx[:,:,2:end-1] = ζxx
	nζxx[:,:,1] = nζxx[:,:,2]
	nζxx[:,:,end] = nζxx[:,:,end-1]
	
	nζyy   = Array{Float64, 3}(undef, nx, ny, nz+2)
	nζyy[:,:,2:end-1] = ζyy
	nζyy[:,:,1] = nζyy[:,:,2]
	nζyy[:,:,end] = nζyy[:,:,end-1]
	
	nψ̄xπxx  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nψ̄xπxx[:,:,2:end-1] = ψ̄xπxx
	nψ̄xπxx[:,:,1] = nψ̄xπxx[:,:,2]
	nψ̄xπxx[:,:,end] = nψ̄xπxx[:,:,end-1]
	
	nψ̄xπyy  = Array{Float64, 3}(undef, nx, ny, nz+2) 
	nψ̄xπyy[:,:,2:end-1] = ψ̄xπyy
	nψ̄xπyy[:,:,1] = nψ̄xπyy[:,:,2]
	nψ̄xπyy[:,:,end] = nψ̄xπyy[:,:,end-1]
	
	nψ̄yπxx  = Array{Float64, 3}(undef, nx, ny, nz+2) 
	nψ̄yπxx[:,:,2:end-1] = ψ̄yπxx
	nψ̄yπxx[:,:,1] = nψ̄yπxx[:,:,2]
	nψ̄yπxx[:,:,end] = nψ̄yπxx[:,:,end-1]
		
	nψ̄yπyy  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nψ̄yπyy[:,:,2:end-1] =ψ̄yπyy
	nψ̄yπyy[:,:,1] = ψ̄yπyy[:,:,2]
	nψ̄yπyy[:,:,end] = ψ̄yπyy[:,:,end-1]
	
	nψ  = Array{Float64, 3}(undef, nx, ny, nz+2)
	nψ[:,:,2:end-1] = ψ5
	nψ[:,:,1] = nψ[:,:,2]
	nψ[:,:,end] = nψ[:,:,end-1]
end

# ╔═╡ dc939159-6cec-4fa2-b3ef-8d43a53b8bf7
# npzwrite("varphi5lev_inv_consts.npz", Dict("ϕ̄ππ"=>ϕ̄ππ, "ϕ̄xπ"=>ϕ̄xπ, "ϕ̄yπ"=>ϕ̄yπ, "∇²ϕ̄ππ"=>∇²ϕ̄ππ, "∇²ϕ̄xπ"=>∇²ϕ̄xπ, "∇²ϕ̄yπ"=>∇²ϕ̄yπ, "ζ"=>ζ, "ψ̄xπ"=>ψ̄xπ, "ψ̄yπ"=>ψ̄yπ, "ψ̄xx"=>ψ̄xx, "ψ̄yy"=>ψ̄yy, "∇f∇ζ"=>∇f∇ζ, "∇f∇ψ̄xπ"=>∇f∇ψ̄xπ, "∇f∇ψ̄yπ"=>∇f∇ψ̄yπ, "ζxx"=>ζxx, "ζyy"=>ζyy, "ψ̄xπxx"=>ψ̄xπxx, "ψ̄xπyy"=>ψ̄xπyy, "ψ̄yπxx"=>ψ̄yπxx, "ψ̄yπyy"=>ψ̄yπyy, "ψ"=>ψ5))

# ╔═╡ 8f125b3f-fb96-4829-98b0-1aeb237201fa
# npzwrite("varphi5extlev_inv_consts.npz", Dict("ϕ̄ππ"=>nϕ̄ππ, "ϕ̄xπ"=>nϕ̄xπ, "ϕ̄yπ"=>nϕ̄yπ, "∇²ϕ̄ππ"=>n∇²ϕ̄ππ, "∇²ϕ̄xπ"=>n∇²ϕ̄xπ, "∇²ϕ̄yπ"=>n∇²ϕ̄yπ, "ζ"=>nζ, "ψ̄xπ"=>nψ̄xπ, "ψ̄yπ"=>nψ̄yπ, "ψ̄xx"=>nψ̄xx, "ψ̄yy"=>nψ̄yy, "∇f∇ζ"=>n∇f∇ζ, "∇f∇ψ̄xπ"=>n∇f∇ψ̄xπ, "∇f∇ψ̄yπ"=>n∇f∇ψ̄yπ, "ζxx"=>nζxx, "ζyy"=>nζyy, "ψ̄xπxx"=>nψ̄xπxx, "ψ̄xπyy"=>nψ̄xπyy, "ψ̄yπxx"=>nψ̄yπxx, "ψ̄yπyy"=>nψ̄yπyy, "ψ"=>nψ))

# ╔═╡ 76986a6f-6c90-4133-a368-3b24c2520843
# npzwrite("general5lev_consts.npz", Dict("fᵢ"=>fᵢ[:,:,1], "fⱼ"=>fⱼ[:,:,1], "A"=>A5, "rdx"=>rdx, "rdx2"=>rdx^2, "rdx3"=>rdx^3, "rdy"=>rdy, "rdy2"=>rdy^2, "rdy3"=>rdy^3, "m"=>msfm_2d, "m2"=>msfm_2d.^2, "m3"=>msfm_2d.^3, "dπ"=>dπ5, "f"=> f_2d))

# ╔═╡ 028992f1-a88a-43ca-b2f4-ddb5ed8ff2e4
# begin
# 	M = gconst(fᵢ[1,:,:], fⱼ[1,:,:], A_py, rdx, rdx^2, rdx^3, rdy, rdy^2, rdy^3, msfm_2d, msfm_2d.^2, msfm_2d.^3, dπ_py, f_2d)
# 	Cφ = φPconst(ϕ̄ππ, ϕ̄xπ, ϕ̄yπ, ∇²ϕ̄ππ, ∇²ϕ̄xπ, ∇²ϕ̄yπ, ζ, ψ̄xπ, ψ̄yπ, ψ̄xx, ψ̄yy, ∇f∇ζ, ∇f∇ψ̄xπ, ∇f∇ψ̄yπ, ζxx, ζyy, ψ̄xπxx, ψ̄xπyy, ψ̄yπxx, ψ̄yπyy)	
# end

# ╔═╡ 0eb8b998-89c1-4591-9e8f-0c391ac7dd32

## ϕ̄ππ = zeros(nx, ny, nzz)
## ϕ̄xπ = zeros(nx, ny, nzz)
## ϕ̄yπ = zeros(nx, ny, nzz)
## ∇²ϕ̄ππ = zeros(nx, ny, nzz)
## ∇²ϕ̄xπ = zeros(nx, ny, nzz)
## ∇²ϕ̄yπ = zeros(nx, ny, nzz)

## ζ = ∇²ψ̄ .+ f_xyz
# not sure about this one ∇f∇ζ = ∇f∇const(ζ, rdx, rdy, msfm)

## ψ̄xπ = zeros(nx, ny, nzz)
## ψ̄yπ = zeros(nx, ny, nzz)
## ψ̄xx = zeros(nx, ny, nzz)
## ψ̄yy = zeros(nx, ny, nzz)
# do we need this? ∇²ψ̄ = zeros(nx, ny, nzz)

## ∇f∇ζ = zeros(nx, ny, nzz)
## ∇f∇ψ̄xπ = zeros(nx, ny, nzz)
## ∇f∇ψ̄yπ = zeros(nx, ny, nzz)
## ∇f∇ψ̄xπ = ∇f∇const(ψ̄xπ, rdx, rdy, msfm)
## ∇f∇ψ̄yπ = ∇f∇const(ψ̄yπ, rdx, rdy, msfm)

## ζxx = zeros(nx, ny, nzz) 
## ζyy = zeros(nx, ny, nzz)
## ψ̄xπxx = zeros(nx, ny, nzz)
## ψ̄xπyy = zeros(nx, ny, nzz)
## ψ̄yπxx = zeros(nx, ny, nzz)
## ψ̄yπyy = zeros(nx, ny, nzz)



# ∇²target = zeros(nx, ny, nzz)
# ∇⁴target = zeros(nx, ny, nzz)
# target_xπ= zeros(nx, ny, nzz)
# target_yπ = zeros(nx, ny, nzz)
# target_ππ = zeros(nx, ny, nzz)
# target_x2π2 = zeros(nx, ny, nzz)
# target_y2π2 = zeros(nx, ny, nzz)
# target_∇²xπ = zeros(nx, ny, nzz)
# target_∇²yπ = zeros(nx, ny, nzz)
# ∇f∇target = zeros(nx, ny, nzz)
# ∇f∇∂yπtarget = zeros(nx, ny, nzz)
# ∇f∇∂xπtarget = zeros(nx, ny, nzz)
# ∇f∇∂ππtarget = zeros(nx, ny, nzz)
# target_x2π2[i,j,k] = zeros(nx, ny, nzz)
# target_y2π2[i,j,k] = zeros(nx, ny, nzz)
# target_xy2π[i,j,k] = zeros(nx, ny, nzz)
# target_x2yπ[i,j,k] = zeros(nx, ny, nzz)
# target_x3π[i,j,k] = zeros(nx, ny, nzz)
# target_y3π[i,j,k] = zeros(nx, ny, nzz)


# m2 = msfm.^2
# m3 = msfm.^3
# rdx2 = rdx^2
# rdy2 = rdy^2
# fx = zeros(nx, ny)
# fy = zeros(nx, ny)




# TϕRψ = -A_xyz.* (∇²ϕ̄ππ .* ∇²target + ϕ̄ππ .* ∇⁴target - ∇²ϕ̄xπ.*target_xπ - ϕ̄xπ.* target_∇²xπ - ∇²ϕ̄yπ.*target_yπ - ϕ̄yπ.* target_∇²yπ);

# TψRϕ1 = ∇f∇ζ*target_ππ + ζ*∇f∇∂ππtarget - ∇f∇ψ̄xπ*target_xπ - ψ̄xπ*target_∇²xπ - ∇f∇ψ̄yπ*target_yπ - ψ̄yπ*target_∇²yπ
# TψRϕ2 = ψ̄xx*ζyy*target_ππ + ψ̄xx*ζ*target_y2π2 - ψ̄xx*ψ̄xπyy*target_xπ - ψ̄xx*ψ̄xπ*target_xy2π - ψ̄xx*ψ̄yπyy*target_yπ - ψ̄xx*ψ̄yπ*target_y3π
# TψRϕ3 = ψ̄yy*ζxx*target_ππ + ψ̄yy*ζ*target_x2π2 - ψ̄yy*ψ̄yπxx*target_xπ - ψ̄yy*ψ̄xπ*target_x3π - ψ̄yy*ψ̄yπxx*target_yπ - ψ̄yy*ψ̄yπ*target_x2yπ
# TψRϕ4 = -2*ψ_bs_xyz*(ζ*target_ππ - ψ̄xπ*target_xπ - ψ̄yπ*target_yπ)
# TψRϕ = A_xyz.*(TψRϕ1+TψRϕ2+TψRϕ3+TψRϕ4)

# φP = TϕRψ - TψRϕ

# ╔═╡ 00408950-2605-4456-b05b-2b6955d46a1b


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
NPZ = "15e1cf62-19b3-5cfa-8e77-841668bca605"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"

[compat]
NPZ = "~0.4.1"
Parameters = "~0.12.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NPZ]]
deps = ["Compat", "ZipFile"]
git-tree-sha1 = "fbfb3c151b0308236d854c555b43cdd84c1e5ebf"
uuid = "15e1cf62-19b3-5cfa-8e77-841668bca605"
version = "0.4.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

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

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

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
# ╠═cb296012-bf4a-11ec-3b05-432705488eaa
# ╠═61b500b3-9a14-423f-9e4c-8592c855f368
# ╠═09954537-cb30-4094-96b8-01c6d02efd44
# ╠═945e1189-c20c-4f8a-b78d-4e9272c346fe
# ╠═c8ef101d-520f-409e-87d9-68cfd55d01f3
# ╠═856e1504-827d-47a3-b031-cf89ec346829
# ╠═befba63a-994f-4f28-8514-b88ac6d91ddb
# ╠═68ccb535-595e-4309-bda1-cea3df0cc9d2
# ╠═1d39f3e0-3b9b-402e-8cec-7c6cb4bd7b3b
# ╠═dc939159-6cec-4fa2-b3ef-8d43a53b8bf7
# ╠═8f125b3f-fb96-4829-98b0-1aeb237201fa
# ╠═76986a6f-6c90-4133-a368-3b24c2520843
# ╠═028992f1-a88a-43ca-b2f4-ddb5ed8ff2e4
# ╠═0eb8b998-89c1-4591-9e8f-0c391ac7dd32
# ╠═00408950-2605-4456-b05b-2b6955d46a1b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
