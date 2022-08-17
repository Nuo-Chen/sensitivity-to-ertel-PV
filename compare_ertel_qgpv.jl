### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ f7aff232-1db8-11ed-232f-ffe9bb113be2
begin
    import Pkg
   	Pkg.activate(".")
	using NPZ
	using Plots
end

# ╔═╡ 3e27fd70-99a1-403b-8d66-318aae5823f5
begin
	qg = npzread("../2021041515z_inverted_20lev.npz")
	aqgpv = reshape(qg["aq"], (20,143,209))
	aqgpv_xyz = permutedims(aqgpv,[3,2,1])

	aepv = npzread("sens2EPV.npy")
	String("load data")
end

# ╔═╡ f119d432-be44-46f7-8e7d-e56d9a759e56
size(aepv)

# ╔═╡ b11aa029-b3de-42ad-a86d-3b034e5742d6
begin
	# # 5 levels at 1000, 900, 500, 200, 50
	idx = [1,2,3, 11, 17, 19,20]
	itv = 2
	slicex = 45:itv:130
	slicey = 40:itv:110

	k = 5
	ly = @layout [a b]
	p1 = contour(aqgpv_xyz[slicex, slicey, idx[k]], title = "q̂₉")
	p2 = contour(aepv[:,:,k], title = "q̂ₑ")
	
	plot(p1, p2, layout=ly)
	plot!(size=(600,300))
end

# ╔═╡ Cell order:
# ╠═f7aff232-1db8-11ed-232f-ffe9bb113be2
# ╠═3e27fd70-99a1-403b-8d66-318aae5823f5
# ╠═f119d432-be44-46f7-8e7d-e56d9a759e56
# ╠═b11aa029-b3de-42ad-a86d-3b034e5742d6
