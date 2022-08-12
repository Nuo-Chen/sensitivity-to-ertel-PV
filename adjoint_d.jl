include("gradient_d.jl")
function adj_Rϕ(psi_bs, target, dx, dy, dx2, dy2,
                h, mf, mf2, A, f, i, j, k) 
    rphi1 = (f[i,j] + ∇²(psi_bs,dx2,dy2,mf2,i,j,k)) * ∂²π(target,h, i, j, k)
    rphi2 = ∂xπ(psi_bs, dx, h, mf, i, j, k) * ∂xπ(target, dx, h, mf, i, j, k)
    rphi3 = ∂yπ(psi_bs, dy, h, mf, i, j, k) * ∂yπ(target, dy, h, mf, i, j, k)
    rphi = A[k] * (rphi1 - rphi2 - rphi3)
end

function adj_Rψ(target, phi_bs, dx, dy, dx2, dy2, 
                h, mf, mf2, A, i, j, k) 
    rpsi1 = ∂²π(phi_bs, h, i, j, k) * ∇²(target, dx2, dy2, mf2, i, j, k)
    rpsi2 = ∂xπ(phi_bs, dx, h, mf, i, j, k) * ∂xπ(target, dx, h, mf, i, j, k)
    rpsi3 = ∂yπ(phi_bs, dy, h, mf, i, j, k) * ∂yπ(target, dy, h, mf, i, j, k)
    rpsi = A[k] * (rpsi1 - rpsi2 - rpsi3)
end

function adj_Tϕ(target, dx2, dy2, mf2, i, j, k) 
    tphi = -∇²(target, dx2, dy2, mf2, i, j, k)
end

function adj_Tψ(psi_bs, phi_bs, target, dx, dy, dx2, dy2,
                mf, mf2, f, i, j, k) 
    tpsi1a = f[i,j] * ∇²(psi_bs, dx2, dy2, mf2, i, j, k)
    tpsi1b = diff2nd(f[i,j+1], f[i,j-1], dy, mf[i,j]) * ∂y(psi_bs, dy, mf, i, j, k)
    tpsi1 = tpsi1a + tpsi1b   
    tpsi2 = ∂²x(psi_bs, dx2, mf2, i, j, k) * ∂²y(target, dy2, mf2, i, j, k)
    tpsi3 = ∂²y(phi_bs, dx2, mf2, i, j, k) * ∂²x(target, dy2, mf2, i, j, k)
    tpsi4 = 2 * psi_bs[i,j,k]
    tpsi = tpsi1 + tpsi2 + tpsi3 - tpsi4 
end

"""target is q̂"""
function φP(target, i, j, k, A, dx, dx2, dx3, dy, dy2, dy3, mf, mf2, mf3, h, f, fx, fy,
	ϕ̄ππ, ϕ̄xπ , ϕ̄yπ , ∇²ϕ̄ππ , ∇²ϕ̄xπ, ∇²ϕ̄yπ, ζ, ψ̄xπ, ψ̄yπ, ψ̄xx, ψ̄yy, ∇f∇ζ, ∇f∇ψ̄xπ , ∇f∇ψ̄yπ, ζxx , ζyy, ψ̄xπxx, ψ̄xπyy, ψ̄yπxx , ψ̄yπyy, ψ)
	# mf is 2d (nx,ny)
	# A is 1d (nz)
	TϕRψ1 = ∇²ϕ̄ππ[i,j,k] * ∇²(target, dx2, dy2, mf2, i,j,k) + ϕ̄ππ[i,j,k] * ∇⁴(target, dx2, dy2, mf2, i,j,k) 
	
	target_∇f∇∂xπ, target_∇²xπ = ∇f∇∂xπ(target, dx, dx2, dy, h, mf, mf2, f, fx, fy, i,j,k)
	target_∇f∇∂yπ, target_∇²yπ = ∇f∇∂yπ(target, dx, dy, dy2, h, mf, mf2, f, fx, fy, i,j,k)
	target_xπ = ∂xπ(target, dx, h, mf, i,j,k) 
	target_yπ = ∂yπ(target, dy, h, mf, i,j,k) 
	target_ππ = ∂²π(target, h, i,j,k) 
	
	TϕRψ2 = - ∇²ϕ̄xπ[i,j,k]* target_xπ - ϕ̄xπ[i,j,k]* target_∇²xπ 
	TϕRψ3 = - ∇²ϕ̄yπ[i,j,k]* target_yπ - ϕ̄yπ[i,j,k]* target_∇²yπ

	TϕRψ = - A[k] * (TϕRψ1+TϕRψ2+TϕRψ3)

	# #################
	TψRϕ1 = (∇f∇ζ[i,j,k] * target_ππ
		+ ζ[i,j,k] * ∇f∇∂ππ(target,dx, dx2, dy, h, mf, mf2, f, fx, fy, i,j,k) 
		- ∇f∇ψ̄xπ[i,j,k] * target_xπ
		- ψ̄xπ[i,j,k] * target_∇f∇∂xπ 
		- ∇f∇ψ̄yπ[i,j,k] * target_yπ
		- ψ̄yπ[i,j,k] * target_∇f∇∂yπ)
	
	TψRϕ2 = (ψ̄xx[i,j,k] * ζyy[i,j,k] * target_ππ
			+ ψ̄xx[i,j,k] * ζ[i,j,k] * ∂⁴y2π2(target, dy2, h, mf2, i,j,k)
			- ψ̄xx[i,j,k] * ψ̄xπyy[i,j,k] * target_xπ
			- ψ̄xx[i,j,k] * ψ̄xπ[i,j,k] * ∂⁴xy2π(target, dx, dy2, h, mf3, i,j,k)
			- ψ̄xx[i,j,k] * ψ̄yπyy[i,j,k] * target_yπ
			- ψ̄xx[i,j,k] * ψ̄yπ[i,j,k] * ∂⁴y3π(target, dy3, h, mf3, i,j,k))
	
	TψRϕ3 = (ψ̄yy[i,j,k] * ζxx[i,j,k] * target_ππ
			+ ψ̄yy[i,j,k] * ζ[i,j,k] * ∂⁴x2π2(target, dx2, h, mf2, i,j,k)
			- ψ̄yy[i,j,k] * ψ̄xπxx[i,j,k] * target_xπ
			- ψ̄yy[i,j,k] * ψ̄xπ[i,j,k] * ∂⁴x3π(target, dx3, h, mf3, i,j,k)
			- ψ̄yy[i,j,k] * ψ̄yπxx[i,j,k] * target_yπ
			- ψ̄yy[i,j,k] * ψ̄yπ[i,j,k] * ∂⁴x2yπ(target, dx2, dy, h, mf3, i,j,k))
	
	TψRϕ4 = -2*ψ[i,j,k]*(ζ[i,j,k]*target_ππ 
					- ψ̄xπ[i,j,k]*target_xπ 
					- ψ̄yπ[i,j,k]*target_yπ)
	
	TψRϕ = A[k] * (TψRϕ1+TψRϕ2+TψRϕ3+TψRϕ4)

	ret = TϕRψ - TψRϕ
end

# """target is aq̂"""
# function φB(target, psi_bs, phi_bs, h,
#             dx, dy, mf, A, f ) 
#     varphib1 = adj_Rϕ(psi_bs, target, dx, dy, h, mf, A,f) * adj_Tψ(psi_bs, phi_bs, target, dx, dy, mf, f)
#     varphib2 = adj_Rψ(phi_bs, target, dx, dy, h, mf, A) * adj_Tϕ(target, dx, dy, mf)
#     varphib = varphib1 - varphib2
# end

"""target is ψ̂"""
function UP(target, psi_bs, phi_bs, dx, dy, dx2, dy2, mf, mf2, f , i, j, k) 
    up = adj_Tϕ(target, dx, dy, mf, i, j, k) - adj_Tψ(psi_bs, phi_bs, target, dx, dy, dx2, dy2, mf, mf2, f, i, j, k)
end

# """target is ϕ̂"""
# function UB(target, psi_bs, phi_bs, h,
#             dx, dy, mf, A, f ) 
#     ub = adj_Rϕ(psi_bs, target, dx, dy, h, mf, A, f) - adj_Rψ(phi_bs, target, dx, dy, h, mf, A)
# end
