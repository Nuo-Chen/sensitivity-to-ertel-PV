include("gradient1.jl")
function adj_Rϕ(psi_bs, target, dx, dy, 
                h, mf, A, f) 
    rphi1 = (f .+ ∇²(psi_bs,dx,dy,mf)) .* ∂²π(target,h)
    rphi2 = ∂xπ(psi_bs, dx, h, mf) .* ∂xπ(target, dx, h, mf)
    rphi3 = ∂yπ(psi_bs, dy, h, mf) .* ∂yπ(target, dy, h, mf)
    rphi = A .* (rphi1 .- rphi2 .- rphi3)
end

function adj_Rψ(phi_bs, target, dx, dy, 
                h, mf, A) 
    rpsi1 = ∂²π(phi_bs, h) .* ∇²(target, dx, dy, mf)
    rpsi2 = ∂xπ(phi_bs, dx, h, mf) .* ∂xπ(target, dx, h, mf)
    rpsi3 = ∂yπ(phi_bs, dy, h, mf) .* ∂yπ(target, dy, h, mf)
    rpsi = A .* (rpsi1 .- rpsi2 .- rpsi3)
end

function adj_Tϕ(target, dx, dy, mf) 
    tphi = -∇²(target, dx, dy, mf)
end

function adj_Tψ(psi_bs, phi_bs, target, dx, dy, 
                mf, f) 
    ∂ψbs∂x, ∂ψbs∂y = ∇(psi_bs, dx, dy, mf)
    tpsi1 = ∇div(f .* ∂ψbs∂x, f .* ∂ψbs∂y, dx, dy, mf)
    tpsi2 = ∂²x(psi_bs, dx, mf) .* ∂²y(target, dy, mf)
    tpsi3 = ∂²y(phi_bs, dx, mf) .* ∂²x(target, dy, mf)
    tpsi4 = 2 * psi_bs
    tpsi = tpsi1 .+ tpsi2 .+ tpsi3 .- tpsi4 
end

"""target is aq̂"""
function φB(target, psi_bs, phi_bs, h,
            dx, dy, mf, A, f ) 
    varphib1 = adj_Rϕ(psi_bs, target, dx, dy, h, mf, A,f) .* adj_Tψ(psi_bs, phi_bs, target, dx, dy, mf, f)
    varphib2 = adj_Rψ(phi_bs, target, dx, dy, h, mf, A) .* adj_Tϕ(target, dx, dy, mf)
    varphib = varphib1 .- varphib2
end

"""target is ψ̂"""
function UP(target, psi_bs, phi_bs, 
            dx, dy, mf, f ) 
    up = adj_Tϕ(target, dx, dy, mf) .- adj_Tψ(psi_bs, phi_bs, target, dx, dy, mf, f)
end

"""target is ϕ̂"""
function UB(target, psi_bs, phi_bs, h,
            dx, dy, mf, A, f ) 
    ub = adj_Rϕ(psi_bs, target, dx, dy, h, mf, A, f) .- adj_Rψ(phi_bs, target, dx, dy, h, mf, A)
end
