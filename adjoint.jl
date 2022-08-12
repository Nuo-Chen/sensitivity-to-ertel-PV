include("gradient.jl")
function adj_Rϕ(psi_bs::Array{T,3}, target::Array{T,3}, dx::T, dy::T, 
                h::Array{T,3}, mf::Array{T,3}, A::Array{T,3}, f::Array{T,3}) where {T<:AbstractFloat}
    rphi1 = (f .+ ∇²(psi_bs,dx,dy,mf)) .* ∂²π(target,h)
    rphi2 = ∂xπ(psi_bs, dx, h, mf) .* ∂xπ(target, dx, h, mf)
    rphi3 = ∂yπ(psi_bs, dy, h, mf) .* ∂yπ(target, dy, h, mf)
    rphi = A .* (rphi1 .- rphi2 .- rphi3)
end

function adj_Rψ(phi_bs::Array{T,3}, target::Array{T,3}, dx::T, dy::T, 
                h::Array{T,3}, mf::Array{T,3}, A::Array{T,3}) where {T<:AbstractFloat}
    rpsi1 = ∂²π(phi_bs, h) .* ∇²(target, dx, dy, mf)
    rpsi2 = ∂xπ(phi_bs, dx, h, mf) .* ∂xπ(target, dx, h, mf)
    rpsi3 = ∂yπ(phi_bs, dy, h, mf) .* ∂yπ(target, dy, h, mf)
    rpsi = A .* (rpsi1 .- rpsi2 .- rpsi3)
end

function adj_Tϕ(target, dx::T, dy::T, mf::Array{T,3}) where {T<:AbstractFloat}
    tphi = -∇²(target, dx, dy, mf)
end

function adj_Tψ(psi_bs::Array{T,3}, phi_bs::Array{T,3}, target::Array{T,3}, dx::T, dy::T, 
                mf::Array{T,3}, f::Array{T,3}) where {T<:AbstractFloat}
    ∂ψbs∂x, ∂ψbs∂y = ∇(psi_bs, dx, dy, mf)
    tpsi1 = ∇div(f .* ∂ψbs∂x, f .* ∂ψbs∂y, dx, dy, mf)
    tpsi2 = ∂²x(psi_bs, dx, mf) .* ∂²y(target, dy, mf)
    tpsi3 = ∂²y(phi_bs, dy, mf) .* ∂²x(target, dx, mf)
    tpsi4 = 2 * psi_bs
    tpsi = tpsi1 .+ tpsi2 .+ tpsi3 .- tpsi4 
end

"""target is q̂"""
function φP(target, psi_bs::Array{T,3}, phi_bs::Array{T,3}, h::Array{T,3},
            dx::T, dy::T, mf::Array{T,3}, A::Array{T,3}, f::Array{T,3} ) where {T<:AbstractFloat}
    varphip1 = adj_Tϕ(target, dx, dy, mf) .* adj_Rψ(phi_bs, target, dx, dy, h, mf, A) 
    varphip2 = adj_Tψ(psi_bs, phi_bs, target, dx, dy, mf, f) .* adj_Rϕ(psi_bs, target, dx, dy, h, mf, A, f)
    varphip = varphip1 .- varphip2
end

"""target is aq̂"""
function φB(target::Array{T,3}, psi_bs::Array{T,3}, phi_bs::Array{T,3}, h::Array{T,3},
            dx::T, dy::T, mf::Array{T,3}, A::Array{T,3}, f::Array{T,3} ) where {T<:AbstractFloat}
    varphib1 = adj_Rϕ(psi_bs, target, dx, dy, h, mf, A,f) .* adj_Tψ(psi_bs, phi_bs, target, dx, dy, mf, f)
    varphib2 = adj_Rψ(phi_bs, target, dx, dy, h, mf, A) .* adj_Tϕ(target, dx, dy, mf)
    varphib = varphib1 .- varphib2
end

"""target is ψ̂"""
function UP(target::Array{T,3}, psi_bs::Array{T,3}, phi_bs::Array{T,3}, 
            dx::T, dy::T, mf::Array{T,3}, f::Array{T,3} ) where {T<:AbstractFloat}
    up = adj_Tϕ(target, dx, dy, mf) .- adj_Tψ(psi_bs, phi_bs, target, dx, dy, mf, f)
end

"""target is ϕ̂"""
function UB(target::Array{T,3}, psi_bs::Array{T,3}, phi_bs::Array{T,3}, h::Array{T,3},
            dx::T, dy::T, mf::Array{T,3}, A::Array{T,3}, f::Array{T,3} ) where {T<:AbstractFloat}
    ub = adj_Rϕ(psi_bs, target, dx, dy, h, mf, A, f) .- adj_Rψ(phi_bs, target, dx, dy, h, mf, A)
end
