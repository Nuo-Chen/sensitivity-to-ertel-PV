struct φPconst{T<:AbstractFloat}
	ϕ̄ππ::Array{T,3}
	ϕ̄xπ::Array{T,3}
	ϕ̄yπ::Array{T,3}
	∇²ϕ̄ππ::Array{T,3}
	∇²ϕ̄xπ::Array{T,3}
	∇²ϕ̄yπ::Array{T,3}

	ζ::Array{T,3}
	ψ̄xπ::Array{T,3}
	ψ̄yπ::Array{T,3}
	ψ̄xx::Array{T,3}
	ψ̄yy::Array{T,3}

	∇f∇ζ::Array{T,3}
	∇f∇ψ̄xπ::Array{T,3}
	∇f∇ψ̄yπ::Array{T,3}
	ζxx::Array{T,3}
 	ζyy::Array{T,3}

	ψ̄xπxx::Array{T,3}
	ψ̄xπyy::Array{T,3}
	ψ̄yπxx::Array{T,3}
	ψ̄yπyy::Array{T,3}
    ψ::Array{T,3}
end

struct gconst{T}
	fᵢ::Array{T,2}
	fⱼ::Array{T,2}

	A::Array{T,1}
	rdx::T
	rdx2::T
	rdx3::T
	rdy::T
	rdy2::T
	rdy3::T
	m::Array{T,2}
	m2::Array{T,2}
	m3::Array{T,2}
	dπ::Array{T,1}
	f::Array{T,2}
end