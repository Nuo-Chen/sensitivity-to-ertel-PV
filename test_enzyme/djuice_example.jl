using gradient
using Enzyme

function ∇²(s, dx::T, dy::T, mf::Array{T,2}) where {T<:AbstractFloat}
    dx2 = dx^2
    mf2 = mf.^2
    ret = zero(mf)
    ret[2:end-1,2:end-1] = (s[2:end-1,3:end]+s[2:end-1,1:end-2]+s[3:end,2:end-1]+s[1:end-2,2:end-1]- 4*s[2:end-1,2:end-1])/dx2
    ret = ret .* mf2
end

rhs = forcing

function cost(rhs, x)
    lhs = ∇²(x, 1, 1, ones(Int8, 1,100,100))
    return sum(sqrt.(lhs - rhs).^2)
end

∂J_∂x = zero(x)
@time autodiff(cost, Active, rhs, Duplicated(x, ∂J_∂x))
