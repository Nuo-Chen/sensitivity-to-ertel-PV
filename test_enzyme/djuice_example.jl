# using Enzyme
using Plots
using RegularizedLeastSquares
using LinearOperators

# https://tknopp.github.io/RegularizedLeastSquares.jl/latest/gettingStarted/

n = 25
tn = 2n+1

x = -n:n
y = -n:n
xx = x' .* ones(tn)
yy = ones(tn)' .* y

d = sqrt.(xx.^2 + yy.^2)

nx = tn; ny = tn; nz=1; nxyz=tn*tn


function logop(row, m, n, l, op, a)
	col = ((l-1)*ny + n -1)*nx + m
	a[col, row] = op
	return a
end

global ∇² = zeros(Float32, nxyz, nxyz)
for k in 1:nz, j in 2:ny-1, i in 2:nx-1
    tijk = ((k-1)*ny +j-1)*nx + i
    jm1 = max(j-1,1)
    im1 = max(i-1,1)
    jp1 = min(j+1,ny)
    ip1 = min(i+1,nx)

    global ∇² = logop(tijk, i, jm1, k, 1., ∇²)
    global ∇² = logop(tijk, im1, j, k, 1., ∇²)
    global ∇² = logop(tijk, i, j, k, -4., ∇²)
    global ∇² = logop(tijk, i, jp1, k, 1., ∇²)
    global ∇² = logop(tijk, im1, j, k, 1., ∇²)
end

rhs = exp.(-(d.-4).^2/2) 
# A = LinearOperator(∇²)
vrhs= Float32.(vec(rhs))


reg = Regularization("L2", 0.01; shape=(nxyz))
solver = createLinearSolver("cgnr", ∇², iterations=100, regMatrix=reg)

x_approx = solve(solver,vrhs)
ima = reshape(x_approx, tn, tn)

ly = @layout [a b]
p1 = contour(rhs, title = "ζ")
p2 = contour(ima, title = "ψ")

plot(p1, p2, layout=ly)
plot!(size=(600,300))
