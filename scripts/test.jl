using Convex, MosekTools

# Generate random problem data
m = 400;  n = 20
A = randn(m, n)

# Create a (column vector) variable of size n x 1.
x = Variable(n, n)

# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# This can be done by: minimize(objective, constraints)
A
problem = minimize(norm(A * x - A, 2), [diag(x) == 0])
sparsity = []
for i in 1:10
    # Solve the problem by calling solve!
    solve!(problem, Mosek.Optimizer, warmstart = i > 1 ? true : false);
    x̂ = x.value
    idxs = abs.(x̂) .< 1e-2
    push!(sparsity, sum(idxs))
    x̂[idxs] .= 0.0
end
x.value
sparsity
# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimal value
problem.optval
x.value

norm(A*x.value - A, 2) + norm(x.value, 1)
