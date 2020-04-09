using GeometricRobotDynID

# All the solvers
using SCS, ECOS, MosekTools

# And additional packages
using LinearAlgebra
using Plots
gr()

# Define the measurement uncertainty
σ = Diagonal([0.119; 0.216; 1].^2)
Σ = kron(I(Int64(990/3)), σ)

# Load the given data
φ₀, A, b, Q = load_data("./data/Train_data_human.mat", Σ = Σ)
γ = 1e-2 * tr(A'*A)

Atrain, btrain, Atest, btest = split_data(A, b, percentage = 0.5)

# Initial error
rmse(A*φ₀, b)
rmse(Atrain*φ₀, btrain)
rmse(Atest*φ₀, btest)

# Set the solver
solver = Mosek.Optimizer # Works also with SCS.Optimizer, but very slow
solver = SCS.Optimizer
φ_e = estimate_parameters(Atrain, btrain, φ₀, γ, metric = EuclideanDistance(), solver = solver, Q = Q)
rmse(A*φ_e, b)
rmse(Atrain*φ_e, btrain)
rmse(Atest*φ_e, btest)

φ_p = estimate_parameters(Atrain, btrain, φ₀, γ, metric = PullbackDistance(), solver = solver, Q = Q)
rmse(A*φ_p, b)
rmse(Atrain*φ_p, btrain)
rmse(Atest*φ_p, btest)


φ_en = estimate_parameters(Atrain, btrain, φ₀, γ, metric = EntropicDistance(), solver = solver, Q = Q)
rmse(A*φ_en, b)
rmse(Atrain*φ_en, btrain)
rmse(Atest*φ_en, btest)

scatter(φ₀, label = "Initial")
scatter!(φ_e, label = "Euclidean")
scatter!(φ_p, label = "Pullback")
scatter!(φ_en, label = "Entropic")

# Estimate the torques
plot(b, xlim = (0, 300))
plot(A*φ_e-b, xlim = (0, 300))
plot!(A*φ_p-b, xlim = (0, 300))
plot!(A*φ_en-b, xlim = (0, 300))
