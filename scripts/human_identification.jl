using GeometricRobotDynID

# Add the Mosek solver
using MosekTools

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

Atrain, btrain, Atest, btest = split_data(A, b, percentage = 0.1)

# Initial error
rmse(A*φ₀, b) # ≈ 8.15
rmse(Atrain*φ₀, btrain) # ≈ 8.34
rmse(Atest*φ₀, btest) # ≈ 6.23

# Set the solver
solver = Mosek.Optimizer # Works also with SCS.Optimizer, but slower
φ_e = estimate_parameters(Atrain, btrain, φ₀, γ, metric = EuclideanDistance(), solver = solver, Q = Q)
rmse(A*φ_e, b) # ≈ 6.33
rmse(Atrain*φ_e, btrain) # ≈ 6-48
rmse(Atest*φ_e, btest) # 4.69

φ_p = estimate_parameters(Atrain, btrain, φ₀, γ, metric = PullbackDistance(), solver = solver, Q = Q)
rmse(A*φ_p, b) # ≈ 6.35
rmse(Atrain*φ_p, btrain) # ≈ 6.48
rmse(Atest*φ_p, btest) # ≈ 4.96

# I am not 100 % sure this works as it should
φ_en = estimate_parameters(Atrain, btrain, φ₀, γ, metric = EntropicDistance(), solver = solver, Q = Q)
rmse(A*φ_en, b) # ≈ 6.58
rmse(Atrain*φ_en, btrain) # ≈ 6.48
rmse(Atest*φ_en, btest) # ≈ 7.37

scatter(φ₀, label = "Initial")
scatter!(φ_e, label = "Euclidean")
scatter!(φ_p, label = "Pullback")
scatter!(φ_en, label = "Entropic")
