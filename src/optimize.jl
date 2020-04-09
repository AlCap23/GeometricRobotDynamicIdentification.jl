
function estimate_parameters(A, b, φ₀, γ; metric::T = PullbackDistance(), solver = SCS.Optimizer(), Σ = nothing, Q = nothing, kwargs...) where T <: GeometricRobotDynID.abstractDistanceMetric
    # Assertions
    @assert size(A, 1) == length(b) "Provide equal measurements"
    @assert size(A, 2) % 10 == 0 "Regressor dimension has to be compatible to parameter size."
    @assert γ > zero(eltype(γ)) "Gamma has to be positive definite!"

    # Initialize the varibales
    n_bodies = Int64(size(A,2)/10)

    # Create the optimization prob
    φ = Variable(size(A,2))
    P = [Semidefinite(4) for i in 1:n_bodies]
    lsqs = sumsquares(A*φ-b)
    dmetric = distance_metric(EuclideanDistance(), φ, φ₀)
    # Add some constraints
    positive_constraints = positive_definite_constraints(φ, P)
    density_constraints = isnothing(Q) ? Constraint[] : GeometricRobotDynID.density_realizability(φ, Q)


    # Solve the initial estimate
    problem = minimize(lsqs+γ*dmetric, [density_constraints...; positive_constraints...])
    solve!(problem, solver(kwargs...));

    # Compute the upper bound by evaluating the lsqs function
    ϵ = evaluate(lsqs)

    # Solve the setpoint problem
    smetric = distance_metric(metric, φ, φ₀)
    problem2 = minimize(smetric, [lsqs <= ϵ; density_constraints...; positive_constraints...])
    solve!(problem2, solver(kwargs...));

    println("Optimization Summary")
    println("Initial Least Squares with euclidean metric regularization")
    println("   Status : $(problem.status)")
    println("   Objective: $(evaluate(problem.objective))")
    println("   Least-Squares Error: $ϵ \n")
    println("Set Point Optimization with $metric regularization")
    println("   Status : $(problem2.status)")
    println("   Objective: $(evaluate(problem2.objective))")
    println("   Least-Squares Error: $(evaluate(lsqs))\n")

    return φ.value
end
