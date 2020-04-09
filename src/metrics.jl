
struct EuclideanDistance <: abstractDistanceMetric end;
struct PullbackDistance <: abstractDistanceMetric end;
struct EntropicDistance <: abstractDistanceMetric end;

function distance_metric(reg::EuclideanDistance, φ::Variable, ϕ::AbstractArray{T}) where T <: Real
    return sumsquares(φ-ϕ)
end

function distance_metric(reg::PullbackDistance, φ::Variable, ϕ::AbstractArray{T}) where T <: Real
    # Generate the const pullback metric
    pullback = zeros(T, 10,10)
    reg = []
    for (j,i) in enumerate(1:10:length(φ))
        lower = i
        upper = (i-1)+10
        pullback .= real(sqrt(pullback_metric(ϕ[lower:upper])))
        push!(reg, sumsquares(pullback*(φ[lower:upper]-ϕ[lower:upper])))
    end
    return sum(reg)
end

function distance_metric(reg::EntropicDistance, φ::Variable, ϕ::AbstractArray{T}) where T <: Real
    # Generate the const pullback metric
    p_inv = zeros(T, 4,4)
    reg = []
    for (j,i) in enumerate(1:10:length(φ))
        lower = i
        upper = (i-1)+10
        p_inv .= inv(pseudo_inertia(ϕ[lower:upper]))

        push!(reg, -logdet(pseudo_inertia(φ[lower:upper]))+tr(p_inv*pseudo_inertia(φ[lower:upper])))
        #push!(reg, entropy(p_inv*pseudo_inertia(φ[lower:upper])))
    end
    return sum(reg)
end

function density_realizability(φ::Variable, Q::AbstractArray)
    constraints = Constraint[]
    for (j, i) in enumerate(1:10:length(φ))
        lower = i
        upper = (i-1)+10
        push!(constraints, tr(Q[j]*pseudo_inertia(φ[lower:upper]))  >= 0)
    end
    return constraints
end

function positive_definite_constraints(φ::Variable, P::AbstractArray{Variable})
    constraints = Constraint[]
    for (j, i) in enumerate(1:10:length(φ))
        lower = i
        upper = (i-1)+10
        push!(constraints, P[j] == pseudo_inertia(φ[lower:upper]))
    end
    return constraints
end


rmse(x::AbstractArray, y::AbstractArray) = sqrt(1/length(x)*sum((x-y).^2) )
