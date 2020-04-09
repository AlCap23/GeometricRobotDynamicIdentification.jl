using LinearAlgebra

function pseudo_inertia(p::AbstractArray{T}) where T <: Real
    m = p[1]
    h = p[2:4]
    J = [p[5] p[8] p[10]; p[8] p[6] p[9]; p[10] p[9] p[7]]

    return [convert(eltype(p), 0.5)*tr(J)*I(3)-J h; h' m]
end

function pseudo_inertia(p::AbstractArray{T}) where T <: Variable
    J = [p[5] p[8] p[10]; p[8] p[6] p[9]; p[10] p[9] p[7]]
    Js = 0.5*tr(J)*I(3)-J
    P = [Js[1,1] Js[1,2] Js[1,3] p[2]; Js[2,1] Js[2,2] Js[2,3] p[3]; Js[3,1] Js[3,2] Js[3,3] p[4]; p[2] p[3] p[4] p[1]]
    return P
end

function pseudo_inertia(p::Convex.IndexAtom)
    J = [p[5] p[8] p[10]; p[8] p[6] p[9]; p[10] p[9] p[7]]
    Js = 0.5*tr(J)*I(3)-J
    P = [Js[1,1] Js[1,2] Js[1,3] p[2]; Js[2,1] Js[2,2] Js[2,3] p[3]; Js[3,1] Js[3,2] Js[3,3] p[4]; p[2] p[3] p[4] p[1]]
    return P
end

function pullback_metric(p::AbstractArray{T}) where T <: Real
    g = zeros(10,10)
    P_inv = inv(pseudo_inertia(p))
    for i in 1:10
        for j in 1:10
            ei = zeros(10)
            ei[i] = one(eltype(p))
            ej = zeros(10)
            ej[j] = one(eltype(p))
            Vi = pseudo_inertia(ei)
            Vj = pseudo_inertia(ej)
            g[i,j] = tr(P_inv* Vi * P_inv*Vj)
        end
    end
    return g
end

function get_com(p::AbstractArray{T,1}) where T <: Real
    @assert length(p) == 10 "Provide proper inertial parameters"
    return p[2:4]
end

function get_mass(p::AbstractArray{T,1}) where T <: Real
    @assert length(p) == 10 "Provide proper inertial parameters"
    return p[1]
end

function get_inertia(p::AbstractArray{T,1}) where T <: Real
    @assert length(p) == 10 "Provide proper inertial parameters"
    return p[5:end]
end

function split_parameter(p::AbstractArray{T}) where T <: Real
    @assert length(p)%10 == 0 "Provide proper inertial parameters"
    n_bodies = Int64(length(p)/10)
    masses = zeros(T, n_bodies)
    coms = Array{T,1}[]
    inertias = Array{T,1}[]

    for i in 1:n_bodies
        lower = (i-1)*10+1
        upper = i*10
        p_tmp = p[lower:upper]
        push!(masses, get_mass(p_tmp))
        push!(coms, get_com(p_tmp))
        push!(inertias, get_inertia(p_tmp))
    end
    return masses, coms, inertias
end
