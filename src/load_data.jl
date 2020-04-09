function load_data(f::String; Σ::AbstractArray = [])
    # Load the data
    data = matopen(f)
    Φ₀ = read(data, "Phi_prior")
    A = read(data, "A")
    b = read(data, "b")
    Q = read(data, "Q")
    close(data)

    # Preprocessing
    if !isempty(Σ)
        @assert (size(Σ, 2) == size(A, 1) && size(Σ, 1) == size(A, 1)) "Provide consistent measurement uncertainty"
        A .= sqrt(Σ)*A
        b .= sqrt(Σ)*b
    end

    Qs = Array{Float64, 2}[]
    for i in Q
        push!(Qs, i)
    end
    φ₀ = vcat(Φ₀...)

    return φ₀, A, b, Qs
end


function split_idxs(N::Int64, target_percentage=0.10)
  splitindex = round(Integer, target_percentage * N)
  return splitindex+1:N, 1:splitindex
end

function split_data(A::AbstractArray{T,2}, b::AbstractArray{T, 2}; percentage::T = 0.1) where T <: Real
    @assert size(A, 1) == size(b,1)
    @assert zero(T) < percentage < one(T)
    train_idx, test_idx = split_idxs(size(A, 1), percentage)
    return view(A, train_idx, :), view(b, train_idx), view(A, test_idx, :), view(b, test_idx)
end
