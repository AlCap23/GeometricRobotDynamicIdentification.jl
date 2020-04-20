function inertia_parameter(x::SpatialInertia)
    y = Array{eltype(x)}(undef, 10)
    y[1] = x.mass
    y[2:4] = x.cross_part
    y[5] = x.moment[1, 1]
    y[6] = x.moment[2, 2]
    y[7] = x.moment[3, 3]
    y[8] = x.moment[1, 2]
    y[9] = x.moment[2, 3]
    y[10] = x.moment[1, 3]
    return y
end

function inertia_parameter(x::RigidBody)
    @assert has_defined_inertia(x) "RigidBody $(x.name) has no defined inertia!"
    return to_array(spatial_inertia(x))
end

function inertia_parameter(m::Mechanism)
    inertias = Array{SpatialInertia{Float64}, 1}()
    for body in bodies(m)
        if has_defined_inertia(body)
            push!(inertias, spatial_inertia(body))
        end
    end
    Φ = vcat(to_array.(inertias)...)
end
