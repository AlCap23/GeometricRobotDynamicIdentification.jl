using GeometricRobotDynID
using RigidBodyDynamics
using MeshCatMechanisms

using LinearAlgebra

srcdir = dirname(pathof(RigidBodyDynamics))
urdf = joinpath(srcdir, "..", "test", "urdf", "Acrobot.urdf")

mechanism = parse_urdf(urdf)
state = MechanismState(mechanism)



inertia_parameter(mechanism)



function regressor(mechanism, state)
    # Get the body ids
    bds = Int64[]
    for body in bodies(mechanism)
        if has_defined_inertia(body)
            push!(bds, body.id)
        end
    end
    bds_keys = state.twists_wrt_world.keys[bds]
    # Get the twist and accelerations
    ẋ = state.twists_wrt_world
    ẍ = state.bias_accelerations_wrt_world
    Ẋ, Ẍ = Twist{Float64}[], SpatialAcceleration{Float64}[]
    for k in bds_keys
        push!(Ẋ, ẋ[k])
        push!(Ẍ, ẍ[k])
    end
    g = mechanism.gravitational_acceleration.v
    js = joints(mechanism)
    A = zeros(Float64, length(bds_keys), 10)
    for i in 1:length(bds_keys)
        T = transform_to_root(state, bds_keys[i]).mat
        sx = [0 -T[3, 4] -T[2, 4]; T[3, 4] 0 T[1, 4]; -T[2, 4] T[1, 4] 0]
        Ad = [T[1:3, 1:3] zeros(3,3); sx*T[1:3, 1:3] T[1:3, 1:3]]
        p̈ = Ẋ[i].linear - g
        dω = Ẍ[i].angular
        ω = Ẋ[i].angular
        px = [0 -p̈[3] p̈[2]; p̈[3] 0 -p̈[1]; -p̈[2] p̈[1] 0]
        ωx = [0 -ω[3] ω[2]; ω[3] 0 -ω[1]; -ω[2] ω[1] 0]
        dωx = [0 -dω[3] dω[2]; dω[3] 0 -dω[1]; -dω[2] dω[1] 0]
        pω = [ω[1] ω[2] ω[3] 0 0 0; 0 ω[1] 0 ω[2] ω[3] 0; 0 0 ω[1] 0 ω[2] ω[3]]
        pdω = [dω[1] dω[2] dω[3] 0 0 0; 0 dω[1] 0 dω[2] dω[3] 0; 0 0 dω[1] 0 dω[2] dω[3]]
        t = motion_subspace(js[1],configuration(state, js[1]))
        A[i, 1:10] = [t.linear... t.angular...] *Ad* [p̈ dωx+ωx*ωx zeros(eltype(g), 3, 6); zeros(3, 1) px pdω+ωx*pω]
    end
    return A
end

rand_configuration!(state)
rand_velocity!(state)
regressor(mechanism, state)

res = DynamicsResult(mechanism)
v̇ = velocity(state)
rand_velocity!(state)

res.v̇ = v̇
spatial_accelerations!(res, state)
res.accelerations
state.twists_wrt_world[BodyID(2)]



accelerations = res.accelerations
twists = state.joint_twists

# Compute the spatial accelerations of each link
# and updates twists
spatial_accelerations!(accelerations, state, v̇)
# Transform into parent frame

twists = state.joint_twists
j = joints(mechanism)
elbow = j[1]
elbow.
a = zeros(Float64, 2*10)
for j in joints(mechanism)
    println(get_p)
end

a = accelerations[BodyID(2)]

transform(a,)

function dotw(x)
    a = zeros(eltype(x), 3, 6)
    a[1, 1:3] = x
    a[2, 2] = x[1]
    a[2,4:5] = x[2:3]
    a[3,3] = x[1]
    a[3, 5:6] = x[2:3]
    return a
end

function regressor!(bodyid, jointid, accelerations, twists)
    ddp = linear(accelerations[bodyid])
    dω = angular(accelerations[bodyid])
    ω = angular(twists[jointid])
    ωx = Spatial.vector_to_skew_symmetric(angular(twists[jointid]))
    ωx2 = Spatial.vector_to_skew_symmetric_squared(angular(twists[jointid]))
    dωx = Spatial.vector_to_skew_symmetric(dω)
    ddpx = Spatial.vector_to_skew_symmetric(ddp)
    X = dotw(dω)+dωx*dotw(ω)
    a = [ddp dωx+ωx2 zeros(3,6); zero(ddp) ddpx X]
end

A = regressor!(BodyID(2), JointID(1), acceleartions, twists)

X = Array(motion_subspace(elbow, configuration(state, elbow)))

X'*A
