module GeometricRobotDynID

using Convex
using MAT

# Load the data if not present
function __init__()
    if !ispath(joinpath(pwd(), "data"))
        mkdir(joinpath(pwd(), "data"))
    end

    # Download the data from git
    data_urls = [
        "https://github.com/alex07143/Geometric-Robot-DynID/raw/master/Data/Body_frames_human.mat",
        "https://github.com/alex07143/Geometric-Robot-DynID/raw/master/Data/Train_data_human.mat"
    ]
    for dfile in data_urls
        f_name = splitpath(dfile)[end]
        if !isfile(joinpath(pwd(), "data", f_name))
            download(dfile, joinpath(pwd(), "data", f_name))
        end
    end
end

include("./inertia_processing.jl")
export pseudo_inertia, pullback_metric

include("./load_data.jl")
export load_data, split_data

abstract type abstractDistanceMetric end;

include("./metrics.jl")
export EuclideanDistance, PullbackDistance, EntropicDistance
export distance_metric, rmse #, density_realizability, positive_definite_constraints

include("./optimize.jl")
export estimate_parameters

end # module
