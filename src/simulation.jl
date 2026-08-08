function SpatiotemporalData2D(xs, ys, times, data)
    axes = (
        range(first(xs), last(xs); length=length(xs)),
        range(first(ys), last(ys); length=length(ys)),
        range(first(times), last(times); length=length(times)),
    )
    interpolation = cubic_spline_interpolation(axes, data)
    return SpatiotemporalData2D(collect(xs), collect(ys), collect(times), data, interpolation)
end

function generate_temporal_process(rng, count, step, kernel::Matern12)
    model = state_space_model(kernel, step)
    state = cholesky(initial_covariance(kernel)).L * randn(rng, dims(model))
    noise_root = cholesky(model.W).L
    values = zeros(count)
    values[1] = only(model.C * state)
    for index in 2:count
        state = model.Φ * state + noise_root * randn(rng, dims(model))
        values[index] = only(model.C * state)
    end
    return values
end

function generate_spatiotemporal_process(rng, xs, ys, step, horizon, spatial_kernel, temporal_kernel)
    points = vec([@SVector[x, y] for x in xs, y in ys])
    problem = STGPKFProblem(points, spatial_kernel, temporal_kernel, step)
    count = Int(ceil(horizon / step))
    independent = reduce(hcat, [
        generate_temporal_process(rng, count, step, temporal_kernel)
        for _ in points
    ])'
    field = problem.sqrt_K_gg * independent
    times = collect(step .* (0:(count - 1)))
    cube = reshape(field, length(xs), length(ys), count)
    return SpatiotemporalData2D(xs, ys, times, cube)
end

function measure(rng, data::SpatiotemporalData2D, x, y, time, noise_std=0.1)
    return data.itp(x, y, time) + noise_std * randn(rng)
end

function rand_point(rng, xs, ys)
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    return @SVector[
        xmin + (xmax - xmin) * rand(rng),
        ymin + (ymax - ymin) * rand(rng),
    ]
end

function deep_merge_config(base::AbstractDict, override::AbstractDict)
    merged = deepcopy(base)
    for (key, value) in override
        merged[key] = if haskey(merged, key) && merged[key] isa AbstractDict && value isa AbstractDict
            deep_merge_config(merged[key], value)
        else
            deepcopy(value)
        end
    end
    return merged
end

function load_experiment_config(
    config_root::AbstractString, experiment::AbstractString, profile::AbstractString,
)
    isdir(config_root) || throw(ArgumentError("Configuration directory not found: $config_root"))
    basename(experiment) == experiment || throw(ArgumentError(
        "Experiment must be a file stem, not a path: $experiment",
    ))
    common_path = joinpath(config_root, "common.toml")
    experiment_path = joinpath(config_root, "experiments", "$experiment.toml")
    isfile(common_path) || throw(ArgumentError("Common configuration not found: $common_path"))
    isfile(experiment_path) || throw(ArgumentError(
        "Experiment configuration not found: $experiment_path",
    ))

    common = TOML.parsefile(common_path)
    profiles = TOML.parsefile(experiment_path)
    haskey(profiles, profile) || throw(ArgumentError(
        "Profile '$profile' not found in $experiment_path",
    ))
    return deep_merge_config(common, profiles[profile])
end

"""
    measurement_std_at_step(fixed_std, sigma_c_squared, step, fix_sigma_c)

Return the per-measurement noise standard deviation at `step`. When
`fix_sigma_c` is true, use `sigma_m^2 = sigma_c_squared / step`;
otherwise preserve the per-measurement standard deviation `fixed_std`.
"""
function measurement_std_at_step(
    fixed_std::Real, sigma_c_squared::Real, step::Real, fix_sigma_c::Bool,
)
    fixed_std > 0 || throw(ArgumentError(
        "Fixed measurement noise must be positive.",
    ))
    sigma_c_squared > 0 || throw(ArgumentError(
        "Continuous-time noise intensity must be positive.",
    ))
    step > 0 || throw(ArgumentError("Sampling interval must be positive."))
    return fix_sigma_c ? sqrt(sigma_c_squared / step) : fixed_std
end

function build_problem(
    config; dx=config["domain"]["dx"], dt=config["simulation"]["dt"], continuous=false,
)
    xs = collect(0.0:dx:config["domain"]["max"])
    ys = collect(0.0:dx:config["domain"]["max"])
    points = vec([@SVector[x, y] for x in xs, y in ys])
    temporal = Matern(
        1 / 2, config["temporal_kernel"]["sigma"],
        config["temporal_kernel"]["length_scale"],
    )
    spatial = Matern(
        1 / 2, config["spatial_kernel"]["sigma"],
        config["spatial_kernel"]["length_scale"],
    )
    problem = continuous ? STGPKFProblemContinuous(points, spatial, temporal) :
                           STGPKFProblem(points, spatial, temporal, dt)
    return problem, xs, ys, spatial, temporal
end
