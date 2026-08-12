using Pkg

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(REPOSITORY_ROOT)

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using JLD2
using Kronecker
using LinearAlgebra
using Random
using SpatiotemporalGPs
using StaticArrays
using Statistics
using TOML

const EXPERIMENT_DETAILS = Dict(
    "01_field_reconstruction" => (
        title="Field reconstruction",
        purpose="Reconstruct one field snapshot and compute its clarity map.",
    ),
    "02_clarity_vs_time" => (
        title="Mean clarity over time",
        purpose="Compare mean clarity trajectories for the configured sensor counts.",
    ),
    "03_discrete_continuous_bounds" => (
        title="Discrete and continuous covariance bounds",
        purpose="Compare the steady discrete-time and continuous-time covariance bounds.",
    ),
    "04_grid_convergence" => (
        title="Grid convergence",
        purpose="Evaluate clarity-bound convergence over grid sizes and sensor counts.",
    ),
    "05_sensor_number_table" => (
        title="Minimum sensor-number table",
        purpose="Find the minimum sensor count for each target clarity and compare with simulation.",
    ),
    "06_clarity_vs_sensors" => (
        title="Clarity versus sensor count",
        purpose="Compare analytic and empirical clarity across sensor counts.",
    ),
    "07_noise_rate_tradeoff" => (
        title="Noise-model and sampling-interval sweep",
        purpose="Evaluate the clarity bound over the configured noise parameter and sampling interval.",
    ),
)

function option_value(arguments, name, default)
    index = findfirst(==(name), arguments)
    isnothing(index) && return default
    index < length(arguments) || error("Missing value after $name")
    return arguments[index + 1]
end

function experiment_context(arguments=ARGS)
    profile = option_value(arguments, "--profile", "paper")
    config_root = abspath(option_value(
        arguments, "--config", joinpath(REPOSITORY_ROOT, "config"),
    ))
    requested_output = option_value(arguments, "--output", nothing)
    experiment = splitext(basename(PROGRAM_FILE))[1]
    output_root = isnothing(requested_output) ?
        joinpath(REPOSITORY_ROOT, "results", experiment) :
        requested_output
    output = abspath(output_root)
    data_dir = joinpath(output, "data")
    figure_dir = joinpath(output, "figures")
    config_dir = joinpath(output, "configs")
    config = load_experiment_config(config_root, experiment, profile)
    mkpath(data_dir)
    mkpath(figure_dir)
    return (
        experiment=experiment,
        profile=profile,
        config=config,
        output=output,
        data_dir=data_dir,
        figure_dir=figure_dir,
        config_dir=config_dir,
        config_root=config_root,
    )
end

function write_resolved_config(context)
    path = joinpath(context.config_dir, "$(context.experiment).toml")
    mkpath(dirname(path))
    open(path, "w") do stream
        TOML.print(stream, context.config; sorted=true)
    end
    absolute_path = abspath(path)
    println("Saved resolved configuration (TOML): $absolute_path")
    return absolute_path
end

function experiment_measurement_std(config, step=config["simulation"]["dt"])
    measurement = config["measurement"]
    return measurement_std_at_step(
        measurement["std"],
        measurement["sigma_c_squared"],
        step,
        measurement["fix_sigma_c"],
    )
end

function measurement_noise_description(context)
    measurement = context.config["measurement"]
    if context.experiment == "07_noise_rate_tradeoff"
        return measurement["fix_sigma_c"] ?
            "fixed sigma_c^2 within each sweep case (revised paper)" :
            "fixed sigma_m within each sweep case (root_v2 paper)"
    end
    return measurement["fix_sigma_c"] ?
        "fixed sigma_c^2 = $(measurement["sigma_c_squared"]) (revised paper)" :
        "fixed sigma_m = $(measurement["std"]) (root_v2 paper)"
end

function begin_experiment(context)
    details = get(EXPERIMENT_DETAILS, context.experiment, (
        title=context.experiment,
        purpose="Run the selected paper experiment.",
    ))
    println()
    println("=== $(details.title) ===")
    println("Experiment: $(context.experiment)")
    println("Profile: $(context.profile)")
    println("Purpose: $(details.purpose)")
    println("Output directory: $(context.output)")
    write_resolved_config(context)
    println("Measurement-noise model: $(measurement_noise_description(context))")
    return details
end

function complete_experiment(context, outputs)
    paths = abspath.(collect(outputs))
    details = EXPERIMENT_DETAILS[context.experiment]
    println("Completed: $(details.title)")
    println("Saved $(length(paths)) result file(s):")
    for path in paths
        println("  - $path")
    end
    return paths
end

function write_csv(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do stream
        println(stream, join(header, ','))
        for row in rows
            println(stream, join(row, ','))
        end
    end
    absolute_path = abspath(path)
    println("Saved data (CSV): $absolute_path")
    return absolute_path
end

function write_csv(path, rows)
    isempty(rows) && error("Cannot write an empty CSV result.")
    columns = propertynames(first(rows))
    all(row -> propertynames(row) == columns, rows) || error(
        "CSV rows do not have consistent columns.",
    )
    # Preserve the existing CSV representation for mixed integer/float rows.
    values = (
        promote((getproperty(row, column) for column in columns)...)
        for row in rows
    )
    return write_csv(path, string.(columns), values)
end

function write_jld2(path; kwargs...)
    mkpath(dirname(path))
    jldsave(path; kwargs...)
    absolute_path = abspath(path)
    println("Saved data (JLD2): $absolute_path")
    return absolute_path
end

function make_field(
    config;
    dx=config["domain"]["dx"],
    dt=config["simulation"]["dt"],
    horizon=config["simulation"]["horizon"],
    seed=config["seed"],
)
    problem, xs, ys, spatial, temporal = build_problem(config; dx=dx, dt=dt)
    rng = MersenneTwister(seed)
    data = generate_spatiotemporal_process(rng, xs, ys, dt, horizon, spatial, temporal)
    return problem, xs, ys, data
end

function run_filter(problem, xs, ys, data, sensor_count, measurement_std, seed; capture_index=nothing)
    rng = MersenneTwister(seed)
    state = stgpkf_initialize(problem)
    captured = nothing
    clarity = zeros(length(data.ts))
    covariance = measurement_std^2 * I(sensor_count)
    for (index, time) in enumerate(data.ts)
        points = [rand_point(rng, xs, ys) for _ in 1:sensor_count]
        observations = [measure(rng, data, point[1], point[2], time, measurement_std) for point in points]
        corrected = stgpkf_correct(problem, state, points, observations, covariance)
        clarity[index] = mean(get_estimate_clarity(problem, corrected))
        index == capture_index && (captured = corrected)
        state = stgpkf_predict(problem, corrected)
    end
    return clarity, captured
end
