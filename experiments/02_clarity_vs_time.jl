include("common.jl")

function run_experiment(context)
    begin_experiment(context)
    config = context.config
    curves = Dict{String, Vector{Float64}}()
    times = nothing
    data_paths = String[]
    for sensor_count in config["estimate_sensors"]
        problem, xs, ys, data = make_field(config; seed=config["seed"])
        clarity, _ = run_filter(
            problem, xs, ys, data, sensor_count, config["measurement"]["std"],
            config["seed"] + sensor_count,
        )
        label = sensor_count == 1 ? "1 sensor" : "$sensor_count sensors"
        curves[label] = clarity
        times = data.ts
        path = joinpath(context.output, "data", "mean_clarity_Nr$(sensor_count).jld2")
        path = write_jld2(
            path; times=data.ts, mean_clarity=clarity, sensor_count=sensor_count,
            seed=config["seed"] + sensor_count,
        )
        push!(data_paths, path)
    end
    figures = plot_mean_clarity(times, curves, joinpath(context.output, "figures"))
    return complete_experiment(context, [data_paths; figures.svg; figures.pdf])
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiment(experiment_context())
