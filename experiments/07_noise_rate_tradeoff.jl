include("common.jl")

function run_experiment(context)
    begin_experiment(context)
    config = context.config
    problem, _, _, _, _ = build_problem(config; continuous=true)
    noise_variances = collect(range(
        config["heatmap_sigma_min"]^2, config["heatmap_sigma_max"]^2;
        length=config["heatmap_nx"],
    ))
    steps = collect(range(
        config["heatmap_dt_min"], config["heatmap_dt_max"];
        length=config["heatmap_ny"],
    ))
    sensor_count = config["sensor_count"]
    clarity = Matrix{Float64}(undef, length(steps), length(noise_variances))
    for (row, step) in enumerate(steps), (column, variance) in enumerate(noise_variances)
        metrics = clarity_metrics(problem, sensor_count, sqrt(variance), step)
        clarity[row, column] = metrics.mean_state_clarity
    end

    data_path = joinpath(context.output, "data", "clarity_grid_Nr$(sensor_count).jld2")
    data_path = write_jld2(
        data_path; noise_variances=noise_variances, steps=steps, clarity=clarity,
        sensor_count=sensor_count, seed=config["seed"],
    )
    figures = plot_noise_rate_tradeoff(
        noise_variances, steps, clarity, sensor_count,
        joinpath(context.output, "figures"),
    )
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiment(experiment_context())
