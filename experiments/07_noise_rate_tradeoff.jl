include("common.jl")

function compute_tradeoff_grid(config)
    config["measurement"]["fix_sigma_c"] === false || error(
        "Experiment 07 must sweep sigma_m^2 and Delta t / N_r rather than fix sigma_c^2.",
    )
    problem, _, _, _, _ = build_problem(config; continuous=true)
    cache = analytic_clarity_cache(problem)
    measurement_variances = collect(range(
        config["measurement_variance_min"],
        config["measurement_variance_max"];
        length=config["measurement_variance_points"],
    ))
    normalized_intervals = collect(range(
        config["normalized_interval_min"],
        config["normalized_interval_max"];
        length=config["normalized_interval_points"],
    ))
    sensor_count = config["sensor_count"]
    clarity_lower_bound = Matrix{Float64}(
        undef,
        length(normalized_intervals),
        length(measurement_variances),
    )

    for (row, normalized_interval) in enumerate(normalized_intervals),
        (column, measurement_variance) in enumerate(measurement_variances)
        step = normalized_interval * sensor_count
        metrics = analytic_clarity_metrics(
            cache,
            sensor_count,
            sqrt(measurement_variance),
            step,
        )
        clarity_lower_bound[row, column] = metrics.mean_field_clarity
    end

    return (;
        measurement_variances,
        normalized_intervals,
        clarity_lower_bound,
        sensor_count,
    )
end

function save_tradeoff_data(context, result)
    return write_jld2(
        joinpath(
            context.data_dir,
            "clarity_grid_Nr$(result.sensor_count).jld2",
        );
        measurement_variances=result.measurement_variances,
        normalized_intervals=result.normalized_intervals,
        clarity_lower_bound=result.clarity_lower_bound,
        sensor_count=result.sensor_count,
        seed=context.config["seed"],
    )
end

function create_tradeoff_figures(context, result)
    return plot_noise_rate_tradeoff(
        result.measurement_variances,
        result.normalized_intervals,
        result.clarity_lower_bound,
        context.figure_dir;
        contour_targets=context.config["contour_clarity_levels"],
    )
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Compute the noise-rate tradeoff grid.
    result = compute_tradeoff_grid(context.config)

    # Step 3: Save numerical data.
    data_path = save_tradeoff_data(context, result)

    # Step 4: Create figures.
    figure_path = create_tradeoff_figures(context, result)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_path, figure_path])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
