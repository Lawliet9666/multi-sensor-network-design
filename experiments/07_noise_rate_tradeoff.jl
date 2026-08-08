include("common.jl")

function compute_tradeoff_grid(config)
    problem, _, _, _, _ = build_problem(config; continuous=true)
    noise_variances = collect(range(
        config["heatmap_sigma_min"]^2,
        config["heatmap_sigma_max"]^2;
        length=config["heatmap_nx"],
    ))
    steps = collect(range(
        config["heatmap_dt_min"],
        config["heatmap_dt_max"];
        length=config["heatmap_ny"],
    ))
    sensor_count = config["sensor_count"]
    fix_sigma_c = config["measurement"]["fix_sigma_c"]
    clarity = Matrix{Float64}(undef, length(steps), length(noise_variances))

    for (row, step) in enumerate(steps), (column, variance) in enumerate(noise_variances)
        measurement_std = measurement_std_at_step(
            sqrt(variance), variance, step, fix_sigma_c,
        )
        metrics = clarity_metrics(problem, sensor_count, measurement_std, step)
        clarity[row, column] = metrics.mean_state_clarity
    end

    noise_axis = fix_sigma_c ? "sigma_c_squared" : "sigma_m_squared"
    return (; noise_variances, steps, clarity, sensor_count, fix_sigma_c, noise_axis)
end

function save_tradeoff_data(context, result)
    return write_jld2(
        joinpath(
            context.data_dir,
            "clarity_grid_Nr$(result.sensor_count).jld2",
        );
        noise_variances=result.noise_variances,
        steps=result.steps,
        clarity=result.clarity,
        sensor_count=result.sensor_count,
        fix_sigma_c=result.fix_sigma_c,
        noise_axis=result.noise_axis,
        seed=context.config["seed"],
    )
end

function create_tradeoff_figures(context, result)
    return plot_noise_rate_tradeoff(
        result.noise_variances,
        result.steps,
        result.clarity,
        result.sensor_count,
        context.figure_dir;
        fix_sigma_c=result.fix_sigma_c,
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
    figures = create_tradeoff_figures(context, result)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
