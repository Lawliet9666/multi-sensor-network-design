include("common.jl")

function compute_sensor_case(
    config,
    continuous_problem,
    problem,
    data,
    sensor_count,
    index,
)
    measurement_std = experiment_measurement_std(config)
    metrics = clarity_metrics(
        continuous_problem,
        sensor_count,
        measurement_std,
        config["simulation"]["dt"],
    )
    simulation_seed = config["seed"] + 20_000 * index
    point_sets = sample_point_sets(
        MersenneTwister(simulation_seed),
        problem.pts,
        sensor_count;
        configuration_count=config["configuration_count"],
    )
    covariance = measurement_std^2 * I(sensor_count)
    clarity_statistics = simulate_expected_clarity(
        problem,
        data,
        point_sets,
        measurement_std,
        covariance;
        trials=config["trials"],
        seed=simulation_seed,
    )
    return (
        N_robots=sensor_count,
        mean_field_variance=metrics.mean_field_variance,
        max_field_variance=metrics.max_field_variance,
        bound_clarity=metrics.mean_field_clarity,
        expected_spatial_mean_clarity_estimate=
            clarity_statistics.expected_spatial_mean_clarity_estimate,
        spatial_mean_clarity_of_expected_covariance_estimate=
            clarity_statistics.spatial_mean_clarity_of_expected_covariance_estimate,
        spatial_mean_clarity_sample_std=
            clarity_statistics.spatial_mean_clarity_sample_std,
        spatial_mean_clarity_standard_error=
            clarity_statistics.spatial_mean_clarity_standard_error,
        monte_carlo_trials=clarity_statistics.monte_carlo_trials,
    )
end

function compute_sensor_curve(config)
    continuous_problem, _, _, _, _ = build_problem(config; continuous=true)
    problem, _, _, data = make_field(config; seed=config["seed"])
    return [
        compute_sensor_case(
            config,
            continuous_problem,
            problem,
            data,
            sensor_count,
            index,
        )
        for (index, sensor_count) in enumerate(config["sensor_curve"])
    ]
end

function save_sensor_curve(context, rows)
    return write_csv(joinpath(context.data_dir, "clarity_vs_sensors.csv"), rows)
end

function create_sensor_curve_figures(context, rows)
    return plot_sensor_table_curve(
        rows,
        context.figure_dir;
        target_clarity=context.config["plot_target_clarity"],
    )
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Compute clarity for every sensor count.
    rows = compute_sensor_curve(context.config)

    # Step 3: Save numerical data.
    data_path = save_sensor_curve(context, rows)

    # Step 4: Create figures.
    figure_path = create_sensor_curve_figures(context, rows)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_path, figure_path])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
