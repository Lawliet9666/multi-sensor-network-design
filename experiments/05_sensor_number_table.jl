include("common.jl")

function simulate_expected_spatial_clarity(config, sensor_count, seed_offset)
    problem, _, _, data = make_field(config; seed=config["seed"])
    measurement_std = experiment_measurement_std(config)
    simulation_seed = config["seed"] + seed_offset
    point_sets = sample_point_sets(
        MersenneTwister(simulation_seed),
        problem.pts,
        sensor_count;
        configuration_count=config["configuration_count"],
    )
    covariance = measurement_std^2 * I(sensor_count)
    return simulate_expected_clarity(
        problem,
        data,
        point_sets,
        measurement_std,
        covariance;
        trials=config["trials"],
        seed=simulation_seed,
    )
end

function compute_sensor_table_row(config, continuous_problem, target, index)
    measurement_std = experiment_measurement_std(config)
    metrics = minimum_sensor_count(
        continuous_problem,
        target,
        measurement_std,
        config["simulation"]["dt"],
    )
    clarity_statistics = simulate_expected_spatial_clarity(
        config, metrics.sensor_count, 10_000 * index,
    )
    return (
        target=target,
        N_robots=metrics.sensor_count,
        sensing_intensity=metrics.sensing_intensity,
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

function compute_sensor_table(config)
    continuous_problem, _, _, _, _ = build_problem(config; continuous=true)
    return [
        compute_sensor_table_row(config, continuous_problem, target, index)
        for (index, target) in enumerate(config["sensor_targets"])
    ]
end

function save_sensor_table(context, rows)
    return write_csv(joinpath(context.data_dir, "sensor_number_table.csv"), rows)
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Compute the minimum-sensor table.
    rows = compute_sensor_table(context.config)

    # Step 3: Save the table.
    data_path = save_sensor_table(context, rows)

    # Step 4: Report the saved output.
    return complete_experiment(context, [data_path])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
