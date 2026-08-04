include("common.jl")

function empirical_clarity(config, sensor_count, seed_offset)
    problem, _, _, data = make_field(config; seed=config["seed"])
    point_sets = sample_point_sets(
        MersenneTwister(config["seed"] + seed_offset), problem.pts, sensor_count;
        configuration_count=config["configuration_count"],
    )
    covariance = config["measurement"]["std"]^2 * I(sensor_count)
    expected = simulate_expected_covariance(
        problem, data, point_sets, config["measurement"]["std"], covariance;
        trials=config["trials"], seed=config["seed"] + seed_offset,
    )
    return mean(get_clarity(problem, last(expected)))
end

function run_experiment(context)
    begin_experiment(context)
    config = context.config
    continuous_problem, _, _, _, _ = build_problem(config; continuous=true)
    rows = NamedTuple[]
    for (index, target) in enumerate(config["sensor_targets"])
        metrics = minimum_sensor_count(
            continuous_problem, target, config["measurement"]["std"],
            config["simulation"]["dt"],
        )
        empirical = empirical_clarity(config, metrics.sensor_count, 10_000 * index)
        push!(rows, (
            target=target, N_robots=metrics.sensor_count,
            sensing_intensity=metrics.sensing_intensity,
            mean_field_variance=metrics.mean_field_variance,
            max_field_variance=metrics.max_field_variance,
            bound_clarity=metrics.mean_field_clarity,
            empirical_clarity=empirical,
        ))
    end

    data_path = joinpath(context.output, "data", "sensor_number_table.csv")
    write_csv(
        data_path,
        ["target", "N_robots", "sensing_intensity", "mean_field_variance",
         "max_field_variance", "bound_clarity", "empirical_clarity"],
        [[row.target, row.N_robots, row.sensing_intensity, row.mean_field_variance,
          row.max_field_variance, row.bound_clarity, row.empirical_clarity] for row in rows],
    )
    return complete_experiment(context, [data_path])
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiment(experiment_context())
