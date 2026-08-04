include("common.jl")

function run_experiment(context)
    begin_experiment(context)
    config = context.config
    continuous_problem, _, _, _, _ = build_problem(config; continuous=true)
    problem, _, _, data = make_field(config; seed=config["seed"])
    rows = NamedTuple[]

    for (index, sensor_count) in enumerate(config["sensor_curve"])
        metrics = clarity_metrics(
            continuous_problem, sensor_count, config["measurement"]["std"],
            config["simulation"]["dt"],
        )
        point_sets = sample_point_sets(
            MersenneTwister(config["seed"] + 20_000 * index), problem.pts, sensor_count;
            configuration_count=config["configuration_count"],
        )
        covariance = config["measurement"]["std"]^2 * I(sensor_count)
        expected = simulate_expected_covariance(
            problem, data, point_sets, config["measurement"]["std"], covariance;
            trials=config["trials"], seed=config["seed"] + 20_000 * index,
        )
        push!(rows, (
            N_robots=sensor_count,
            mean_field_variance=metrics.mean_field_variance,
            max_field_variance=metrics.max_field_variance,
            bound_clarity=metrics.mean_field_clarity,
            empirical_clarity=mean(get_clarity(problem, last(expected))),
        ))
    end

    data_path = joinpath(context.output, "data", "clarity_vs_sensors.csv")
    write_csv(
        data_path,
        ["N_robots", "mean_field_variance", "max_field_variance",
         "bound_clarity", "empirical_clarity"],
        [[row.N_robots, row.mean_field_variance, row.max_field_variance,
          row.bound_clarity, row.empirical_clarity] for row in rows],
    )
    figures = plot_sensor_table_curve(rows, joinpath(context.output, "figures"))
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiment(experiment_context())
