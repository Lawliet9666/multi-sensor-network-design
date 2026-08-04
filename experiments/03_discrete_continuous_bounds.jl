include("common.jl")

function run_experiment(context)
    begin_experiment(context)
    config = context.config
    sensor_count = config["sensor_count"]
    measurement_std = config["measurement"]["std"]
    rows = NamedTuple[]

    for step in config["bound_steps"]
        problem, _, _, _, _ = build_problem(config; dt=step)
        continuous_problem, _, _, _, _ = build_problem(config; dt=step, continuous=true)
        measurement_covariance = measurement_std^2 * I(sensor_count)
        sensor_models = precompute_sensor_configs(
            problem, sensor_count, measurement_covariance,
        )
        times = step .* (0:(Int(ceil(config["simulation"]["horizon"] / step)) - 1))
        discrete_history = discrete_covariance_bound(
            problem, times, sensor_models;
            sample_count=config["bound_sample_count"],
            seed=config["seed"], beta=config["beta"],
        )

        fixed_noise_rate = config["simulation"]["dt"] * measurement_std^2
        continuous_std = sqrt(fixed_noise_rate / step)
        continuous_covariance = continuous_std^2 * I(sensor_count)
        continuous_models = precompute_sensor_configs(
            problem, sensor_count, continuous_covariance,
        )
        grid_count = length(problem.pts)
        A = I(grid_count) ⊗ continuous_problem.ss_model.A
        B = I(grid_count) ⊗ continuous_problem.ss_model.B
        G = config["beta"] .* G_from_samples(continuous_models, grid_count, step)
        continuous_bound = continuous_covariance_bound(A, G, B * B')

        push!(rows, (
            Ng=grid_count,
            dx=config["domain"]["dx"],
            N_robots=sensor_count,
            dt=step,
            measurement_std=measurement_std,
            continuous=linear_operator(continuous_bound),
            discrete=linear_operator(last(discrete_history)),
        ))
    end

    data_path = joinpath(context.output, "data", "discrete_continuous_bounds.csv")
    write_csv(
        data_path,
        ["Ng", "dx", "N_robots", "dt", "measurement_std", "continuous", "discrete"],
        [[row.Ng, row.dx, row.N_robots, row.dt, row.measurement_std,
          row.continuous, row.discrete] for row in rows],
    )
    figures = plot_bound_comparison(
        [row.dt for row in rows], [row.continuous for row in rows],
        [row.discrete for row in rows], joinpath(context.output, "figures"),
    )
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiment(experiment_context())
