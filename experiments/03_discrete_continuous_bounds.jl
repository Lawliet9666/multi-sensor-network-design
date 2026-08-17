include("common.jl")

function compute_finite_horizon_at_step(config, step)
    sensor_count = config["sensor_count"]
    measurement_std = experiment_measurement_std(config, step)
    problem, _, _, _, _ = build_problem(config; dt=step)
    continuous_problem, _, _, _, _ = build_problem(
        config; dt=step, continuous=true,
    )

    measurement_covariance = measurement_std^2 * I(sensor_count)
    sensor_models = precompute_sensor_configs(
        problem, sensor_count, measurement_covariance,
    )
    comparison = finite_horizon_covariance_comparison(
        problem,
        continuous_problem,
        sensor_models;
        horizon=config["simulation"]["horizon"],
        trials=config["monte_carlo_trials"],
        seed=config["seed"],
    )

    grid_count = length(problem.pts)
    summary = (
        Ng=grid_count,
        dx=config["domain"]["dx"],
        N_robots=sensor_count,
        dt=step,
        horizon=config["simulation"]["horizon"],
        measurement_std=measurement_std,
        monte_carlo_trials=config["monte_carlo_trials"],
        empirical_max_frobenius_error=comparison.empirical_max_frobenius_error,
        discrete_max_frobenius_error=comparison.discrete_max_frobenius_error,
        empirical_max_linear_operator_error=
            comparison.empirical_max_linear_operator_error,
        discrete_max_linear_operator_error=
            comparison.discrete_max_linear_operator_error,
        empirical_max_relative_linear_operator_error_percent=
            comparison.empirical_max_relative_linear_operator_error_percent,
        discrete_max_relative_linear_operator_error_percent=
            comparison.discrete_max_relative_linear_operator_error_percent,
    )
    return summary, comparison
end

function compute_finite_horizon_results(config)
    trajectory_step = config["trajectory_step"]
    any(step -> isapprox(step, trajectory_step), config["bound_steps"]) || error(
        "trajectory_step must be one of bound_steps.",
    )
    summaries = NamedTuple[]
    trajectory = nothing
    for step in config["bound_steps"]
        summary, comparison = compute_finite_horizon_at_step(config, step)
        push!(summaries, summary)
        isapprox(step, trajectory_step) && (trajectory = comparison)
    end
    isnothing(trajectory) && error("No trajectory was computed at trajectory_step.")
    return summaries, trajectory
end

function save_bound_table(context, rows)
    return write_csv(
        joinpath(context.data_dir, "discrete_continuous_bounds.csv"),
        rows,
    )
end

function save_trajectory_table(context, trajectory)
    rows = [
        (
            time=trajectory.times[index],
            empirical_mean_covariance=trajectory.empirical_mean_covariance[index],
            discrete_bound_mean_covariance=trajectory.discrete_bound_mean_covariance[index],
            continuous_bound_mean_covariance=trajectory.continuous_bound_mean_covariance[index],
        )
        for index in eachindex(trajectory.times)
    ]
    return write_csv(
        joinpath(context.data_dir, "discrete_continuous_trajectory.csv"),
        rows,
    )
end

function create_bound_figures(context, summaries, trajectory)
    return plot_bound_comparison(
        trajectory,
        summaries,
        context.figure_dir,
    )
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Compute the finite-horizon covariance comparisons.
    summaries, trajectory = compute_finite_horizon_results(context.config)

    # Step 3: Save the convergence errors and representative trajectory.
    data_path = save_bound_table(context, summaries)
    trajectory_path = save_trajectory_table(context, trajectory)

    # Step 4: Create figures.
    figure_path = create_bound_figures(context, summaries, trajectory)

    # Step 5: Report saved outputs.
    return complete_experiment(
        context, [data_path, trajectory_path, figure_path],
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
