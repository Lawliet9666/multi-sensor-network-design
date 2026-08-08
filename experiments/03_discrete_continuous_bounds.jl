include("common.jl")

function compute_bound_at_step(config, step)
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
    times = step .* (0:(Int(ceil(config["simulation"]["horizon"] / step)) - 1))
    discrete_history = discrete_covariance_bound(
        problem,
        times,
        sensor_models;
        sample_count=config["bound_sample_count"],
        seed=config["seed"],
        beta=config["beta"],
    )

    grid_count = length(problem.pts)
    A = I(grid_count) ⊗ continuous_problem.ss_model.A
    B = I(grid_count) ⊗ continuous_problem.ss_model.B
    G = config["beta"] .* G_from_samples(sensor_models, grid_count, step)
    continuous_bound = continuous_covariance_bound(A, G, B * B')

    return (
        Ng=grid_count,
        dx=config["domain"]["dx"],
        N_robots=sensor_count,
        dt=step,
        measurement_std=measurement_std,
        continuous=linear_operator(continuous_bound),
        discrete=linear_operator(last(discrete_history)),
    )
end

function compute_bound_rows(config)
    return [
        compute_bound_at_step(config, step)
        for step in config["bound_steps"]
    ]
end

function save_bound_table(context, rows)
    return write_csv(
        joinpath(context.data_dir, "discrete_continuous_bounds.csv"),
        rows,
    )
end

function create_bound_figures(context, rows)
    return plot_bound_comparison(
        getproperty.(rows, :dt),
        getproperty.(rows, :continuous),
        getproperty.(rows, :discrete),
        context.figure_dir,
    )
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Compute discrete and continuous bounds.
    rows = compute_bound_rows(context.config)

    # Step 3: Save the comparison table.
    data_path = save_bound_table(context, rows)

    # Step 4: Create figures.
    figures = create_bound_figures(context, rows)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
