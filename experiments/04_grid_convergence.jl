include("common.jl")

function compute_grid_case(config, dx, sensor_count)
    problem, _, _, _, _ = build_problem(
        config; dx=dx, dt=config["grid_dt"], continuous=true,
    )
    metrics = clarity_metrics(
        problem,
        sensor_count,
        config["measurement"]["std"],
        config["grid_dt"],
    )
    return (
        Ng=length(problem.pts),
        dx=dx,
        N_robots=sensor_count,
        dt=config["grid_dt"],
        measurement_std=config["measurement"]["std"],
        mean_state_variance=metrics.mean_state_variance,
        max_state_variance=metrics.max_state_variance,
        clarity=metrics.mean_state_clarity,
    )
end

function compute_grid_rows(config)
    return [
        compute_grid_case(config, dx, sensor_count)
        for dx in config["grid_dx"]
        for sensor_count in config["grid_sensors"]
    ]
end

function save_grid_table(context, rows)
    return write_csv(joinpath(context.data_dir, "grid_convergence.csv"), rows)
end

function create_grid_figures(context, rows)
    return plot_grid_convergence(rows, context.figure_dir)
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Evaluate every grid and sensor-count combination.
    rows = compute_grid_rows(context.config)

    # Step 3: Save the convergence table.
    data_path = save_grid_table(context, rows)

    # Step 4: Create figures.
    figures = create_grid_figures(context, rows)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
