include("common.jl")

function compute_grid_case(config, cache, grid_count, dx, sensor_count)
    measurement_std = experiment_measurement_std(config, config["grid_dt"])
    metrics = analytic_clarity_metrics(
        cache,
        sensor_count,
        measurement_std,
        config["grid_dt"],
    )
    return (
        Ng=grid_count,
        dx=dx,
        N_robots=sensor_count,
        dt=config["grid_dt"],
        measurement_std=measurement_std,
        mean_state_variance=metrics.mean_state_variance,
        max_state_variance=metrics.max_state_variance,
        clarity=metrics.mean_field_clarity,
    )
end

function compute_grid_rows(config)
    rows = NamedTuple[]
    for dx in config["grid_dx"]
        problem, _, _, _, _ = build_problem(
            config; dx=dx, dt=config["grid_dt"], continuous=true,
        )
        cache = analytic_clarity_cache(problem)
        for sensor_count in config["grid_sensors"]
            push!(rows, compute_grid_case(
                config, cache, length(problem.pts), dx, sensor_count,
            ))
        end
    end
    return rows
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
    figure_path = create_grid_figures(context, rows)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_path, figure_path])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
