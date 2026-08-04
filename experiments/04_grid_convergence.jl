include("common.jl")

function run_experiment(context)
    begin_experiment(context)
    config = context.config
    rows = NamedTuple[]
    for dx in config["grid_dx"], sensor_count in config["grid_sensors"]
        problem, _, _, _, _ = build_problem(
            config; dx=dx, dt=config["grid_dt"], continuous=true,
        )
        metrics = clarity_metrics(
            problem, sensor_count, config["measurement"]["std"], config["grid_dt"],
        )
        push!(rows, (
            Ng=length(problem.pts), dx=dx, N_robots=sensor_count,
            dt=config["grid_dt"], measurement_std=config["measurement"]["std"],
            mean_state_variance=metrics.mean_state_variance,
            max_state_variance=metrics.max_state_variance,
            clarity=metrics.mean_state_clarity,
        ))
    end

    data_path = joinpath(context.output, "data", "grid_convergence.csv")
    write_csv(
        data_path,
        ["Ng", "dx", "N_robots", "dt", "measurement_std",
         "mean_state_variance", "max_state_variance", "clarity"],
        [[row.Ng, row.dx, row.N_robots, row.dt, row.measurement_std,
          row.mean_state_variance, row.max_state_variance, row.clarity] for row in rows],
    )
    figures = plot_grid_convergence(rows, joinpath(context.output, "figures"))
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiment(experiment_context())
