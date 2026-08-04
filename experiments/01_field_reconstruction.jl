include("common.jl")

function run_experiment(context)
    begin_experiment(context)
    config = context.config
    dx = config["reconstruction_dx"]
    problem, xs, ys, data = make_field(config; dx=dx, seed=config["seed"])
    capture_index = max(1, length(data.ts) - config["snapshot_offset"])
    clarity, captured = run_filter(
        problem, xs, ys, data, config["reconstruction_sensors"],
        config["measurement"]["std"], config["seed"] + 100;
        capture_index=capture_index,
    )
    isnothing(captured) && error("The requested reconstruction state was not captured.")
    data_path = joinpath(context.output, "data", "field_reconstruction.jld2")
    data_path = write_jld2(
        data_path; times=data.ts, clarity=clarity, capture_index=capture_index,
        seed=config["seed"], sensors=config["reconstruction_sensors"],
    )
    figures = plot_field_reconstruction(
        data, problem, captured, capture_index, joinpath(context.output, "figures"),
    )
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiment(experiment_context())
