include("common.jl")

function compute_reconstruction(config)
    dx = config["reconstruction_dx"]
    problem, xs, ys, data = make_field(config; dx=dx, seed=config["seed"])
    capture_index = max(1, length(data.ts) - config["snapshot_offset"])
    clarity, captured = run_filter(
        problem, xs, ys, data, config["reconstruction_sensors"],
        experiment_measurement_std(config), config["seed"] + 100;
        capture_index=capture_index,
    )
    isnothing(captured) && error("The requested reconstruction state was not captured.")
    return (; problem, data, clarity, captured, capture_index)
end

function save_reconstruction_data(context, result)
    config = context.config
    return write_jld2(
        joinpath(context.data_dir, "field_reconstruction.jld2");
        times=result.data.ts,
        clarity=result.clarity,
        capture_index=result.capture_index,
        seed=config["seed"],
        sensors=config["reconstruction_sensors"],
    )
end

function create_reconstruction_figures(context, result)
    return plot_field_reconstruction(
        result.data,
        result.problem,
        result.captured,
        result.capture_index,
        context.figure_dir,
    )
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Reconstruct the field.
    result = compute_reconstruction(context.config)

    # Step 3: Save numerical data.
    data_path = save_reconstruction_data(context, result)

    # Step 4: Create figures.
    figures = create_reconstruction_figures(context, result)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_path, figures.svg, figures.pdf])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
