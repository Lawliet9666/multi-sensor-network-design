include("common.jl")

function compute_clarity_curves(config)
    problem, xs, ys, data = make_field(config; seed=config["seed"])
    curves = Dict{Int, Vector{Float64}}()
    sensor_results = NamedTuple[]

    for sensor_count in config["estimate_sensors"]
        run_seed = config["seed"] + sensor_count
        clarity, _ = run_filter(
            problem,
            xs,
            ys,
            data,
            sensor_count,
            experiment_measurement_std(config),
            run_seed,
        )
        curves[sensor_count] = clarity
        push!(sensor_results, (; sensor_count, clarity, run_seed))
    end

    return (; times=data.ts, curves, sensor_results)
end

function save_clarity_curves(context, result)
    paths = String[]
    for sensor_result in result.sensor_results
        path = write_jld2(
            joinpath(
                context.data_dir,
                "mean_clarity_Nr$(sensor_result.sensor_count).jld2",
            );
            times=result.times,
            mean_clarity=sensor_result.clarity,
            sensor_count=sensor_result.sensor_count,
            seed=sensor_result.run_seed,
        )
        push!(paths, path)
    end
    return paths
end

function create_clarity_figures(context, result)
    return plot_mean_clarity(result.times, result.curves, context.figure_dir)
end

function main(arguments=ARGS)
    # Step 1: Prepare the experiment.
    context = experiment_context(arguments)
    begin_experiment(context)

    # Step 2: Compute one clarity curve per sensor count.
    result = compute_clarity_curves(context.config)

    # Step 3: Save numerical data.
    data_paths = save_clarity_curves(context, result)

    # Step 4: Create figures.
    figure_path = create_clarity_figures(context, result)

    # Step 5: Report saved outputs.
    return complete_experiment(context, [data_paths; figure_path])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
