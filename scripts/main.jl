using Pkg

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(REPOSITORY_ROOT)

function option_value(arguments, name, default)
    index = findfirst(==(name), arguments)
    isnothing(index) && return default
    index < length(arguments) || error("Missing value after $name")
    return arguments[index + 1]
end

function prepare_output(arguments)
    requested = option_value(arguments, "--output", nothing)
    output = isnothing(requested) ?
        joinpath(REPOSITORY_ROOT, "results") :
        abspath(requested)
    mkpath(output)
    return output
end

function required_outputs(output)
    return [
        joinpath(output, "01_field_reconstruction", "data", "field_reconstruction.jld2"),
        joinpath(output, "02_clarity_vs_time", "data", "mean_clarity_Nr1.jld2"),
        joinpath(output, "02_clarity_vs_time", "data", "mean_clarity_Nr6.jld2"),
        joinpath(output, "02_clarity_vs_time", "data", "mean_clarity_Nr20.jld2"),
        joinpath(output, "03_discrete_continuous_bounds", "data", "discrete_continuous_bounds.csv"),
        joinpath(output, "03_discrete_continuous_bounds", "data", "discrete_continuous_trajectory.csv"),
        joinpath(output, "04_grid_convergence", "data", "grid_convergence.csv"),
        joinpath(output, "05_sensor_number_table", "data", "sensor_number_table.csv"),
        joinpath(output, "06_clarity_vs_sensors", "data", "clarity_vs_sensors.csv"),
        joinpath(output, "07_noise_rate_tradeoff", "data", "clarity_grid_Nr1.jld2"),
        joinpath(output, "01_field_reconstruction", "figures", "estimate_1.pdf"),
        joinpath(output, "02_clarity_vs_time", "figures", "estimate_2.pdf"),
        joinpath(output, "03_discrete_continuous_bounds", "figures", "continuous_vs_discrete2.pdf"),
        joinpath(output, "04_grid_convergence", "figures", "ng_converge_clarity.pdf"),
        joinpath(output, "06_clarity_vs_sensors", "figures", "clarity_vs_nr_compact.pdf"),
        joinpath(output, "07_noise_rate_tradeoff", "figures", "clarity_heatmap_Nr1_compact.pdf"),
    ]
end

function generated_files(output)
    files = String[]
    for (directory, _, names) in walkdir(output), name in names
        path = normpath(joinpath(directory, name))
        relpath(path, output) == "README.md" && continue
        push!(files, path)
    end
    return sort(files)
end

function main(arguments=ARGS)
    profile = option_value(arguments, "--profile", "paper")
    config_root = abspath(option_value(
        arguments, "--config", joinpath(REPOSITORY_ROOT, "config"),
    ))
    isdir(config_root) || error("Configuration directory not found: $config_root")
    output = prepare_output(arguments)
    experiment_files = sort(filter(
        path -> occursin(r"/\d\d_.*\.jl$", path),
        readdir(joinpath(REPOSITORY_ROOT, "experiments"); join=true),
    ))
    length(experiment_files) == 7 || error(
        "Expected seven paper experiments; found $(length(experiment_files)).",
    )

    println()
    println("=== STGPKF paper reproduction ===")
    println("Profile: $profile")
    println("Experiments: $(length(experiment_files))")
    println("Output directory: $output")

    julia = Base.julia_cmd()
    for (index, experiment) in enumerate(experiment_files)
        experiment_name = splitext(basename(experiment))[1]
        experiment_output = joinpath(output, experiment_name)
        command = `$julia --startup-file=no --project=$REPOSITORY_ROOT $experiment --profile $profile --config $config_root --output $experiment_output`
        println()
        println("[$index/$(length(experiment_files))] Launching $(basename(experiment))")
        withenv("GKSwstype" => "100") do
            run(command)
        end
    end

    resolved_configs = [
        joinpath(
            output,
            splitext(basename(path))[1],
            "configs",
            "$(splitext(basename(path))[1]).toml",
        )
        for path in experiment_files
    ]
    expected = sort(vcat(required_outputs(output), resolved_configs))
    missing = filter(path -> !isfile(path) || filesize(path) == 0, expected)
    isempty(missing) || error(
        "Reproduction finished with missing or empty outputs: $(join(missing, ", "))",
    )
    actual = generated_files(output)
    unexpected = setdiff(actual, expected)
    isempty(unexpected) || error(
        "Reproduction finished with unexpected outputs: $(join(unexpected, ", "))",
    )
    println()
    println("=== Reproduction complete ===")
    println("Profile: $profile")
    println("Experiments completed: $(length(experiment_files))/$(length(experiment_files))")
    println("Verified result files: $(length(expected))")
    println("Resolved configurations: $(length(resolved_configs))")
    println("Output directory: $output")
    return output
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
