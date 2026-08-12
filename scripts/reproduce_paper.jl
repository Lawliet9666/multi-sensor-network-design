using Pkg

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(REPOSITORY_ROOT)

using Dates
using TOML

function option_value(arguments, name, default)
    index = findfirst(==(name), arguments)
    isnothing(index) && return default
    index < length(arguments) || error("Missing value after $name")
    return arguments[index + 1]
end

function repository_revision(root)
    try
        git_root = readchomp(pipeline(
            `git -C $root rev-parse --show-toplevel`; stderr=devnull,
        ))
        realpath(git_root) == realpath(root) || return "unversioned"
        return readchomp(pipeline(`git -C $root rev-parse HEAD`; stderr=devnull))
    catch
        return "unversioned"
    end
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
        joinpath(output, "03_discrete_continuous_bounds", "data", "discrete_continuous_bounds.csv"),
        joinpath(output, "04_grid_convergence", "data", "grid_convergence.csv"),
        joinpath(output, "05_sensor_number_table", "data", "sensor_number_table.csv"),
        joinpath(output, "06_clarity_vs_sensors", "data", "clarity_vs_sensors.csv"),
        joinpath(output, "07_noise_rate_tradeoff", "data", "clarity_grid_Nr1.jld2"),
        joinpath(output, "01_field_reconstruction", "figures", "estimate_1.pdf"),
        joinpath(output, "02_clarity_vs_time", "figures", "estimate_2.pdf"),
        joinpath(output, "03_discrete_continuous_bounds", "figures", "continuous_vs_discrete2.pdf"),
        joinpath(output, "04_grid_convergence", "figures", "ng_converge_clarity.pdf"),
        joinpath(output, "06_clarity_vs_sensors", "figures", "clarity_vs_nr.pdf"),
        joinpath(output, "07_noise_rate_tradeoff", "figures", "clarity_heatmap_Nr1.pdf"),
    ]
end

function write_metadata(output, profile, config_root, commands)
    metadata = Dict(
        "profile" => profile,
        "seed" => TOML.parsefile(joinpath(config_root, "common.toml"))["seed"],
        "config_root" => abspath(config_root),
        "resolved_config_directory" => "configs",
        "generated_at_utc" => string(now(UTC)),
        "julia_version" => string(VERSION),
        "repository_revision" => repository_revision(REPOSITORY_ROOT),
        "commands" => commands,
    )
    metadata_path = joinpath(output, "metadata.toml")
    open(metadata_path, "w") do stream
        TOML.print(stream, metadata; sorted=true)
    end
    println("Saved run metadata (TOML): $(abspath(metadata_path))")
    return abspath(metadata_path)
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

    metadata_paths = String[]
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
        push!(metadata_paths, write_metadata(
            experiment_output, profile, config_root, [string(command)],
        ))
    end

    missing = filter(path -> !isfile(path) || filesize(path) == 0, required_outputs(output))
    isempty(missing) || error("Reproduction finished with missing outputs: $(join(missing, ", "))")
    resolved_configs = [
        joinpath(
            output,
            splitext(basename(path))[1],
            "configs",
            "$(splitext(basename(path))[1]).toml",
        )
        for path in experiment_files
    ]
    missing_configs = filter(path -> !isfile(path) || filesize(path) == 0, resolved_configs)
    isempty(missing_configs) || error(
        "Reproduction finished with missing resolved configurations: $(join(missing_configs, ", "))",
    )
    println()
    println("=== Reproduction complete ===")
    println("Profile: $profile")
    println("Experiments completed: $(length(experiment_files))/$(length(experiment_files))")
    println("Verified result files: $(length(required_outputs(output)))")
    println("Resolved configurations: $(length(resolved_configs))")
    println("Metadata files: $(length(metadata_paths))")
    println("Output directory: $output")
    return output
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
