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

function prepare_output(arguments, profile)
    requested = option_value(arguments, "--output", nothing)
    if isnothing(requested)
        stamp = Dates.format(now(UTC), "yyyymmdd-HHMMSS")
        output = joinpath(REPOSITORY_ROOT, "results", "$profile-$stamp")
    else
        output = abspath(requested)
        isdir(output) && !isempty(readdir(output)) && error(
            "Output directory is not empty: $output",
        )
    end
    mkpath(output)
    return output
end

function required_outputs(output)
    return [
        joinpath(output, "data", "field_reconstruction.jld2"),
        joinpath(output, "data", "mean_clarity_Nr1.jld2"),
        joinpath(output, "data", "discrete_continuous_bounds.csv"),
        joinpath(output, "data", "grid_convergence.csv"),
        joinpath(output, "data", "sensor_number_table.csv"),
        joinpath(output, "data", "clarity_vs_sensors.csv"),
        joinpath(output, "data", "clarity_grid_Nr1.jld2"),
        joinpath(output, "figures", "estimate_1.pdf"),
        joinpath(output, "figures", "estimate_2.pdf"),
        joinpath(output, "figures", "continuous_vs_discrete2.pdf"),
        joinpath(output, "figures", "ng_converge_clarity.pdf"),
        joinpath(output, "figures", "clarity_vs_nr.pdf"),
        joinpath(output, "figures", "clarity_heatmap_Nr1.pdf"),
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
    output = prepare_output(arguments, profile)
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

    commands = String[]
    julia = Base.julia_cmd()
    for (index, experiment) in enumerate(experiment_files)
        command = `$julia --startup-file=no --project=$REPOSITORY_ROOT $experiment --profile $profile --config $config_root --output $output`
        push!(commands, string(command))
        println()
        println("[$index/$(length(experiment_files))] Launching $(basename(experiment))")
        withenv("GKSwstype" => "100") do
            run(command)
        end
    end

    missing = filter(path -> !isfile(path) || filesize(path) == 0, required_outputs(output))
    isempty(missing) || error("Reproduction finished with missing outputs: $(join(missing, ", "))")
    resolved_configs = [
        joinpath(output, "configs", "$(splitext(basename(path))[1]).toml")
        for path in experiment_files
    ]
    missing_configs = filter(path -> !isfile(path) || filesize(path) == 0, resolved_configs)
    isempty(missing_configs) || error(
        "Reproduction finished with missing resolved configurations: $(join(missing_configs, ", "))",
    )
    metadata_path = write_metadata(output, profile, config_root, commands)
    println()
    println("=== Reproduction complete ===")
    println("Profile: $profile")
    println("Experiments completed: $(length(experiment_files))/$(length(experiment_files))")
    println("Verified result files: $(length(required_outputs(output)))")
    println("Resolved configurations: $(length(resolved_configs))")
    println("Metadata files: 1")
    println("Output directory: $output")
    return output
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
