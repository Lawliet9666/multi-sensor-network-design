using SpatiotemporalGPs
using LinearAlgebra
using Plots
using Random
using StaticArrays
using Statistics
using Test

@testset "Matérn-1/2 model" begin
    kernel = Matern(1 / 2, 2.0, 60.0)
    @test kernel(0.0, 0.0) ≈ 4.0
    @test kernel(0.0, 60.0) ≈ 4 / exp(1)
    @test_throws ArgumentError Matern(3 / 2, 1.0, 1.0)

    model = state_space_model(kernel, 0.05)
    steady_state = only(model.W) / (1 - only(model.Φ)^2)
    @test steady_state ≈ 30.0 atol=1e-10
    @test only(model.C)^2 * steady_state ≈ 4.0 atol=1e-10
end

@testset "CPU square-root filter" begin
    state = KFState(; μ=[0.0], Σ=reshape([1.0], 1, 1))
    corrected = correct(state, [1.0], reshape([1.0], 1, 1), reshape([1.0], 1, 1))
    @test only(get_μ(corrected)) ≈ 0.5 atol=1e-12
    @test only(get_Σ(corrected)) ≈ 0.5 atol=1e-12
end

@testset "Paper STGPKF path" begin
    points = vec([@SVector[x, y] for x in 0.0:1.0:2.0, y in 0.0:1.0:2.0])
    spatial = Matern(1 / 2, 1.0, 2.0)
    temporal = Matern(1 / 2, 2.0, 60.0)
    problem = STGPKFProblem(points, spatial, temporal, 0.5)
    state = stgpkf_initialize(problem)
    point = [@SVector[0.5, 0.5]]
    corrected = stgpkf_correct(problem, state, point, [0.2], reshape([4.0], 1, 1))
    predicted = stgpkf_predict(problem, corrected)

    @test length(get_estimate(problem, predicted)) == length(points)
    @test all(value -> 0 < value <= 1, get_estimate_clarity(problem, predicted))
    @test isposdef(Matrix(get_Σ(predicted)))
end

@testset "Expected spatial-mean clarity" begin
    points = vec([@SVector[x, 0.0] for x in 0.0:1.0:1.0])
    spatial = Matern(1 / 2, 1.0, 2.0)
    temporal = Matern(1 / 2, 2.0, 60.0)
    problem = STGPKFProblem(points, spatial, temporal, 0.05)
    trial_covariances = [
        [0.4 0.0; 0.0 2.5],
        [3.0 0.0; 0.0 0.7],
    ]
    statistics = SpatiotemporalGPs._expected_clarity_statistics(
        problem, trial_covariances,
    )
    direct_trial_mean = mean(
        spatial_mean_clarity.(Ref(problem), trial_covariances),
    )

    @test statistics.expected_spatial_mean_clarity_estimate ≈ direct_trial_mean
    @test statistics.expected_spatial_mean_clarity_estimate >
          statistics.spatial_mean_clarity_of_expected_covariance_estimate
    @test statistics.spatial_mean_clarity_sample_std > 0
    @test statistics.spatial_mean_clarity_standard_error ≈
          statistics.spatial_mean_clarity_sample_std / sqrt(2)
    @test statistics.monte_carlo_trials == 2
    @test_throws ArgumentError SpatiotemporalGPs._expected_clarity_statistics(
        problem, Matrix{Float64}[],
    )

    xs = collect(0.0:1.0:3.0)
    ys = collect(0.0:1.0:3.0)
    simulation_problem = STGPKFProblem(
        vec([@SVector[x, y] for x in xs, y in ys]),
        spatial,
        temporal,
        0.05,
    )
    data = generate_spatiotemporal_process(
        MersenneTwister(3), xs, ys, 0.05, 0.2, spatial, temporal,
    )
    point_sets = sample_point_sets(
        MersenneTwister(5), simulation_problem.pts, 1; configuration_count=8,
    )
    first_result = simulate_expected_clarity(
        simulation_problem,
        data,
        point_sets,
        2.0,
        reshape([4.0], 1, 1);
        trials=3,
        seed=7,
    )
    repeated_result = simulate_expected_clarity(
        simulation_problem,
        data,
        point_sets,
        2.0,
        reshape([4.0], 1, 1);
        trials=3,
        seed=7,
    )
    @test first_result == repeated_result
    @test 0 < first_result.expected_spatial_mean_clarity_estimate <= 1
    @test first_result.expected_spatial_mean_clarity_estimate >=
          first_result.spatial_mean_clarity_of_expected_covariance_estimate
    @test first_result.monte_carlo_trials == 3
end

@testset "I.i.d. grid-point sampling" begin
    point = @SVector[0.0, 0.0]
    point_sets = sample_point_sets(
        MersenneTwister(1), [point], 2; configuration_count=3,
    )

    @test length(point_sets) == 3
    @test all(points -> points == [point, point], point_sets)
    @test_throws ArgumentError sample_point_sets(
        MersenneTwister(1), [point], 0; configuration_count=1,
    )
end

@testset "Analytic network-design bound" begin
    config_root = joinpath(@__DIR__, "..", "config")
    config = load_experiment_config(
        config_root, "04_grid_convergence", "paper",
    )
    measurement_std = measurement_std_at_step(
        config["measurement"]["std"],
        config["measurement"]["sigma_c_squared"],
        config["grid_dt"],
        config["measurement"]["fix_sigma_c"],
    )
    problem, _, _, _, _ = build_problem(config; continuous=true)
    one_sensor = clarity_metrics(
        problem, 1, measurement_std, config["grid_dt"],
    )
    two_sensors = clarity_metrics(
        problem, 2, measurement_std, config["grid_dt"],
    )
    @test two_sensors.mean_field_clarity > one_sensor.mean_field_clarity
    @test minimum_sensor_count(
        problem, 0.5, measurement_std, config["grid_dt"],
    ).sensor_count >= 1
end

@testset "Experiment configuration composition" begin
    config_root = joinpath(@__DIR__, "..", "config")
    paper = load_experiment_config(config_root, "02_clarity_vs_time", "paper")

    @test paper["temporal_kernel"]["length_scale"] == 60.0
    @test paper["domain"]["max"] == 5.0
    @test paper["domain"]["dx"] == 0.5
    @test paper["simulation"]["dt"] == 0.02
    @test paper["estimate_sensors"] == [1, 6, 20]
    bound_config = load_experiment_config(
        config_root, "03_discrete_continuous_bounds", "paper",
    )
    @test bound_config["measurement"]["fix_sigma_c"] === true
    @test bound_config["measurement"]["sigma_c_squared"] == 0.2
    @test bound_config["trajectory_step"] == 0.02
    @test bound_config["monte_carlo_trials"] == 30
    @test bound_config["trajectory_step"] in bound_config["bound_steps"]
    @test !haskey(bound_config, "steady_state_tolerance")
    @test !haskey(bound_config, "steady_state_maximum_iterations")
    bound_measurement_std = measurement_std_at_step(
        bound_config["measurement"]["std"],
        bound_config["measurement"]["sigma_c_squared"],
        bound_config["trajectory_step"],
        bound_config["measurement"]["fix_sigma_c"],
    )
    @test bound_measurement_std^2 * bound_config["trajectory_step"] ≈
          bound_config["measurement"]["sigma_c_squared"]
    sensor_config = load_experiment_config(
        config_root, "06_clarity_vs_sensors", "paper",
    )
    @test sensor_config["trials"] == 30
    @test 21 in sensor_config["sensor_curve"]
    @test sensor_config["plot_target_clarity"] == 0.8
    fixed_sigma_c_experiments = [
        "01_field_reconstruction",
        "02_clarity_vs_time",
        "03_discrete_continuous_bounds",
        "04_grid_convergence",
        "05_sensor_number_table",
        "06_clarity_vs_sensors",
    ]
    @test all(fixed_sigma_c_experiments) do experiment
        config = load_experiment_config(config_root, experiment, "paper")
        config["measurement"]["fix_sigma_c"] === true
    end
    tradeoff_config = load_experiment_config(
        config_root, "07_noise_rate_tradeoff", "paper",
    )
    @test tradeoff_config["measurement"]["fix_sigma_c"] === false
    @test tradeoff_config["measurement_variance_min"] == 0.01
    @test tradeoff_config["measurement_variance_max"] == 4.0
    @test tradeoff_config["normalized_interval_min"] == 0.1
    @test tradeoff_config["normalized_interval_max"] == 2.0
    @test tradeoff_config["contour_clarity_levels"] == [0.4, 0.5, 0.6]
    @test !haskey(tradeoff_config, "heatmap_sigma_min")
    @test !haskey(tradeoff_config, "heatmap_dt_min")
    @test_throws ArgumentError load_experiment_config(
        config_root, "02_clarity_vs_time", "unknown",
    )
    @test_throws ArgumentError load_experiment_config(
        config_root, "missing_experiment", "paper",
    )
end

@testset "Noise-rate tradeoff parameterization" begin
    config = load_experiment_config(
        joinpath(@__DIR__, "..", "config"),
        "07_noise_rate_tradeoff",
        "paper",
    )
    problem, _, _, _, _ = build_problem(config; continuous=true)
    sensor_count = config["sensor_count"]

    first = clarity_metrics(
        problem,
        sensor_count,
        sqrt(0.5),
        0.4 * sensor_count,
    )
    same_theta = clarity_metrics(
        problem,
        sensor_count,
        sqrt(1.0),
        0.2 * sensor_count,
    )
    lower_theta = clarity_metrics(
        problem,
        sensor_count,
        sqrt(0.5),
        0.8 * sensor_count,
    )

    @test first.sensing_intensity ≈ same_theta.sensing_intensity
    @test first.mean_field_clarity ≈ same_theta.mean_field_clarity
    @test lower_theta.mean_field_clarity < first.mean_field_clarity
end

@testset "Paper plots save PDF only" begin
    mktempdir() do directory
        figure = plot([0.0, 1.0], [0.0, 1.0]; label=false)
        path = SpatiotemporalGPs.save_paper_plot(figure, directory, "paper_plot")
        @test path == joinpath(directory, "paper_plot.pdf")
        @test isfile(path)
        @test filesize(path) > 0
        @test !isfile(joinpath(directory, "paper_plot.svg"))
    end
end

@testset "Finite-horizon covariance comparison" begin
    points = vec([@SVector[x, y] for x in 0.0:1.0:1.0, y in 0.0:1.0:1.0])
    spatial = Matern(1 / 2, 1.0, 2.0)
    temporal = Matern(1 / 2, 2.0, 60.0)
    step = 0.05
    sigma_c_squared = 0.2
    problem = STGPKFProblem(points, spatial, temporal, step)
    continuous_problem = STGPKFProblemContinuous(points, spatial, temporal)
    sensor_models = precompute_sensor_configs(
        problem, 1, (sigma_c_squared / step) * I(1),
    )

    first_result = finite_horizon_covariance_comparison(
        problem,
        continuous_problem,
        sensor_models;
        horizon=0.2,
        trials=4,
        seed=7,
    )
    repeated_result = finite_horizon_covariance_comparison(
        problem,
        continuous_problem,
        sensor_models;
        horizon=0.2,
        trials=4,
        seed=7,
    )

    @test length(first_result.times) == 5
    @test first_result.times == collect(0.0:step:0.2)
    initial_mean = linear_operator(Matrix(get_Σ(stgpkf_initialize(problem))))
    @test first_result.empirical_mean_covariance[1] == initial_mean
    @test first_result.discrete_bound_mean_covariance[1] == initial_mean
    @test first_result.continuous_bound_mean_covariance[1] == initial_mean
    @test first_result == repeated_result
    @test all(isfinite, first_result.empirical_mean_covariance)
    @test all(isfinite, first_result.discrete_bound_mean_covariance)
    @test all(isfinite, first_result.continuous_bound_mean_covariance)
    @test first_result.empirical_max_frobenius_error >= 0
    @test first_result.discrete_max_frobenius_error >= 0
    @test first_result.empirical_max_linear_operator_error == maximum(abs.(
        first_result.empirical_mean_covariance .-
        first_result.continuous_bound_mean_covariance
    ))
    @test first_result.discrete_max_linear_operator_error == maximum(abs.(
        first_result.discrete_bound_mean_covariance .-
        first_result.continuous_bound_mean_covariance
    ))
    @test first_result.empirical_max_relative_linear_operator_error_percent ≈
          100 * maximum(abs.((
              first_result.empirical_mean_covariance .-
              first_result.continuous_bound_mean_covariance
          ) ./ first_result.empirical_mean_covariance))
    @test first_result.discrete_max_relative_linear_operator_error_percent ≈
          100 * maximum(abs.((
              first_result.discrete_bound_mean_covariance .-
              first_result.continuous_bound_mean_covariance
          ) ./ first_result.discrete_bound_mean_covariance))

    trial_sequences = [
        rand(MersenneTwister(7 + 1000 * trial), 1:length(sensor_models), 20)
        for trial in 1:4
    ]
    @test allunique(trial_sequences)
end

@testset "Exact discrete configuration average" begin
    points = vec([@SVector[x, y] for x in 0.0:1.0:1.0, y in 0.0:1.0:1.0])
    problem = STGPKFProblem(
        points,
        Matern(1 / 2, 1.0, 2.0),
        Matern(1 / 2, 2.0, 60.0),
        0.05,
    )
    sensor_models = precompute_sensor_configs(problem, 1, 4.0I(1))
    covariance = Matrix(get_Σ(stgpkf_initialize(problem)))
    observation_vectors, measurement_variances =
        SpatiotemporalGPs._single_sensor_parameters(sensor_models, size(covariance, 1))
    transition_squared = abs2(only(problem.ss_model.Φ))
    process_variance = only(problem.ss_model.W)
    actual = copy(covariance)
    SpatiotemporalGPs._single_sensor_bound_step!(
        actual,
        observation_vectors,
        measurement_variances,
        transition_squared,
        process_variance,
        zeros(size(covariance, 1)),
        zeros(size(covariance)),
    )

    expected = zeros(size(covariance))
    for (observation, measurement_covariance) in sensor_models
        projected = covariance * observation'
        corrected = covariance - projected * (
            (measurement_covariance + observation * projected) \ projected'
        )
        expected .+= transition_squared .* corrected
        expected .+= process_variance .* I(size(covariance, 1))
    end
    expected ./= length(sensor_models)
    @test actual ≈ expected atol=1e-10 rtol=1e-10
    @test issymmetric(actual)
    @test isposdef(Symmetric(actual))

    selected = copy(covariance)
    SpatiotemporalGPs._single_sensor_covariance_step!(
        selected,
        first(observation_vectors),
        first(measurement_variances),
        transition_squared,
        process_variance,
        zeros(size(covariance, 1)),
    )
    @test issymmetric(selected)
    @test isposdef(Symmetric(selected))
end

@testset "Continuous Riccati RK4" begin
    initial = 5.0
    drift = -0.1
    process_variance = 1.0
    information = 0.2
    duration = 0.2
    modes = [initial]
    workspaces = [zeros(1) for _ in 1:5]
    SpatiotemporalGPs._rk4_riccati_modes_step!(
        modes,
        [information],
        drift,
        process_variance,
        duration,
        workspaces...,
    )

    rate = sqrt(drift^2 + process_variance * information)
    positive_root = (drift + rate) / information
    negative_root = (drift - rate) / information
    ratio = (
        (initial - positive_root) / (initial - negative_root)
    ) * exp(-2 * rate * duration)
    analytic = (positive_root - ratio * negative_root) / (1 - ratio)
    @test only(modes) ≈ analytic atol=1e-8 rtol=1e-8
end

@testset "Finite-horizon discrete-bound convergence" begin
    points = vec([@SVector[x, y] for x in 0.0:1.0:1.0, y in 0.0:1.0:1.0])
    spatial = Matern(1 / 2, 1.0, 2.0)
    temporal = Matern(1 / 2, 2.0, 60.0)
    continuous_problem = STGPKFProblemContinuous(points, spatial, temporal)
    errors = Float64[]
    for step in (0.1, 0.05, 0.025)
        problem = STGPKFProblem(points, spatial, temporal, step)
        sensor_models = precompute_sensor_configs(
            problem, 1, (0.2 / step) * I(1),
        )
        result = finite_horizon_covariance_comparison(
            problem,
            continuous_problem,
            sensor_models;
            horizon=0.2,
            trials=1,
            seed=1,
        )
        push!(errors, result.discrete_max_frobenius_error)
    end
    @test errors[1] > errors[2] > errors[3]
end

@testset "Discrete covariance-bound steady state" begin
    points = vec([@SVector[x, y] for x in 0.0:1.0:2.0, y in 0.0:1.0:2.0])
    problem = STGPKFProblem(
        points,
        Matern(1 / 2, 1.0, 2.0),
        Matern(1 / 2, 2.0, 60.0),
        0.5,
    )
    sensor_models = precompute_sensor_configs(problem, 1, reshape([0.4], 1, 1))
    result = discrete_covariance_bound_steady_state(
        problem,
        sensor_models;
        tolerance=1e-6,
        maximum_iterations=20_000,
        sample_count=length(sensor_models),
        seed=1,
        beta=1.0,
    )
    @test result.residual < 1e-6
    @test result.iterations <= 20_000
    @test isposdef(Matrix(result.covariance))
end

@testset "Measurement-noise scaling" begin
    @test measurement_std_at_step(2.0, 0.2, 0.05, true) ≈ 2.0
    @test measurement_std_at_step(2.0, 0.2, 0.20, true) ≈ 1.0
    @test measurement_std_at_step(2.0, 0.2, 0.20, false) ≈ 2.0
    @test_throws ArgumentError measurement_std_at_step(0.0, 0.2, 0.20, true)
    @test_throws ArgumentError measurement_std_at_step(2.0, 0.0, 0.20, true)
    @test_throws ArgumentError measurement_std_at_step(2.0, 0.2, 0.0, true)
end
