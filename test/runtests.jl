using SpatiotemporalGPs
using LinearAlgebra
using Random
using StaticArrays
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
    @test paper["simulation"]["dt"] == 0.05
    @test paper["estimate_sensors"] == [1, 6, 20]
    bound_config = load_experiment_config(
        config_root, "03_discrete_continuous_bounds", "paper",
    )
    @test bound_config["measurement"]["fix_sigma_c"] === true
    @test bound_config["measurement"]["sigma_c_squared"] == 0.2
    @test bound_config["steady_state_tolerance"] == 1e-6
    @test bound_config["steady_state_maximum_iterations"] == 20_000
    @test bound_config["measurement"]["std"]^2 *
          bound_config["simulation"]["dt"] ≈
          bound_config["measurement"]["sigma_c_squared"]
    experiment_names = [
        "01_field_reconstruction",
        "02_clarity_vs_time",
        "03_discrete_continuous_bounds",
        "04_grid_convergence",
        "05_sensor_number_table",
        "06_clarity_vs_sensors",
        "07_noise_rate_tradeoff",
    ]
    @test all(experiment_names) do experiment
        config = load_experiment_config(config_root, experiment, "paper")
        config["measurement"]["fix_sigma_c"] === true
    end
    @test_throws ArgumentError load_experiment_config(
        config_root, "02_clarity_vs_time", "unknown",
    )
    @test_throws ArgumentError load_experiment_config(
        config_root, "missing_experiment", "paper",
    )
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
