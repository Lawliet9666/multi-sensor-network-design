module SpatiotemporalGPs

using Combinatorics
using Interpolations
using Kronecker
using LinearAlgebra
using MatrixEquations
using Plots
using Random
using StaticArrays
using Statistics
using StatsBase
using TOML

include("types.jl")
include("kernels.jl")
include("kalman_filter.jl")
include("stgpkf.jl")
include("simulation.jl")
include("bounds.jl")
include("plotting.jl")

export Matern, Matern12, kernel_matrix, state_space_model
export STGPKFProblem, STGPKFProblemContinuous, SpatiotemporalData2D
export KFState, get_μ, get_Σ, predict, correct
export stgpkf_initialize, stgpkf_predict, stgpkf_correct
export get_estimate, get_estimate_clarity, get_clarity, get_covariance
export generate_spatiotemporal_process, measure, rand_point
export discrete_covariance_bound, discrete_covariance_bound_steady_state
export continuous_covariance_bound, analytic_covariance_bound
export G_from_samples, G_analytic
export precompute_sensor_configs, sample_point_sets
export simulate_expected_covariance, expected_covariance, linear_operator
export clarity_metrics, minimum_sensor_count
export load_experiment_config, measurement_std_at_step, build_problem
export plot_field_reconstruction, plot_mean_clarity, plot_bound_comparison
export plot_grid_convergence, plot_sensor_table_curve, plot_noise_rate_tradeoff

end
