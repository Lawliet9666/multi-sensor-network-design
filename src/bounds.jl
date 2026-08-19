linear_operator(covariance::AbstractMatrix) = mean(diag(covariance))

function expected_covariance(trajectories)
    isempty(trajectories) && throw(ArgumentError("At least one trajectory is required."))
    steps = length(first(trajectories))
    return [mean([trajectory[index] for trajectory in trajectories]) for index in 1:steps]
end

function precompute_sensor_configs(problem, sensor_count::Int, measurement_covariance)
    grid_count = length(problem.pts)
    1 <= sensor_count <= grid_count || throw(ArgumentError(
        "Sensor count must be between 1 and $grid_count.",
    ))
    indices = collect(combinations(1:grid_count, sensor_count))
    matrices = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, length(indices))
    for (position, selection) in enumerate(indices)
        points = [problem.pts[index] for index in selection]
        matrices[position] = sensor_matrices(problem, points, measurement_covariance)
    end
    return matrices
end

function sample_point_sets(rng, points, sensor_count::Int; configuration_count::Int)
    grid_count = length(points)
    grid_count > 0 || throw(ArgumentError("At least one grid point is required."))
    sensor_count > 0 || throw(ArgumentError("Sensor count must be positive."))
    configuration_count > 0 || throw(ArgumentError("Configuration count must be positive."))
    return [
        [points[index] for index in sample(rng, 1:grid_count, sensor_count; replace=true)]
        for _ in 1:configuration_count
    ]
end

function discrete_bound_step(
    problem, covariance, sensor_models, selected_indices::AbstractVector{<:Integer};
                             beta::Real)
    grid_count = length(problem.pts)
    transition = I(grid_count) ⊗ problem.ss_model.Φ
    process_noise = I(grid_count) ⊗ problem.ss_model.W
    predicted = transition * covariance * transition'
    cross_left = transition * covariance
    cross_right = covariance * transition'

    correction = zeros(eltype(covariance), size(covariance))
    for index in selected_indices
        H, V = sensor_models[index]
        correction .+= cross_left * H' * ((V + H * covariance * H') \ (H * cross_right))
    end
    correction ./= length(selected_indices)
    next_covariance = process_noise + predicted - beta .* correction
    next_covariance = 0.5 .* (next_covariance + next_covariance')
    next_covariance += 1e-12I
    return Matrix(next_covariance)
end

function discrete_bound_step(problem, covariance, sensor_models, rng::AbstractRNG;
                             sample_count::Int, beta::Real)
    count = min(sample_count, length(sensor_models))
    selected_indices = sample(rng, 1:length(sensor_models), count; replace=false)
    return discrete_bound_step(
        problem, covariance, sensor_models, selected_indices; beta=beta,
    )
end

function discrete_covariance_bound(problem::STGPKFProblem, times, sensor_models;
                                   sample_count::Int=500, seed::Int=1, beta::Real=1.0)
    covariance = Matrix(get_Σ(stgpkf_initialize(problem)))
    history = Matrix{Float64}[]
    for step_index in eachindex(times)
        rng = MersenneTwister(seed + step_index)
        covariance = discrete_bound_step(
            problem, covariance, sensor_models, rng;
            sample_count=sample_count, beta=beta,
        )
        push!(history, copy(covariance))
    end
    return history
end

function discrete_covariance_bound_steady_state(
    problem::STGPKFProblem,
    sensor_models;
    tolerance::Real=1e-6,
    maximum_iterations::Int=20_000,
    sample_count::Int=500,
    seed::Int=1,
    beta::Real=1.0,
)
    tolerance > 0 || throw(ArgumentError("Steady-state tolerance must be positive."))
    maximum_iterations > 0 || throw(ArgumentError(
        "Maximum steady-state iterations must be positive.",
    ))
    sample_count > 0 || throw(ArgumentError("Sample count must be positive."))
    isempty(sensor_models) && throw(ArgumentError("At least one sensor model is required."))

    count = min(sample_count, length(sensor_models))
    selected_indices = sample(
        MersenneTwister(seed), 1:length(sensor_models), count; replace=false,
    )
    covariance = Matrix(get_Σ(stgpkf_initialize(problem)))

    for iteration in 1:maximum_iterations
        next_covariance = discrete_bound_step(
            problem, covariance, sensor_models, selected_indices; beta=beta,
        )
        residual = norm(next_covariance - covariance) / max(1.0, norm(covariance))
        if residual < tolerance
            return (
                covariance=next_covariance,
                iterations=iteration,
                residual=residual,
            )
        end
        covariance = next_covariance
    end

    error(
        "Discrete covariance bound did not converge within " *
        "$maximum_iterations iterations at tolerance $tolerance.",
    )
end

function _finite_horizon_step_count(horizon::Real, step::Real)
    horizon > 0 || throw(ArgumentError("Finite horizon must be positive."))
    step > 0 || throw(ArgumentError("Sampling interval must be positive."))
    step_count = round(Int, horizon / step)
    isapprox(step_count * step, horizon; rtol=1e-10, atol=1e-12) || throw(
        ArgumentError("Finite horizon $horizon must be an integer multiple of step $step."),
    )
    return step_count
end

function _single_sensor_parameters(sensor_models, state_dimension::Int)
    isempty(sensor_models) && throw(ArgumentError("At least one sensor model is required."))
    observation_vectors = Vector{Vector{Float64}}(undef, length(sensor_models))
    measurement_variances = Vector{Float64}(undef, length(sensor_models))
    for (index, (observation, measurement_covariance)) in enumerate(sensor_models)
        size(observation) == (1, state_dimension) || throw(ArgumentError(
            "Finite-horizon comparison currently requires one sensor per configuration.",
        ))
        size(measurement_covariance) == (1, 1) || throw(ArgumentError(
            "Finite-horizon comparison requires scalar measurement covariance.",
        ))
        variance = only(measurement_covariance)
        isfinite(variance) && variance > 0 || throw(ArgumentError(
            "Measurement covariance must be finite and positive.",
        ))
        observation_vectors[index] = vec(copy(observation))
        measurement_variances[index] = variance
    end
    return observation_vectors, measurement_variances
end

function _check_innovation_variance(value::Real)
    isfinite(value) && value > 0 || error(
        "Encountered a non-positive or non-finite innovation variance: $value",
    )
    return value
end

function _single_sensor_covariance_step!(
    covariance::Matrix{Float64}, observation::Vector{Float64},
    measurement_variance::Float64, transition_squared::Float64,
    process_variance::Float64, projected_covariance::Vector{Float64},
)
    mul!(projected_covariance, covariance, observation)
    innovation_variance = _check_innovation_variance(
        measurement_variance + dot(observation, projected_covariance),
    )
    covariance .*= transition_squared
    BLAS.ger!(
        -transition_squared / innovation_variance,
        projected_covariance,
        projected_covariance,
        covariance,
    )
    @inbounds for index in axes(covariance, 1)
        covariance[index, index] += process_variance
    end
    LinearAlgebra.copytri!(covariance, 'U')
    return covariance
end

function _single_sensor_bound_step!(
    covariance::Matrix{Float64}, observation_vectors,
    measurement_variances, transition_squared::Float64,
    process_variance::Float64, projected_covariance::Vector{Float64},
    mean_correction::Matrix{Float64},
)
    fill!(mean_correction, 0.0)
    for index in eachindex(observation_vectors)
        observation = observation_vectors[index]
        mul!(projected_covariance, covariance, observation)
        innovation_variance = _check_innovation_variance(
            measurement_variances[index] + dot(observation, projected_covariance),
        )
        BLAS.ger!(
            inv(innovation_variance),
            projected_covariance,
            projected_covariance,
            mean_correction,
        )
    end
    covariance .*= transition_squared
    covariance .-= (
        transition_squared / length(observation_vectors)
    ) .* mean_correction
    @inbounds for index in axes(covariance, 1)
        covariance[index, index] += process_variance
    end
    LinearAlgebra.copytri!(covariance, 'U')
    return covariance
end

function _validate_isotropic_continuous_model(problem, initial_covariance)
    size(problem.ss_model.A) == (1, 1) || throw(ArgumentError(
        "Finite-horizon comparison currently supports the paper's scalar temporal model.",
    ))
    size(problem.ss_model.B) == (1, 1) || throw(ArgumentError(
        "Finite-horizon comparison currently supports scalar process excitation.",
    ))
    state_dimension = size(initial_covariance, 1)
    initial_scale = tr(initial_covariance) / state_dimension
    isotropic_initial = Matrix(initial_covariance - initial_scale * I)
    norm(isotropic_initial) <= 1e-10 * max(1.0, norm(initial_covariance)) || throw(
        ArgumentError("Initial covariance must be isotropic for the finite-horizon solver."),
    )
    return (
        drift=Float64(only(problem.ss_model.A)),
        process_variance=Float64(abs2(only(problem.ss_model.B))),
        initial_scale=Float64(initial_scale),
    )
end

function _riccati_mode_derivative(value, information, drift, process_variance)
    return 2 * drift * value + process_variance - information * value^2
end

function _rk4_riccati_modes_step!(
    modes, information_eigenvalues, drift, process_variance, step,
    k1, k2, k3, k4, intermediate,
)
    substep_count = max(1, ceil(Int, step / 0.01))
    substep = step / substep_count
    for _ in 1:substep_count
        @inbounds for index in eachindex(modes)
            k1[index] = _riccati_mode_derivative(
                modes[index], information_eigenvalues[index], drift, process_variance,
            )
            intermediate[index] = modes[index] + substep * k1[index] / 2
        end
        @inbounds for index in eachindex(modes)
            k2[index] = _riccati_mode_derivative(
                intermediate[index], information_eigenvalues[index], drift, process_variance,
            )
            intermediate[index] = modes[index] + substep * k2[index] / 2
        end
        @inbounds for index in eachindex(modes)
            k3[index] = _riccati_mode_derivative(
                intermediate[index], information_eigenvalues[index], drift, process_variance,
            )
            intermediate[index] = modes[index] + substep * k3[index]
        end
        @inbounds for index in eachindex(modes)
            k4[index] = _riccati_mode_derivative(
                intermediate[index], information_eigenvalues[index], drift, process_variance,
            )
            modes[index] += substep * (
                k1[index] + 2k2[index] + 2k3[index] + k4[index]
            ) / 6
            isfinite(modes[index]) && modes[index] > 0 || error(
                "Continuous Riccati integration produced an invalid covariance mode.",
            )
        end
    end
    return modes
end

function _covariance_from_modes!(covariance, scaled_eigenvectors, eigenvectors, modes)
    @inbounds for column in axes(eigenvectors, 2), row in axes(eigenvectors, 1)
        scaled_eigenvectors[row, column] = eigenvectors[row, column] * modes[column]
    end
    mul!(covariance, scaled_eigenvectors, eigenvectors')
    LinearAlgebra.copytri!(covariance, 'U')
    return covariance
end

"""
    finite_horizon_covariance_comparison(
        problem, continuous_problem, sensor_models; horizon, trials, seed,
    )

Compare the Monte Carlo estimate of the expected discrete covariance, its exact
uniform-configuration discrete upper bound, and the continuous Riccati solution
on a shared finite time grid. The returned scalar histories include the common
initial covariance at `t = 0`; both maximum Frobenius errors and maximum errors
after applying `linear_operator` are evaluated over the complete horizon. The
relative linear-operator errors are reported as percentages relative to the
empirical covariance and discrete bound, respectively.

This paper experiment uses one grid-point sensor and the scalar Matérn-1/2
temporal model. Unsupported sensing or temporal models fail explicitly.
"""
function finite_horizon_covariance_comparison(
    problem::STGPKFProblem,
    continuous_problem::STGPKFProblemContinuous,
    sensor_models;
    horizon::Real,
    trials::Int,
    seed::Int,
)
    trials > 0 || throw(ArgumentError("Monte Carlo trial count must be positive."))
    step = problem.ΔT
    step_count = _finite_horizon_step_count(horizon, step)
    initial_covariance = Matrix(get_Σ(stgpkf_initialize(problem)))
    state_dimension = size(initial_covariance, 1)
    observation_vectors, measurement_variances = _single_sensor_parameters(
        sensor_models, state_dimension,
    )
    size(problem.ss_model.Φ) == (1, 1) || throw(ArgumentError(
        "Finite-horizon comparison currently supports the paper's scalar temporal model.",
    ))
    size(problem.ss_model.W) == (1, 1) || throw(ArgumentError(
        "Finite-horizon comparison currently supports scalar process covariance.",
    ))
    transition_squared = Float64(abs2(only(problem.ss_model.Φ)))
    discrete_process_variance = Float64(only(problem.ss_model.W))

    continuous_parameters = _validate_isotropic_continuous_model(
        continuous_problem, initial_covariance,
    )
    information = Matrix(G_from_samples(sensor_models, state_dimension, step))
    information_decomposition = eigen(Symmetric(information))
    information_eigenvalues = information_decomposition.values
    minimum(information_eigenvalues) >= -1e-10 || error(
        "The averaged information matrix must be positive semidefinite.",
    )
    information_eigenvalues = max.(information_eigenvalues, 0.0)
    eigenvectors = Matrix(information_decomposition.vectors)

    trial_covariances = [copy(initial_covariance) for _ in 1:trials]
    trial_rngs = [MersenneTwister(seed + 1000 * trial) for trial in 1:trials]
    bound_covariance = copy(initial_covariance)
    empirical_covariance = similar(initial_covariance)
    continuous_covariance = copy(initial_covariance)
    difference = similar(initial_covariance)
    projected_covariances = [zeros(state_dimension) for _ in 1:trials]
    bound_projection = zeros(state_dimension)
    mean_correction = zeros(state_dimension, state_dimension)

    modes = fill(continuous_parameters.initial_scale, state_dimension)
    k1, k2, k3, k4, intermediate = (zeros(state_dimension) for _ in 1:5)
    scaled_eigenvectors = similar(eigenvectors)

    times = collect(range(0.0; step=step, length=step_count + 1))
    empirical_mean_covariance = Vector{Float64}(undef, step_count + 1)
    discrete_bound_mean_covariance = similar(empirical_mean_covariance)
    continuous_bound_mean_covariance = similar(empirical_mean_covariance)
    initial_mean = linear_operator(initial_covariance)
    empirical_mean_covariance[1] = initial_mean
    discrete_bound_mean_covariance[1] = initial_mean
    continuous_bound_mean_covariance[1] = initial_mean
    empirical_max_frobenius_error = 0.0
    discrete_max_frobenius_error = 0.0
    empirical_max_linear_operator_error = 0.0
    discrete_max_linear_operator_error = 0.0
    empirical_max_relative_linear_operator_error_percent = 0.0
    discrete_max_relative_linear_operator_error_percent = 0.0

    for time_index in 2:(step_count + 1)
        fill!(empirical_covariance, 0.0)
        for trial in 1:trials
            configuration = rand(trial_rngs[trial], eachindex(observation_vectors))
            _single_sensor_covariance_step!(
                trial_covariances[trial],
                observation_vectors[configuration],
                measurement_variances[configuration],
                transition_squared,
                discrete_process_variance,
                projected_covariances[trial],
            )
            empirical_covariance .+= trial_covariances[trial]
        end
        empirical_covariance ./= trials
        all(isfinite, empirical_covariance) || error(
            "Monte Carlo covariance average contains non-finite values.",
        )

        _single_sensor_bound_step!(
            bound_covariance,
            observation_vectors,
            measurement_variances,
            transition_squared,
            discrete_process_variance,
            bound_projection,
            mean_correction,
        )
        all(isfinite, bound_covariance) || error(
            "Discrete covariance bound contains non-finite values.",
        )
        _rk4_riccati_modes_step!(
            modes,
            information_eigenvalues,
            continuous_parameters.drift,
            continuous_parameters.process_variance,
            step,
            k1,
            k2,
            k3,
            k4,
            intermediate,
        )
        _covariance_from_modes!(
            continuous_covariance, scaled_eigenvectors, eigenvectors, modes,
        )

        empirical_mean_covariance[time_index] = linear_operator(empirical_covariance)
        discrete_bound_mean_covariance[time_index] = linear_operator(bound_covariance)
        continuous_bound_mean_covariance[time_index] = mean(modes)
        empirical_reference = empirical_mean_covariance[time_index]
        discrete_reference = discrete_bound_mean_covariance[time_index]
        isfinite(empirical_reference) && empirical_reference > 0 || error(
            "Empirical mean covariance must be finite and positive.",
        )
        isfinite(discrete_reference) && discrete_reference > 0 || error(
            "Discrete-bound mean covariance must be finite and positive.",
        )
        empirical_linear_operator_error = abs(
            empirical_reference - continuous_bound_mean_covariance[time_index],
        )
        discrete_linear_operator_error = abs(
            discrete_reference - continuous_bound_mean_covariance[time_index],
        )
        empirical_max_linear_operator_error = max(
            empirical_max_linear_operator_error,
            empirical_linear_operator_error,
        )
        discrete_max_linear_operator_error = max(
            discrete_max_linear_operator_error,
            discrete_linear_operator_error,
        )
        empirical_max_relative_linear_operator_error_percent = max(
            empirical_max_relative_linear_operator_error_percent,
            100 * empirical_linear_operator_error / empirical_reference,
        )
        discrete_max_relative_linear_operator_error_percent = max(
            discrete_max_relative_linear_operator_error_percent,
            100 * discrete_linear_operator_error / discrete_reference,
        )
        difference .= empirical_covariance .- continuous_covariance
        empirical_max_frobenius_error = max(
            empirical_max_frobenius_error, norm(difference),
        )
        difference .= bound_covariance .- continuous_covariance
        discrete_max_frobenius_error = max(
            discrete_max_frobenius_error, norm(difference),
        )
    end

    return (
        times=times,
        empirical_mean_covariance=empirical_mean_covariance,
        discrete_bound_mean_covariance=discrete_bound_mean_covariance,
        continuous_bound_mean_covariance=continuous_bound_mean_covariance,
        empirical_max_frobenius_error=empirical_max_frobenius_error,
        discrete_max_frobenius_error=discrete_max_frobenius_error,
        empirical_max_linear_operator_error=empirical_max_linear_operator_error,
        discrete_max_linear_operator_error=discrete_max_linear_operator_error,
        empirical_max_relative_linear_operator_error_percent=
            empirical_max_relative_linear_operator_error_percent,
        discrete_max_relative_linear_operator_error_percent=
            discrete_max_relative_linear_operator_error_percent,
    )
end

function continuous_covariance_bound(A::AbstractMatrix, G::AbstractMatrix, Q::AbstractMatrix)
    covariance, _, _ = arec(A', G, Q)
    return Symmetric(covariance)
end

function analytic_covariance_bound(A::AbstractMatrix, grid_count::Int, G::AbstractMatrix)
    a = only(A)
    inverse_information = inv(G)
    covariance = a * inverse_information + sqrt(a^2 * I(grid_count) + G) * inverse_information
    return Symmetric(covariance)
end

function G_from_samples(sensor_models, grid_count::Int, step::Real)
    information = zeros(grid_count, grid_count)
    for (H, V) in sensor_models
        information .+= H' * ((V * step) \ H)
    end
    return Symmetric(information ./ length(sensor_models))
end

function G_analytic(problem, sensor_count::Int, measurement_std::Real, step::Real)
    measurement_std > 0 || throw(ArgumentError("Measurement noise must be positive."))
    step > 0 || throw(ArgumentError("Sampling interval must be positive."))
    grid_count = length(problem.pts)
    kernel = kernel_matrix(problem.ks, problem.pts)
    observation = I(grid_count) ⊗ problem.ss_model.C
    scale = sensor_count / (measurement_std^2 * step * grid_count)
    return Symmetric(scale .* (observation' * kernel * observation))
end

"""Cached spatial spectrum for the Theorem 20 analytical design criterion."""
struct AnalyticClarityCache
    normalized_kernel_eigenvalues::Vector{Float64}
    grid_count::Int
    drift::Float64
    output_scale::Float64
end

function analytic_clarity_cache(problem::STGPKFProblemContinuous)
    size(problem.ss_model.A) == (1, 1) || throw(ArgumentError(
        "The spectral clarity criterion requires a scalar temporal state.",
    ))
    size(problem.ss_model.C) == (1, 1) || throw(ArgumentError(
        "The spectral clarity criterion requires C = C0 I.",
    ))
    drift = Float64(only(problem.ss_model.A))
    output_scale = Float64(only(problem.ss_model.C))
    drift < 0 || throw(ArgumentError(
        "The spectral clarity criterion requires a stable scalar drift.",
    ))
    output_scale != 0 || throw(ArgumentError(
        "The spectral clarity criterion requires a nonzero output scale.",
    ))

    grid_count = length(problem.pts)
    kernel = kernel_matrix(problem.ks, problem.pts)
    eigenvalues = eigvals(Symmetric(kernel ./ grid_count))
    tolerance = 1e-10 * max(1.0, maximum(abs, eigenvalues))
    minimum(eigenvalues) >= -tolerance || throw(ArgumentError(
        "The spatial kernel must have nonnegative eigenvalues.",
    ))
    eigenvalues = max.(Float64.(eigenvalues), 0.0)
    return AnalyticClarityCache(
        eigenvalues,
        grid_count,
        drift,
        output_scale,
    )
end

function analytic_clarity_metrics(
    cache::AnalyticClarityCache,
    sensor_count::Int,
    measurement_std::Real,
    step::Real,
)
    sensor_count > 0 || throw(ArgumentError("Sensor count must be positive."))
    measurement_std > 0 || throw(ArgumentError("Measurement noise must be positive."))
    step > 0 || throw(ArgumentError("Sampling interval must be positive."))

    sensing_intensity = sensor_count / (measurement_std^2 * step)
    output_scale_squared = cache.output_scale^2
    information_modes = sensing_intensity * output_scale_squared .*
                        cache.normalized_kernel_eigenvalues
    covariance_modes = @. inv(
        sqrt(cache.drift^2 + information_modes) - cache.drift
    )
    normalized_field_modes = output_scale_squared .*
                             cache.normalized_kernel_eigenvalues .* covariance_modes
    mean_field_variance = sum(normalized_field_modes)

    return (
        sensor_count=sensor_count,
        sensing_intensity=sensing_intensity,
        mean_state_variance=mean(covariance_modes),
        max_state_variance=maximum(covariance_modes),
        mean_field_variance=mean_field_variance,
        max_field_variance=cache.grid_count * maximum(normalized_field_modes),
        mean_field_clarity=inv(1 + mean_field_variance),
    )
end

function _simulate_covariance_trial(
    problem::STGPKFProblem,
    data,
    point_sets,
    measurement_std::Real,
    measurement_covariance,
    seed::Int;
    record_trajectory::Bool,
)
    isempty(point_sets) && throw(ArgumentError("At least one sensor configuration is required."))
    rng = MersenneTwister(seed)
    state = stgpkf_initialize(problem)
    trajectory = record_trajectory ? Matrix{Float64}[] : nothing
    for time in data.ts
        points = point_sets[rand(rng, eachindex(point_sets))]
        observations = [
            measure(rng, data, point[1], point[2], time, measurement_std)
            for point in points
        ]
        corrected = stgpkf_correct(
            problem, state, points, observations, measurement_covariance,
        )
        state = stgpkf_predict(problem, corrected)
        record_trajectory && push!(trajectory, Matrix(get_Σ(state)))
    end
    return record_trajectory ? trajectory : Matrix(get_Σ(state))
end

function simulate_expected_covariance(problem::STGPKFProblem, data, point_sets,
                                      measurement_std::Real, measurement_covariance;
                                      trials::Int, seed::Int)
    trials > 0 || throw(ArgumentError("Trial count must be positive."))
    trajectories = [
        _simulate_covariance_trial(
            problem,
            data,
            point_sets,
            measurement_std,
            measurement_covariance,
            seed + 1000 * trial;
            record_trajectory=true,
        )
        for trial in 1:trials
    ]
    return expected_covariance(trajectories)
end

spatial_mean_clarity(problem, covariance) = mean(get_clarity(problem, covariance))

function _expected_clarity_statistics(problem, terminal_covariances)
    isempty(terminal_covariances) && throw(ArgumentError(
        "At least one terminal covariance is required.",
    ))
    trial_clarities = spatial_mean_clarity.(Ref(problem), terminal_covariances)
    monte_carlo_trials = length(trial_clarities)
    sample_std = monte_carlo_trials == 1 ? 0.0 : std(trial_clarities)
    expected_covariance_estimate = mean(terminal_covariances)
    return (
        expected_spatial_mean_clarity_estimate=mean(trial_clarities),
        spatial_mean_clarity_of_expected_covariance_estimate=
            spatial_mean_clarity(problem, expected_covariance_estimate),
        spatial_mean_clarity_sample_std=sample_std,
        spatial_mean_clarity_standard_error=sample_std / sqrt(monte_carlo_trials),
        monte_carlo_trials=monte_carlo_trials,
    )
end

"""
    simulate_expected_clarity(
        problem, data, point_sets, measurement_std, measurement_covariance;
        trials, seed,
    )

Estimate the expected spatially averaged clarity at the terminal simulation
time. Clarity is evaluated separately for every randomized sensing trial before
the trial values are averaged. The returned
`spatial_mean_clarity_of_expected_covariance_estimate` is the distinct surrogate
obtained by applying clarity to the Monte Carlo mean covariance.
"""
function simulate_expected_clarity(
    problem::STGPKFProblem,
    data,
    point_sets,
    measurement_std::Real,
    measurement_covariance;
    trials::Int,
    seed::Int,
)
    trials > 0 || throw(ArgumentError("Trial count must be positive."))
    terminal_covariances = [
        _simulate_covariance_trial(
            problem,
            data,
            point_sets,
            measurement_std,
            measurement_covariance,
            seed + 1000 * trial;
            record_trajectory=false,
        )
        for trial in 1:trials
    ]
    return _expected_clarity_statistics(problem, terminal_covariances)
end

function clarity_metrics(problem::STGPKFProblemContinuous, sensor_count::Int,
                         measurement_std::Real, step::Real)
    grid_count = length(problem.pts)
    information = G_analytic(problem, sensor_count, measurement_std, step)
    covariance = analytic_covariance_bound(problem.ss_model.A, grid_count, information)
    field_covariance = get_covariance(problem, covariance)
    mean_variance = tr(field_covariance) / grid_count
    return (
        sensor_count=sensor_count,
        sensing_intensity=sensor_count / (measurement_std^2 * step),
        mean_state_clarity=mean(get_clarity(problem, covariance)),
        mean_state_variance=tr(covariance) / grid_count,
        max_state_variance=eigmax(covariance),
        mean_field_clarity=inv(1 + mean_variance),
        mean_field_variance=mean_variance,
        max_field_variance=eigmax(field_covariance),
    )
end

function minimum_sensor_count(problem::STGPKFProblemContinuous, target::Real,
                              measurement_std::Real, step::Real; maximum::Int=200)
    0 < target < 1 || throw(ArgumentError("Target clarity must be in (0, 1)."))
    low, high = 1, maximum
    answer = nothing
    while low <= high
        middle = (low + high) ÷ 2
        metrics = clarity_metrics(problem, middle, measurement_std, step)
        if metrics.mean_field_clarity >= target
            answer = metrics
            high = middle - 1
        else
            low = middle + 1
        end
    end
    isnothing(answer) && error("Target clarity $target is infeasible with at most $maximum sensors.")
    return answer
end

function minimum_sensor_count(cache::AnalyticClarityCache, target::Real,
                              measurement_std::Real, step::Real; maximum::Int=200)
    0 < target < 1 || throw(ArgumentError("Target clarity must be in (0, 1)."))
    low, high = 1, maximum
    answer = nothing
    while low <= high
        middle = (low + high) ÷ 2
        metrics = analytic_clarity_metrics(
            cache, middle, measurement_std, step,
        )
        if metrics.mean_field_clarity >= target
            answer = metrics
            high = middle - 1
        else
            low = middle + 1
        end
    end
    isnothing(answer) && error(
        "Target clarity $target is infeasible with at most $maximum sensors.",
    )
    return answer
end
