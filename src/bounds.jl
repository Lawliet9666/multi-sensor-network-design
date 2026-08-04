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
    1 <= sensor_count <= grid_count || throw(ArgumentError(
        "Sensor count must be between 1 and $grid_count.",
    ))
    configuration_count > 0 || throw(ArgumentError("Configuration count must be positive."))
    return [
        [points[index] for index in sample(rng, 1:grid_count, sensor_count; replace=false)]
        for _ in 1:configuration_count
    ]
end

function discrete_bound_step(problem, covariance, sensor_models, rng;
                             sample_count::Int, beta::Real)
    grid_count = length(problem.pts)
    transition = I(grid_count) ⊗ problem.ss_model.Φ
    process_noise = I(grid_count) ⊗ problem.ss_model.W
    predicted = transition * covariance * transition'
    cross_left = transition * covariance
    cross_right = covariance * transition'

    count = min(sample_count, length(sensor_models))
    selected = sample(rng, 1:length(sensor_models), count; replace=false)
    correction = zeros(eltype(covariance), size(covariance))
    for index in selected
        H, V = sensor_models[index]
        correction .+= cross_left * H' * ((V + H * covariance * H') \ (H * cross_right))
    end
    correction ./= count
    next_covariance = process_noise + predicted - beta .* correction
    next_covariance = 0.5 .* (next_covariance + next_covariance')
    next_covariance += 1e-12I
    return Matrix(next_covariance)
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

function simulate_expected_covariance(problem::STGPKFProblem, data, point_sets,
                                      measurement_std::Real, measurement_covariance;
                                      trials::Int, seed::Int)
    trials > 0 || throw(ArgumentError("Trial count must be positive."))
    trajectories = Vector{Vector{Matrix{Float64}}}(undef, trials)
    for trial in 1:trials
        rng = MersenneTwister(seed + 1000 * trial)
        state = stgpkf_initialize(problem)
        trajectory = Matrix{Float64}[]
        for time in data.ts
            points = point_sets[rand(rng, eachindex(point_sets))]
            observations = [
                measure(rng, data, point[1], point[2], time, measurement_std)
                for point in points
            ]
            corrected = stgpkf_correct(problem, state, points, observations, measurement_covariance)
            state = stgpkf_predict(problem, corrected)
            push!(trajectory, Matrix(get_Σ(state)))
        end
        trajectories[trial] = trajectory
    end
    return expected_covariance(trajectories)
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
