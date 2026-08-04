function STGPKFProblem(points, spatial_kernel, temporal_kernel, step)
    isempty(points) && throw(ArgumentError("Grid points must be non-empty."))
    model = state_space_model(temporal_kernel, step)
    kernel = kernel_matrix(spatial_kernel, points)
    square_root = Symmetric(sqrt(kernel))
    inverse_square_root = Symmetric(inv(square_root))
    return STGPKFProblem(
        points, spatial_kernel, temporal_kernel, Float64(step), model,
        square_root, inverse_square_root,
    )
end

function STGPKFProblemContinuous(points, spatial_kernel, temporal_kernel)
    isempty(points) && throw(ArgumentError("Grid points must be non-empty."))
    model = state_space_model(temporal_kernel)
    kernel = kernel_matrix(spatial_kernel, points)
    square_root = Symmetric(sqrt(kernel))
    inverse_square_root = Symmetric(inv(square_root))
    return STGPKFProblemContinuous(
        points, spatial_kernel, temporal_kernel, model,
        square_root, inverse_square_root,
    )
end

function check_dimensions(problem::STGPKFProblem, state::KFState)
    expected = length(problem.pts) * dims(problem.ss_model)
    length(state) == expected || throw(DimensionMismatch(
        "Filter state has length $(length(state)); expected $expected.",
    ))
    return nothing
end

function get_estimate(problem::STGPKFProblem, state::KFState)
    n = length(problem.pts)
    return problem.sqrt_K_gg * ((I(n) ⊗ problem.ss_model.C) * get_μ(state))
end

function field_standard_deviation(problem, covariance_root::AbstractMatrix)
    n = length(problem.pts)
    observation = problem.sqrt_K_gg * (I(n) ⊗ problem.ss_model.C)
    transformed = covariance_root * observation'
    return sqrt.(vec(sum(abs2, transformed; dims=1)))
end

function get_estimate_clarity(problem::STGPKFProblem, state::KFState)
    σ = field_standard_deviation(problem, state.U)
    return @. inv(1 + σ^2)
end

function get_clarity(problem::Union{STGPKFProblem, STGPKFProblemContinuous}, covariance)
    root = cholesky(Symmetric(covariance)).U
    σ = field_standard_deviation(problem, root)
    return @. inv(1 + σ^2)
end

function get_covariance(problem::Union{STGPKFProblem, STGPKFProblemContinuous}, covariance)
    root = cholesky(Symmetric(covariance)).U
    n = length(problem.pts)
    transformed = problem.sqrt_K_gg * ((I(n) ⊗ problem.ss_model.C) * root')
    return Symmetric(transformed * transformed')
end

function stgpkf_initialize(problem::STGPKFProblem)
    n = length(problem.pts)
    state_dimension = dims(problem.ss_model)
    mean = zeros(n * state_dimension)
    covariance = (2.0 * I(n)) ⊗ initial_covariance(problem.kt)
    return KFState(; μ=mean, Σ=covariance)
end

function stgpkf_predict(problem::STGPKFProblem, state::KFState)
    check_dimensions(problem, state)
    n = length(problem.pts)
    transition = I(n) ⊗ problem.ss_model.Φ
    process_noise = I(n) ⊗ problem.ss_model.W
    return predict(state, transition, process_noise)
end

function sensor_matrices(problem::STGPKFProblem, points, measurement_covariance)
    Kmm = kernel_matrix(problem.ks, points)
    Kmg = kernel_matrix(problem.ks, points, problem.pts)
    interpolation = Kmg * problem.inv_sqrt_K_gg
    H = interpolation * (I(length(problem.pts)) ⊗ problem.ss_model.C)
    V = Symmetric(measurement_covariance) + Symmetric(Kmm) -
        Symmetric(interpolation * interpolation')
    return Matrix(H), Matrix(V)
end

function stgpkf_correct(problem::STGPKFProblem, state::KFState,
                        points::AbstractVector, observations::AbstractVector,
                        measurement_covariance::AbstractMatrix)
    check_dimensions(problem, state)
    length(points) == length(observations) || throw(DimensionMismatch(
        "Measurement points and observations must have the same length.",
    ))
    size(measurement_covariance) == (length(points), length(points)) ||
        throw(DimensionMismatch("Measurement covariance has the wrong size."))

    H, V = sensor_matrices(problem, points, measurement_covariance)
    return correct(state, observations, H, V)
end
