"Construct the Matérn-1/2 kernel used by every paper experiment."
function Matern(order, σ, length_scale)
    order == 1 / 2 || throw(ArgumentError("This paper repository supports only Matérn-1/2."))
    σ > 0 || throw(ArgumentError("Kernel standard deviation must be positive."))
    length_scale > 0 || throw(ArgumentError("Kernel length scale must be positive."))
    return Matern12(σ^2, inv(length_scale))
end

dims(model::ContinuousTimeStateSpaceModel) = size(model.A, 1)
dims(model::DiscreteTimeStateSpaceModel) = size(model.Φ, 1)

function (kernel::Matern12)(x, y)
    return kernel.σsq * exp(-kernel.λ * norm(x - y))
end

function kernel_matrix(kernel::AbstractKernel, points::AbstractVector)
    n = length(points)
    matrix = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in i:n
        matrix[i, j] = kernel(points[i], points[j])
        matrix[j, i] = matrix[i, j]
    end
    return Symmetric(matrix)
end

function kernel_matrix(kernel::AbstractKernel, x::AbstractVector, y::AbstractVector)
    matrix = Matrix{Float64}(undef, length(x), length(y))
    for i in eachindex(x), j in eachindex(y)
        matrix[i, j] = kernel(x[i], y[j])
    end
    return matrix
end

function state_space_model(kernel::Matern12{F}) where {F}
    λ = kernel.λ
    σ = sqrt(kernel.σsq)
    return ContinuousTimeStateSpaceModel(
        @SMatrix([-λ;;]),
        @SMatrix([one(F);;]),
        @SMatrix([σ * sqrt(2λ);;]),
    )
end

function state_space_model(kernel::Matern12{F}, step::Real) where {F}
    step > 0 || throw(ArgumentError("Sampling interval must be positive."))
    λ = kernel.λ
    σ = sqrt(kernel.σsq)
    Φ = @SMatrix([exp(-step * λ);;])
    W = @SMatrix([(1 - exp(-2 * step * λ)) / (2 * λ);;])
    C = @SMatrix([σ * sqrt(2λ);;])
    return DiscreteTimeStateSpaceModel(Φ, W, C, F(step))
end

initial_covariance(kernel::Matern12) = @SMatrix([inv(2 * kernel.λ);;])
