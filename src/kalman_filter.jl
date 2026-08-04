"Square-root Kalman-filter state with covariance `U' * U`."
struct KFState{T, V <: AbstractVector{T}, M}
    μ::V
    U::UpperTriangular{T, M}
end

function KFState(; μ, Σ)
    size(Σ) == (length(μ), length(μ)) || throw(DimensionMismatch(
        "Mean length $(length(μ)) does not match covariance size $(size(Σ)).",
    ))
    return KFState(μ, chol_sqrt(Symmetric(Σ)))
end

get_μ(state::KFState) = state.μ
get_Σ(state::KFState) = Symmetric(state.U' * state.U)
Base.length(state::KFState) = length(state.μ)

function predict(state::KFState, A, W)
    μ = A * state.μ
    U = qrr(state.U * A', chol_sqrt(W))
    return KFState(μ, U)
end

function correct(state::KFState, observation, H, V)
    Γv = chol_sqrt(V)
    innovation = observation - H * state.μ
    gain = kalman_gain(state, H, Γv)
    μ = state.μ + gain * innovation
    U = qrr(state.U - (state.U * H') * gain', Γv * gain')
    return KFState(μ, U)
end

chol_sqrt(matrix) = cholesky(matrix).U

function qrr(a::AbstractMatrix, b::AbstractMatrix)
    stacked = Matrix([a; b])
    LinearAlgebra.LAPACK.geqrf!(stacked)
    n = minimum(size(stacked))
    return UpperTriangular(stacked[1:n, 1:n])
end

function kalman_gain(state::KFState, H, Γv)
    innovation_root = qrr(state.U * H', Γv)
    return ((state.U' * state.U * H') / innovation_root) / innovation_root'
end
