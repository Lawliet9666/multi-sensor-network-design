abstract type AbstractKernel end
abstract type AbstractMaternKernel <: AbstractKernel end

"Matérn-1/2 covariance kernel parameterized by variance and inverse length scale."
struct Matern12{F} <: AbstractMaternKernel
    σsq::F
    λ::F
end

abstract type AbstractStateSpaceModel end

struct ContinuousTimeStateSpaceModel{MA, MB, MC} <: AbstractStateSpaceModel
    A::MA
    B::MB
    C::MC
end

struct DiscreteTimeStateSpaceModel{MΦ, MW, MC, F} <: AbstractStateSpaceModel
    Φ::MΦ
    W::MW
    C::MC
    dt::F
end

struct STGPKFProblem{P, F, VP <: AbstractVector{P}, KS <: AbstractKernel,
                     KT <: AbstractKernel, SS <: DiscreteTimeStateSpaceModel,
                     M1 <: AbstractMatrix{F}, M2 <: AbstractMatrix{F}}
    pts::VP
    ks::KS
    kt::KT
    ΔT::F
    ss_model::SS
    sqrt_K_gg::M1
    inv_sqrt_K_gg::M2
end

struct STGPKFProblemContinuous{P, F, VP <: AbstractVector{P}, KS <: AbstractKernel,
                               KT <: AbstractKernel, SS <: ContinuousTimeStateSpaceModel,
                               M1 <: AbstractMatrix{F}, M2 <: AbstractMatrix{F}}
    pts::VP
    ks::KS
    kt::KT
    ss_model::SS
    sqrt_K_gg::M1
    inv_sqrt_K_gg::M2
end

struct SpatiotemporalData2D{F, VF <: AbstractVector{F}, AF <: AbstractArray{F}, TI}
    xs::VF
    ys::VF
    ts::VF
    data::AF
    itp::TI
end
