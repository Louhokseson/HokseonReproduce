using DrWatson
@quickactivate "HokseonReproduce"

include(srcdir("distribution.jl"))

"""
    BrandbygeAborbateModel is given at https://doi.org/10.1103/PhysRevB.52.6042 by Brandbyge et al. 1995

    The DOS of the impurity follows the equation (22 a,b)
"""
struct BrandbygeAborbateModel
    Δ₀::AbstractFloat
    β::AbstractFloat
    ε∞::AbstractFloat
    C::AbstractFloat
    α::AbstractFloat
end

function BrandbygeAborbateModel(;
        Δ₀= 0.2,
        β= 1.0,
        ε∞= 5.0,
        C= 3.0,
        α = 0.5,
    )
    return BrandbygeAborbateModel(Δ₀, β, ε∞, C, α)
end


function DOS(r::Real ,Δϵ::Real, m::BrandbygeAborbateModel)
    """
    DOS of the impurity follows the equation (22 a,b)

    we assume the absorbate is a Cauchy/Lorentz distribution

    r: distance from the impurity to the reservoir

    Δϵ: energy range for plotting the absorbate DOS
    
    m: BrandbygeAborbateModel
    """
    Δ = m.Δ₀ * exp(-m.β * r)
    ϵₐ = m.ε∞ - m.C * exp(-m.α * r)

    lorentzian(ω) = Lorentzian(ω, ϵₐ, Δ)

    dos_r = collect(range(ϵₐ - Δϵ, ϵₐ + Δϵ, length = 100))

    dos = lorentzian.(dos_r)
    return dos_r, dos
end