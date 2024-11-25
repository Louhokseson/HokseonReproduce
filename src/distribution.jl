"""
Lorentzian/ Cauchy distribution function
"""
function Lorentzian(ω::Real, ω0::Real, Γ::Real) 
    return 1/π * Γ / ((ω - ω0)^2 + Γ^2)
end

"""
We need to use a lot of times the Fermi-Dirac distribution function
"""
function FermiDirac(ϵ::Real, ϵf::Real, T::Real)
    return 1 / (exp((ϵ - ϵf) / T) + 1)
end
