#
#  AdaptTimestep.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 08/10/2020.
#
#

module AdaptTimestep

using LinearAlgebra
#using Threads

@inline @views function adaptTimestep!(N,F,τ,ξr,ξΩ,σ,kT,L)

    Fmax_sq = maximum(sum(F[:,:,1].*F[:,:,1],dims=2))
    τmax_sq = maximum(sum(τ[:,:,1].*τ[:,:,1],dims=2))
    ξrmax_sq = maximum(ξr.⋅ξr)
    ξΩmax_sq = maximum(ξΩ.⋅ξΩ)

    Δt = min(σ^2/(4.0*ξrmax_sq),((σ*L)^2)/(4.0*ξΩmax_sq),σ/(2.0*sqrt(Fmax_sq)),σ*L/(2.0*sqrt(τmax_sq)))

    return Δt

end

export adaptTimestep!

end
