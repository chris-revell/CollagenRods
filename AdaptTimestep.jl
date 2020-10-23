#
#  AdaptTimestep.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 08/10/2020.
#
#

module AdaptTimestep

using LinearAlgebra
#using .Threads

@inline function adaptTimestep!(N,F,τ,σ,D,kT)

    Fmax_sq = maximum(sum(F[:,:,1].*F[:,:,1],dims=2))
    τmax_sq = maximum(sum(τ[:,:,1].*τ[:,:,1],dims=2))

    dt = min(σ^2/(32*D),σ/sqrt(Fmax_sq),σ/sqrt(τmax_sq))
    #dt = min(0.00001,kT*σ/(2.0*D*sqrt(Fmax_sq)))

    return dt

end

export adaptTimestep!

end
