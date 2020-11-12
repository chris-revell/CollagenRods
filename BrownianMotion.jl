#
#  BrownianMotion.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module BrownianMotion

using LinearAlgebra
using Random
using Distributions
using Base.Threads

# Update rod positions and orientations according to established Brownian rod theory (Löwen Phys Rev E 1994)
@inline function brownianMotion!(N,Ω,ξr,ξΩ,E,DParallel,DPerpendicular,DRotation,RNG)
    @threads for i=1:N
        # Translational component of brownian motion
        ξr[i,:] .= DParallel*rand(RNG[threadid()],Normal(0.0,1.0)).*view(Ω,i,:) .+ DPerpendicular*rand(RNG[threadid()],Normal(0.0,1.0)).*view(E,i,:,1) .+ DPerpendicular*rand(RNG[threadid()],Normal(0.0,1.0)).*view(E,i,:,2)
        # Rotational component of brownian motion
        ξΩ[i,:] .= DRotation*rand(RNG[threadid()],Normal(0.0,1.0)).*view(E,i,:,1) .+ DRotation*rand(RNG[threadid()],Normal(0.0,1.0)).*view(E,i,:,2)
    end
    return nothing
end
export brownianMotion!

end
