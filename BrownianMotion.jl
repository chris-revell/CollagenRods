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

# Update rod positions and orientations according to established Brownian rod theory (Löwen Phys Rev E 1994)
@inline function brownianMotion!(N,Ω,ξr,ξΩ,E,stds,RNG)
    for i=1:N
        # Translational component of brownian motion
        ξr[i,:] .= rand(RNG,Normal(0.0,stds[1])).*Ω[i,:] .+ rand(RNG,Normal(0.0,stds[2])).*E[i,:,1] .+ rand(RNG,Normal(0.0,stds[2])).*E[i,:,2]
        # Rotational component of brownian motion
        ξΩ[i,:] .= rand(RNG,Normal(0.0,stds[3])).*E[i,:,1] .+ rand(RNG,Normal(0.0,stds[3])).*E[i,:,2]
    end
end
export brownianMotion!

end
