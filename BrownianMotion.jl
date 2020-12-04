#
#  BrownianMotion.jl
#  collagen-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#
# Function to calculate the translational and rotational stochastic term for each particle

module BrownianMotion

using LinearAlgebra
using Random
using Distributions
using Base.Threads

# Update rod positions and orientations according to established Brownian rod theory (Löwen Phys Rev E 1994)
@inline @views function brownianMotion!(N,Ω,ξr,ξΩ,E,DParallel,DPerpendicular,DRotation,RNG,allRands)
    for i=1:N
        tID = threadid()
        allRands = rand(RNG[tID],Normal(0.0,1.0),5)
        # Translational component of brownian motion
        ξr[i,:] .= sqrt(2.0*DParallel)*allRands[1].*Ω[i,:] .+ sqrt(2.0*DPerpendicular).*(allRands[2].*E[i,1,:] .+ allRands[3].*E[i,2,:])
        # Rotational component of brownian motion
        ξΩ[i,:] .= sqrt(2.0*DRotation).*(allRands[4].*E[i,1,:] .+ allRands[5].*E[i,2,:])
    end
    return nothing
end
export brownianMotion!

end
