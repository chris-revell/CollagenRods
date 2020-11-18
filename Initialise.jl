#
#  Initialise.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 17/11/2020.
#
#

module Initialise

using LinearAlgebra
using StaticArrays
using Random
using OutputData

@inline @views function initialise!(N,r,Ω,boxSize,threadRNG,nthreads,outputToggle)

    # Create random number generators for each thread
    for i in 1:nthreads
        threadRNG[i] = Random.MersenneTwister()
    end
    
    # Initialise random rod positions
    for i=1:N
        r[i] = SVector{3}((rand(threadRNG[1],Float64,3).-0.5).*boxSize)
    end
    # Initialise random rod orientations
    for i=1:N
        Ω[i] = SVector{3}(rand(threadRNG[1],Float64,3).-0.5)
        # Normalise magnitude
        Ω[i] = Ω[i]./sqrt(Ω[i]⋅Ω[i])
    end

    return nothing
end

export initialise!

end
