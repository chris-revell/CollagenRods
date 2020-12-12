#
#  Initialise.jl
#  collagen-rods
#
#  Created by Christopher Revell on 03/12/2020.
#
#
# Function to calculate repulsive forces between given rods

module Initialise

using LinearAlgebra
using Random
using Base.Threads
using OutputData
using CreateRunDirectory

@inline @views function initialise!(L,σ,ϵ,kT,Q,η,N,tMax,containerRadius,containerVolume,interactionThresh,r,Ω,outputToggle)

    # Diffusion constants from Löwen Phys Rev E 1994
    p                 = L/σ         # Rod aspect ratio
    D₀                = kT/(η*L)
    DParallel         = D₀*(log(p)+-0.207+0.980/p-0.133/p^2)/2π
    DPerpendicular    = D₀*(log(p)+0.839+0.185/p+0.233/p^2)/4π
    DRotation         = 3*D₀*(log(p)-0.662+0.917/p-0.05/p^2)/(π*L^2)

    containerHeight = containerVolume/(π*containerRadius^2)

    for i=1:N
        foundPosition = 0
        while foundPosition==0
            r[i,:] = (rand(Float64,3).-0.5).*containerHeight
            Ω[i,:] = normalize(rand(Float64,3))
            if norm((r[i,:].+Ω[i,:].*L/2.0)[1:2]) < containerRadius && norm((r[i,:].-Ω[i,:].*L/2.0)[1:2]) < containerRadius
                foundPosition=1
            else
                nothing
            end
        end
    end

    if outputToggle==1
        foldername = createRunDirectory(N,L,σ,ϵ,Q,p,η,kT,tMax,containerRadius,containerVolume,containerHeight,D₀,DParallel,DPerpendicular,DRotation,interactionThresh)
        outfile = open("output/$(foldername)/output.txt","w")
        outputData(r,Ω,outfile,0,tMax)
    end

    # Create random number generators for each thread
    threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())
    for i in 1:nthreads()
        threadRNG[i] = Random.MersenneTwister()
    end

    return D₀, DParallel, DPerpendicular, DRotation, foldername, outfile, threadRNG
end

export initialise!

end
