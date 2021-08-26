#
#  RepulsiveForces.jl
#  CollagenRods
#
#  Created by Christopher Revell on 03/12/2020.
#
#
# Function to calculate repulsive forces between given rods

module RepulsiveForces

using LinearAlgebra
using ShortestDistance
using LennardJones

@inline @views function repulsiveForces!(N,r,Ω,E,F,τ,rᵢⱼ,x,y,L,dummyVectors,tID,DParallel,DPerpendicular,DRotation,kT,ϵ,σ)


    # Find shortest distance between rod x and rod y
    (μ,λ) = shortestRodToRod!(r,Ω,rᵢⱼ[:,tID],x,y,L,dummyVectors[:,:,tID])
    rMag = sqrt(rᵢⱼ[:,tID]⋅rᵢⱼ[:,tID])

    # If shortest distance is less than rod radius, calculate hard core repulsion
    if rMag < σ
        FMag  = lennardJones(rMag,ϵ,σ)
        rᵢⱼ[:,tID] .*= FMag/rMag
        # Linear forces acting on rods x and y using orthonormal basis vectors
        F[x,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[x,1,:]).*E[x,1,:] .+ (rᵢⱼ[:,tID]⋅E[x,2,:]).*E[x,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[x,:]).*Ω[x,:]
        F[y,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[y,1,:]).*E[y,1,:] .+ (rᵢⱼ[:,tID]⋅E[y,2,:]).*E[y,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[y,:]).*Ω[y,:]
        # Moments on rods x and y
        τ[x,:,tID] .+= (DRotation/kT).*((λ.*Ω[x,:])×rᵢⱼ[:,tID])×Ω[x,:]
        τ[y,:,tID] .-= (DRotation/kT).*((μ.*Ω[y,:])×rᵢⱼ[:,tID])×Ω[y,:]
    end

    return nothing
end

export repulsiveForces!

end
