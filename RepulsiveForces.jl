#
#  RepulsiveForces.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module RepulsiveForces

using LinearAlgebra
include("./ShortestDistance.jl")
using .ShortestDistance
using Base.Threads

@inline function repulsiveForces!(pairsList,N,r,Ω,F,τ,E,ϵ,σ,DParallel,DPerpendicular,DRotation,kT)

    @threads for (i,j) in pairsList

        rᵢⱼ = zeros(Float64,3)  # Vector between rod centrepoints

        # Find shortest distance between rod i and rod j
        (μ,λ) = shortestDistance!(r,Ω,rᵢⱼ,i,j)
        minSeparation = sqrt(rᵢⱼ⋅rᵢⱼ)

        # If shortest distance is less than rod radius, calculate hard core repulsion
        if minSeparation < σ
            FMag  = 24.0*ϵ*((σ^6.0)/(minSeparation^7.0) - 2.0*(σ^12.0)/(minSeparation^13.0))
            rᵢⱼ .*= FMag/minSeparation

            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
            F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]

            # Moments on rods i and j
            τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].-λ.*Ω[i,:])×rᵢⱼ)×Ω[i,:]
            τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].-μ.*Ω[j,:])×rᵢⱼ)×Ω[j,:]
        end
    end
    return nothing
end

export repulsiveForces!

end
