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

@inline function repulsiveForces!(N,r,Ω,F,τ,E,ϵ,σ,DParallel,DPerpendicular,DRotation,kT)

    @threads for i=1:N
        rᵢⱼ = zeros(Float64,3)
        for j=i+1:N

            # Check for parallel rods
            #if Ω[i,:].-Ω[j,:]

            # Find shortest distance between rod i and rod j

            (μ,λ) = shortestDistance!(r,Ω,rᵢⱼ,i,j)
            minSeparation = sqrt(rᵢⱼ⋅rᵢⱼ)

            # If shortest distance is less than rod radius, calculate hard core repulsion
            if minSeparation < σ/2.0
                FMag  = 24.0*ϵ*((σ^6.0)/(minSeparation^7.0) - 2.0*(σ^12.0)/(minSeparation^13.0))
                rᵢⱼ .*= FMag/minSeparation

                # Linear forces acting on rods i and j using orthonormal basis vectors
                F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
                F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]

                # Moment on rod i
                τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].-λ.*Ω[i,:])×rᵢⱼ)×Ω[i,:]
                # Moment on rod j
                τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].-μ.*Ω[j,:])×rᵢⱼ)×Ω[j,:]
            end
        end
    end
end

export repulsiveForces!

end
