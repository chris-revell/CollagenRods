#
#  AttractiveForces.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 08/10/2020.
#
#

module AttractiveForces

using LinearAlgebra
include("./ShortestDistance.jl")
using .ShortestDistance
using Base.Threads

@inline function attractiveForces!(N,r,Ω,F,τ,E,DParallel,DPerpendicular,DRotation,kT,L)

    k₀ = 1.0
    ϵ₀ = 1.0
    Q  = 1.0

    @threads for i=1:N
        rᵢⱼ = zeros(Float64,3)
        for j=i+1:N

            # N end
            rᵢⱼ .= (r[j,:] .+ Ω[j,:].*L/2.0) .- (r[i,:] .+ Ω[i,:].*L/2.0)
            rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
            FMag  = (1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
            F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].+Ω[i,:].*L/2.0)×rᵢⱼ)×Ω[i,:]
            τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].+Ω[j,:].*L/2.0)×rᵢⱼ)×Ω[j,:]

            # C end
            rᵢⱼ .= (r[i,:] .- Ω[i,:].*L/2.0) .- (r[j,:] .- Ω[j,:].*L/2.0)
            rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
            FMag  = (1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
            F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].-Ω[i,:].*L/2.0)×rᵢⱼ)×Ω[i,:]
            τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].-Ω[j,:].*L/2.0)×rᵢⱼ)×Ω[j,:]

        end
    end
end

export attractiveForces!

end
