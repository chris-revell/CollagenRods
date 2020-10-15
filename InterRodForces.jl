#
#  InterRodForces.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 15/10/2020.
#
#

module InterRodForces

using LinearAlgebra
include("./ShortestDistance.jl")
using .ShortestDistance
using Base.Threads

@inline function interRodForces!(pairsList,N,r,Ω,F,τ,E,A,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ)

    k₀ = 1.0
    ϵ₀ = 1.0
    Q  = 20.0

    @threads for (x,y) in pairsList

        rᵢⱼ = zeros(Float64,3)

        # ---- RepulsiveForces ----

        # Find shortest distance between rod x and rod y
        (μ,λ) = shortestDistance!(r,Ω,rᵢⱼ,x,y)
        minSeparation = sqrt(rᵢⱼ⋅rᵢⱼ)

        # If shortest distance is less than rod radius, calculate hard core repulsion
        if minSeparation < σ
            FMag  = 24.0*ϵ*((σ^6.0)/(minSeparation^7.0) - 2.0*(σ^12.0)/(minSeparation^13.0))
            rᵢⱼ .*= FMag/minSeparation

            # Linear forces acting on rods x and y using orthonormal basis vectors
            F[x,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[x,:,1]).*E[x,:,1] .+ (rᵢⱼ⋅E[x,:,2]).*E[x,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[x,:]).*Ω[x,:]
            F[y,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[y,:,1]).*E[y,:,1] .+ (rᵢⱼ⋅E[y,:,2]).*E[y,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[y,:]).*Ω[y,:]

            # Moments on rods x and y
            τ[x,:,threadid()] .+= (DRotation/kT).*((r[x,:].-λ.*Ω[x,:])×rᵢⱼ)×Ω[x,:]
            τ[y,:,threadid()] .-= (DRotation/kT).*((r[y,:].-μ.*Ω[y,:])×rᵢⱼ)×Ω[y,:]
        end

        # ---- AttractiveForces ----

        # Repeat forces in both directions within pair
        for (i,j) in [(x,y),(y,x)]
            # N end to N adjacent internal point
            rᵢⱼ .= (r[j,:] .+ Ω[j,:].*L/2.0) .- (r[i,:] .+ Ω[i,:].*L/3.0)
            rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
            FMag  = (1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
            F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].+Ω[i,:].*L/3.0)×rᵢⱼ)×Ω[i,:]
            τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].+Ω[j,:].*L/2.0)×rᵢⱼ)×Ω[j,:]

            # N end to C adjacent internal point
            rᵢⱼ .= (r[j,:] .+ Ω[j,:].*L/2.0) .- (r[i,:] .- Ω[i,:].*L/3.0)
            rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
            FMag  = (1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
            F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].-Ω[i,:].*L/3.0)×rᵢⱼ)×Ω[i,:]
            τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].+Ω[j,:].*L/2.0)×rᵢⱼ)×Ω[j,:]

            # C end to C adjacent internal point
            rᵢⱼ .= (r[j,:] .- Ω[j,:].*L/2.0) .- (r[i,:] .- Ω[i,:].*L/3.0)
            rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
            FMag  = (1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
            F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].-Ω[i,:].*L/3.0)×rᵢⱼ)×Ω[i,:]
            τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].-Ω[j,:].*L/2.0)×rᵢⱼ)×Ω[j,:]

            # C end to N adjacent internal point
            rᵢⱼ .= (r[j,:] .- Ω[j,:].*L/2.0) .- (r[i,:] .+ Ω[i,:].*L/3.0)
            rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
            FMag  = (1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]
            F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,threadid()] .+= (DRotation/kT).*((r[i,:].+Ω[i,:].*L/3.0)×rᵢⱼ)×Ω[i,:]
            τ[j,:,threadid()] .-= (DRotation/kT).*((r[j,:].-Ω[j,:].*L/2.0)×rᵢⱼ)×Ω[j,:]
        end
    end
    return nothing
end

export interRodForces!

end
