#
#  CovalentForces.jl
#  collagen-rods
#
#  Created by Christopher Revell on 03/12/2020.
#
#
# Function to calculate covalent forces between given rods at N-C terminal ends

module CovalentForces

using LinearAlgebra
using LennardJones

@inline @views function covalentForces!(N,r,Ω,E,F,τ,rᵢⱼ,x,y,L,dummyVectors,tID,DParallel,DPerpendicular,DRotation,kT,ϵ,σ,Δu,covalentThresh,Q)

    # Repeat forces in both directions within pair
    for (i,j) in [(x,y),(y,x)]
        # N end to C adjacent internal point
        rᵢⱼ[:,tID] .= (r[j,:] .- Δu.*Ω[j,:]) .- (r[i,:] .+ 0.5*L.*Ω[i,:])
        rMag = sqrt(rᵢⱼ[:,tID]⋅rᵢⱼ[:,tID])
        if σ < rMag < covalentThresh
            FMag  = 10.0*Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ[:,tID] .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[i,1,:]).*E[i,1,:] .+ (rᵢⱼ[:,tID]⋅E[i,2,:]).*E[i,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[i,:]).*Ω[i,:]
            F[j,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[j,1,:]).*E[j,1,:] .+ (rᵢⱼ[:,tID]⋅E[j,2,:]).*E[j,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,tID] .+= (DRotation/kT).*((Ω[i,:].*L/2.0)×rᵢⱼ[:,tID])×Ω[i,:]
            τ[j,:,tID] .-= (DRotation/kT).*((-Ω[j,:].*Δu)×rᵢⱼ[:,tID])×Ω[j,:]
        end
        # C end to N adjacent internal point
        rᵢⱼ[:,tID] .= (r[j,:] .+ Δu.*Ω[j,:]) .- (r[i,:] .- 0.5*L.*Ω[i,:])
        rMag = sqrt(rᵢⱼ[:,tID]⋅rᵢⱼ[:,tID])
        if σ < rMag < covalentThresh
            FMag  = 10.0*Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
            rᵢⱼ[:,tID] .*= FMag/rMag
            # Linear forces acting on rods i and j using orthonormal basis vectors
            F[i,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[i,1,:]).*E[i,1,:] .+ (rᵢⱼ[:,tID]⋅E[i,2,:]).*E[i,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[i,:]).*Ω[i,:]
            F[j,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[j,1,:]).*E[j,1,:] .+ (rᵢⱼ[:,tID]⋅E[j,2,:]).*E[j,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[j,:]).*Ω[j,:]
            # Moment on rod i
            τ[i,:,tID] .+= (DRotation/kT).*((-Ω[i,:].*L/2.0)×rᵢⱼ[:,tID])×Ω[i,:]
            τ[j,:,tID] .-= (DRotation/kT).*((Ω[j,:].*Δu)×rᵢⱼ[:,tID])×Ω[j,:]
        end
    end

    return nothing
end

export covalentForces!

end
