#
#  BoundaryForces.jl
#  collagen-rods
#
#  Created by Christopher Revell on 11/12/2020.
#
#
# Function to calculate boundary forces both tips of each rod

module BoundaryForces

using LinearAlgebra
using LennardJones

@inline @views function boundaryForces!(N,L,r,Ω,F,τ,E,containerRadius,σ,ϵ,DParallel,DPerpendicular,DRotation,kT)

    for i=1:N

        # + tip interacting with cylindrical boundary
        dR = r[i,:] .+ Ω[i,:].*L/2.0
        rMag = norm(dR[1:2])
        overlap = containerRadius-rMag
        if overlap < σ/2.0
            FMag  = lennardJones(overlap,ϵ,σ/2.0)
            dR[1:2] .*= FMag/rMag
            dR[3] = 0.0
            F[i,:,1] .+= (DPerpendicular/kT)*((dR⋅E[i,1,:]).*E[i,1,:] .+ (dR⋅E[i,2,:]).*E[i,2,:]) .+ (DParallel/kT)*(dR⋅Ω[i,:]).*Ω[i,:]
            τ[i,:,1] .+= (DRotation/kT).*((Ω[i,:].*L/2.0)×dR)×Ω[i,:]
        end

        # - tip interacting with cylindrical boundary
        dR .= r[i,:] .- Ω[i,:].*L/2.0
        rMag = norm(dR[1:2])
        overlap = containerRadius-rMag
        if overlap < σ/2.0
            FMag  = lennardJones(overlap,ϵ,σ/2.0)
            dR[1:2] .*= FMag/rMag
            dR[3] = 0.0
            F[i,:,1] .+= (DPerpendicular/kT)*((dR⋅E[i,1,:]).*E[i,1,:] .+ (dR⋅E[i,2,:]).*E[i,2,:]) .+ (DParallel/kT)*(dR⋅Ω[i,:]).*Ω[i,:]
            τ[i,:,1] .+= (DRotation/kT).*((-Ω[i,:].*L/2.0)×dR)×Ω[i,:]
        end

    end
    return nothing
end

export boundaryForces!

end
