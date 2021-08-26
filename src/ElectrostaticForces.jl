#
#  ElectrostaticForces.jl
#  CollagenRods
#
#  Created by Christopher Revell on 03/12/2020.
#
#
# Function to calculate electrostatic forces between given rods

module ElectrostaticForces

using LinearAlgebra
using LennardJones

@inline @views function electrostaticForces!(N,r,Ω,E,F,τ,rᵢⱼ,x,y,dummyVectors,tID,DParallel,DPerpendicular,DRotation,kT,ϵ,σ,electrostaticThresh,electrostaticPairs,blackWidth,colouredWidth,Q)

    # for (i,j) in electrostaticPairs
    #     μ = (i-5)*colouredWidth/2.0 + (i-5)*blackWidth/2.0
    #     λ = (j-5)*colouredWidth/2.0 + (j-5)*blackWidth/2.0
    #
    #     rᵢⱼ[:,tID] .= r[y,:] .- r[x,:] .+ μ.*Ω[y,:] .- λ.*Ω[x,:]
    #     rMag = norm(rᵢⱼ[:,tID])
    #
    #     if rMag > σ
    #         FMag  = Q/(4*π*rMag^2)
    #         rᵢⱼ[:,tID] .*= FMag/rMag
    #         # Linear forces acting on rods x and y using orthonormal basis vectors
    #         F[x,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[x,1,:]).*E[x,1,:] .+ (rᵢⱼ[:,tID]⋅E[x,2,:]).*E[x,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[x,:]).*Ω[x,:]
    #         F[y,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[y,1,:]).*E[y,1,:] .+ (rᵢⱼ[:,tID]⋅E[y,2,:]).*E[y,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[y,:]).*Ω[y,:]
    #         # Moments on rods x and y
    #         τ[x,:,tID] .+= (DRotation/kT).*((λ.*Ω[x,:])×rᵢⱼ[:,tID])×Ω[x,:]
    #         τ[y,:,tID] .-= (DRotation/kT).*((μ.*Ω[y,:])×rᵢⱼ[:,tID])×Ω[y,:]
    #     end
    # end

    for i=1:9
        μ = (i-5)*colouredWidth/2.0 + (i-5)*blackWidth/2.0
        for j=1:9
            if i==j
                #skip
            else
                λ = (j-5)*colouredWidth/2.0 + (j-5)*blackWidth/2.0

                rᵢⱼ[:,tID] .= r[y,:] .- r[x,:] .+ μ.*Ω[y,:] .- λ.*Ω[x,:]
                rMag = norm(rᵢⱼ[:,tID])

                if σ < rMag < electrostaticThresh

                    iseven(i-j) ? FMag = Q/(4*π*rMag^2) : FMag = -Q/(4*π*rMag^2)

                    rᵢⱼ[:,tID] .*= FMag/rMag
                    # Linear forces acting on rods x and y using orthonormal basis vectors
                    F[x,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[x,1,:]).*E[x,1,:] .+ (rᵢⱼ[:,tID]⋅E[x,2,:]).*E[x,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[x,:]).*Ω[x,:]
                    F[y,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[y,1,:]).*E[y,1,:] .+ (rᵢⱼ[:,tID]⋅E[y,2,:]).*E[y,2,:]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[y,:]).*Ω[y,:]
                    # Moments on rods x and y
                    τ[x,:,tID] .+= (DRotation/kT).*((λ.*Ω[x,:])×rᵢⱼ[:,tID])×Ω[x,:]
                    τ[y,:,tID] .-= (DRotation/kT).*((μ.*Ω[y,:])×rᵢⱼ[:,tID])×Ω[y,:]
                end
            end
        end
    end

    return nothing
end

export electrostaticForces!

end
