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
include("./LennardJones.jl")
using .LennardJones
using Base.Threads

@inline @views function interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,dummyVectors)

    #Δu = (1.8/4.4)*L

    @threads for (x,y) in pairsList

        tID = threadid()

        # ---- RepulsiveForces ----

        # Find shortest distance between rod x and rod y
        (μ,λ) = shortestRodToRod!(r,Ω,rᵢⱼ[:,tID],x,y,L,dummyVectors[tID,:])
        rᵢⱼ[:,tID] .= (μ.*Ω[y] .- λ.*Ω[x])
        rMag = sqrt(rᵢⱼ[:,tID]⋅rᵢⱼ[:,tID])

        # If shortest distance is less than rod radius, calculate hard core repulsion
        if rMag < σ
            FMag  = lennardJones(rMag,ϵ,σ)
            rᵢⱼ[:,tID] .*= FMag/rMag
            # Linear forces acting on rods x and y using orthonormal basis vectors
            F[x,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[x][:,1]).*E[x][:,1] .+ (rᵢⱼ[:,tID]⋅E[x][:,2]).*E[x][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[x]).*Ω[x]
            F[y,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[y][:,1]).*E[y][:,1] .+ (rᵢⱼ[:,tID]⋅E[y][:,2]).*E[y][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[y]).*Ω[y]
            # Moments on rods x and y
            τ[x,:,tID] .+= (DRotation/kT).*((λ.*Ω[x])×rᵢⱼ[:,tID])×Ω[x]
            τ[y,:,tID] .-= (DRotation/kT).*((μ.*Ω[y])×rᵢⱼ[:,tID])×Ω[y]
        end


        # ---- Electrostatic ----
        colouredWidth = L/10
        blackWidth    = L/(2*8)

        electrostaticPairs = [(1,3),(2,4),(3,5),(4,6),(5,7),(6,8),(7,9),(9,1)]

        for (i,j) in electrostaticPairs
            μ = (i-5)*colouredWidth/2.0 + (i-5)*blackWidth/2.0
            λ = (j-5)*colouredWidth/2.0 + (j-5)*blackWidth/2.0

            rᵢⱼ[:,tID] .= r[y] .- r[x] .+ μ.*Ω[y] .- λ.*Ω[x]
            rMag = sqrt(rᵢⱼ[:,tID]⋅rᵢⱼ[:,tID])

            if rMag > σ
                FMag  = Q/(4*π*rMag^2)
                rᵢⱼ[:,tID] .*= FMag/rMag
                # Linear forces acting on rods x and y using orthonormal basis vectors
                F[x,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[x][:,1]).*E[x][:,1] .+ (rᵢⱼ[:,tID]⋅E[x][:,2]).*E[x][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[x]).*Ω[x]
                F[y,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[y][:,1]).*E[y][:,1] .+ (rᵢⱼ[:,tID]⋅E[y][:,2]).*E[y][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[y]).*Ω[y]
                # Moments on rods x and y
                τ[x,:,tID] .+= (DRotation/kT).*((λ.*Ω[x])×rᵢⱼ[:,tID])×Ω[x]
                τ[y,:,tID] .-= (DRotation/kT).*((μ.*Ω[y])×rᵢⱼ[:,tID])×Ω[y]
            end
        end


        # ---- Covalent Forces ----

        # # Repeat forces in both directions within pair
        # for (i,j) in [(x,y),(y,x)]
        #     # N end to C adjacent internal point
        #     rᵢⱼ[:,tID] .= (r[j,:] .- Δu.*Ω[j]) .- (r[i,:] .+ 0.5*L.*Ω[i])
        #     rMag = sqrt(rᵢⱼ[:,tID]⋅rᵢⱼ[:,tID])
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         rᵢⱼ[:,tID] .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[i][:,1]).*E[i][:,1] .+ (rᵢⱼ[:,tID]⋅E[i][:,2]).*E[i][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[i]).*Ω[i]
        #         F[j,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[j][:,1]).*E[j][:,1] .+ (rᵢⱼ[:,tID]⋅E[j][:,2]).*E[j][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[j]).*Ω[j]
        #         # Moment on rod i
        #         τ[i,:,tID] .+= (DRotation/kT).*((Ω[i].*L/2.0)×rᵢⱼ[:,tID])×Ω[i]
        #         τ[j,:,tID] .-= (DRotation/kT).*((-Ω[j].*Δu)×rᵢⱼ[:,tID])×Ω[j]
        #     end
        #     # C end to N adjacent internal point
        #     rᵢⱼ[:,tID] .= (r[j,:] .+ Δu.*Ω[j]) .- (r[i,:] .- 0.5*L.*Ω[i])
        #     rMag = sqrt(rᵢⱼ[:,tID]⋅rᵢⱼ[:,tID])
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         rᵢⱼ[:,tID] .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,tID] .+= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[i][:,1]).*E[i][:,1] .+ (rᵢⱼ[:,tID]⋅E[i][:,2]).*E[i][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[i]).*Ω[i]
        #         F[j,:,tID] .-= (DPerpendicular/kT)*((rᵢⱼ[:,tID]⋅E[j][:,1]).*E[j][:,1] .+ (rᵢⱼ[:,tID]⋅E[j][:,2]).*E[j][:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,tID]⋅Ω[j]).*Ω[j]
        #         # Moment on rod i
        #         τ[i,:,tID] .+= (DRotation/kT).*((-Ω[i].*L/2.0)×rᵢⱼ[:,tID])×Ω[i]
        #         τ[j,:,tID] .-= (DRotation/kT).*((Ω[j].*Δu)×rᵢⱼ[:,tID])×Ω[j]
        #     end
        # end
    end

    # Sum forces calculated in each thread
    F[:,:,1] = sum(F,dims=3)
    τ[:,:,1] = sum(τ,dims=3)
    #println(τ[:,:,:])

    return nothing
end

export interRodForces!

end
