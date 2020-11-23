#
#  InterRodForces.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 15/10/2020.
#
#

module InterRodForces

using LinearAlgebra
using ShortestDistance
using LennardJones
using Base.Threads

@inline @views function interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,dummyVectors)

    #k₀ = 1.0
    #ϵ₀ = 1.0
    #Q  = 10.0

    @threads for (x,y) in pairsList

        Δu = (1.8/4.4)*L

        # ---- RepulsiveForces ----

        # Find shortest distance between rod x and rod y
        (μ,λ) = shortestRodToRod!(r,Ω,rᵢⱼ[:,threadid()],x,y,L,dummyVectors[:,:,threadid()])
        rMag = sqrt(rᵢⱼ[:,threadid()]⋅rᵢⱼ[:,threadid()])

        # If shortest distance is less than rod radius, calculate hard core repulsion
        if rMag < σ
            FMag  = lennardJones(rMag,ϵ,σ)
            rᵢⱼ[:,threadid()] .*= FMag/rMag
            # Linear forces acting on rods x and y using orthonormal basis vectors
            F[x,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[x,:,1]).*E[x,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[x,:,2]).*E[x,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[x,:]).*Ω[x,:]
            F[y,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[y,:,1]).*E[y,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[y,:,2]).*E[y,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[y,:]).*Ω[y,:]
            # Moments on rods x and y
            τ[x,:,threadid()] .+= (DRotation/kT).*((λ.*Ω[x,:])×rᵢⱼ[:,threadid()])×Ω[x,:]
            τ[y,:,threadid()] .-= (DRotation/kT).*((μ.*Ω[y,:])×rᵢⱼ[:,threadid()])×Ω[y,:]
        end


        # ---- Electrostatic ----
        colouredWidth = L/10
        blackWidth    = L/(2*8)

        electrostaticPairs = [(1,3),(2,4),(3,5),(4,6),(5,7),(6,8),(7,9),(9,1)]

        for (i,j) in electrostaticPairs
            μ = (i-5)*colouredWidth/2.0 + (i-5)*blackWidth/2.0
            λ = (j-5)*colouredWidth/2.0 + (j-5)*blackWidth/2.0

            rᵢⱼ[:,threadid()] .= r[y,:] .- r[x,:] .+ μ.*Ω[y,:] .- λ.*Ω[x,:]
            rMag = sqrt(rᵢⱼ[:,threadid()]⋅rᵢⱼ[:,threadid()])

            if rMag > σ
                FMag  = Q/(4*π*rMag^2)
                rᵢⱼ[:,threadid()] .*= FMag/rMag
                # Linear forces acting on rods x and y using orthonormal basis vectors
                F[x,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[x,:,1]).*E[x,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[x,:,2]).*E[x,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[x,:]).*Ω[x,:]
                F[y,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[y,:,1]).*E[y,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[y,:,2]).*E[y,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[y,:]).*Ω[y,:]
                # Moments on rods x and y
                τ[x,:,threadid()] .+= (DRotation/kT).*((λ.*Ω[x,:])×rᵢⱼ[:,threadid()])×Ω[x,:]
                τ[y,:,threadid()] .-= (DRotation/kT).*((μ.*Ω[y,:])×rᵢⱼ[:,threadid()])×Ω[y,:]
            end
        end


        # ---- Covalent Forces ----

        # # Repeat forces in both directions within pair
        # for (i,j) in [(x,y),(y,x)]
        #     # N end to C adjacent internal point
        #     rᵢⱼ[:,threadid()] .= (r[j,:] .- Δu.*Ω[j,:]) .- (r[i,:] .+ 0.5*L.*Ω[i,:])
        #     rMag = sqrt(rᵢⱼ[:,threadid()]⋅rᵢⱼ[:,threadid()])
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         rᵢⱼ[:,threadid()] .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[i,:]).*Ω[i,:]
        #         F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[j,:]).*Ω[j,:]
        #         # Moment on rod i
        #         τ[i,:,threadid()] .+= (DRotation/kT).*((Ω[i,:].*L/2.0)×rᵢⱼ[:,threadid()])×Ω[i,:]
        #         τ[j,:,threadid()] .-= (DRotation/kT).*((-Ω[j,:].*Δu)×rᵢⱼ[:,threadid()])×Ω[j,:]
        #     end
        #     # C end to N adjacent internal point
        #     rᵢⱼ[:,threadid()] .= (r[j,:] .+ Δu.*Ω[j,:]) .- (r[i,:] .- 0.5*L.*Ω[i,:])
        #     rMag = sqrt(rᵢⱼ[:,threadid()]⋅rᵢⱼ[:,threadid()])
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         rᵢⱼ[:,threadid()] .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[i,:,1]).*E[i,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[i,:,2]).*E[i,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[i,:]).*Ω[i,:]
        #         F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ[:,threadid()]⋅E[j,:,1]).*E[j,:,1] .+ (rᵢⱼ[:,threadid()]⋅E[j,:,2]).*E[j,:,2]) .+ (DParallel/kT)*(rᵢⱼ[:,threadid()]⋅Ω[j,:]).*Ω[j,:]
        #         # Moment on rod i
        #         τ[i,:,threadid()] .+= (DRotation/kT).*((-Ω[i,:].*L/2.0)×rᵢⱼ[:,threadid()])×Ω[i,:]
        #         τ[j,:,threadid()] .-= (DRotation/kT).*((Ω[j,:].*Δu)×rᵢⱼ[:,threadid()])×Ω[j,:]
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
