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

@inline function interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,dummyVectors)

    #k₀ = 1.0
    #ϵ₀ = 1.0
    #Q  = 10.0

    @threads for (x,y) in pairsList

        Δu = (1.8/4.4)*L

        # ---- RepulsiveForces ----

        # Find shortest distance between rod x and rod y
        (μ,λ) = shortestRodToRod!(r,Ω,view(rᵢⱼ,:,threadid()),x,y,L,dummyVectors[:,:,threadid()])
        rMag = sqrt(view(rᵢⱼ,:,threadid())⋅view(rᵢⱼ,:,threadid()))

        # If shortest distance is less than rod radius, calculate hard core repulsion
        if rMag < σ
            FMag  = lennardJones(rMag,ϵ,σ)
            rᵢⱼ[:,threadid()] .*= FMag/rMag
            # Linear forces acting on rods x and y using orthonormal basis vectors
            F[x,:,threadid()] .+= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,x,:,1)).*view(E,x,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,x,:,2)).*view(E,x,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,x,:)).*view(Ω,x,:)
            F[y,:,threadid()] .-= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,y,:,1)).*view(E,y,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,y,:,2)).*view(E,y,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,y,:)).*view(Ω,y,:)
            # Moments on rods x and y
            τ[x,:,threadid()] .+= (DRotation/kT).*((λ.*view(Ω,x,:))×view(rᵢⱼ,:,threadid()))×view(Ω,x,:)
            τ[y,:,threadid()] .-= (DRotation/kT).*((μ.*view(Ω,y,:))×view(rᵢⱼ,:,threadid()))×view(Ω,y,:)
        end


        # ---- Electrostatic ----
        colouredWidth = L/10
        blackWidth    = L/(2*8)

        electrostaticPairs = [(1,3),(2,4),(3,5),(4,6),(5,7),(6,8),(7,9),(9,1)]

        for (i,j) in electrostaticPairs
            μ = (i-5)*colouredWidth/2.0 + (i-5)*blackWidth/2.0
            λ = (j-5)*colouredWidth/2.0 + (j-5)*blackWidth/2.0

            rᵢⱼ[:,threadid()] .= view(r,y,:) .- view(r,x,:) .+ μ.*view(Ω,y,:) .- λ.*view(Ω,x,:)
            rMag = sqrt(view(rᵢⱼ,:,threadid())⋅view(rᵢⱼ,:,threadid()))

            if rMag > σ
                FMag  = Q/(4*π*rMag^2)
                rᵢⱼ[:,threadid()] .*= FMag/rMag
                # Linear forces acting on rods x and y using orthonormal basis vectors
                F[x,:,threadid()] .+= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,x,:,1)).*view(E,x,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,x,:,2)).*view(E,x,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,x,:)).*view(Ω,x,:)
                F[y,:,threadid()] .-= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,y,:,1)).*view(E,y,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,y,:,2)).*view(E,y,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,y,:)).*view(Ω,y,:)
                # Moments on rods x and y
                τ[x,:,threadid()] .+= (DRotation/kT).*((λ.*view(Ω,x,:))×view(rᵢⱼ,:,threadid()))×view(Ω,x,:)
                τ[y,:,threadid()] .-= (DRotation/kT).*((μ.*view(Ω,y,:))×view(rᵢⱼ,:,threadid()))×view(Ω,y,:)
            end
        end


        # ---- Covalent Forces ----

        # # Repeat forces in both directions within pair
        # for (i,j) in [(x,y),(y,x)]
        #     # N end to C adjacent internal point
        #     view(rᵢⱼ,:,threadid()) .= (r[j,:] .- Δu.*view(Ω,j,:)) .- (r[i,:] .+ 0.5*L.*view(Ω,i,:))
        #     rMag = sqrt(view(rᵢⱼ,:,threadid())⋅view(rᵢⱼ,:,threadid()))
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         view(rᵢⱼ,:,threadid()) .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,threadid()] .+= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,i,:,1)).*view(E,i,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,i,:,2)).*view(E,i,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,i,:)).*view(Ω,i,:)
        #         F[j,:,threadid()] .-= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,j,:,1)).*view(E,j,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,j,:,2)).*view(E,j,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,j,:)).*view(Ω,j,:)
        #         # Moment on rod i
        #         τ[i,:,threadid()] .+= (DRotation/kT).*((view(Ω,i,:).*L/2.0)×view(rᵢⱼ,:,threadid()))×view(Ω,i,:)
        #         τ[j,:,threadid()] .-= (DRotation/kT).*((-view(Ω,j,:).*Δu)×view(rᵢⱼ,:,threadid()))×view(Ω,j,:)
        #     end
        #     # C end to N adjacent internal point
        #     view(rᵢⱼ,:,threadid()) .= (r[j,:] .+ Δu.*view(Ω,j,:)) .- (r[i,:] .- 0.5*L.*view(Ω,i,:))
        #     rMag = sqrt(view(rᵢⱼ,:,threadid())⋅view(rᵢⱼ,:,threadid()))
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         view(rᵢⱼ,:,threadid()) .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,threadid()] .+= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,i,:,1)).*view(E,i,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,i,:,2)).*view(E,i,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,i,:)).*view(Ω,i,:)
        #         F[j,:,threadid()] .-= (DPerpendicular/kT)*((view(rᵢⱼ,:,threadid())⋅view(E,j,:,1)).*view(E,j,:,1) .+ (view(rᵢⱼ,:,threadid())⋅view(E,j,:,2)).*view(E,j,:,2)) .+ (DParallel/kT)*(view(rᵢⱼ,:,threadid())⋅view(Ω,j,:)).*view(Ω,j,:)
        #         # Moment on rod i
        #         τ[i,:,threadid()] .+= (DRotation/kT).*((-view(Ω,i,:).*L/2.0)×view(rᵢⱼ,:,threadid()))×view(Ω,i,:)
        #         τ[j,:,threadid()] .-= (DRotation/kT).*((view(Ω,j,:).*Δu)×view(rᵢⱼ,:,threadid()))×view(Ω,j,:)
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
