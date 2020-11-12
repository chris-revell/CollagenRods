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

@inline function interRodForces!(pairsList,N,r,Ω,F,τ,E,A,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q)

    #k₀ = 1.0
    #ϵ₀ = 1.0
    #Q  = 10.0

    @threads for (x,y) in pairsList

        rᵢⱼ = zeros(Float64,3)

        Δu = (1.8/4.4)*L

        # ---- RepulsiveForces ----

        # Find shortest distance between rod x and rod y
        (μ,λ) = shortestRodToRod!(r,Ω,rᵢⱼ,x,y,L)
        rMag = sqrt(rᵢⱼ⋅rᵢⱼ)

        # If shortest distance is less than rod radius, calculate hard core repulsion
        if rMag < σ
            FMag  = lennardJones(rMag,ϵ,σ)
            rᵢⱼ .*= FMag/rMag
            # Linear forces acting on rods x and y using orthonormal basis vectors
            F[x,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,x,:,1)).*view(E,x,:,1) .+ (rᵢⱼ⋅view(E,x,:,2)).*view(E,x,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,x,:)).*view(Ω,x,:)
            F[y,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,y,:,1)).*view(E,y,:,1) .+ (rᵢⱼ⋅view(E,y,:,2)).*view(E,y,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,y,:)).*view(Ω,y,:)
            # Moments on rods x and y
            τ[x,:,threadid()] .+= (DRotation/kT).*((λ.*view(Ω,x,:))×rᵢⱼ)×view(Ω,x,:)
            τ[y,:,threadid()] .-= (DRotation/kT).*((μ.*view(Ω,y,:))×rᵢⱼ)×view(Ω,y,:)
        end


        # ---- Electrostatic ----
        for i=-5:5
            μ = i*L/10.0
            for j=-5:5
                λ = i*L/10.0
                if rMag > σ
                    rᵢⱼ .= view(r,y,:) .- view(r,x,:) .+ μ.*view(Ω,y,:) .- λ.*view(Ω,x,:)
                    rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
                    FMag  = Q/(4*π*rMag^2)
                    rᵢⱼ .*= FMag/rMag
                    # Linear forces acting on rods x and y using orthonormal basis vectors
                    F[x,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,x,:,1)).*view(E,x,:,1) .+ (rᵢⱼ⋅view(E,x,:,2)).*view(E,x,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,x,:)).*view(Ω,x,:)
                    F[y,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,y,:,1)).*view(E,y,:,1) .+ (rᵢⱼ⋅view(E,y,:,2)).*view(E,y,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,y,:)).*view(Ω,y,:)
                    # Moments on rods x and y
                    τ[x,:,threadid()] .+= (DRotation/kT).*((λ.*view(Ω,x,:))×rᵢⱼ)×view(Ω,x,:)
                    τ[y,:,threadid()] .-= (DRotation/kT).*((μ.*view(Ω,y,:))×rᵢⱼ)×view(Ω,y,:)
                end
            end
        end



        # ---- Covalent Forces ----

        # # Repeat forces in both directions within pair
        # for (i,j) in [(x,y),(y,x)]
        #     # N end to C adjacent internal point
        #     rᵢⱼ .= (r[j,:] .- Δu.*view(Ω,j,:)) .- (r[i,:] .+ 0.5*L.*view(Ω,i,:))
        #     rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         rᵢⱼ .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,i,:,1)).*view(E,i,:,1) .+ (rᵢⱼ⋅view(E,i,:,2)).*view(E,i,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,i,:)).*view(Ω,i,:)
        #         F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,j,:,1)).*view(E,j,:,1) .+ (rᵢⱼ⋅view(E,j,:,2)).*view(E,j,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,j,:)).*view(Ω,j,:)
        #         # Moment on rod i
        #         τ[i,:,threadid()] .+= (DRotation/kT).*((view(Ω,i,:).*L/2.0)×rᵢⱼ)×view(Ω,i,:)
        #         τ[j,:,threadid()] .-= (DRotation/kT).*((-view(Ω,j,:).*Δu)×rᵢⱼ)×view(Ω,j,:)
        #     end
        #     # C end to N adjacent internal point
        #     rᵢⱼ .= (r[j,:] .+ Δu.*view(Ω,j,:)) .- (r[i,:] .- 0.5*L.*view(Ω,i,:))
        #     rMag = sqrt(rᵢⱼ⋅rᵢⱼ)
        #     if σ < rMag < 0.5*L
        #         FMag  = Q/(4.0*π*rMag^2) #(1.0+k₀*rMag)*(Q*exp(-k₀*rMag))/(4.0*π*ϵ₀*rMag^2)
        #         rᵢⱼ .*= FMag/rMag
        #         # Linear forces acting on rods i and j using orthonormal basis vectors
        #         F[i,:,threadid()] .+= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,i,:,1)).*view(E,i,:,1) .+ (rᵢⱼ⋅view(E,i,:,2)).*view(E,i,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,i,:)).*view(Ω,i,:)
        #         F[j,:,threadid()] .-= (DPerpendicular/kT)*((rᵢⱼ⋅view(E,j,:,1)).*view(E,j,:,1) .+ (rᵢⱼ⋅view(E,j,:,2)).*view(E,j,:,2)) .+ (DParallel/kT)*(rᵢⱼ⋅view(Ω,j,:)).*view(Ω,j,:)
        #         # Moment on rod i
        #         τ[i,:,threadid()] .+= (DRotation/kT).*((-view(Ω,i,:).*L/2.0)×rᵢⱼ)×view(Ω,i,:)
        #         τ[j,:,threadid()] .-= (DRotation/kT).*((view(Ω,j,:).*Δu)×rᵢⱼ)×view(Ω,j,:)
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
