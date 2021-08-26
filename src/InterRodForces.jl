#
#  InterRodForces.jl
#  CollagenRods
#
#  Created by Christopher Revell on 15/10/2020.
#
#
# Function to calculate forces and moments between all particles across multiple threads.

module InterRodForces

using LinearAlgebra
using ShortestDistance
using LennardJones
using RepulsiveForces
using CovalentForces
using ElectrostaticForces
using Base.Threads

@inline @views function interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,dummyVectors,electrostaticPairs)

    Δu            = (1.8/4.4)*L
    colouredWidth = L/10
    blackWidth    = L/(2*8)

    covalentThresh      = 3.0*σ
    electrostaticThresh = 100.0*σ

    @threads for (x,y) in pairsList

        tID = threadid()

        # ---- RepulsiveForces ----
        repulsiveForces!(N,r,Ω,E,F,τ,rᵢⱼ,x,y,L,dummyVectors,tID,DParallel,DPerpendicular,DRotation,kT,ϵ,σ)


        # ---- Electrostatic ----
        electrostaticForces!(N,r,Ω,E,F,τ,rᵢⱼ,x,y,dummyVectors,tID,DParallel,DPerpendicular,DRotation,kT,ϵ,σ,electrostaticThresh,electrostaticPairs,blackWidth,colouredWidth,Q)


        # ---- Covalent Forces ----
        covalentForces!(N,r,Ω,E,F,τ,rᵢⱼ,x,y,L,dummyVectors,tID,DParallel,DPerpendicular,DRotation,kT,ϵ,σ,Δu,covalentThresh,Q)

    end

    # Sum forces calculated in each thread
    sum!(F[:,:,1],F)
    sum!(τ[:,:,1],τ)

    return nothing
end

export interRodForces!

end
