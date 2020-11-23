#
#  mail.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 05/10/2020.
#
#

module Simulate

# Load official packages
using LinearAlgebra
using Random
using Distributions
using Base.Threads
using StaticArrays

using TimerOutputs

# Load local modules
using OrthonormalBases
using InterRodForces
using BrownianMotion
using AdaptTimestep
using CellListFunctions
using OutputData
using Initialise
using CreateRunDirectory
using Visualise

@inline @views function simulate(N::Int64,L::Float64,σ::Float64,ϵ::Float64,Q::Float64,tMax::Float64,boxSize::Float64,outputToggle::Int64)

    to = TimerOutput()

    η               = 1.0   # Solvent shear viscocity
    kT              = 1.0   # Boltzman constant*Temperature
    # Diffusion constants from Löwen Phys Rev E 1994
    p               = L/σ         # Rod aspect ratio
    D₀              = kT/(η*L)
    DParallel       = D₀*(log(p)+-0.207+0.980/p-0.133/p^2)/2π
    DPerpendicular  = D₀*(log(p)+0.839+0.185/p+0.233/p^2)/4π
    DRotation       = 3*D₀*(log(p)-0.662+0.917/p-0.05/p^2)/(π*L^2)
    interactionThresh = 1.3*L

    r  = SizedVector{N}(fill(SVector{3}(zeros(Float64,3)),N))                 # Centrepoint positions of all rods as a vector of static 3-vectors
    Ω  = SizedVector{N}(fill(SVector{3}(zeros(Float64,3)),N))                 # Orientations of all rods as a vector of static 3-vectors
    τ  = SizedArray{Tuple{N,3,nthreads()}}(zeros(Float64,N,3,nthreads()))     # Moments on each rod
    F  = SizedArray{Tuple{N,3,nthreads()}}(zeros(Float64,N,3,nthreads()))     # Forces on each rod
    ξr = SizedVector{N}(fill(SVector{3}(zeros(Float64,3)),N))                 # Translational stochastic component
    ξΩ = SizedVector{N}(fill(SVector{3}(zeros(Float64,3)),N))                 # Rotational stochastic component
    E  = SizedArray{Tuple{N,3,2}}(zeros(Float64,N,3,2))                       # Matrix for orthonormal bases
    pairsList = Tuple{Int64, Int64}[]                                         # Array storing tuple of particle interaction pairs eg pairsList[2]=(1,5) => 2nd element of array shows that particles 1 and 5 are in interaction range
    neighbourCells = MVector{13}(Vector{Tuple{Int64,Int64,Int64}}(undef, 13)) # Vector storing 13 neighbouring cells for a given cell
    threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())             # Vector of random number generators for each thread
    rᵢⱼ = MMatrix{3,nthreads()}(zeros(Float64,3,nthreads()))                  # Matric of dummy vectors for later calculations
    dummyVectors = SizedMatrix{nthreads(),5}(fill(SVector{3}(zeros(Float64,3)),nthreads(),5))

    initialise!(N,r,Ω,boxSize,threadRNG,nthreads(),outputToggle)

    # Set up output folder
    if outputToggle==1
        foldername = createRunDirectory(N,L,σ,ϵ,p,η,kT,tMax,boxSize,D₀,DParallel,DPerpendicular,DRotation,interactionThresh)
        outfile = open("output/$(foldername)/output.txt","w")
        outputData(r,Ω,outfile,0,tMax)
    end

    # Iterate until max run time reached
    t = 0.0 # Initialise system time
    while t < tMax

        # Create list of particle interaction pairs based on cell lists algorithm
        @timeit to "findPairs" pairsList = findPairs!(N,r,interactionThresh,neighbourCells)

        # Find perpendicular basis vectors for each rod
        @timeit to "orthonormalBases" orthonormalBases!(N,Ω,E)

        # Calculate attractive and repulsive forces between rods
        @timeit to "interRodForces" interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,dummyVectors)

        # Calculate stochastic component of Langevin equation
        @timeit to "brownianMotion" brownianMotion!(N,Ω,ξr,ξΩ,E,DParallel,DPerpendicular,DRotation,threadRNG)

        # Adapt timestep according to force magnitudes
        @timeit to "adaptTimestep" Δt = adaptTimestep!(N,F,τ,ξr,ξΩ,σ,kT,L)

        # Forward Euler integration of overdamped Langevin equation for position and orientation, given drift and stochastic terms.
        @timeit to "integrate" for i=1:N # @threads?
            Ω[i] = Ω[i] .+ τ[i,:,1].*Δt .+ Ω[i].*sqrt(Δt)
            Ω[i] = Ω[i]./sqrt(Ω[i]⋅Ω[i])
            r[i] = r[i] .+ F[i,:,1].*Δt .+ ξr[i].*sqrt(Δt)
        end
        t += Δt

        if t%(tMax/100.0) < Δt && outputToggle==1
            outputData(r,Ω,outfile,t,tMax)
        end

        # Refresh force and moment arrays
        fill!(τ,0.0)
        fill!(F,0.0)

    end
    if outputToggle==1
        visualise("output/"*foldername,N,L,σ)
    end

    display(to)

end

export simulate

end