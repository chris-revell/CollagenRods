#
#  Simulate.jl
#  collagen-rods
#
#  Created by Christopher Revell on 23/11/2020.
#
#
# Function to combine all other functions and run a full simulation given some input parameters.

module Simulate

# Julia packages
using LinearAlgebra
using Random
using Distributions
using Base.Threads
using StaticArrays

# Local modules
using OrthonormalBases
using InterRodForces
using BrownianMotion
using AdaptTimestep
using CellListFunctions
using CreateRunDirectory
using OutputData
using Visualise
using Initialise
using Integrate

@inline @views function simulate(N,L,σ,ϵ,Q,tMax,containerRadius,containerVolume,outputToggle,renderToggle)

    η                 = 1.0   # Solvent shear viscocity
    kT                = 1.0   # Boltzman constant*Temperature

    # Interaction threshold for cell list grid
    interactionThresh = 1.3*L

    r         = SizedArray{Tuple{N,3}}(zeros(Float64,N,3))                      # Random initial centrepoint positions of all rods
    Ω         = SizedArray{Tuple{N,3}}(zeros(Float64,N,3))                      # Random initial orientations of all rods
    τ         = SizedArray{Tuple{N,3,nthreads()}}(zeros(Float64,N,3,nthreads()))# Moments on each rod
    F         = SizedArray{Tuple{N,3,nthreads()}}(zeros(Float64,N,3,nthreads()))# Forces on each rod
    ξr        = SizedArray{Tuple{N,3}}(zeros(Float64,N,3))                      # Translational stochastic component
    ξΩ        = SizedArray{Tuple{N,3}}(zeros(Float64,N,3))                      # Rotational stochastic component
    E         = SizedArray{Tuple{N,2,3}}(zeros(Float64,N,2,3))                  # Matrix for orthonormal bases
    rᵢⱼ       = SizedArray{Tuple{3,nthreads()}}(zeros(Float64,3,nthreads()))    # Matrix of dummy vectors for later calculations
    # TODO Set fixed size for pairsList
    pairsList          = Tuple{Int64, Int64}[]                                           # Array storing tuple of particle interaction pairs eg pairsList[2]=(1,5) => 2nd element of array shows that particles 1 and 5 are in interaction range
    neighbourCells     = MVector{13}(Vector{Tuple{Int64,Int64,Int64}}(undef, 13))        # Vector storing 13 neighbouring cells for a given cell
    dummyVectors       = SizedArray{Tuple{3,3,nthreads()}}(zeros(Float64,3,3,nthreads()))# Array of vectors to reuse in calculations and avoid allocations
    xVector            = SVector{3}([1.0,0.0,0.0])                                       # Vector along x axis for orthonormal basis calculations
    allRands           = MVector{5}(zeros(Float64,5))
    electrostaticPairs = SVector{8}([(1,3),(2,4),(3,5),(4,6),(5,7),(6,8),(7,9),(9,1)])

    p,D₀,DParallel,DPerpendicular,DRotation,foldername,outfile,threadRNG = initialise!(L,σ,ϵ,kT,Q,η,N,tMax,containerRadius,containerVolume,interactionThresh,r,Ω,outputToggle)

    # Iterate until max run time reached
    t = 0.0 # Initialise system time
    while t < tMax

        # Create list of particle interaction pairs based on cell lists algorithm
        pairsList = findPairs!(N,r,interactionThresh,neighbourCells,rᵢⱼ[:,1])

        # Find perpendicular basis vectors for each rod
        orthonormalBases!(N,Ω,E,xVector)

        # Calculate attractive and repulsive forces between rods
        interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,dummyVectors,electrostaticPairs)

        # Calculate stochastic component of Langevin equation
        brownianMotion!(N,Ω,ξr,ξΩ,E,DParallel,DPerpendicular,DRotation,threadRNG,allRands)

        # Adapt timestep according to force magnitudes
        Δt = adaptTimestep!(N,F,τ,ξr,ξΩ,σ,kT,L)

        # Forward Euler integration of overdamped Langevin equation for position and orientation, given drift and stochastic terms.
        t = integrate!(r,Ω,F,τ,ξr,ξΩ,Δt)

        if t%(tMax/100.0) < Δt && outputToggle==1
            outputData(r,Ω,outfile,t,tMax)
        end

        # Refresh force and moment arrays
        fill!(τ,0.0)
        fill!(F,0.0)
    end

    outputToggle==1 && renderToggle==1 ? visualise("output/"*foldername) : nothing

end

export simulate

end
