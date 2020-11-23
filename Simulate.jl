
module Simulate

# Julia packages
using LinearAlgebra
using Random
using Distributions
using Base.Threads
using StaticArrays
using TimerOutputs

# Local modules
using OrthonormalBases
using InterRodForces
using BrownianMotion
using AdaptTimestep
using CellListFunctions
using CreateRunDirectory
using OutputData
using Visualise

@inline @views function simulate(N,L,σ,ϵ,Q,tMax,boxSize,outputToggle)

    to = TimerOutput()

    η                 = 1.0   # Solvent shear viscocity
    kT                = 1.0   # Boltzman constant*Temperature
    # Diffusion constants from Löwen Phys Rev E 1994
    p                 = L/σ         # Rod aspect ratio
    D₀                = kT/(η*L)
    DParallel         = D₀*(log(p)+-0.207+0.980/p-0.133/p^2)/2π
    DPerpendicular    = D₀*(log(p)+0.839+0.185/p+0.233/p^2)/4π
    DRotation         = 3*D₀*(log(p)-0.662+0.917/p-0.05/p^2)/(π*L^2)
    interactionThresh = 1.3*L

    # Create random number generators for each thread
    threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())
    for i in 1:nthreads()
        threadRNG[i] = Random.MersenneTwister()
    end

    r  = SizedArray{Tuple{N,3}}((rand(Float64,N,3).-0.5).*boxSize)            # Random initial centrepoint positions of all rods
    Ω  = SizedArray{Tuple{N,3}}(rand(Float64,N,3).*2.0.-1.0)                  # Random initial orientations of all rods
    normalize!.(eachslice(Ω,dims=1))
    #Ω .= Ω./sqrt.(sum(Ω.^2,dims=2))                                           # Normalise magnitude
    τ  = SizedArray{Tuple{N,3,nthreads()}}(zeros(Float64,N,3,nthreads()))     # Moments on each rod
    F  = SizedArray{Tuple{N,3,nthreads()}}(zeros(Float64,N,3,nthreads()))     # Forces on each rod
    ξr = SizedArray{Tuple{N,3}}(zeros(Float64,N,3))                           # Translational stochastic component
    ξΩ = SizedArray{Tuple{N,3}}(zeros(Float64,N,3))                           # Rotational stochastic component
    E  = SizedArray{Tuple{N,3,2}}(zeros(Float64,N,3,2))                       # Matrix for orthonormal bases
    rᵢⱼ= SizedArray{Tuple{3,nthreads()}}(zeros(Float64,3,nthreads()))         # Matric of dummy vectors for later calculations
    pairsList = Tuple{Int64, Int64}[]                                         # Array storing tuple of particle interaction pairs eg pairsList[2]=(1,5) => 2nd element of array shows that particles 1 and 5 are in interaction range
    neighbourCells = MVector{13}(Vector{Tuple{Int64,Int64,Int64}}(undef, 13)) # Vector storing 13 neighbouring cells for a given cell

    dummyVectors = SizedArray{Tuple{3,3,nthreads()}}(zeros(Float64,3,3,nthreads()))

    if outputToggle==1
        foldername = createRunDirectory(N,L,σ,ϵ,p,η,kT,tMax,boxSize,D₀,DParallel,DPerpendicular,DRotation,interactionThresh)
        outfile = open("output/$(foldername)/output.txt","w")
        outputData(r,Ω,outfile,0,tMax)
    end

    # Iterate until max run time reached
    t = 0.0 # Initialise system time
    while t < tMax

        # Create list of particle interaction pairs based on cell lists algorithm
        #pairsList = findPairs!(N,r,interactionThresh,neighbourCells)
        @timeit to "findPairs" pairsList = findPairs!(N,r,interactionThresh,neighbourCells)

        # Find perpendicular basis vectors for each rod
        #orthonormalBases!(N,Ω,E)
        @timeit to "orthonormalBases" orthonormalBases!(N,Ω,E)

        # Calculate attractive and repulsive forces between rods
        #interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q)
        @timeit to "interRodForces" interRodForces!(pairsList,N,r,Ω,F,τ,E,rᵢⱼ,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,dummyVectors)

        # Calculate stochastic component of Langevin equation
        #brownianMotion!(N,Ω,ξr,ξΩ,E,DParallel,DPerpendicular,DRotation,threadRNG)
        @timeit to "brownianMotion" brownianMotion!(N,Ω,ξr,ξΩ,E,DParallel,DPerpendicular,DRotation,threadRNG)

        # Adapt timestep according to force magnitudes
        #Δt = adaptTimestep!(N,F,τ,ξr,ξΩ,σ,kT,L)
        @timeit to "adaptTimestep" Δt = adaptTimestep!(N,F,τ,ξr,ξΩ,σ,kT,L)

        # Forward Euler integration of overdamped Langevin equation for position and orientation, given drift and stochastic terms.
        #Ω .+= τ[:,:,1].*Δt .+ ξΩ.*sqrt(Δt)
        @timeit to "integrate" Ω .+= τ[:,:,1].*Δt .+ ξΩ.*sqrt(Δt)
            #Ω .= Ω./sqrt.(sum(Ω.^2,dims=2)) # Normalise magnitude
            #@timeit to "integrate" Ω .= Ω./sqrt.(sum(Ω.^2,dims=2)) # Normalise magnitude
        #normalize!.(eachslice(Ω,dims=1))
        @timeit to "integrate" normalize!.(eachslice(Ω,dims=1))

        #r .+= F[:,:,1].*Δt .+ ξr.*sqrt(Δt)
        @timeit to "integrate" r .+= F[:,:,1].*Δt .+ ξr.*sqrt(Δt)

        t += Δt

        if t%(tMax/100.0) < Δt && outputToggle==1
            outputData(r,Ω,outfile,t,tMax)
        end

        # Refresh force and moment arrays
        fill!(τ,0.0)
        fill!(F,0.0)
    end

    if outputToggle==1
        #run(`python3 visualise.py $("output/"*foldername)`)
        visualise("output/"*foldername,N,L,σ)
    end

    display(to)

end

export simulate

end
