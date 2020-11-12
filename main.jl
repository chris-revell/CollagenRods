using LinearAlgebra
using Random
using Distributions
#using Plots
using Base.Threads

include("./OrthonormalBases.jl")
using .OrthonormalBases
include("./InterRodForces.jl")
using .InterRodForces
include("./BrownianMotion.jl")
using .BrownianMotion
include("./VisualiseStep.jl")
using .VisualiseStep
include("./AdaptTimestep.jl")
using .AdaptTimestep
include("./CellListFunctions.jl")
using .CellListFunctions
include("./CreateRunDirectory.jl")
using .CreateRunDirectory
include("./OutputData.jl")
using .OutputData

#%%

function main(N,L,σ,ϵ,Q,tMax,boxSize,outputToggle)

    η               = 1.0   # Solvent shear viscocity
    kT              = 1.0   # Boltzman constant*Temperature
    # Diffusion constants from Löwen Phys Rev E 1994
    p               = L/σ         # Rod aspect ratio
    D₀              = kT/(η*L)
    DParallel       = D₀*(log(p)+-0.207+0.980/p-0.133/p^2)/2π
    DPerpendicular  = D₀*(log(p)+0.839+0.185/p+0.233/p^2)/4π
    DRotation       = 3*D₀*(log(p)-0.662+0.917/p-0.05/p^2)/(π*L^2)
    interactionThresh = 1.3*L

    realTimePlotToggle = 1

    # Create random number generators for each thread
    threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())
    for i in 1:nthreads()
        threadRNG[i] = Random.MersenneTwister()
    end

    r  = (rand(Float64,N,3).-0.5).*boxSize # Random initial centrepoint positions of all rods
    Ω  = rand(Float64,N,3).*2.0.-1.0       # Random initial orientations of all rods
    Ω .= Ω./sqrt.(sum(Ω.^2,dims=2))        # Normalise magnitude
    τ  = zeros(Float64,N,3,nthreads())     # Moments on each rod
    F  = zeros(Float64,N,3,nthreads())     # Forces on each rod
    ξr = zeros(Float64,N,3)                # Translational stochastic component
    ξΩ = zeros(Float64,N,3)                # Rotational stochastic component
    E  = zeros(Float64,N,3,2)              # Matrix for orthonormal bases
    A  = zeros(Float64,3,nthreads())       # Matric of dummy vectors for later calculations
    rᵢⱼ = zeros(Float64,3,nthreads())
    pairsList = Tuple{Int64, Int64}[]      # Array storing tuple of particle interaction pairs eg pairsList[2]=(1,5) => 2nd element of array shows that particles 1 and 5 are in interaction range
    neighbourCells = Vector{Tuple{Int64,Int64,Int64}}(undef, 13) # Vector storing 13 neighbouring cells for a given cell

    if outputToggle==1
        foldername = createRunDirectory(N,L,σ,ϵ,p,η,kT,tMax,boxSize,D₀,DParallel,DPerpendicular,DRotation,interactionThresh)
        outfile = open("output/$(foldername)/output.txt","w")
        outputData(r,Ω,outfile,0,tMax)
        if realTimePlotToggle = 1
            run(`python3 visualiseSingle.py $("output/"*foldername)`)
        end
    end

    # Iterate until max run time reached
    t = 0.0 # Initialise system time
    while t < tMax

        # Create list of particle interaction pairs based on cell lists algorithm
        pairsList = findPairs!(N,r,interactionThresh,neighbourCells)

        # Find perpendicular basis vectors for each rod
        orthonormalBases!(N,Ω,E)

        # Calculate attractive and repulsive forces between rods
        interRodForces!(pairsList,N,r,Ω,F,τ,E,A,DParallel,DPerpendicular,DRotation,kT,L,ϵ,σ,Q,rᵢⱼ)

        # Calculate stochastic component of Langevin equation
        brownianMotion!(N,Ω,ξr,ξΩ,E,DParallel,DPerpendicular,DRotation,threadRNG)

        # Adapt timestep according to force magnitudes
        Δt = adaptTimestep!(N,F,τ,ξr,ξΩ,σ,kT,L)

        # Forward Euler integration of overdamped Langevin equation for position and orientation, given drift and stochastic terms.
        Ω .+= view(τ,:,:,1).*Δt .+ ξΩ.*sqrt(Δt)
        Ω .= Ω./sqrt.(sum(Ω.^2,dims=2)) # Normalise magnitude
        r .+= view(F,:,:,1).*Δt .+ ξr.*sqrt(Δt)

        t += Δt

        if t%(tMax/100.0) < Δt && outputToggle==1
            outputData(r,Ω,outfile,t,tMax)
            if realTimePlotToggle = 1
                run(`python3 visualiseSingle.py $("output/"*foldername)`)
            end 
        end

        # Refresh force and moment arrays
        τ .= 0.0
        F .= 0.0
    end
    if outputToggle==1
        if realTimePlotToggle != 1
            run(`python3 visualise.py $("output/"*foldername)`)
        end
    end
end

#%%

const N       = 10     # Number of rods
const L       = 0.5   # Rod length
const σ       = 0.005  # Rod diameter
const ϵ       = 1000.0   # Hard core repulsion L-J potential depth
const Q       = 10.0
const tMax    = 0.0005  # Simulation duration
const boxSize = 1.0   # Dimensions of box in which rods are initialised

main(2,L,L/5,ϵ/100.0,Q,tMax/100.0,boxSize,0)

main(N,L,σ,ϵ,Q,tMax,boxSize,1)
