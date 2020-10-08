using LinearAlgebra
using Random
using Distributions
using Plots
using Base.Threads

include("./OrthonormalBases.jl")
using .OrthonormalBases
include("./RepulsiveForces.jl")
using .RepulsiveForces
include("./AttractiveForces.jl")
using .AttractiveForces
include("./BrownianMotion.jl")
using .BrownianMotion
include("./VisualiseStep.jl")
using .VisualiseStep
include("./AdaptTimestep.jl")
using .AdaptTimestep

#%%

function main(N,L,σ,ϵ,p,η,kT,tMax,boxSize)

    # Diffusion constants from Löwen Phys Rev E 1994
    D₀              = kT/(η*L)
    DParallel       = D₀*(log(p)+-0.207+0.980/p-0.133/p^2)/2π
    DPerpendicular  = D₀*(log(p)+0.839+0.185/p+0.233/p^2)/4π
    DRotation       = 3*D₀*(log(p)-0.662+0.917/p-0.05/p^2)/(π*L^2)
    
    # Create random number generators for each thread
    threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())
    for i in 1:nthreads()
        threadRNG[i] = Random.MersenneTwister()
    end

    r = (rand(Float64,N,3).-0.5).*boxSize # Random initial centrepoint positions of all rods
    Ω = rand(Float64,N,3).*2.0.-1.0       # Random initial orientations of all rods
    Ω .= Ω./sqrt.(sum(Ω.^2,dims=2))       # Normalise magnitude
    τ = zeros(Float64,N,3,nthreads())     # Moments on each rod
    F = zeros(Float64,N,3,nthreads())     # Forces on each rod
    ξr = zeros(Float64,N,3)               # Translational stochastic component
    ξΩ = zeros(Float64,N,3)               # Rotational stochastic component
    E = zeros(Float64,N,3,2)              # Matrix for orthonormal bases
    #A = zeros(Float64,3,nthreads())      # Dummy vector for later calculations

    anim = Animation() # Animation object to which plots are appended

    t = 0.0 # Initialise system time

    # Iterate until reaching max run time.
    while t < tMax

        orthonormalBases!(N,Ω,E)
        repulsiveForces!(N,r,Ω,F,τ,E,ϵ,σ,DParallel,DPerpendicular,DRotation,kT)
        attractiveForces!(N,r,Ω,F,τ,E,DParallel,DPerpendicular,DRotation,kT,L)
        brownianMotion!(N,Ω,ξr,ξΩ,E,stds,threadRNG)

        F[:,:,1] = sum(F,dims=3)
        τ[:,:,1] = sum(τ,dims=3)
        Δt = adaptTimestep!(N,F,σ,D₀,kT)

        Ω .+= τ[:,:,1].*Δt .+ ξΩ
        Ω .= Ω./sqrt.(sum(Ω.^2,dims=2)) # Normalise magnitude
        r .+= F[:,:,1].*Δt .+ ξr

        t += Δt
        if t%(tMax/100.0) < Δt
            visualiseStep(N,r,Ω,L)
            frame(anim)
            println(t)
        end

        τ .= 0.0
        F .= 0.0

    end
    gif(anim, "test.gif", fps = 15);
end

#%%

N       = 10          # Number of rods
L       = 1.0         # Rod length
σ       = 0.005       # Rod diameter
ϵ       = 0.000000001 # Hard core repulsion L-J potential depth
p       = L/σ         # Rod aspect ratio
η       = 1.0         # Solvent shear viscocity
kT      = 1.0         # Boltzman constant*Temperature
#Δt      = 0.000000001 # Time step
tMax    = 0.1         # Simulation duration
boxSize = 0.5         # Dimensions of box in which rods are initialised

main(N,L,σ,ϵ,p,η,kT,tMax,boxSize)
