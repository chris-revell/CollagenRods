using LinearAlgebra
using Random
using Distributions
using Plots

include("./OrthonormalBases.jl")
using .OrthonormalBases
include("./RepulsiveForces.jl")
using .RepulsiveForces
include("./BrownianMotion.jl")
using .BrownianMotion
include("./VisualiseStep.jl")
using .VisualiseStep


function main(N,L,σ,ϵ,p,η,kT,Δt,tMax,boxSize)

    # Diffusion constants from Löwen Phys Rev E 1994
    D₀              = kT/(η*L)
    DParallel       = D₀*(log(p)+-0.207+0.980/p-0.133/p^2)/2π
    DPerpendicular  = D₀*(log(p)+0.839+0.185/p+0.233/p^2)/4π
    DRotation       = 3*D₀*(log(p)-0.662+0.917/p-0.05/p^2)/(π*L^2)
    # Standard deviations for Wiener process noise in parallel, perpendicular, and rotational directions
    stds = [sqrt(2.0*DParallel*Δt),sqrt(2.0*DPerpendicular*Δt),sqrt(2.0*DRotation*Δt)]

    # Random number seeding
    RNG = Random.MersenneTwister()

    r = (rand(Float64,N,3).-0.5).*boxSize # Random initial centrepoint positions of all rods
    Ω = rand(Float64,N,3).*2.0.-1.0       # Random initial orientations of all rods
    Ω .= Ω./sqrt.(sum(Ω.^2,dims=2))       # Normalise magnitude
    τ = zeros(Float64,N,3)                # Moments on each rod
    F = zeros(Float64,N,3)                # Forces on each rod
    ξr = zeros(Float64,N,3)               # Translational stochastic component
    ξΩ = zeros(Float64,N,3)               # Rotational stochastic component
    E = zeros(Float64,N,3,2)              # Matrix for orthonormal bases
    A = zeros(Float64,3)                  # Dummy vector for later calculations



    @gif for t=0:Δt:tMax

        orthonormalBases!(N,Ω,E)
        repulsiveForces!(N,r,Ω,F,τ,E,A,ϵ,σ,DParallel,DPerpendicular,DRotation,kT)
        brownianMotion!(N,Ω,ξr,ξΩ,E,stds,RNG)


        Ω .+= τ.*Δt .+ ξΩ
        Ω .= Ω./sqrt.(sum(Ω.^2,dims=2)) # Normalise magnitude

        r .+= F.*Δt .+ ξr

        τ .= 0.0
        F .= 0.0
        visualiseStep(N,r,Ω,L)
    end every 100

end

#%%

N       = 4           # Number of rods
L       = 1.0         # Rod length
σ       = 0.005       # Rod diameter
ϵ       = 0.000000001 # Hard core repulsion L-J potential depth
p       = L/σ         # Rod aspect ratio
η       = 1.0         # Solvent shear viscocity
kT      = 1.0         # Boltzman constant*Temperature
Δt      = 0.0000001 # Time step
tMax    = 0.1      # Simulation duration
boxSize = 0.5         # Dimensions of box in which rods are initialised

main(N,L,σ,ϵ,p,η,kT,Δt,tMax,boxSize)
