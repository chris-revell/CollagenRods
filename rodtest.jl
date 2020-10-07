using LinearAlgebra
using Random
using Distributions
using Plots
#using Threads

N       = 4      # Number of rods
L       = 1.0    # Rod length
σ       = 0.005  # Rod diameter
ϵ       = 0.000000001    # Hard core repulsion L-J potential depth
p       = L/σ    # Rod aspect ratio
η       = 1.0    # Solvent shear viscocity
kT      = 1.0    # Boltzman constant*Temperature
Δt      = 0.00000001 # Time step
tMax    = 1.0    # Simulation duration

# Diffusion constants
D₀ = kT/(η*L)
DParallel = D₀*(log(p)+-0.207+0.980/p-0.133/p^2)/2π # From Löwen Phys Rev E 1994
DPerpendicular = D₀*(log(p)+0.839+0.185/p+0.233/p^2)/4π
DRotation = 3*D₀*(log(p)-0.662+0.917/p-0.05/p^2)/(π*L^2)
# Standard deviations for Wiener process noise
stds = [sqrt(2.0*DParallel*Δt),sqrt(2.0*DPerpendicular*Δt),sqrt(2.0*DRotation*Δt)]

# Random number seeding
RNG = Random.MersenneTwister()

# Positions of all rod centre points
r = rand(Float64,N,3).*0.5.-0.25
# Orientations of all rods
Ω = rand(Float64,N,3).*2.0 .- 1.0    # Random initial orientations
Ω .= Ω./sqrt.(sum(Ω.^2,dims=2)) # Normalise magnitude
τ = zeros(Float64,N,3)  # Moments on each rod
F = zeros(Float64,N,3)  # Forces on each rod

# Initialise matrix for orthonormal bases
E = zeros(Float64,3,2)
# Dummy vector for later calculations
A = zeros(Float64,3)
B = zeros(Float64,3)


# Update rod positions and orientations according to established Brownian rod theory
function brownianMotion!(N,Ω,r,E,stds)
    for i=1:N
        # Create orthonormal basis vectors around rod axis
        E  .= nullspace(Matrix(Ω[i,:]'))
        # Update centre point position
        r[i,:] .+= rand(RNG,Normal(0.0,stds[1])).*Ω[i,:] .+ rand(RNG,Normal(0.0,stds[2])).*E[:,1] .+ rand(RNG,Normal(0.0,stds[2])).*E[:,2]
        # Update orientation
        Ω[i,:] .+= rand(RNG,Normal(0.0,stds[3])).*E[:,1] .+ rand(RNG,Normal(0.0,stds[3])).*E[:,2]
    end
    Ω .= Ω./sqrt.(sum(Ω.^2,dims=2)) # Normalise orientation magnitudes
end

# Plot all rods
function visualiseStep(N,r,Ω)
    plot()
    for i=1:N
        plot!([r[i,1]-Ω[i,1]./2.0,r[i,1]+Ω[i,1]./2.0],[r[i,2]-Ω[i,2]./2.0,r[i,2]+Ω[i,2]./2.0],[r[i,3]-Ω[i,3]./2.0,r[i,3]+Ω[i,3]./2.0],xlims=[-1,1],ylims=[-1,1],zlims=[-1,1],linecolor=:green,aspect_ratio=:equal,legend=:false)
    end
end


function shortestDistance!(r,Ω,r₁₂,i,j)
    r₁₂ .= r[i,:]-r[j,:]
    λᵢ = (r₁₂⋅Ω[i,:] - (Ω[i,:]⋅Ω[j,:])*(r₁₂⋅Ω[j,:]))/(1.0-(Ω[i,:]⋅Ω[j,:])^2)
    μⱼ = ((Ω[i,:]⋅Ω[j,:])*(r₁₂⋅Ω[i,:]) - r₁₂⋅Ω[j,:])/(1.0-(Ω[i,:]⋅Ω[j,:])^2)
    r₁₂ .+= μⱼ.*Ω[j,:] .- λᵢ.*Ω[i,:]
    return (μⱼ,λᵢ)
end

function repulsiveForces!(N,r,Ω,F,τ,rᵢⱼ,ϵ)

    for i=1:N
        Eᵢ = nullspace(Matrix(Ω[i,:]'))
        for j=i+1:N
            Eⱼ = nullspace(Matrix(Ω[i,:]'))

            # Find shortest distance between rod i and rod j
            (μ,λ) = shortestDistance!(r,Ω,rᵢⱼ,i,j)
            minSeparation = sqrt(rᵢⱼ⋅rᵢⱼ)

            # If shortest distance is less than rod radius, calculate hard core repulsion
            if minSeparation < σ/2.0
                F_mag  = 24.0*ϵ*((σ^6.0)/(minSeparation^7.0) - 2.0*(σ^12.0)/(minSeparation^13.0))
                rᵢⱼ .= (F_mag/minSeparation).*A

                # Linear forces acting on rods i and j using orthonormal basis vectors
                F[i,:] .+= (DPerpendicular/kT)*((rᵢⱼ⋅Eᵢ[:,1]).*Eᵢ[:,1] .+ (rᵢⱼ⋅Eᵢ[:,2]).*Eᵢ[:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[i,:]).*Ω[i,:]

                F[j,:] .-= (DPerpendicular/kT)*((rᵢⱼ⋅Eⱼ[:,1]).*Eⱼ[:,1] .+ (rᵢⱼ⋅Eⱼ[:,2]).*Eⱼ[:,2]) .+ (DParallel/kT)*(rᵢⱼ⋅Ω[j,:]).*Ω[j,:]

                # Moment on rod i
                τ[i,:] .+= (DRotation/kT).*((r[i,:].-λ.*Ω[i,:])×rᵢⱼ)×Ω[i,:]
                # Moment on rod j
                τ[j,:] .-= (DRotation/kT).*((r[j,:].-μ.*Ω[j,:])×rᵢⱼ)×Ω[j,:]
            end
        end
    end
end



@gif for t=0:Δt:0.01
    repulsiveForces!(N,r,Ω,F,τ,A,ϵ)
    brownianMotion!(N,Ω,r,E,stds)
    Ω .+= τ.*Δt
    r .+= F.*Δt
    if F != zeros(Float64,N,3)
        print(F)
    end
    τ .= 0.0
    F .= 0.0
    visualiseStep(N,r,Ω)
end every 100
