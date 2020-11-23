
push!(LOAD_PATH,"./")
using Simulate

N       = 4     # Number of rods
L       = 0.5   # Rod length
σ       = 0.005  # Rod diameter
ϵ       = 1.0   # Hard core repulsion L-J potential depth
Q       = 10.0
tMax    = 0.0001  # Simulation duration
boxSize = 2.0   # Dimensions of box in which rods are initialised

simulate(N,L,σ,ϵ,Q,tMax,boxSize,0)
