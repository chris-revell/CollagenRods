#
#  main.jl
#  collagen-rods
#
#  Created by Christopher Revell on 23/11/2020.
#
#

# Load simulate function
"./" in LOAD_PATH ? nothing : push!(LOAD_PATH,"./")
using Simulate

# Define run parameters
N       = 96     # Number of rods
L       = 1.0   # Rod length
σ       = 0.005  # Rod diameter
ϵ       = 1.0   # Hard core repulsion L-J potential depth
Q       = 1.0
tMax    = 0.0001  # Simulation duration
boxSize = 2.0   # Dimensions of box in which rods are initialised

# Run simulation
simulate(N,L,σ,ϵ,Q,tMax,boxSize,0)
