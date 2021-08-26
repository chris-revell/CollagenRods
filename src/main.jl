#
#  main.jl
#  CollagenRods
#
#  Created by Christopher Revell on 23/11/2020.
#
#

# Load simulate function
"./" in LOAD_PATH ? nothing : push!(LOAD_PATH,"./")
using Simulate

# Define run parameters
N       = 56    # Number of rods
L       = 1.0   # Rod length
σ       = 0.01  # Rod diameter
ϵ       = 10.0  # Hard core repulsion L-J potential depth
Q       = 10.0  #
tMax    = 1.0   # Simulation duration
boxSize = 4.0   # Dimensions of box in which rods are initialised

# Run simulation
simulate(N,L,σ,ϵ,Q,tMax,boxSize,1)
