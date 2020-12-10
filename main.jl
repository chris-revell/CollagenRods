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
N       = 56    # Number of rods
L       = 1.0   # Rod length
σ       = 0.01  # Rod diameter
ϵ       = 10.0  # Hard core repulsion L-J potential depth
Q       = 10.0  # Electrostatic charge
tMax    = 0.0001# Simulation duration
contVol = 256.0 # Volume of container in which rods are initialised
contRad = 2.0   # Radius of container in which rods are initialised
outTgle = 1     # Toggle controlling whether data are saved to file
renTgle = 1     # Toggle controlling whether data are rendered to images

# Run simulation
simulate(N,L,σ,ϵ,Q,tMax,contRad,contVol,outTgle,renTgle)
