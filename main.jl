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
N               = 20    # Number of rods
L               = 0.5   # Rod length
σ               = 0.01  # Rod diameter
ϵ               = 10.0  # Hard core repulsion L-J potential depth
Q               = 10.0  # Electrostatic charge
tMax            = 0.1   # Simulation duration
containerRadius = 0.3   # Radius of container in which rods are initialised
containerVolume = 1.0   # Volume of container in which rods are initialised
outputToggle    = 1     # Toggle controlling whether data are saved to file
renderToggle    = 1     # Toggle controlling whether data are rendered to images

# Run simulation
simulate(N,L,σ,ϵ,Q,tMax,containerRadius,containerVolume,outputToggle,renderToggle)
