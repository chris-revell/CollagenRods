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
Qₑ              = 10.0  # Electrostatic charge
Qcov            = 100.0 # Charge used for covalent interactions
tMax            = 0.1   # Simulation duration
cellListThresh  = 1.3*L # Grid element size for cell list
electroThresh   = 100.0*σ
covalentThresh  = 3.0*σ
containerRadius = 0.3   # Radius of container in which rods are initialised
containerVolume = 1.0   # Volume of container in which rods are initialised
outputToggle    = 1     # Toggle controlling whether data are saved to file
renderToggle    = 1     # Toggle controlling whether data are rendered to images

# Run simulation
simulate(N,L,σ,ϵ,Qₑ,Qcov,tMax,electroThresh,covalentThresh,cellListThresh,containerRadius,containerVolume,outputToggle,renderToggle)
