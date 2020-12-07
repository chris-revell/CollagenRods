#
#  runVisualise.jl
#  collagen-rods
#
#  Created by Christopher Revell on 01/12/2020.
#
#
# Independent script to run visualise() from scratch

"./" in LOAD_PATH ? nothing : push!(LOAD_PATH,"./")
using Visualise

visualise(ARGS[1])
