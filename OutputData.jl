#
#  OutputData.jl
#  collagen-rods
#
#  Created by Christopher Revell on 15/10/2020.
#
#
# Function to save run data to file as a list of positions and orientations at each timestep.

module OutputData

using DelimitedFiles
using LinearAlgebra

@inline @views function outputData(r,Ω,outfile,t,tMax)

    writedlm(outfile,r,", ")
    writedlm(outfile,Ω,", ")
    flush(outfile)
    println("Simulating: $t/$tMax")

    return nothing

end

export outputData

end
