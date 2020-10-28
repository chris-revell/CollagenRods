#
#  OutputData.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 15/10/2020.
#
#

module OutputData

using DelimitedFiles
using LinearAlgebra

@inline function outputData(r,Ω,outfile,t,tMax)

    writedlm(outfile,r,", ")
    writedlm(outfile,Ω,", ")
    flush(outfile)
    println("Simulating: $t/$tMax")

    return nothing

end

export outputData

end
