#
#  LennardJones.jl
#  collagen-rods
#
#  Created by Christopher Revell on 30/10/2020.
#
#

module LennardJones

using DelimitedFiles
using LinearAlgebra

@inline @views function lennardJones(r,ϵ,σ)

    F = 24.0*ϵ*((σ^6.0)/(r^7.0) - 2.0*(σ^12.0)/(r^13.0))

    return F

end

export lennardJones

end
